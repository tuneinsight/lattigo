package ring

import (
	"math"
	"math/bits"
)

//=========================================
//===== CRT SIMPLE SCAliNG PARAMETERS =====
//=========================================

// SimpleScaler is the structure storing the parameters to reconstruct a polynomial, scale it by t/Q and return the results modulo t.
// Algorithm from https://eprint.iacr.org/2018/117.pdf
type SimpleScaler struct {
	context *Context

	t  uint64     //Plaintext modulus
	wi []uint64   //Integer parts of   ([Q/Qi]^(-1))_{Qi} * t/Qi
	ti []Float128 //Fractional part of ([Q/Qi]^(-1))_{Qi} * t/Qi

	reducealgoMul func(x, y uint64) uint64
	reducealgoAdd func(x uint64) uint64

	reducealgoAddParam uint64
	reducealgoMulParam uint64
}

// NewSimpleScaler creates a new SimpleScaler from t (the modulus under which the reconstruction is returned) and context (the context in which the polynomial
// to reconstruct will be represented).
func NewSimpleScaler(t uint64, context *Context) (newParams *SimpleScaler) {

	newParams = new(SimpleScaler)

	var tmp Float128
	var QiB Int     // Qi
	var QiStar Int  // Q/Qi
	var QiBarre Int // (Q/Qi)^(-1) mod Qi

	newParams.t = t
	newParams.context = context

	// Assigns the correct reduction algorithm depending on the provided t
	// If t is a power of 2
	if (t&(t-1)) == 0 && t != 0 {

		newParams.reducealgoAddParam = t - 1
		newParams.reducealgoMulParam = t - 1

		newParams.reducealgoMul = func(x, y uint64) uint64 {
			return (x * y) & newParams.reducealgoAddParam
		}

		newParams.reducealgoAdd = func(x uint64) uint64 {
			return x & newParams.reducealgoMulParam
		}

		// Else (we can use montgomery reduction)
	} else {
		newParams.reducealgoAddParam = BRedParams(t)[0]
		newParams.reducealgoMulParam = MRedParams(t)

		newParams.reducealgoMul = func(x, y uint64) uint64 {
			ahi, alo := bits.Mul64(x, y)
			R := alo * newParams.reducealgoMulParam
			H, _ := bits.Mul64(R, newParams.t)
			r := ahi - H + newParams.t

			if r >= newParams.t {
				r -= newParams.t
			}

			return r
		}

		newParams.reducealgoAdd = func(x uint64) uint64 {

			s0, _ := bits.Mul64(x, newParams.reducealgoAddParam)

			r := x - s0*newParams.t

			if r >= newParams.t {
				r -= newParams.t
			}

			return r
		}
	}

	// Integer and rational part of QiBarre * t/Qi
	newParams.wi = make([]uint64, len(context.Modulus))
	newParams.ti = make([]Float128, len(context.Modulus))

	for i, qi := range context.Modulus {

		QiB.SetUint(qi)
		QiStar.Div(context.ModulusBigint, &QiB)
		QiBarre.Inv(&QiStar, &QiB)
		QiBarre.Mod(&QiBarre, &QiB)

		tmp = Float128Div(Float128SetUint53(t), Float128SetUint64(qi))

		tmp = Float128Mul(tmp, Float128SetUint64(QiBarre.Uint64()))

		//floor( ([Q/Qi]^(-1))_{Qi} * t/Qi )
		newParams.wi[i] = Float128ToUint53(tmp)

		// If t is not a power of 2 converts in montgomery form
		if (t&(t-1)) != 0 && t != 0 {
			newParams.wi[i] = MForm(newParams.wi[i], t, BRedParams(t))
		}

		QiBarre.Mul(&QiBarre, NewUint(t))
		QiBarre.Mod(&QiBarre, &QiB)

		newParams.ti[i] = Float128Div(Float128SetUint64(QiBarre.Uint64()), Float128SetUint64(qi)) //floor( ([Q/Qi]^(-1))_{Qi} * t/Qi ) - ( ([Q/Qi]^(-1))_{Qi} * t/Qi )
	}

	return
}

// Scale returns the reconstruction of p1 scaled by a factor t/Q and mod t on the reciever p2.
func (parameters *SimpleScaler) Scale(p1, p2 *Poly) {

	var a uint64
	var b Float128

	for i := uint64(0); i < parameters.context.N; i++ {

		a = 0
		b[0], b[1] = 0, 0

		for j := range parameters.context.Modulus {
			// round(xi*wi + xi*ti)%t
			a += parameters.reducealgoMul(parameters.wi[j], p1.Coeffs[j][i])

			b = Float128Add(b, Float128Mul(parameters.ti[j], Float128SetUint64(p1.Coeffs[j][i])))
		}

		a += Float128ToUint64(b)
		a = parameters.reducealgoAdd(a)

		for j := 0; j < len(p2.Coeffs); j++ {
			p2.Coeffs[j][i] = a
		}
	}
}

//==========================================
//===== CRT COMPLEX SCAliNG PARAMETERS =====
//==========================================

// ComplexScaler is the structure holding the parameters for the complex scaling, which is the operation of reducing a polynomial
// in basis Q + P to a polynomial in basis Q while at the same time rescaling it by a factor t/Q.
// Algorithm from https://eprint.iacr.org/2018/117.pdf
type ComplexScaler struct {

	// 1. General Parameters
	contextQ *Context
	contextP *Context

	// 2. Scaling Parameters

	// Scale factore
	t uint64
	// t*P
	tP *Int

	// Integer part of (t * P)/Qi (mod Qi)
	wiiMont [][]uint64
	// Rational part of (t * P)/Qi
	fi []Float128

	// (t * P * r)/Pj (mod Qi)
	sijMont [][]uint64

	// (t * P) mod Qi
	li [][]uint64

	// (Q * P * r)/Qi
	qpiMont []uint64
	// (Q * P * r)/Pj
	qpjMont []uint64

	// Qi and Pi in float128 type
	qiFloat128 []Float128
	pjFloat128 []Float128
}

// NewComplexScaler creates a new ComplexScaler from t and the provided contexts.
func NewComplexScaler(t uint64, contextQ, contextP *Context) (newParams *ComplexScaler) {

	newParams = new(ComplexScaler)

	newParams.contextQ = contextQ
	newParams.contextP = contextP

	newParams.t = t
	newParams.tP = NewInt(0).Mul(NewInt(int64(t)), contextP.ModulusBigint)

	newParams.wiiMont = make([][]uint64, len(contextQ.Modulus))
	newParams.sijMont = make([][]uint64, len(contextP.Modulus))

	newParams.li = make([][]uint64, len(contextQ.Modulus))
	newParams.fi = make([]Float128, len(contextQ.Modulus))

	newParams.qpiMont = make([]uint64, len(contextQ.Modulus))
	newParams.qpjMont = make([]uint64, len(contextP.Modulus))

	newParams.qiFloat128 = make([]Float128, len(contextQ.Modulus))
	newParams.pjFloat128 = make([]Float128, len(contextP.Modulus))

	var QiB Int
	var PjB Int
	var v uint64

	tmp := NewInt(1) // Temporary variable

	QPiStar := make([]Int, len(contextQ.Modulus))
	QPjStar := make([]Int, len(contextP.Modulus))

	// ModulusBigint Q*P
	QP := NewInt(1).Mul(contextQ.ModulusBigint, contextP.ModulusBigint)

	wi := make([]Int, len(contextQ.Modulus))

	for i, qi := range contextQ.Modulus {
		QiB.SetUint(qi)
		wi[i].Div(newParams.tP, &QiB)
	}

	for i, qi := range contextQ.Modulus {
		QiB.SetUint(qi)

		// Integer part of tp/q_i mod each Qi
		newParams.wiiMont[i] = make([]uint64, len(contextQ.Modulus))
		for j := 0; j < len(contextQ.Modulus); j++ {
			newParams.wiiMont[i][j] = MForm(tmp.Mod(&wi[j], &QiB).Uint64(), qi, contextQ.bredParams[i])
		}

		// Rational part of tp/q_i (mod each Qi)
		newParams.fi[i] = Float128Div(Float128SetUint64(tmp.Mod(newParams.tP, &QiB).Uint64()), Float128SetUint64(qi))

		// li = -(t*P) mod Qi
		newParams.li[i] = make([]uint64, len(contextQ.Modulus)+len(contextP.Modulus)+1)
		v = qi - tmp.Mod(newParams.tP, &QiB).Uint64()
		newParams.li[i][0] = 0
		for j := 1; j < (len(contextQ.Modulus) + len(contextP.Modulus) + 1); j++ {
			newParams.li[i][j] = CRed(newParams.li[i][j-1]+v, qi)
		}

		// (Q * P * r)/Qi (in Montgomery form)
		QPiStar[i].Div(QP, &QiB) // QP/Qi

		newParams.qpiMont[i] = MForm(tmp.Inv(&QPiStar[i], &QiB).Uint64(), qi, contextQ.bredParams[i])

		newParams.qiFloat128[i] = Float128SetUint64(qi)
	}

	for j, pj := range contextP.Modulus {

		// (t*P * r)/Pj mod Qi (in Montgomery form)
		newParams.sijMont[j] = make([]uint64, len(contextQ.Modulus))
		PjB.SetUint(pj)
		for i, qi := range contextQ.Modulus {
			QiB.SetUint(qi)
			tmp.Div(newParams.tP, &PjB)
			tmp.Mod(tmp, &QiB)
			newParams.sijMont[j][i] = MForm(tmp.Uint64(), qi, contextQ.bredParams[i])
		}

		// (Q * P * r)/Pj
		QPjStar[j].Div(QP, &PjB)
		newParams.qpjMont[j] = MForm(tmp.Inv(&QPjStar[j], &PjB).Uint64(), pj, contextP.bredParams[j])

		newParams.pjFloat128[j] = Float128SetUint64(pj)
	}

	return newParams
}

// Scale takes a polynomial in basis {Q0,Q1....Qi,P0,P1...Pj}, rescales it by a factor t/Q and returns the result in basis {Q0,Q1....Qi}.
func (parameters *ComplexScaler) Scale(p1, p2 *Poly) {

	var tmp, yjFLoat128 Float128
	var v uint64
	var aInt uint64
	var aFloat uint64

	yi := make([]uint64, len(parameters.contextQ.Modulus))
	yj := make([]uint64, len(parameters.contextP.Modulus))
	yiFloat128 := make([]Float128, len(parameters.contextQ.Modulus))

	// Given a polynomial represented in basis Q0, P1, P2
	//
	// [[a, b, c, d]_Q0
	//  [a, b, c, d]_P1
	//  [a, b, c, d]_P2]
	//
	// Loops over each column of coefficient (first a, then b...), and scales it by t/Q and return the result
	// in base Q0. Then discards the base P.
	for x := uint64(0); x < parameters.contextQ.N; x++ {

		tmp[0], tmp[1] = 0, 0

		for i, qi := range parameters.contextQ.Modulus {

			yi[i] = MRed(p1.Coeffs[i][x], parameters.qpiMont[i], qi, parameters.contextQ.mredParams[i])

			yiFloat128[i][0] = float64(yi[i] >> 12)
			yiFloat128[i][1] = float64(yi[i]&0xfff) / float64(4096)

			tmp = Float128Add(tmp, Float128Div(yiFloat128[i], parameters.qiFloat128[i]))
		}

		for j, pj := range parameters.contextP.Modulus {

			yj[j] = MRed(p1.Coeffs[j+len(parameters.contextQ.Modulus)][x], parameters.qpjMont[j], pj, parameters.contextP.mredParams[j])

			yjFLoat128[0] = float64(yj[j] >> 12)
			yjFLoat128[1] = float64(yj[j]&0xfff) / float64(4096)

			tmp = Float128Add(tmp, Float128Div(yjFLoat128, parameters.pjFloat128[j]))
		}

		v = uint64(math.Round(tmp[0]))

		tmp[0], tmp[1] = 0, 0

		for i := range parameters.contextQ.Modulus {
			tmp = Float128Add(tmp, Float128Mul(yiFloat128[i], parameters.fi[i]))
		}

		// Isolates the integer part of the float128 from the floatting point part, then adds the rounded floating point
		// part to the integer part (this prevents occasional rounding errors when the second floating element is negative)
		aFloat = uint64(tmp[0]*4096) + uint64(math.Round((tmp[0]*4096)-float64(uint64(tmp[0]*4096))+tmp[1]*4096))

		var reduce uint64

		for u, qi := range parameters.contextQ.Modulus {

			aInt = aFloat

			reduce = 1

			for i := range parameters.contextQ.Modulus {

				aInt += MRed(yi[i], parameters.wiiMont[u][i], qi, parameters.contextQ.mredParams[u])

				if reduce&7 == 7 {
					aInt = BRedAdd(aInt, qi, parameters.contextQ.bredParams[u])
				}
			}

			for j := range parameters.contextP.Modulus {

				aInt += MRed(yj[j], parameters.sijMont[j][u], qi, parameters.contextQ.mredParams[u])

				if reduce&7 == 7 {
					aInt = BRedAdd(aInt, qi, parameters.contextQ.bredParams[u])
				}
			}

			p2.Coeffs[u][x] = BRedAdd(aInt+parameters.li[u][v], qi, parameters.contextQ.bredParams[u])
		}
	}

	p2.Coeffs = p2.Coeffs[:len(parameters.contextQ.Modulus)]
}
