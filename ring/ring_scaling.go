package ring

import (
	"math/big"
	"math/bits"
)

func (context *Context) DivFloorByLastModulusNTT(p0 *Poly) {

	level := len(p0.Coeffs) - 1

	p_tmp := make([]uint64, context.N)

	InvNTT(p0.Coeffs[level], p0.Coeffs[level], context.N, context.nttPsiInv[level], context.nttNInv[level], context.Modulus[level], context.mredParams[level])

	for i := 0; i < level; i++ {

		NTT(p0.Coeffs[level], p_tmp, context.N, context.nttPsi[i], context.Modulus[i], context.mredParams[i], context.bredParams[i])

		p0tmp := p0.Coeffs[i]

		qi := context.Modulus[i]
		mredParams := context.mredParams[i]
		rescalParams := context.rescaleParams[level-1][i]

		// (x[i] - x[-1]) * InvQ
		for j := uint64(0); j < context.N; j++ {
			p0tmp[j] = MRed(p0tmp[j]+(qi-p_tmp[j]), rescalParams, qi, mredParams)
		}
	}

	p0.Coeffs = p0.Coeffs[:level]
}

func (context *Context) DivFloorByLastModulus(p0 *Poly) {

	level := len(p0.Coeffs) - 1

	for i := 0; i < level; i++ {
		p0tmp := p0.Coeffs[level]
		p1tmp := p0.Coeffs[i]
		qi := context.Modulus[i]
		bredParams := context.bredParams[i]
		mredParams := context.mredParams[i]
		rescaleParams := context.rescaleParams[level-1][i]
		// (x[i] - x[-1]) * InvQ
		for j := uint64(0); j < context.N; j++ {
			p1tmp[j] = MRed(p1tmp[j]+(qi-BRedAdd(p0tmp[j], qi, bredParams)), rescaleParams, qi, mredParams)
		}
	}

	p0.Coeffs = p0.Coeffs[:level]
}

func (context *Context) DivFloorByLastModulusManyNTT(p0 *Poly, nbRescales uint64) {
	context.InvNTTLvl(uint64(len(p0.Coeffs)-1), p0, p0)
	context.DivFloorByLastModulusMany(p0, nbRescales)
	context.NTTLvl(uint64(len(p0.Coeffs)-1), p0, p0)
}

func (context *Context) DivFloorByLastModulusMany(p0 *Poly, nbRescales uint64) {
	for k := uint64(0); k < nbRescales; k++ {
		context.DivFloorByLastModulus(p0)
	}
}

func (context *Context) DivRoundByLastModulusNTT(p0 *Poly) {

	var pHalf, pHalfNegQi uint64

	level := len(p0.Coeffs) - 1

	p_tmp := make([]uint64, context.N)

	InvNTT(p0.Coeffs[level], p0.Coeffs[level], context.N, context.nttPsiInv[level], context.nttNInv[level], context.Modulus[level], context.mredParams[level])

	// Centers by (p-1)/2
	pHalf = (context.Modulus[level] - 1) >> 1
	p0tmp := p0.Coeffs[level]
	pj := context.Modulus[level]
	for i := uint64(0); i < context.N; i++ {
		p0tmp[i] = CRed(p0tmp[i]+pHalf, pj)
	}

	for i := 0; i < level; i++ {

		p1tmp := p0.Coeffs[i]

		qi := context.Modulus[i]
		bredParams := context.bredParams[i]
		mredParams := context.mredParams[i]
		rescaleParams := context.rescaleParams[level-1][i]

		pHalfNegQi = context.Modulus[i] - BRedAdd(pHalf, qi, bredParams)

		for j := uint64(0); j < context.N; j++ {
			p_tmp[j] = p0tmp[j] + pHalfNegQi
		}

		NTT(p_tmp, p_tmp, context.N, context.nttPsi[i], qi, mredParams, bredParams)

		// (x[i] - x[-1]) * InvQ
		for j := uint64(0); j < context.N; j++ {
			p1tmp[j] = MRed(p1tmp[j]+(qi-p_tmp[j]), rescaleParams, qi, mredParams)
		}
	}

	p0.Coeffs = p0.Coeffs[:level]
}

func (context *Context) DivRoundByLastModulus(p0 *Poly) {

	var pHalf, pHalfNegQi uint64

	level := len(p0.Coeffs) - 1

	// Centers by (p-1)/2
	pHalf = (context.Modulus[level] - 1) >> 1
	p0tmp := p0.Coeffs[level]
	pj := context.Modulus[level]
	pHalf = (context.Modulus[level] - 1) >> 1
	for i := uint64(0); i < context.N; i++ {
		p0tmp[i] = CRed(p0tmp[i]+pHalf, pj)
	}

	for i := 0; i < level; i++ {

		p1tmp := p0.Coeffs[i]

		qi := context.Modulus[i]
		bredParams := context.bredParams[i]
		mredParams := context.mredParams[i]
		rescaleParams := context.rescaleParams[level-1][i]

		pHalfNegQi = context.Modulus[i] - BRedAdd(pHalf, qi, bredParams)

		// (x[i] - x[-1]) * InvQ
		for j := uint64(0); j < context.N; j++ {
			p1tmp[j] = MRed(p1tmp[j]+(qi-BRedAdd(p0tmp[j]+pHalfNegQi, qi, bredParams)), rescaleParams, qi, mredParams)
		}
	}

	p0.Coeffs = p0.Coeffs[:level]
}

func (context *Context) DivRoundByLastModulusManyNTT(p0 *Poly, nbRescales uint64) {
	context.InvNTTLvl(uint64(len(p0.Coeffs)-1), p0, p0)
	context.DivRoundByLastModulusMany(p0, nbRescales)
	context.NTTLvl(uint64(len(p0.Coeffs)-1), p0, p0)
}

func (context *Context) DivRoundByLastModulusMany(p0 *Poly, nbRescales uint64) {
	for k := uint64(0); k < nbRescales; k++ {
		context.DivRoundByLastModulus(p0)
	}
}

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
	QiB := new(big.Int)     // Qi
	QiStar := new(big.Int)  // Q/Qi
	QiBarre := new(big.Int) // (Q/Qi)^(-1) mod Qi

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

		QiB.SetUint64(qi)
		QiStar.Quo(context.ModulusBigint, QiB)
		QiBarre.ModInverse(QiStar, QiB)
		QiBarre.Mod(QiBarre, QiB)

		tmp = Float128Div(Float128SetUint53(t), Float128SetUint64(qi))

		tmp = Float128Mul(tmp, Float128SetUint64(QiBarre.Uint64()))

		//floor( ([Q/Qi]^(-1))_{Qi} * t/Qi )
		newParams.wi[i] = Float128ToUint53(tmp)

		// If t is not a power of 2 converts in montgomery form
		if (t&(t-1)) != 0 && t != 0 {
			newParams.wi[i] = MForm(newParams.wi[i], t, BRedParams(t))
		}

		QiBarre.Mul(QiBarre, NewUint(t))
		QiBarre.Mod(QiBarre, QiB)

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
