package ring

import (
	"math/big"
	"math/bits"
)

// DivFloorByLastModulusNTT divides floor the polynomial by its last modulus. Input must be in the NTT domain.
func (context *Context) DivFloorByLastModulusNTT(p0 *Poly) {

	level := len(p0.Coeffs) - 1

	pTmp := make([]uint64, context.N)

	InvNTT(p0.Coeffs[level], p0.Coeffs[level], context.N, context.nttPsiInv[level], context.nttNInv[level], context.Modulus[level], context.mredParams[level])

	for i := 0; i < level; i++ {

		NTT(p0.Coeffs[level], pTmp, context.N, context.nttPsi[i], context.Modulus[i], context.mredParams[i], context.bredParams[i])

		p0tmp := p0.Coeffs[i]

		qi := context.Modulus[i]
		mredParams := context.mredParams[i]
		rescalParams := context.rescaleParams[level-1][i]

		// (x[i] - x[-1]) * InvQ
		for j := uint64(0); j < context.N; j++ {
			p0tmp[j] = MRed(p0tmp[j]+(qi-pTmp[j]), rescalParams, qi, mredParams)
		}
	}

	p0.Coeffs = p0.Coeffs[:level]
}

// DivFloorByLastModulus divides floor the polynomial by its last modulus.
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

// DivFloorByLastModulusManyNTT divides floor sequentially nbRescales times the polynmial by its last modulus. Input must be in the NTT domain.
func (context *Context) DivFloorByLastModulusManyNTT(p0 *Poly, nbRescales uint64) {
	context.InvNTTLvl(uint64(len(p0.Coeffs)-1), p0, p0)
	context.DivFloorByLastModulusMany(p0, nbRescales)
	context.NTTLvl(uint64(len(p0.Coeffs)-1), p0, p0)
}

// DivFloorByLastModulusMany divides floor sequentially nbRescales times the polynmial by its last modulus.
func (context *Context) DivFloorByLastModulusMany(p0 *Poly, nbRescales uint64) {
	for k := uint64(0); k < nbRescales; k++ {
		context.DivFloorByLastModulus(p0)
	}
}

// DivRoundByLastModulusNTT divides round the polynomial by its last modulus. Input must be in the NTT domain.
func (context *Context) DivRoundByLastModulusNTT(p0 *Poly) {

	var pHalf, pHalfNegQi uint64

	level := len(p0.Coeffs) - 1

	pTmp := make([]uint64, context.N)

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
			pTmp[j] = p0tmp[j] + pHalfNegQi
		}

		NTT(pTmp, pTmp, context.N, context.nttPsi[i], qi, mredParams, bredParams)

		// (x[i] - x[-1]) * InvQ
		for j := uint64(0); j < context.N; j++ {
			p1tmp[j] = MRed(p1tmp[j]+(qi-pTmp[j]), rescaleParams, qi, mredParams)
		}
	}

	p0.Coeffs = p0.Coeffs[:level]
}

// DivRoundByLastModulus divides round the polynomial by its last modulus. Input must be in the NTT domain.
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

// DivRoundByLastModulusManyNTT divides round sequentially nbRescales times the polynmial by its last modulus. Input must be in the NTT domain.
func (context *Context) DivRoundByLastModulusManyNTT(p0 *Poly, nbRescales uint64) {
	context.InvNTTLvl(uint64(len(p0.Coeffs)-1), p0, p0)
	context.DivRoundByLastModulusMany(p0, nbRescales)
	context.NTTLvl(uint64(len(p0.Coeffs)-1), p0, p0)
}

// DivRoundByLastModulusMany divides round sequentially nbRescales times the polynmial by its last modulus.
func (context *Context) DivRoundByLastModulusMany(p0 *Poly, nbRescales uint64) {
	for k := uint64(0); k < nbRescales; k++ {
		context.DivRoundByLastModulus(p0)
	}
}

// SimpleScaler is the structure storing the parameters to reconstruct a polynomial, scale it by t/Q and return the results modulo t.
// Used during the BFV decoding.
// Algorithm from https://eprint.iacr.org/2018/117.pdf.
type SimpleScaler struct {
	context *Context

	t  uint64   //Plaintext modulus
	wi []uint64 //Integer parts of   ([Q/Qi]^(-1))_{Qi} * t/Qi
	ti []*big.Float

	reducealgoMul func(x, y uint64) uint64
	reducealgoAdd func(x uint64) uint64

	reducealgoAddParam uint64
	reducealgoMulParam uint64
}

// SimpleScalerFloatPrecision is the precision in bits for the big.Float in the scaling by t/Q.
const SimpleScalerFloatPrecision = 80

// NewSimpleScaler creates a new SimpleScaler from t (the modulus under which the reconstruction is returned) and
// context (the context in which the polynomial to reconstruct will be represented).
// Used during the BFV decoding.
// Algorithm from https://eprint.iacr.org/2018/117.pdf.
func NewSimpleScaler(t uint64, context *Context) (newParams *SimpleScaler) {

	newParams = new(SimpleScaler)

	var tmp *big.Float
	QiB := new(big.Int)         // Qi
	QiBF := new(big.Float)      // Qi
	QiStar := new(big.Int)      // Q/Qi
	QiBarre := new(big.Int)     // (Q/Qi)^(-1) mod Qi
	QiBarreBF := new(big.Float) // (Q/Qi)^(-1) mod Qi
	tBF := new(big.Float).SetUint64(t)

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

		// Else (we can use Montgomery reduction)
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
	newParams.ti = make([]*big.Float, len(context.Modulus))

	for i, qi := range context.Modulus {

		QiB.SetUint64(qi)
		QiBF.SetUint64(qi)
		QiStar.Quo(context.ModulusBigint, QiB)
		QiBarre.ModInverse(QiStar, QiB)
		QiBarre.Mod(QiBarre, QiB)
		QiBarreBF.SetInt(QiBarre)

		tmp = new(big.Float).Quo(tBF, QiBF)
		tmp.Mul(tmp, QiBarreBF)

		//floor( ([Q/Qi]^(-1))_{Qi} * t/Qi )
		newParams.wi[i], _ = tmp.Uint64()

		// If t is not a power of 2 converts in Montgomery form
		if (t&(t-1)) != 0 && t != 0 {
			newParams.wi[i] = MForm(newParams.wi[i], t, BRedParams(t))
		}

		QiBarre.Mul(QiBarre, NewUint(t))
		QiBarre.Mod(QiBarre, QiB)
		QiBarreBF.SetInt(QiBarre)

		newParams.ti[i] = new(big.Float).SetPrec(SimpleScalerFloatPrecision).Quo(QiBarreBF, QiBF)
	}

	return
}

// Scale returns the reconstruction of p1 scaled by a factor t/Q and mod t on the reciever p2.
func (parameters *SimpleScaler) Scale(p1, p2 *Poly) {

	for i := uint64(0); i < parameters.context.N; i++ {

		var a uint64
		var bBF big.Float
		var p1j big.Float
		for j := range parameters.context.Modulus {
			// round(xi*wi + xi*ti)%t
			a += parameters.reducealgoMul(parameters.wi[j], p1.Coeffs[j][i])

			p1j.SetPrec(SimpleScalerFloatPrecision).SetUint64(p1.Coeffs[j][i])
			bBF.SetPrec(SimpleScalerFloatPrecision).Add(&bBF, p1j.Mul(&p1j, parameters.ti[j]))
		}

		bBF.Add(&bBF, new(big.Float).SetFloat64(0.5))
		buint64, _ := bBF.Uint64()

		abf := a + buint64

		a = parameters.reducealgoAdd(abf)

		for j := 0; j < len(p2.Coeffs); j++ {
			p2.Coeffs[j][i] = a
		}
	}
}
