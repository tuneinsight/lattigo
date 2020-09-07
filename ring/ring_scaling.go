package ring

import (
	"math/big"
	"math/bits"
)

// Scaler is an interface that rescales polynomial coefficients by a fraction t/Q.
type Scaler interface {
	// DivByQOverTRounded returns p1 scaled by a factor t/Q and mod t on the reciever p2.
	DivByQOverTRounded(p1, p2 *Poly)
}

// RNSScaler implements the Scaler interface by operating scaling by t/Q in the RNS domain.
// This implementation of the Scaler interface should be favored instead of the SimpleScaler
// implementation.
type RNSScaler struct {
	ringQ           *Ring
	paramsQT        *modupParams
	modDownParamsQT uint64
	polypoolT       *Poly

	qHalf     *big.Int // (q-1)/2
	qHalfModT uint64   // (q-1)/2 mod t

	t    uint64
	qInv uint64 //(q mod t)^-1 mod t

	bredParamsT []uint64
	mredParamsT uint64

	paramsQP *modupParams
}

// NewRNSScaler creates a new SimpleScaler from t, the modulus under which the reconstruction is returned, and
// the ring in which the polynomial to reconstruct is represented.
func NewRNSScaler(t uint64, ringQ *Ring) (rnss *RNSScaler) {

	rnss = new(RNSScaler)

	rnss.ringQ = ringQ

	rnss.mredParamsT = MRedParams(t)

	rnss.polypoolT = NewPoly(ringQ.N, 1)

	rnss.t = t
	rnss.qHalf = new(big.Int)
	rnss.qInv = rnss.qHalf.Mod(ringQ.ModulusBigint, NewUint(t)).Uint64()
	rnss.qInv = ModExp(rnss.qInv, t-2, t)
	rnss.qInv = MForm(rnss.qInv, t, BRedParams(t))

	rnss.qHalf.Set(ringQ.ModulusBigint)
	rnss.qHalf.Rsh(rnss.qHalf, 1)
	rnss.qHalfModT = rnss.qHalf.Mod(rnss.qHalf, NewUint(t)).Uint64()

	rnss.qHalf.Set(ringQ.ModulusBigint)
	rnss.qHalf.Rsh(rnss.qHalf, 1)

	rnss.paramsQP = basisextenderparameters(ringQ.Modulus, []uint64{t})

	return
}

// DivByQOverTRounded returns p1 scaled by a factor t/Q and mod t on the reciever p2.
func (rnss *RNSScaler) DivByQOverTRounded(p1Q, p2T *Poly) {

	ringQ := rnss.ringQ

	T := rnss.t
	p2tmp := p2T.Coeffs[0]
	p3tmp := rnss.polypoolT.Coeffs[0]
	mredParams := rnss.mredParamsT
	qInv := rnss.qInv
	qHalfModT := rnss.qHalfModT

	// Multiplies P_{Q} by t and extends the basis from P_{Q} to t*(P_{Q}||P_{t})
	// Since the coefficients of P_{t} are multiplied by t, they are all zero,
	// hence the basis extension can be omited
	ringQ.MulScalar(p1Q, T, p1Q)

	// Centers  t*P_{Q} around (Q-1)/2 to round instead of floor during the division
	ringQ.AddScalarBigint(p1Q, rnss.qHalf, p1Q)

	// Extends the basis of (t*P_{Q} + (Q-1)/2) to (t*P_{t} + (Q-1)/2)
	modUpExact(p1Q.Coeffs, rnss.polypoolT.Coeffs, rnss.paramsQP)

	// Computes [Q^{-1} * (t*P_{t} -   (t*P_{Q} - ((Q-1)/2 mod t)))] mod t which returns round(t/Q * P_{Q}) mod t
	for j := uint64(0); j < ringQ.N; j++ {
		p2tmp[j] = MRed(qHalfModT+T-p3tmp[j], qInv, T, mredParams)
	}
}

// SimpleScaler implements the Scaler interface by operating RNS reconstruction and scaling by t/Q.
// This implementation of the Scaler interface is less efficient than the RNSScaler, but uses simple
// multiprecision arithmetic of the math/big package. The BFV implementation should use the RNSScaler.
type SimpleScaler struct {
	ring *Ring

	r, one *big.Int

	t   uint64 //Plaintext modulus
	tBI *big.Int
	wi  []uint64 //Integer parts of   ([Q/Qi]^(-1))_{Qi} * t/Qi
	ti  []*big.Float

	p1BI []*big.Int

	reducealgoMul func(x, y uint64) uint64
	reducealgoAdd func(x uint64) uint64

	reducealgoAddParam uint64
	reducealgoMulParam uint64
}

// SimpleScalerFloatPrecision is the precision in bits for the big.Float in the scaling by t/Q.
const SimpleScalerFloatPrecision = 80

// NewSimpleScaler creates a new SimpleScaler from t, the modulus under which the reconstruction is returned, and
// the ring in which the polynomial to reconstruct is represented.
func NewSimpleScaler(t uint64, r *Ring) (ss *SimpleScaler) {

	ss = new(SimpleScaler)

	ss.r = new(big.Int)
	ss.one = NewInt(1)

	var tmp *big.Float
	QiB := new(big.Int)         // Qi
	QiBF := new(big.Float)      // Qi
	QiStar := new(big.Int)      // Q/Qi
	QiBarre := new(big.Int)     // (Q/Qi)^(-1) mod Qi
	QiBarreBF := new(big.Float) // (Q/Qi)^(-1) mod Qi
	tBF := new(big.Float).SetUint64(t)

	ss.t = t
	ss.tBI = new(big.Int).SetUint64(t)
	ss.ring = r

	// Assigns the correct reduction algorithm depending on the provided t
	// If t is a power of 2
	if (t&(t-1)) == 0 && t != 0 {

		ss.reducealgoAddParam = t - 1
		ss.reducealgoMulParam = t - 1

		ss.reducealgoMul = func(x, y uint64) uint64 {
			return (x * y) & ss.reducealgoAddParam
		}

		ss.reducealgoAdd = func(x uint64) uint64 {
			return x & ss.reducealgoMulParam
		}

		// Else (we can use Montgomery reduction)
	} else {
		ss.reducealgoAddParam = BRedParams(t)[0]
		ss.reducealgoMulParam = MRedParams(t)

		ss.reducealgoMul = func(x, y uint64) uint64 {
			ahi, alo := bits.Mul64(x, y)
			R := alo * ss.reducealgoMulParam
			H, _ := bits.Mul64(R, ss.t)
			r := ahi - H + ss.t

			if r >= ss.t {
				r -= ss.t
			}

			return r
		}

		ss.reducealgoAdd = func(x uint64) uint64 {

			s0, _ := bits.Mul64(x, ss.reducealgoAddParam)

			r := x - s0*ss.t

			if r >= ss.t {
				r -= ss.t
			}

			return r
		}
	}

	// Integer and rational part of QiBarre * t/Qi
	ss.wi = make([]uint64, len(r.Modulus))
	ss.ti = make([]*big.Float, len(r.Modulus))

	ss.p1BI = make([]*big.Int, r.N)
	for i := range ss.p1BI {
		ss.p1BI[i] = new(big.Int)
		ss.p1BI[i].Mul(ss.ring.ModulusBigint, ss.ring.ModulusBigint) // Extends to Q^2
	}

	for i, qi := range r.Modulus {

		QiB.SetUint64(qi)
		QiBF.SetUint64(qi)
		QiStar.Quo(r.ModulusBigint, QiB)
		QiBarre.ModInverse(QiStar, QiB)
		QiBarre.Mod(QiBarre, QiB)
		QiBarreBF.SetInt(QiBarre)

		tmp = new(big.Float).Quo(tBF, QiBF)
		tmp.Mul(tmp, QiBarreBF)

		//floor( ([Q/Qi]^(-1))_{Qi} * t/Qi )
		ss.wi[i], _ = tmp.Uint64()

		// If t is not a power of 2 converts in Montgomery form
		if (t&(t-1)) != 0 && t != 0 {
			ss.wi[i] = MForm(ss.wi[i], t, BRedParams(t))
		}

		QiBarre.Mul(QiBarre, NewUint(t))
		QiBarre.Mod(QiBarre, QiB)
		QiBarreBF.SetInt(QiBarre)

		ss.ti[i] = new(big.Float).SetPrec(SimpleScalerFloatPrecision).Quo(QiBarreBF, QiBF)
	}

	return
}

// DivByQOverTRounded returns p1 scaled by a factor t/Q and mod t on the reciever p2.
func (ss *SimpleScaler) DivByQOverTRounded(p1, p2 *Poly) {
	ss.reconstructThenScale(p1, p2)
}

// reconstructThenScale operates the RNS reconstruction and scaling sequentially.
func (ss *SimpleScaler) reconstructThenScale(p1, p2 *Poly) {

	// reconstruction
	ss.ring.PolyToBigintNoAlloc(p1, ss.p1BI)

	// scaling
	for i, coeff := range ss.p1BI {
		coeff.Mul(coeff, ss.tBI)
		coeff.QuoRem(coeff, ss.ring.ModulusBigint, ss.r)
		if ss.r.Lsh(ss.r, 1).CmpAbs(ss.ring.ModulusBigint) != -1 {
			coeff.Add(coeff, ss.one)
		}

		for j := range p2.Coeffs {
			p2.Coeffs[j][i] = coeff.Uint64()
		}
	}
}

// reconstructAndScale operates the RNS reconstruction and scaling operation in one single pass over the coefficients.
// Algorithm from https://eprint.iacr.org/2018/117.pdf.
func (ss *SimpleScaler) reconstructAndScale(p1, p2 *Poly) {

	for i := uint64(0); i < ss.ring.N; i++ {

		var a uint64
		var bBF big.Float
		var p1j big.Float
		for j := range ss.ring.Modulus {
			// round(xi*wi + xi*ti)%t
			a += ss.reducealgoMul(ss.wi[j], p1.Coeffs[j][i])

			p1j.SetPrec(SimpleScalerFloatPrecision).SetUint64(p1.Coeffs[j][i])
			bBF.SetPrec(SimpleScalerFloatPrecision).Add(&bBF, p1j.Mul(&p1j, ss.ti[j]))
		}

		bBF.Add(&bBF, new(big.Float).SetFloat64(0.5))
		buint64, _ := bBF.Uint64()

		abf := a + buint64

		a = ss.reducealgoAdd(abf)

		for j := 0; j < len(p2.Coeffs); j++ {
			p2.Coeffs[j][i] = a
		}
	}
}

// ============== Scaling-related methods ==============

// DivFloorByLastModulusNTT divides floor the polynomial by its last modulus. Input must be in the NTT domain.
func (r *Ring) DivFloorByLastModulusNTT(p0 *Poly) {

	level := len(p0.Coeffs) - 1

	pTmp := make([]uint64, r.N)

	InvNTT(p0.Coeffs[level], p0.Coeffs[level], r.N, r.NttPsiInv[level], r.NttNInv[level], r.Modulus[level], r.MredParams[level])

	for i := 0; i < level; i++ {

		NTT(p0.Coeffs[level], pTmp, r.N, r.NttPsi[i], r.Modulus[i], r.MredParams[i], r.BredParams[i])

		p0tmp := p0.Coeffs[i]

		qi := r.Modulus[i]
		mredParams := r.MredParams[i]
		rescalParams := r.RescaleParams[level-1][i]

		// (x[i] - x[-1]) * InvQ
		for j := uint64(0); j < r.N; j++ {
			p0tmp[j] = MRed(p0tmp[j]+(qi-pTmp[j]), rescalParams, qi, mredParams)
		}
	}

	p0.Coeffs = p0.Coeffs[:level]
}

// DivFloorByLastModulus divides floor the polynomial by its last modulus.
func (r *Ring) DivFloorByLastModulus(p0 *Poly) {

	level := len(p0.Coeffs) - 1

	for i := 0; i < level; i++ {
		p0tmp := p0.Coeffs[level]
		p1tmp := p0.Coeffs[i]
		qi := r.Modulus[i]
		bredParams := r.BredParams[i]
		mredParams := r.MredParams[i]
		rescaleParams := r.RescaleParams[level-1][i]
		// (x[i] - x[-1]) * InvQ
		for j := uint64(0); j < r.N; j++ {
			p1tmp[j] = MRed(p1tmp[j]+(qi-BRedAdd(p0tmp[j], qi, bredParams)), rescaleParams, qi, mredParams)
		}
	}

	p0.Coeffs = p0.Coeffs[:level]
}

// DivFloorByLastModulusManyNTT divides floor sequentially nbRescales times the polynmial by its last modulus. Input must be in the NTT domain.
func (r *Ring) DivFloorByLastModulusManyNTT(p0 *Poly, nbRescales uint64) {
	r.InvNTTLvl(uint64(len(p0.Coeffs)-1), p0, p0)
	r.DivFloorByLastModulusMany(p0, nbRescales)
	r.NTTLvl(uint64(len(p0.Coeffs)-1), p0, p0)
}

// DivFloorByLastModulusMany divides floor sequentially nbRescales times the polynmial by its last modulus.
func (r *Ring) DivFloorByLastModulusMany(p0 *Poly, nbRescales uint64) {
	for k := uint64(0); k < nbRescales; k++ {
		r.DivFloorByLastModulus(p0)
	}
}

// DivRoundByLastModulusNTT divides round the polynomial by its last modulus. Input must be in the NTT domain.
func (r *Ring) DivRoundByLastModulusNTT(p0 *Poly) {

	var pHalf, pHalfNegQi uint64

	level := len(p0.Coeffs) - 1

	pTmp := make([]uint64, r.N)

	InvNTT(p0.Coeffs[level], p0.Coeffs[level], r.N, r.NttPsiInv[level], r.NttNInv[level], r.Modulus[level], r.MredParams[level])

	// Centers by (p-1)/2
	pHalf = (r.Modulus[level] - 1) >> 1
	p0tmp := p0.Coeffs[level]
	pj := r.Modulus[level]

	for i := uint64(0); i < r.N; i++ {
		p0tmp[i] = CRed(p0tmp[i]+pHalf, pj)
	}

	for i := 0; i < level; i++ {

		p1tmp := p0.Coeffs[i]

		qi := r.Modulus[i]
		bredParams := r.BredParams[i]
		mredParams := r.MredParams[i]
		rescaleParams := r.RescaleParams[level-1][i]

		pHalfNegQi = r.Modulus[i] - BRedAdd(pHalf, qi, bredParams)

		for j := uint64(0); j < r.N; j++ {
			pTmp[j] = p0tmp[j] + pHalfNegQi
		}

		NTT(pTmp, pTmp, r.N, r.NttPsi[i], qi, mredParams, bredParams)

		// (x[i] - x[-1]) * InvQ
		for j := uint64(0); j < r.N; j++ {
			p1tmp[j] = MRed(p1tmp[j]+(qi-pTmp[j]), rescaleParams, qi, mredParams)
		}
	}

	p0.Coeffs = p0.Coeffs[:level]
}

// DivRoundByLastModulus divides round the polynomial by its last modulus. Input must be in the NTT domain.
func (r *Ring) DivRoundByLastModulus(p0 *Poly) {

	var pHalf, pHalfNegQi uint64

	level := len(p0.Coeffs) - 1

	// Centers by (p-1)/2
	pHalf = (r.Modulus[level] - 1) >> 1
	p0tmp := p0.Coeffs[level]
	pj := r.Modulus[level]
	pHalf = (r.Modulus[level] - 1) >> 1

	for i := uint64(0); i < r.N; i++ {
		p0tmp[i] = CRed(p0tmp[i]+pHalf, pj)
	}

	for i := 0; i < level; i++ {

		p1tmp := p0.Coeffs[i]

		qi := r.Modulus[i]
		bredParams := r.BredParams[i]
		mredParams := r.MredParams[i]
		rescaleParams := r.RescaleParams[level-1][i]

		pHalfNegQi = r.Modulus[i] - BRedAdd(pHalf, qi, bredParams)

		// (x[i] - x[-1]) * InvQ
		for j := uint64(0); j < r.N; j++ {
			p1tmp[j] = MRed(p1tmp[j]+(qi-BRedAdd(p0tmp[j]+pHalfNegQi, qi, bredParams)), rescaleParams, qi, mredParams)
		}
	}

	p0.Coeffs = p0.Coeffs[:level]
}

// DivRoundByLastModulusManyNTT divides round sequentially nbRescales times the polynmial by its last modulus. Input must be in the NTT domain.
func (r *Ring) DivRoundByLastModulusManyNTT(p0 *Poly, nbRescales uint64) {
	r.InvNTTLvl(uint64(len(p0.Coeffs)-1), p0, p0)
	r.DivRoundByLastModulusMany(p0, nbRescales)
	r.NTTLvl(uint64(len(p0.Coeffs)-1), p0, p0)
}

// DivRoundByLastModulusMany divides round sequentially nbRescales times the polynmial by its last modulus.
func (r *Ring) DivRoundByLastModulusMany(p0 *Poly, nbRescales uint64) {
	for k := uint64(0); k < nbRescales; k++ {
		r.DivRoundByLastModulus(p0)
	}
}
