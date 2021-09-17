package ring

import (
	"math/big"
	"math/bits"
	"unsafe"
)

// Scaler is an interface that rescales polynomial coefficients by a fraction t/Q.
type Scaler interface {
	// DivByQOverTRounded returns p1 scaled by a factor t/Q and mod t on the receiver p2.
	DivByQOverTRounded(p1, p2 *Poly)
}

// RNSScaler implements the Scaler interface by performing a scaling by t/Q in the RNS domain.
// This implementation of the Scaler interface is preferred over the SimpleScaler implementation.
type RNSScaler struct {
	ringQ, ringT *Ring
	polypoolQ    *Poly
	polypoolT    *Poly

	qHalf     *big.Int // (q-1)/2
	qHalfModT uint64   // (q-1)/2 mod t
	qInv      uint64   //(q mod t)^-1 mod t

	paramsQP modupParams
}

// NewRNSScaler creates a new SimpleScaler from t, the modulus under which the reconstruction is returned, the Ring in which the polynomial to reconstruct is represented.
func NewRNSScaler(ringQ, ringT *Ring) (rnss *RNSScaler) {

	rnss = new(RNSScaler)

	rnss.ringQ = ringQ
	rnss.ringT = ringT

	rnss.polypoolQ = ringQ.NewPoly()
	rnss.polypoolT = ringT.NewPoly()

	t := ringT.Modulus[0]

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

// DivByQOverTRounded returns p1 scaled by a factor t/Q and mod t on the receiver p2.
func (rnss *RNSScaler) DivByQOverTRounded(p1Q, p2T *Poly) {

	ringQ := rnss.ringQ
	ringT := rnss.ringT

	T := ringT.Modulus[0]
	p2tmp := p2T.Coeffs[0]
	p3tmp := rnss.polypoolT.Coeffs[0]
	mredParams := rnss.ringT.MredParams[0]
	qInv := T - rnss.qInv
	qHalfModT := T - rnss.qHalfModT

	// Multiply P_{Q} by t and extend the basis from P_{Q} to t*(P_{Q}||P_{t})
	// Since the coefficients of P_{t} are multiplied by t, they are all zero,
	// hence the basis extension can be omitted
	ringQ.MulScalar(p1Q, T, rnss.polypoolQ)

	// Center t*P_{Q} around (Q-1)/2 to round instead of floor during the division
	ringQ.AddScalarBigint(rnss.polypoolQ, rnss.qHalf, rnss.polypoolQ)

	// Extend the basis of (t*P_{Q} + (Q-1)/2) to (t*P_{t} + (Q-1)/2)
	modUpExact(rnss.polypoolQ.Coeffs, rnss.polypoolT.Coeffs, ringQ, ringT, rnss.paramsQP)

	// Compute [Q^{-1} * (t*P_{t} -   (t*P_{Q} - ((Q-1)/2 mod t)))] mod t which returns round(t/Q * P_{Q}) mod t
	for j := 0; j < ringQ.N; j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

		z[0] = MRed(qHalfModT+x[0], qInv, T, mredParams)
		z[1] = MRed(qHalfModT+x[1], qInv, T, mredParams)
		z[2] = MRed(qHalfModT+x[2], qInv, T, mredParams)
		z[3] = MRed(qHalfModT+x[3], qInv, T, mredParams)
		z[4] = MRed(qHalfModT+x[4], qInv, T, mredParams)
		z[5] = MRed(qHalfModT+x[5], qInv, T, mredParams)
		z[6] = MRed(qHalfModT+x[6], qInv, T, mredParams)
		z[7] = MRed(qHalfModT+x[7], qInv, T, mredParams)
	}
}

// SimpleScaler implements the Scaler interface by performing an RNS reconstruction and scaling by t/Q.
// This implementation of the Scaler interface is less efficient than the RNSScaler, but uses simple
// multi-precision arithmetic of the math/big package.
type SimpleScaler struct {
	ringQ *Ring

	r, one *big.Int

	t   uint64 // Plaintext modulus
	tBI *big.Int
	wi  []uint64 // Integer parts of   ([Q/Qi]^(-1))_{Qi} * t/Qi
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
// ringQ, the Ring in which the polynomial to reconstruct is represented.
func NewSimpleScaler(t uint64, ringQ *Ring) (ss *SimpleScaler) {

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
	ss.ringQ = ringQ

	// Assign the correct reduction algorithm depending on the provided t
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

		// Otherwise, we can use Montgomery reduction
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

	// Integer and rational parts of QiBarre * t/Qi
	ss.wi = make([]uint64, len(ringQ.Modulus))
	ss.ti = make([]*big.Float, len(ringQ.Modulus))

	ss.p1BI = make([]*big.Int, ringQ.N)
	for i := range ss.p1BI {
		ss.p1BI[i] = new(big.Int)
		ss.p1BI[i].Mul(ss.ringQ.ModulusBigint, ss.ringQ.ModulusBigint) // Extend to Q^2
	}

	for i, qi := range ringQ.Modulus {

		QiB.SetUint64(qi)
		QiBF.SetUint64(qi)
		QiStar.Quo(ringQ.ModulusBigint, QiB)
		QiBarre.ModInverse(QiStar, QiB)
		QiBarre.Mod(QiBarre, QiB)
		QiBarreBF.SetInt(QiBarre)

		tmp = new(big.Float).Quo(tBF, QiBF)
		tmp.Mul(tmp, QiBarreBF)

		// floor( ([Q/Qi]^(-1))_{Qi} * t/Qi )
		ss.wi[i], _ = tmp.Uint64()

		// If t is not a power of 2, convert to Montgomery form
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

// DivByQOverTRounded returns p1 scaled by a factor t/Q and mod t on the receiver p2.
func (ss *SimpleScaler) DivByQOverTRounded(p1, p2 *Poly) {
	ss.reconstructThenScale(p1, p2)
}

// reconstructThenScale performs the RNS reconstruction and scaling sequentially.
func (ss *SimpleScaler) reconstructThenScale(p1, p2 *Poly) {

	// reconstruction
	ss.ringQ.PolyToBigint(p1, ss.p1BI)

	// scaling
	for i, coeff := range ss.p1BI {
		coeff.Mul(coeff, ss.tBI)
		coeff.QuoRem(coeff, ss.ringQ.ModulusBigint, ss.r)
		if ss.r.Lsh(ss.r, 1).CmpAbs(ss.ringQ.ModulusBigint) != -1 {
			coeff.Add(coeff, ss.one)
		}

		for j := range p2.Coeffs {
			p2.Coeffs[j][i] = coeff.Uint64()
		}
	}
}

// reconstructAndScale performs the RNS reconstruction and scaling operations in one single pass over the coefficients.
// Algorithm from https://eprint.iacr.org/2018/117.pdf.
func (ss *SimpleScaler) reconstructAndScale(p1, p2 *Poly) {

	for i := 0; i < ss.ringQ.N; i++ {

		var a uint64
		var bBF big.Float
		var p1j big.Float
		for j := range ss.ringQ.Modulus {
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

// DivFloorByLastModulusNTTLvl divides (floored) the polynomial by its last modulus. The input must be in the NTT domain.
// Output poly level must be equal or one less than input level.
func (r *Ring) DivFloorByLastModulusNTTLvl(level int, p0, pool, p1 *Poly) {

	pool0 := pool.Coeffs[0]
	pool1 := pool.Coeffs[1]

	InvNTTLazy(p0.Coeffs[level], pool0, r.N, r.NttPsiInv[level], r.NttNInv[level], r.Modulus[level], r.MredParams[level])

	for i := 0; i < level; i++ {
		NTTLazy(pool0, pool1, r.N, r.NttPsi[i], r.Modulus[i], r.MredParams[i], r.BredParams[i])
		// (-x[i] + x[-1]) * -InvQ
		SubVecAndMulScalarMontgomeryTwoQiVec(pool1, p0.Coeffs[i], p1.Coeffs[i], r.RescaleParams[level-1][i], r.Modulus[i], r.MredParams[i])
	}
}

// DivFloorByLastModulusLvl divides (floored) the polynomial by its last modulus.
// Output poly level must be equal or one less than input level.
func (r *Ring) DivFloorByLastModulusLvl(level int, p0, p1 *Poly) {
	for i := 0; i < level; i++ {
		// (x[i] - x[-1]) * InvQ
		SubVecAndMulScalarMontgomeryTwoQiVec(p0.Coeffs[level], p0.Coeffs[i], p1.Coeffs[i], r.RescaleParams[level-1][i], r.Modulus[i], r.MredParams[i])
	}
}

// DivFloorByLastModulusManyNTTLvl divides (floored) sequentially nbRescales times the polynomial by its last modulus. Input must be in the NTT domain.
// Output poly level must be equal or nbRescales less than input level.
func (r *Ring) DivFloorByLastModulusManyNTTLvl(level, nbRescales int, p0, pool, p1 *Poly) {

	if nbRescales == 0 {

		if p0 != p1 {
			CopyValuesLvl(level, p0, p1)
		}

	} else {

		r.InvNTTLvl(level, p0, pool)

		for i := 0; i < nbRescales; i++ {
			r.DivFloorByLastModulusLvl(level-i, pool, pool)
		}

		r.NTTLvl(level-nbRescales, pool, p1)
	}
}

// DivFloorByLastModulusManyLvl divides (floored) sequentially nbRescales times the polynomial by its last modulus.
// Output poly level must be equal or nbRescales less than input level.
func (r *Ring) DivFloorByLastModulusManyLvl(level, nbRescales int, p0, pool, p1 *Poly) {

	if nbRescales == 0 {

		if p0 != p1 {
			CopyValuesLvl(level, p0, p1)
		}

	} else {

		if nbRescales > 1 {
			r.DivFloorByLastModulusLvl(level, p0, pool)

			for i := 1; i < nbRescales; i++ {

				if i == nbRescales-1 {
					r.DivFloorByLastModulusLvl(level-i, pool, p1)
				} else {
					r.DivFloorByLastModulusLvl(level-i, pool, pool)
				}
			}

		} else {
			r.DivFloorByLastModulusLvl(level, p0, p1)
		}
	}
}

// DivRoundByLastModulusNTTLvl divides (rounded) the polynomial by its last modulus. The input must be in the NTT domain.
// Output poly level must be equal or one less than input level.
func (r *Ring) DivRoundByLastModulusNTTLvl(level int, p0, pool, p1 *Poly) {

	pool0 := pool.Coeffs[0]
	pool1 := pool.Coeffs[1]

	InvNTT(p0.Coeffs[level], pool0, r.N, r.NttPsiInv[level], r.NttNInv[level], r.Modulus[level], r.MredParams[level])

	// Center by (p-1)/2
	pj := r.Modulus[level]
	pHalf := (pj - 1) >> 1

	AddScalarVec(pool0, pool0, pHalf, pj)

	for i := 0; i < level; i++ {
		qi := r.Modulus[i]
		qInv := r.MredParams[i]
		bredParams := r.BredParams[i]

		AddScalarNoModVec(pool0, pool1, r.Modulus[i]-BRedAdd(pHalf, qi, bredParams))
		NTTLazy(pool1, pool1, r.N, r.NttPsi[i], qi, qInv, bredParams)
		SubVecAndMulScalarMontgomeryTwoQiVec(pool1, p0.Coeffs[i], p1.Coeffs[i], r.RescaleParams[level-1][i], r.Modulus[i], r.MredParams[i])
	}
}

// DivRoundByLastModulusLvl divides (rounded) the polynomial by its last modulus. The input must be in the NTT domain.
// Output poly level must be equal or one less than input level.
func (r *Ring) DivRoundByLastModulusLvl(level int, p0, p1 *Poly) {
	// Center by (p-1)/2
	pHalf := (r.Modulus[level] - 1) >> 1
	pj := r.Modulus[level]
	AddScalarVec(p0.Coeffs[level], p0.Coeffs[level], pHalf, pj)
	for i := 0; i < level; i++ {
		qi := r.Modulus[i]
		bredParams := r.BredParams[i]
		AddScalarNoModAndNegTwoQiNoModVec(p0.Coeffs[i], p0.Coeffs[i], r.Modulus[i]-BRedAdd(pHalf, qi, bredParams), qi)
		AddVecNoModAndMulScalarMontgomeryVec(p0.Coeffs[level], p0.Coeffs[i], p1.Coeffs[i], r.RescaleParams[level-1][i], qi, r.MredParams[i])
	}
}

// DivRoundByLastModulusManyNTTLvl divides (rounded) sequentially nbRescales times the polynomial by its last modulus. The input must be in the NTT domain.
// Output poly level must be equal or nbRescales less than input level.
func (r *Ring) DivRoundByLastModulusManyNTTLvl(level, nbRescales int, p0, pool, p1 *Poly) {

	if nbRescales == 0 {

		if p0 != p1 {
			CopyValuesLvl(level, p0, p1)
		}

	} else {

		if nbRescales > 1 {

			r.InvNTTLvl(level, p0, pool)

			for i := 0; i < nbRescales; i++ {
				r.DivRoundByLastModulusLvl(level-i, pool, pool)
			}

			r.NTTLvl(p1.Level(), pool, p1)

		} else {
			r.DivRoundByLastModulusNTTLvl(level, p0, pool, p1)
		}
	}
}

// DivRoundByLastModulusManyLvl divides (rounded) sequentially nbRescales times the polynomial by its last modulus.
// Output poly level must be equal or nbRescales less than input level.
func (r *Ring) DivRoundByLastModulusManyLvl(level, nbRescales int, p0, pool, p1 *Poly) {

	if nbRescales == 0 {

		if p0 != p1 {
			CopyValuesLvl(level, p0, p1)
		}

	} else {

		if nbRescales > 1 {

			r.DivRoundByLastModulusLvl(level, p0, pool)

			for i := 1; i < nbRescales; i++ {

				if i == nbRescales-1 {
					r.DivRoundByLastModulusLvl(level-i, pool, p1)
				} else {
					r.DivRoundByLastModulusLvl(level-i, pool, pool)
				}
			}

		} else {
			r.DivRoundByLastModulusLvl(level, p0, p1)
		}
	}
}
