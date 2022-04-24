package bfv

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// Scaler is an interface that rescales polynomial coefficients by a fraction t/Q.
type Scaler interface {
	// DivByQOverTRoundedLvl returns p1 scaled by a factor t/Q and mod t on the receiver p2.
	DivByQOverTRoundedLvl(level int, p1, p2 *ring.Poly)
	ScaleUpByQOverTLvl(level int, p1, p2 *ring.Poly)
}

// RNSScaler implements the Scaler interface by performing a scaling by t/Q in the RNS domain.
type RNSScaler struct {
	ringQ, ringT *ring.Ring

	buffQ *ring.Poly
	buffP *ring.Poly

	qHalf     []*big.Int // (q-1)/2
	qHalfModT []uint64   // (q-1)/2 mod t
	qInv      []uint64   //(q mod t)^-1 mod t
	tInvModQi []uint64   // t^-1 mod qi

	paramsQP []ring.ModupParams

	tDividesQ bool
}

// NewRNSScaler creates a new RNSScaler from t, the modulus under which the reconstruction is returned, the Ring in which the polynomial to reconstruct is represented.
func NewRNSScaler(ringQ *ring.Ring, T uint64) (rnss *RNSScaler) {

	if utils.IsInSliceUint64(T, ringQ.Modulus) && ringQ.Modulus[0] != T {
		panic("cannot NewRNSScaler: T must be Q[0] if T|Q")
	}

	rnss = new(RNSScaler)

	rnss.ringQ = ringQ

	rnss.buffQ = ringQ.NewPoly()

	rnss.ringT = new(ring.Ring)
	rnss.ringT.N = ringQ.N
	rnss.ringT.Modulus = []uint64{T}
	rnss.ringT.BredParams = [][]uint64{ring.BRedParams(T)}
	rnss.ringT.MredParams = []uint64{ring.MRedParams(T)}
	rnss.buffP = rnss.ringT.NewPoly()

	rnss.tDividesQ = T == ringQ.Modulus[0]

	if !rnss.tDividesQ {

		rnss.tInvModQi = make([]uint64, len(ringQ.Modulus))
		for i, qi := range ringQ.Modulus {
			rnss.tInvModQi[i] = ring.MForm(ring.ModExp(T, qi-2, qi), qi, ringQ.BredParams[i])
		}

		rnss.qHalf = make([]*big.Int, len(ringQ.Modulus))
		rnss.qInv = make([]uint64, len(ringQ.Modulus))
		rnss.qHalfModT = make([]uint64, len(ringQ.Modulus))
		rnss.paramsQP = make([]ring.ModupParams, len(ringQ.Modulus))

		bigQ := new(big.Int).SetUint64(1)
		tmp := new(big.Int)
		bredParams := ring.BRedParams(T)
		TBig := ring.NewUint(T)
		for i := range ringQ.Modulus {
			rnss.paramsQP[i] = ring.GenModUpParams(ringQ.Modulus[:i+1], rnss.ringT.Modulus)

			bigQ.Mul(bigQ, ring.NewUint(ringQ.Modulus[i]))

			rnss.qInv[i] = tmp.Mod(bigQ, TBig).Uint64()
			rnss.qInv[i] = ring.ModExp(rnss.qInv[i], T-2, T)
			rnss.qInv[i] = ring.MForm(rnss.qInv[i], T, bredParams)

			rnss.qHalf[i] = new(big.Int).Set(bigQ)
			rnss.qHalf[i].Rsh(rnss.qHalf[i], 1)

			rnss.qHalfModT[i] = tmp.Mod(rnss.qHalf[i], TBig).Uint64()
		}
	}

	return
}

// DivByQOverTRoundedLvl returns p1 scaled by a factor t/Q and mod t on the receiver p2.
func (rnss *RNSScaler) DivByQOverTRoundedLvl(level int, p1Q, p2T *ring.Poly) {

	ringQ := rnss.ringQ

	if level > 0 {
		if rnss.tDividesQ {
			ringQ.DivRoundByLastModulusManyLvl(level, level, p1Q, rnss.buffQ, p2T)
		} else {

			ringT := rnss.ringT
			T := ringT.Modulus[0]
			tInv := ringT.MredParams[0]
			p2tmp := p2T.Coeffs[0]
			p3tmp := rnss.buffP.Coeffs[0]
			qInv := T - rnss.qInv[level]
			qHalfModT := T - rnss.qHalfModT[level]

			// Multiplies P_{Q} by t and extend the basis from P_{Q} to t*(P_{Q}||P_{t})
			// Since the coefficients of P_{t} are multiplied by t, they are all zero,
			// hence the basis extension can be omitted
			ringQ.MulScalarLvl(level, p1Q, T, rnss.buffQ)

			// Centers t*P_{Q} around (Q-1)/2 to round instead of floor during the division
			ringQ.AddScalarBigintLvl(level, rnss.buffQ, rnss.qHalf[level], rnss.buffQ)

			// Extends the basis of (t*P_{Q} + (Q-1)/2) to (t*P_{t} + (Q-1)/2)
			ring.ModUpExact(rnss.buffQ.Coeffs[:level+1], rnss.buffP.Coeffs, ringQ, ringT, rnss.paramsQP[level])

			// Computes [Q^{-1} * (t*P_{t} - (t*P_{Q} - ((Q-1)/2 mod t)))] mod t which returns round(t/Q * P_{Q}) mod t
			ring.AddScalarNoModAndMulScalarMontgomeryVec(p3tmp, p2tmp, qHalfModT, qInv, T, tInv)
		}
	} else {
		if rnss.tDividesQ {
			copy(p2T.Coeffs[0], p1Q.Coeffs[0])
		} else {
			// In this case lvl = 0 and T < Q. This step has a maximum precision of 53 bits, however
			// since |Q| < 62 bits, and min(logN) = 10, then |<s, e>| > 10 bits, hence there is no
			// possible case where |T| > 51 bits & lvl = 0 that does not lead to an overflow of
			// the error when decrypting.
			qOverT := float64(ringQ.Modulus[0]) / float64(rnss.ringT.Modulus[0])
			tmp0, tmp1 := p2T.Coeffs[0], p1Q.Coeffs[0]
			for i := 0; i < ringQ.N; i++ {
				tmp0[i] = uint64(float64(tmp1[i])/qOverT + 0.5)
			}
		}
	}
}

// ScaleUpByQOverTLvl takes a Poly pIn in ringT, scales its coefficients up by (Q/T) mod Q, and writes the result on pOut.
func (rnss *RNSScaler) ScaleUpByQOverTLvl(level int, pIn, pOut *ring.Poly) {

	if !rnss.tDividesQ {
		ScaleUpTCoprimeWithQVecLvl(level, rnss.ringQ, rnss.ringT, rnss.tInvModQi, rnss.buffQ.Coeffs[0], pIn, pOut)
	} else {
		ScaleUpTIsQ0VecLvl(level, rnss.ringQ, pIn, pOut)
	}
}

// ScaleUpTCoprimeWithQVecLvl takes a Poly pIn in ringT, scales its coefficients up by (Q/T) mod Q, and writes the result in a
// Poly pOut in ringQ.
func ScaleUpTCoprimeWithQVecLvl(level int, ringQ, ringT *ring.Ring, tInvModQi, buffN []uint64, pIn, pOut *ring.Poly) {

	qModTmontgomery := ring.MForm(new(big.Int).Mod(ringQ.ModulusAtLevel[level], ring.NewUint(ringT.Modulus[0])).Uint64(), ringT.Modulus[0], ringT.BredParams[0])

	t := ringT.Modulus[0]
	tHalf := t >> 1
	tInv := ringT.MredParams[0]

	// (x * Q + T/2) mod T
	ring.MulScalarMontgomeryAndAddScalarVec(pIn.Coeffs[0], buffN, tHalf, qModTmontgomery, t, tInv)

	// (x * T^-1 - T/2) mod Qi
	for i := 0; i < level+1; i++ {
		p0tmp := buffN
		p1tmp := pOut.Coeffs[i]
		qi := ringQ.Modulus[i]
		bredParams := ringQ.BredParams[i]
		mredParams := ringQ.MredParams[i]
		rescaleParams := qi - tInvModQi[i]
		tHalfNegQi := qi - ring.BRedAdd(tHalf, qi, bredParams)

		ring.AddScalarNoModAndMulScalarMontgomeryVec(p0tmp, p1tmp, tHalfNegQi, rescaleParams, qi, mredParams)
	}
}

// ScaleUpTIsQ0VecLvl takes a Poly pIn in ringT, scales its coefficients up by (Q/T) mod Q, and writes the result on a
// Poly pOut in ringQ.
// T is in this case assumed to be the first prime in the moduli chain.
func ScaleUpTIsQ0VecLvl(level int, ringQ *ring.Ring, pIn, pOut *ring.Poly) {

	// Q/T mod T
	tmp := new(big.Int)
	tmp.Quo(ringQ.ModulusAtLevel[level], ringQ.ModulusAtLevel[0])
	QOverTMont := ring.MForm(tmp.Mod(tmp, new(big.Int).SetUint64(ringQ.Modulus[0])).Uint64(), ringQ.Modulus[0], ringQ.BredParams[0])

	// pOut = Q/T * pIn
	ring.MulScalarMontgomeryVec(pIn.Coeffs[0], pOut.Coeffs[0], QOverTMont, ringQ.Modulus[0], ringQ.MredParams[0])

	for i := 1; i < level+1; i++ {
		ring.ZeroVec(pOut.Coeffs[i])
	}
}
