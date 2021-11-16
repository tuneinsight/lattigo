package lwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"math/big"
)

// InitLUT takes a function g, and creates an LUT polynomial for the function between the intervals a, b.
// Inputs to the LUT evaluation are assumed to have been normalized with the change of basis (2*x - a - b)/(b-a).
// Interval a, b should take into account the "drift" of the value x, caused by the change of modulus from Q to 2N.
func InitLUT(g func(x float64) (y float64), scale float64, ringQ *ring.Ring, a, b float64) (F *ring.Poly) {
	F = ringQ.NewPoly()
	Q := ringQ.Modulus

	// Discretization interval
	interval := 2.0 / float64(ringQ.N)

	for j, qi := range Q {

		// Interval [-1, 0] of g(x)
		for i := 0; i < (ringQ.N>>1)+1; i++ {
			F.Coeffs[j][i] = scaleUp(g(normalizeInv(-interval*float64(i), a, b)), scale, qi)
		}

		// Interval ]0, 1[ of g(x)
		for i := (ringQ.N >> 1) + 1; i < ringQ.N; i++ {
			F.Coeffs[j][i] = scaleUp(-g(normalizeInv(interval*float64(ringQ.N-i), a, b)), scale, qi)
		}
	}

	ringQ.NTT(F, F)

	return
}

// ExtractAndEvaluateLUT extracts on the fly LWE samples and evaluate the provided LUT on the LWE.
// ct : a rlwe Ciphertext with coefficient encoded values at level 0
// lutPolyWihtSlotIndex : a map with [slot_index] -> LUT
// lutKey : LUTKey
// Returns a map[slot_index] -> LUT(ct[slot_index])
func (h *Handler) ExtractAndEvaluateLUT(ct *rlwe.Ciphertext, lutPolyWihtSlotIndex map[int]*ring.Poly, lutKey *LUTKey) (res map[int]*rlwe.Ciphertext) {

	ks := h.KeySwitcher

	bRLWEMod2N := h.poolMod2N[0]
	aRLWEMod2N := h.poolMod2N[1]

	acc := h.accumulator

	ringQLUT := h.paramsLUT.RingQ()
	ringQLWE := h.paramsLWE.RingQ()
	ringQPLUT := h.paramsLUT.RingQP()

	mask := uint64(ringQLUT.N<<1) - 1

	level := ct.Level()

	ringQLWE.InvNTTLvl(level, ct.Value[0], acc.Value[0])
	ringQLWE.InvNTTLvl(level, ct.Value[1], acc.Value[1])

	h.ModSwitchRLWETo2N(acc.Value[1], acc.Value[1])

	// Copy coefficients multiplied by X^{N-1} in reverse order:
	// a_{0} -a_{N-1} -a2_{N-2} ... -a_{1}
	tmp0 := aRLWEMod2N.Coeffs[0]
	tmp1 := acc.Value[1].Coeffs[0]
	tmp0[0] = tmp1[0]
	for j := 1; j < ringQLWE.N; j++ {
		tmp0[j] = -tmp1[ringQLWE.N-j] & mask
	}

	h.ModSwitchRLWETo2N(acc.Value[0], bRLWEMod2N)

	res = make(map[int]*rlwe.Ciphertext)

	tmpRGSW := rlwe.NewCiphertextRGSWNTT(h.paramsLUT, h.paramsLUT.MaxLevel())

	var prevIndex int
	for index := 0; index < ringQLWE.N; index++ {

		if lut, ok := lutPolyWihtSlotIndex[index]; ok {

			MulBySmallMonomialMod2N(mask, aRLWEMod2N, index-prevIndex)
			prevIndex = index

			a := aRLWEMod2N.Coeffs[0]
			b := bRLWEMod2N.Coeffs[0][index]

			ringQLUT.MulCoeffsMontgomery(lut, h.xPowMinusOne[b].Q, acc.Value[0])
			ringQLUT.Add(acc.Value[0], lut, acc.Value[0])
			acc.Value[1].Zero() // TODO remove

			for j := 0; j < ringQLWE.N; j++ {
				MulRGSWByXPowAlphaMinusOne(lutKey.SkPos[j], h.xPowMinusOne[a[j]], ringQPLUT, tmpRGSW)
				MulRGSWByXPowAlphaMinusOneAndAdd(lutKey.SkNeg[j], h.xPowMinusOne[-a[j]&mask], ringQPLUT, tmpRGSW)
				AddRGSW(lutKey.EncOne, ringQPLUT, tmpRGSW) // TODO : add 1 in plaintext
				ks.MulRGSWSingleModulus(acc, tmpRGSW, acc)
			}

			// Extracts the first coefficient
			res[index] = acc.CopyNew() //ExtractLWEFromRLWESingle(acc, ringQRLWE)
		}
	}

	return
}

//MulBySmallMonomial multiplies pol by x^n, with 0 <= n < N
func MulBySmallMonomialMod2N(mask uint64, pol *ring.Poly, n int) {
	if n != 0 {
		N := len(pol.Coeffs[0])
		pol.Coeffs[0] = append(pol.Coeffs[0][N-n:], pol.Coeffs[0][:N-n]...)
		tmp := pol.Coeffs[0]
		for j := 0; j < n; j++ {
			tmp[j] = -tmp[j] & mask
		}
	}
}

// ModSwitchRLWETo2N applys round(x * 2N / Q) to the coefficients of polQ and returns the
// result on pol2N.
func (h *Handler) ModSwitchRLWETo2N(polQ *ring.Poly, pol2N *ring.Poly) {
	coeffsBigint := make([]*big.Int, len(polQ.Coeffs[0]))

	level := polQ.Level()

	ringQ := h.paramsLWE.RingQ()

	ringQ.PolyToBigintLvl(polQ.Level(), polQ, coeffsBigint)

	QBig := ring.NewUint(1)
	for i := 0; i < level+1; i++ {
		QBig.Mul(QBig, ring.NewUint(ringQ.Modulus[i]))
	}

	twoN := uint64(h.paramsLUT.N() << 1)
	twoNBig := ring.NewUint(twoN)
	tmp := pol2N.Coeffs[0]
	for i := 0; i < ringQ.N; i++ {
		coeffsBigint[i].Mul(coeffsBigint[i], twoNBig)
		ring.DivRound(coeffsBigint[i], QBig, coeffsBigint[i])
		tmp[i] = coeffsBigint[i].Uint64() & (twoN - 1)
	}
}
