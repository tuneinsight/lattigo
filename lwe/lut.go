package lwe

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/rlwe/rgsw"
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

// ExtractAndEvaluateLUTAndRepack extracts on the fly LWE samples and evaluate the provided LUT on the LWE and repacks everything into a single rlwe.Ciphertext.
// ct : a rlwe Ciphertext with coefficient encoded values at level 0
// lutPolyWihtSlotIndex : a map with [slot_index] -> LUT
// repackIndex : a map with [slot_index_have] -> slot_index_want
// lutKey : LUTKey
// Returns a *rlwe.Ciphertext
func (h *Handler) ExtractAndEvaluateLUTAndRepack(ct *rlwe.Ciphertext, lutPolyWihtSlotIndex map[int]*ring.Poly, repackIndex map[int]int, lutKey *LUTKey) (res *rlwe.Ciphertext) {
	cts := h.ExtractAndEvaluateLUT(ct, lutPolyWihtSlotIndex, lutKey)

	ciphertexts := make(map[int]*rlwe.Ciphertext)

	for i := range cts {
		ciphertexts[repackIndex[i]] = cts[i]
	}

	return h.evalRLWE.MergeRLWE(ciphertexts)
}

// EvalGate evaluates the selected binary gate on the input rlwe ciphertext
func (h *Handler) EvalGate(ct *Ciphertext, logNLWE int, gate *ring.Poly, lutKey *LUTKey) {
	eval := h.evalRGSW

	acc := h.accumulator

	ringQLUT := h.paramsLUT.RingQ()
	ringQLWE := h.paramsLWE.RingQ()
	ringQPLUT := h.paramsLUT.RingQP()

	// mod 2N
	mask := uint64(ringQLUT.N<<1) - 1

	tmpRGSW := rgsw.NewCiphertextNTT(h.paramsLUT, h.paramsLUT.MaxLevel())

	a := ct.Value[0][1:]
	b := ct.Value[0][0]

	// LWE = -as + m + e, a
	// LUT = LUT * X^{-as + m + e}
	ringQLUT.MulCoeffsMontgomery(gate, h.xPowMinusOne[b].Q, acc.Value[0])
	ringQLUT.Add(acc.Value[0], gate, acc.Value[0])
	acc.Value[1].Zero() // TODO remove
	for i := 0; i < ringQLWE.N; i++ {
		MulRGSWByXPowAlphaMinusOne(lutKey.SkPos[i], h.xPowMinusOne[a[i]], ringQPLUT, tmpRGSW)
		MulRGSWByXPowAlphaMinusOneAndAdd(lutKey.SkNeg[i], h.xPowMinusOne[-a[i]&mask], ringQPLUT, tmpRGSW)
		AddOneRGSW(lutKey.OneRGSW, ringQLUT, tmpRGSW)
		eval.ExternalProduct(acc, tmpRGSW, acc)
	}

	eval.SwitchKeysInPlace(0, acc.Value[1], nil, eval.Pool[1].Q, eval.Pool[2].Q) // TODO : add RLWE -> LWE Key
	ringQLUT.AddLvl(0, acc.Value[0], eval.Pool[1].Q, acc.Value[0])
	ringQLUT.InvNTT(eval.Pool[2].Q, acc.Value[1])
	ringQLUT.InvNTT(acc.Value[0], acc.Value[0])

	Qflo := float64(ringQLUT.Modulus[0])
	maskLWE := uint64(2<<logNLWE) - 1

	c := acc.Value[0].Coeffs[0][0]
	c = uint64(float64(c<<logNLWE)/Qflo + 0.5)
	ct.Value[0][0] = c & maskLWE

	c = acc.Value[1].Coeffs[0][0]
	c = uint64(float64(c<<logNLWE)/Qflo + 0.5)
	ct.Value[0][1] = c & maskLWE
	for i := 1; i < ringQLWE.N; i++ {
		c = acc.Value[1].Coeffs[0][ringQLUT.N-2*i]
		c = uint64(float64(c<<logNLWE)/Qflo + 0.5)
		ct.Value[0][i+1] = -c & maskLWE
	}

}

// ExtractAndEvaluateLUT extracts on the fly LWE samples and evaluate the provided LUT on the LWE.
// ct : a rlwe Ciphertext with coefficient encoded values at level 0
// lutPolyWihtSlotIndex : a map with [slot_index] -> LUT
// lutKey : LUTKey
// Returns a map[slot_index] -> LUT(ct[slot_index])
func (h *Handler) ExtractAndEvaluateLUT(ct *rlwe.Ciphertext, lutPolyWihtSlotIndex map[int]*ring.Poly, lutKey *LUTKey) (res map[int]*rlwe.Ciphertext) {

	eval := h.evalRGSW

	bRLWEMod2N := h.poolMod2N[0]
	aRLWEMod2N := h.poolMod2N[1]

	acc := h.accumulator

	ringQLUT := h.paramsLUT.RingQ()
	ringQLWE := h.paramsLWE.RingQ()
	ringQPLUT := h.paramsLUT.RingQP()

	// mod 2N
	mask := uint64(ringQLUT.N<<1) - 1

	ringQLWE.InvNTTLvl(ct.Level(), ct.Value[0], acc.Value[0])
	ringQLWE.InvNTTLvl(ct.Level(), ct.Value[1], acc.Value[1])

	// Switch modulus from Q to 2N
	h.ModSwitchRLWETo2NLvl(ct.Level(), acc.Value[1], acc.Value[1])

	// Conversion from Convolution(a, sk) to DotProd(a, sk) for LWE decryption.
	// Copy coefficients multiplied by X^{N-1} in reverse order:
	// a_{0} -a_{N-1} -a2_{N-2} ... -a_{1}
	tmp0 := aRLWEMod2N.Coeffs[0]
	tmp1 := acc.Value[1].Coeffs[0]
	tmp0[0] = tmp1[0]
	for j := 1; j < ringQLWE.N; j++ {
		tmp0[j] = -tmp1[ringQLWE.N-j] & mask
	}

	h.ModSwitchRLWETo2NLvl(ct.Level(), acc.Value[0], bRLWEMod2N)

	res = make(map[int]*rlwe.Ciphertext)

	tmpRGSW := rgsw.NewCiphertextNTT(h.paramsLUT, h.paramsLUT.MaxLevel())

	var prevIndex int
	for index := 0; index < ringQLWE.N; index++ {

		if lut, ok := lutPolyWihtSlotIndex[index]; ok {

			MulBySmallMonomialMod2N(mask, aRLWEMod2N, index-prevIndex)
			prevIndex = index

			a := aRLWEMod2N.Coeffs[0]
			b := bRLWEMod2N.Coeffs[0][index]

			// LWE = -as + m + e, a
			// LUT = LUT * X^{-as + m + e}
			ringQLUT.MulCoeffsMontgomery(lut, h.xPowMinusOne[b].Q, acc.Value[0])
			ringQLUT.Add(acc.Value[0], lut, acc.Value[0])
			acc.Value[1].Zero() // TODO remove

			for j := 0; j < ringQLWE.N; j++ {
				// RGSW[(X^{a} - 1) * sk_{j}[0] + (X^{-a} - 1) * sk_{j}[1] + 1]
				MulRGSWByXPowAlphaMinusOne(lutKey.SkPos[j], h.xPowMinusOne[a[j]], ringQPLUT, tmpRGSW)
				MulRGSWByXPowAlphaMinusOneAndAdd(lutKey.SkNeg[j], h.xPowMinusOne[-a[j]&mask], ringQPLUT, tmpRGSW)
				AddOneRGSW(lutKey.OneRGSW, ringQLUT, tmpRGSW)

				// LUT[RLWE] = LUT[RLWE] x RGSW[(X^{a} - 1) * sk_{j}[0] + (X^{-a} - 1) * sk_{j}[1] + 1]
				eval.ExternalProduct(acc, tmpRGSW, acc)
			}

			res[index] = acc.CopyNew()
		}

		// LUT[RLWE] = LUT[RLWE] * X^{m+e}
	}

	return
}

// AddOneRGSW adds one in plaintext on the output RGSW ciphertext.
func AddOneRGSW(oneRGSW []*ring.Poly, ringQ *ring.Ring, res *rgsw.Ciphertext) {
	nQ := res.LevelQ() + 1
	nP := res.LevelP() + 1

	if nP == 0 {
		nP = 1
	}

	for i := range res.Value[0].Value {
		for j := range res.Value[0].Value[i] {
			start, end := i*nP, (i+1)*nP
			if end > nQ {
				end = nQ
			}
			for k := start; k < end; k++ {
				ring.AddVecNoMod(res.Value[0].Value[i][j][0].Q.Coeffs[k], oneRGSW[j].Coeffs[k], res.Value[0].Value[i][j][0].Q.Coeffs[k])
				ring.AddVecNoMod(res.Value[1].Value[i][j][1].Q.Coeffs[k], oneRGSW[j].Coeffs[k], res.Value[1].Value[i][j][1].Q.Coeffs[k])
			}
		}
	}
}

//MulBySmallMonomialMod2N multiplies pol by x^n, with 0 <= n < N
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

// ModSwitchRLWETo2NLvl applys round(x * 2N / Q) to the coefficients of polQ and returns the result on pol2N.
func (h *Handler) ModSwitchRLWETo2NLvl(level int, polQ *ring.Poly, pol2N *ring.Poly) {
	coeffsBigint := make([]*big.Int, len(polQ.Coeffs[0]))

	ringQ := h.paramsLWE.RingQ()

	ringQ.PolyToBigintLvl(level, polQ, 1, coeffsBigint)

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
