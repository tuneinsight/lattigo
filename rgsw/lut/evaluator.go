package lut

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v3/rgsw"
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
)

// Evaluator is a struct that stores necessary
// data to handle LWE <-> RLWE conversion and
// LUT evaluation.
type Evaluator struct {
	*rgsw.Evaluator
	paramsLUT rlwe.Parameters
	paramsLWE rlwe.Parameters

	xPowMinusOne []ringqp.Poly //X^n - 1 from 0 to 2N LWE

	poolMod2N [2]*ring.Poly

	accumulator *rlwe.Ciphertext
	Sk          *rlwe.SecretKey

	tmpRGSW *rgsw.Ciphertext

	one *rgsw.Plaintext
}

// NewEvaluator creates a new Handler
func NewEvaluator(paramsLUT, paramsLWE rlwe.Parameters, rtks *rlwe.RotationKeySet) (eval *Evaluator) {
	eval = new(Evaluator)
	eval.Evaluator = rgsw.NewEvaluator(paramsLUT, &rlwe.EvaluationKey{Rtks: rtks})
	eval.paramsLUT = paramsLUT
	eval.paramsLWE = paramsLWE

	ringQ := paramsLUT.RingQ()
	ringP := paramsLUT.RingP()

	eval.poolMod2N = [2]*ring.Poly{paramsLWE.RingQ().NewPolyLvl(0), paramsLWE.RingQ().NewPolyLvl(0)}
	eval.accumulator = rlwe.NewCiphertextNTT(paramsLUT, 1, paramsLUT.MaxLevel())

	// Compute X^{n} -  1 from 0 to 2N LWE
	oneNTTMFormQ := ringQ.NewPoly()
	for i := range ringQ.Modulus {
		for j := 0; j < ringQ.N; j++ {
			oneNTTMFormQ.Coeffs[i][j] = ring.MForm(1, ringQ.Modulus[i], ringQ.BredParams[i])
		}
	}

	N := ringQ.N

	eval.xPowMinusOne = make([]ringqp.Poly, 2*N)
	for i := 0; i < N; i++ {
		eval.xPowMinusOne[i].Q = ringQ.NewPoly()
		eval.xPowMinusOne[i+N].Q = ringQ.NewPoly()
		if i == 0 || i == 1 {
			for j := range ringQ.Modulus {
				eval.xPowMinusOne[i].Q.Coeffs[j][i] = ring.MForm(1, ringQ.Modulus[j], ringQ.BredParams[j])
			}

			ringQ.NTT(eval.xPowMinusOne[i].Q, eval.xPowMinusOne[i].Q)

			// Negacyclic wrap-arround for n > N
			ringQ.Neg(eval.xPowMinusOne[i].Q, eval.xPowMinusOne[i+N].Q)

		} else {
			ringQ.MulCoeffsMontgomery(eval.xPowMinusOne[1].Q, eval.xPowMinusOne[i-1].Q, eval.xPowMinusOne[i].Q) // X^{n} = X^{1} * X^{n-1}

			// Negacyclic wrap-arround for n > N
			ringQ.Neg(eval.xPowMinusOne[i].Q, eval.xPowMinusOne[i+N].Q) // X^{2n} = -X^{1} * X^{n-1}
		}
	}

	// Subtract -1 in NTT
	for i := 0; i < 2*N; i++ {
		ringQ.Sub(eval.xPowMinusOne[i].Q, oneNTTMFormQ, eval.xPowMinusOne[i].Q) // X^{n} - 1
	}

	if ringP != nil {
		oneNTTMFormP := ringP.NewPoly()
		for i := range ringP.Modulus {
			for j := 0; j < ringP.N; j++ {
				oneNTTMFormP.Coeffs[i][j] = ring.MForm(1, ringP.Modulus[i], ringP.BredParams[i])
			}
		}

		for i := 0; i < N; i++ {
			eval.xPowMinusOne[i].P = ringP.NewPoly()
			eval.xPowMinusOne[i+N].P = ringP.NewPoly()
			if i == 0 || i == 1 {
				for j := range ringP.Modulus {
					eval.xPowMinusOne[i].P.Coeffs[j][i] = ring.MForm(1, ringP.Modulus[j], ringP.BredParams[j])
				}

				ringP.NTT(eval.xPowMinusOne[i].P, eval.xPowMinusOne[i].P)

				// Negacyclic wrap-arround for n > N
				ringP.Neg(eval.xPowMinusOne[i].P, eval.xPowMinusOne[i+N].P)

			} else {
				// X^{n} = X^{1} * X^{n-1}
				ringP.MulCoeffsMontgomery(eval.xPowMinusOne[1].P, eval.xPowMinusOne[i-1].P, eval.xPowMinusOne[i].P)

				// Negacyclic wrap-arround for n > N
				// X^{2n} = -X^{1} * X^{n-1}
				ringP.Neg(eval.xPowMinusOne[i].P, eval.xPowMinusOne[i+N].P)
			}
		}

		// Subtract -1 in NTT
		for i := 0; i < 2*N; i++ {
			// X^{n} - 1
			ringP.Sub(eval.xPowMinusOne[i].P, oneNTTMFormP, eval.xPowMinusOne[i].P)
		}
	}

	levelQ := paramsLUT.QCount() - 1
	levelP := paramsLUT.PCount() - 1
	decompRNS := paramsLUT.DecompRNS(levelQ, levelP)
	decompPw2 := paramsLUT.DecompPw2(levelQ, levelP)
	ringQP := paramsLUT.RingQP()
	eval.tmpRGSW = rgsw.NewCiphertext(levelQ, levelP, decompRNS, decompPw2, *ringQP)

	eval.one = rgsw.NewPlaintext(uint64(1), levelQ, levelP, paramsLUT.Pow2Base(), decompPw2, *ringQP)

	return
}

// EvaluateAndRepack extracts on the fly LWE samples and evaluate the provided LUT on the LWE and repacks everything into a single rlwe.Ciphertext.
// ct : a rlwe Ciphertext with coefficient encoded values at level 0
// lutPolyWihtSlotIndex : a map with [slot_index] -> LUT
// repackIndex : a map with [slot_index_have] -> slot_index_want
// lutKey : LUTKey
// Returns a *rlwe.Ciphertext
func (eval *Evaluator) EvaluateAndRepack(ct *rlwe.Ciphertext, lutPolyWihtSlotIndex map[int]*ring.Poly, repackIndex map[int]int, key EvaluationKey) (res *rlwe.Ciphertext) {
	cts := eval.Evaluate(ct, lutPolyWihtSlotIndex, key)

	ciphertexts := make(map[int]*rlwe.Ciphertext)

	for i := range cts {
		ciphertexts[repackIndex[i]] = cts[i]
	}

	return eval.MergeRLWE(ciphertexts)
}

// Evaluate extracts on the fly LWE samples and evaluate the provided LUT on the LWE.
// ct : a rlwe Ciphertext with coefficient encoded values at level 0
// lutPolyWihtSlotIndex : a map with [slot_index] -> LUT
// lutKey : lut.Key
// Returns a map[slot_index] -> LUT(ct[slot_index])
func (eval *Evaluator) Evaluate(ct *rlwe.Ciphertext, lutPolyWihtSlotIndex map[int]*ring.Poly, key EvaluationKey) (res map[int]*rlwe.Ciphertext) {

	bRLWEMod2N := eval.poolMod2N[0]
	aRLWEMod2N := eval.poolMod2N[1]

	acc := eval.accumulator

	ringQLUT := eval.paramsLUT.RingQ()
	ringQLWE := eval.paramsLWE.RingQ()
	ringQPLUT := *eval.paramsLUT.RingQP()

	// mod 2N
	mask := uint64(ringQLUT.N<<1) - 1

	ringQLWE.InvNTTLvl(ct.Level(), ct.Value[0], acc.Value[0])
	ringQLWE.InvNTTLvl(ct.Level(), ct.Value[1], acc.Value[1])

	// Switch modulus from Q to 2N
	eval.ModSwitchRLWETo2NLvl(ct.Level(), acc.Value[1], acc.Value[1])

	// Conversion from Convolution(a, sk) to DotProd(a, sk) for LWE decryption.
	// Copy coefficients multiplied by X^{N-1} in reverse order:
	// a_{0} -a_{N-1} -a2_{N-2} ... -a_{1}
	tmp0 := aRLWEMod2N.Coeffs[0]
	tmp1 := acc.Value[1].Coeffs[0]
	tmp0[0] = tmp1[0]
	for j := 1; j < ringQLWE.N; j++ {
		tmp0[j] = -tmp1[ringQLWE.N-j] & mask
	}

	eval.ModSwitchRLWETo2NLvl(ct.Level(), acc.Value[0], bRLWEMod2N)

	levelQ := key.SkPos[0].LevelQ()
	levelP := key.SkPos[0].LevelP()

	res = make(map[int]*rlwe.Ciphertext)

	var prevIndex int
	for index := 0; index < ringQLWE.N; index++ {

		if lut, ok := lutPolyWihtSlotIndex[index]; ok {

			MulBySmallMonomialMod2N(mask, aRLWEMod2N, index-prevIndex)
			prevIndex = index

			a := aRLWEMod2N.Coeffs[0]
			b := bRLWEMod2N.Coeffs[0][index]

			// LWE = -as + m + e, a
			// LUT = LUT * X^{-as + m + e}
			ringQLUT.MulCoeffsMontgomery(lut, eval.xPowMinusOne[b].Q, acc.Value[0])
			ringQLUT.Add(acc.Value[0], lut, acc.Value[0])
			acc.Value[1].Zero() // TODO remove

			for j := 0; j < ringQLWE.N; j++ {
				// RGSW[(X^{a} - 1) * sk_{j}[0] + (X^{-a} - 1) * sk_{j}[1] + 1]
				rgsw.MulByXPowAlphaMinusOneConstantLvl(levelQ, levelP, key.SkPos[j], eval.xPowMinusOne[a[j]], ringQPLUT, eval.tmpRGSW)
				rgsw.MulByXPowAlphaMinusOneAndAddNoModLvl(levelQ, levelP, key.SkNeg[j], eval.xPowMinusOne[-a[j]&mask], ringQPLUT, eval.tmpRGSW)
				rgsw.AddNoModLvl(levelQ, levelP, eval.one, ringQPLUT, eval.tmpRGSW)

				// LUT[RLWE] = LUT[RLWE] x RGSW[(X^{a} - 1) * sk_{j}[0] + (X^{-a} - 1) * sk_{j}[1] + 1]
				eval.ExternalProduct(acc, eval.tmpRGSW, acc)
			}

			res[index] = acc.CopyNew()
		}

		// LUT[RLWE] = LUT[RLWE] * X^{m+e}
	}

	return
}

// ModSwitchRLWETo2NLvl applys round(x * 2N / Q) to the coefficients of polQ and returns the result on pol2N.
func (eval *Evaluator) ModSwitchRLWETo2NLvl(level int, polQ *ring.Poly, pol2N *ring.Poly) {
	coeffsBigint := make([]*big.Int, len(polQ.Coeffs[0]))

	ringQ := eval.paramsLWE.RingQ()

	ringQ.PolyToBigintLvl(level, polQ, 1, coeffsBigint)

	QBig := ring.NewUint(1)
	for i := 0; i < level+1; i++ {
		QBig.Mul(QBig, ring.NewUint(ringQ.Modulus[i]))
	}

	twoN := uint64(eval.paramsLUT.N() << 1)
	twoNBig := ring.NewUint(twoN)
	tmp := pol2N.Coeffs[0]
	for i := 0; i < ringQ.N; i++ {
		coeffsBigint[i].Mul(coeffsBigint[i], twoNBig)
		ring.DivRound(coeffsBigint[i], QBig, coeffsBigint[i])
		tmp[i] = coeffsBigint[i].Uint64() & (twoN - 1)
	}
}
