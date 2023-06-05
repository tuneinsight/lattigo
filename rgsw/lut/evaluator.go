package lut

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v4/rgsw"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// Evaluator is a struct that stores the necessary
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
func NewEvaluator(paramsLUT, paramsLWE rlwe.Parameters, evk rlwe.EvaluationKeySet) (eval *Evaluator) {
	eval = new(Evaluator)
	eval.Evaluator = rgsw.NewEvaluator(paramsLUT, evk)
	eval.paramsLUT = paramsLUT
	eval.paramsLWE = paramsLWE

	ringQ := paramsLUT.RingQ()
	ringP := paramsLUT.RingP()

	eval.poolMod2N = [2]*ring.Poly{paramsLWE.RingQ().NewPoly(), paramsLWE.RingQ().NewPoly()}
	eval.accumulator = rlwe.NewCiphertext(paramsLUT, 1, paramsLUT.MaxLevel())
	eval.accumulator.IsNTT = true // This flag is always true

	N := ringQ.N()

	// Compute X^{n} -  1 from 0 to 2N LWE
	oneNTTMFormQ := ringQ.NewPoly()
	for i := range oneNTTMFormQ.Coeffs {

		coeffs := oneNTTMFormQ.Coeffs[i]

		s := ringQ.SubRings[i]

		for j := 0; j < N; j++ {
			coeffs[j] = ring.MForm(1, s.Modulus, s.BRedConstant)
		}
	}

	eval.xPowMinusOne = make([]ringqp.Poly, 2*N)
	for i := 0; i < N; i++ {
		eval.xPowMinusOne[i].Q = ringQ.NewPoly()
		eval.xPowMinusOne[i+N].Q = ringQ.NewPoly()
		if i == 0 || i == 1 {
			for j, s := range ringQ.SubRings {
				eval.xPowMinusOne[i].Q.Coeffs[j][i] = ring.MForm(1, s.Modulus, s.BRedConstant)
			}

			ringQ.NTT(eval.xPowMinusOne[i].Q, eval.xPowMinusOne[i].Q)

			// Negacyclic wrap-around for n > N
			ringQ.Neg(eval.xPowMinusOne[i].Q, eval.xPowMinusOne[i+N].Q)

		} else {
			ringQ.MulCoeffsMontgomery(eval.xPowMinusOne[1].Q, eval.xPowMinusOne[i-1].Q, eval.xPowMinusOne[i].Q) // X^{n} = X^{1} * X^{n-1}

			// Negacyclic wrap-around for n > N
			ringQ.Neg(eval.xPowMinusOne[i].Q, eval.xPowMinusOne[i+N].Q) // X^{2n} = -X^{1} * X^{n-1}
		}
	}

	// Subtract -1 in NTT
	for i := 0; i < 2*N; i++ {
		ringQ.Sub(eval.xPowMinusOne[i].Q, oneNTTMFormQ, eval.xPowMinusOne[i].Q) // X^{n} - 1
	}

	if ringP != nil {
		oneNTTMFormP := ringP.NewPoly()
		for i := range oneNTTMFormP.Coeffs {

			coeffs := oneNTTMFormP.Coeffs[i]

			table := ringP.SubRings[i]

			for j := 0; j < N; j++ {
				coeffs[j] = ring.MForm(1, table.Modulus, table.BRedConstant)
			}
		}

		for i := 0; i < N; i++ {
			eval.xPowMinusOne[i].P = ringP.NewPoly()
			eval.xPowMinusOne[i+N].P = ringP.NewPoly()
			if i == 0 || i == 1 {
				for j, table := range ringP.SubRings {
					eval.xPowMinusOne[i].P.Coeffs[j][i] = ring.MForm(1, table.Modulus, table.BRedConstant)
				}

				ringP.NTT(eval.xPowMinusOne[i].P, eval.xPowMinusOne[i].P)

				// Negacyclic wrap-around for n > N
				ringP.Neg(eval.xPowMinusOne[i].P, eval.xPowMinusOne[i+N].P)

			} else {
				// X^{n} = X^{1} * X^{n-1}
				ringP.MulCoeffsMontgomery(eval.xPowMinusOne[1].P, eval.xPowMinusOne[i-1].P, eval.xPowMinusOne[i].P)

				// Negacyclic wrap-around for n > N
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

	eval.tmpRGSW = rgsw.NewCiphertext(paramsLUT, levelQ, levelP, decompRNS, decompPw2)
	eval.one = rgsw.NewPlaintext(paramsLUT, uint64(1), levelQ, levelP, paramsLUT.Pow2Base(), decompPw2)

	return
}

// EvaluateAndRepack extracts on the fly LWE samples, evaluates the provided LUT on the LWE and repacks everything into a single rlwe.Ciphertext.
// ct : a rlwe Ciphertext with coefficient encoded values at level 0
// lutPolyWithSlotIndex : a map with [slot_index] -> LUT
// repackIndex : a map with [slot_index_have] -> slot_index_want
// lutKey : LUTKey
// Returns a *rlwe.Ciphertext
func (eval *Evaluator) EvaluateAndRepack(ct *rlwe.Ciphertext, lutPolyWithSlotIndex map[int]*ring.Poly, repackIndex map[int]int, key EvaluationKey) (res *rlwe.Ciphertext) {
	cts := eval.Evaluate(ct, lutPolyWithSlotIndex, key)

	ciphertexts := make(map[int]*rlwe.Ciphertext)

	for i := range cts {
		ciphertexts[repackIndex[i]] = cts[i]
	}

	return eval.Pack(ciphertexts, eval.Parameters().LogN(), true)
}

// Evaluate extracts on the fly LWE samples and evaluates the provided LUT on the LWE.
// ct : a rlwe Ciphertext with coefficient encoded values at level 0
// lutPolyWithSlotIndex : a map with [slot_index] -> LUT
// lutKey : lut.Key
// Returns a map[slot_index] -> LUT(ct[slot_index])
func (eval *Evaluator) Evaluate(ct *rlwe.Ciphertext, lutPolyWithSlotIndex map[int]*ring.Poly, key EvaluationKey) (res map[int]*rlwe.Ciphertext) {

	bRLWEMod2N := eval.poolMod2N[0]
	aRLWEMod2N := eval.poolMod2N[1]

	acc := eval.accumulator

	levelQ := key.SkPos[0].LevelQ()
	levelP := key.SkPos[0].LevelP()

	ringQPLUT := *eval.paramsLUT.RingQP().AtLevel(levelQ, levelP)
	ringQLUT := ringQPLUT.RingQ

	ringQLWE := eval.paramsLWE.RingQ().AtLevel(ct.Level())

	// mod 2N
	mask := uint64(ringQLUT.N()<<1) - 1

	if ct.IsNTT {
		ringQLWE.INTT(ct.Value[0], acc.Value[0])
		ringQLWE.INTT(ct.Value[1], acc.Value[1])
	} else {
		ring.CopyLvl(ct.Level(), ct.Value[0], acc.Value[0])
		ring.CopyLvl(ct.Level(), ct.Value[1], acc.Value[1])
	}

	// Switch modulus from Q to 2N
	eval.ModSwitchRLWETo2NLvl(ct.Level(), acc.Value[1], acc.Value[1])

	// Conversion from Convolution(a, sk) to DotProd(a, sk) for LWE decryption.
	// Copy coefficients multiplied by X^{N-1} in reverse order:
	// a_{0} -a_{N-1} -a2_{N-2} ... -a_{1}
	tmp0 := aRLWEMod2N.Coeffs[0]
	tmp1 := acc.Value[1].Coeffs[0]
	tmp0[0] = tmp1[0]
	NLWE := ringQLWE.N()
	for j := 1; j < NLWE; j++ {
		tmp0[j] = -tmp1[ringQLWE.N()-j] & mask
	}

	eval.ModSwitchRLWETo2NLvl(ct.Level(), acc.Value[0], bRLWEMod2N)

	res = make(map[int]*rlwe.Ciphertext)

	var prevIndex int
	for index := 0; index < NLWE; index++ {

		if lut, ok := lutPolyWithSlotIndex[index]; ok {

			MulBySmallMonomialMod2N(mask, aRLWEMod2N, index-prevIndex)
			prevIndex = index

			a := aRLWEMod2N.Coeffs[0]
			b := bRLWEMod2N.Coeffs[0][index]

			// LWE = -as + m + e, a
			// LUT = LUT * X^{-as + m + e}
			ringQLUT.MulCoeffsMontgomery(lut, eval.xPowMinusOne[b].Q, acc.Value[0])
			ringQLUT.Add(acc.Value[0], lut, acc.Value[0])
			acc.Value[1].Zero()

			for j := 0; j < NLWE; j++ {
				// RGSW[(X^{a} - 1) * sk_{j}[0] + (X^{-a} - 1) * sk_{j}[1] + 1]
				rgsw.MulByXPowAlphaMinusOneLazy(key.SkPos[j], eval.xPowMinusOne[a[j]], ringQPLUT, eval.tmpRGSW)
				rgsw.MulByXPowAlphaMinusOneThenAddLazy(key.SkNeg[j], eval.xPowMinusOne[-a[j]&mask], ringQPLUT, eval.tmpRGSW)
				rgsw.AddLazy(eval.one, ringQPLUT, eval.tmpRGSW)

				// LUT[RLWE] = LUT[RLWE] x RGSW[(X^{a} - 1) * sk_{j}[0] + (X^{-a} - 1) * sk_{j}[1] + 1]
				eval.ExternalProduct(acc, eval.tmpRGSW, acc)

			}
			res[index] = acc.CopyNew()

			if !eval.paramsLUT.NTTFlag() {
				ringQLUT.INTT(res[index].Value[0], res[index].Value[0])
				ringQLUT.INTT(res[index].Value[1], res[index].Value[1])
				res[index].IsNTT = false
			}
		}
		// LUT[RLWE] = LUT[RLWE] * X^{m+e}
	}

	return
}

// ModSwitchRLWETo2NLvl applies round(x * 2N / Q) to the coefficients of polQ and returns the result on pol2N.
func (eval *Evaluator) ModSwitchRLWETo2NLvl(level int, polQ *ring.Poly, pol2N *ring.Poly) {
	coeffsBigint := make([]*big.Int, len(polQ.Coeffs[0]))

	ringQ := eval.paramsLWE.RingQ().AtLevel(level)

	ringQ.PolyToBigint(polQ, 1, coeffsBigint)

	QBig := ringQ.ModulusAtLevel[level]

	twoN := uint64(eval.paramsLUT.N() << 1)
	twoNBig := bignum.NewInt(twoN)
	tmp := pol2N.Coeffs[0]
	N := ringQ.N()
	for i := 0; i < N; i++ {
		coeffsBigint[i].Mul(coeffsBigint[i], twoNBig)
		bignum.DivRound(coeffsBigint[i], QBig, coeffsBigint[i])
		tmp[i] = coeffsBigint[i].Uint64() & (twoN - 1)
	}
}
