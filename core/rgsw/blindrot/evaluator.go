package blindrot

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/core/rgsw"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// Evaluator is a struct that stores the necessary
// data to handle LWE <-> RLWE conversion and
// blind rotations.
type Evaluator struct {
	*rgsw.Evaluator
	paramsBR  rlwe.Parameters
	paramsLWE rlwe.Parameters

	poolMod2N [2]ring.Poly

	accumulator *rlwe.Ciphertext

	galoisGenDiscreteLog map[uint64]int
}

// NewEvaluator instantiates a new [Evaluator].
func NewEvaluator(paramsBR, paramsLWE rlwe.ParameterProvider) (eval *Evaluator) {
	eval = new(Evaluator)
	eval.Evaluator = rgsw.NewEvaluator(paramsBR, nil)
	eval.paramsBR = *paramsBR.GetRLWEParameters()
	eval.paramsLWE = *paramsLWE.GetRLWEParameters()

	eval.poolMod2N = [2]ring.Poly{eval.paramsLWE.RingQ().NewPoly(), eval.paramsLWE.RingQ().NewPoly()}
	eval.accumulator = rlwe.NewCiphertext(paramsBR, 1, eval.paramsBR.MaxLevel())
	eval.accumulator.IsNTT = true // This flag is always true

	// Generates a map for the discrete log of (+/- 1) * GaloisGen^k for 0 <= k < N-1.
	// galoisGenDiscreteLog: map[+/-G^{k} mod 2N] = k
	eval.galoisGenDiscreteLog = getGaloisElementInverseMap(ring.GaloisGen, eval.paramsBR.N())

	return
}

// Evaluate extracts on the fly LWE samples and evaluates the provided blind rotation on the LWE.
// testPolyWithSlotIndex : a map with [slot_index] -> blind rotation
// Returns a map[slot_index] -> BlindRotate(ct[slot_index])
func (eval *Evaluator) Evaluate(ct *rlwe.Ciphertext, testPolyWithSlotIndex map[int]*ring.Poly, BRK BlindRotationEvaluationKeySet) (res map[int]*rlwe.Ciphertext, err error) {

	bRLWEMod2N := eval.poolMod2N[0]
	aRLWEMod2N := eval.poolMod2N[1]

	acc := eval.accumulator

	brk, err := BRK.GetBlindRotationKey(0)

	if err != nil {
		return nil, err
	}

	ringQBR := eval.paramsBR.RingQ().AtLevel(brk.LevelQ())
	ringQLWE := eval.paramsLWE.RingQ().AtLevel(ct.Level())

	if ct.IsNTT {
		ringQLWE.INTT(ct.Value[0], acc.Value[0])
		ringQLWE.INTT(ct.Value[1], acc.Value[1])
	} else {
		acc.Value[0].CopyLvl(ct.Level(), ct.Value[0])
		acc.Value[1].CopyLvl(ct.Level(), ct.Value[1])
	}

	// Switch modulus from Q to 2N and ensure they are odd
	eval.modSwitchRLWETo2NLvl(ct.Level(), acc.Value[1], acc.Value[1], true)

	// Conversion from Convolution(a, sk) to DotProd(a, sk) for LWE decryption.
	// Copy coefficients multiplied by X^{N-1} in reverse order:
	// a_{0} -a_{N-1} -a2_{N-2} ... -a_{1}
	tmp0 := aRLWEMod2N.Coeffs[0]
	tmp1 := acc.Value[1].Coeffs[0]
	tmp0[0] = tmp1[0]
	NLWE := ringQLWE.N()
	mask := uint64(ringQBR.N()<<1) - 1
	for j := 1; j < NLWE; j++ {
		tmp0[j] = -tmp1[ringQLWE.N()-j] & mask
	}

	// Switch modulus from Q to 2N
	eval.modSwitchRLWETo2NLvl(ct.Level(), acc.Value[0], bRLWEMod2N, false)

	res = make(map[int]*rlwe.Ciphertext)

	var prevIndex int
	for index := 0; index < NLWE; index++ {

		if testPoly, ok := testPolyWithSlotIndex[index]; ok {

			mulBySmallMonomialMod2N(mask, aRLWEMod2N, index-prevIndex)
			prevIndex = index

			a := aRLWEMod2N.Coeffs[0]
			b := bRLWEMod2N.Coeffs[0][index]

			// Line 2 of Algorithm 7 of https://eprint.iacr.org/2022/198
			// Acc = (f(X^{-g}) * X^{-g * b}, 0)
			Xb := ringQBR.NewMonomialXi(int(b))
			ringQBR.NTT(Xb, Xb)
			ringQBR.MForm(Xb, Xb)
			ringQBR.MulCoeffsMontgomery(*testPoly, Xb, acc.Value[1]) // use unused buffer because AutomorphismNTT is not in place
			ringQBR.AutomorphismNTT(acc.Value[1], ringQBR.NthRoot()-ring.GaloisGen, acc.Value[0])
			acc.Value[1].Zero()

			// Line 3 of Algorithm 7 https://eprint.iacr.org/2022/198 (Algorithm 3 of https://eprint.iacr.org/2022/198)
			if err = eval.BlindRotateCore(a, acc, BRK); err != nil {
				return nil, fmt.Errorf("BlindRotateCore: %s", err)
			}

			// f(X) * X^{b + <a, s>}
			res[index] = acc.CopyNew()

			if !eval.paramsBR.NTTFlag() {
				ringQBR.INTT(res[index].Value[0], res[index].Value[0])
				ringQBR.INTT(res[index].Value[1], res[index].Value[1])
				res[index].IsNTT = false
			}
		}
	}

	return
}

// BlindRotateCore implements Algorithm 3 of https://eprint.iacr.org/2022/198
func (eval *Evaluator) BlindRotateCore(a []uint64, acc *rlwe.Ciphertext, BRK BlindRotationEvaluationKeySet) (err error) {

	evk, err := BRK.GetEvaluationKeySet()

	if err != nil {
		return err
	}

	eval.Evaluator = eval.Evaluator.WithKey(evk)

	// GaloisElement(k) = GaloisGen^{k} mod 2N
	GaloisElement := eval.paramsBR.GaloisElement

	// Maps a[i] to (+/-) g^{k} mod 2N
	discreteLogSets := eval.getDiscreteLogSets(a)

	Nhalf := eval.paramsBR.N() >> 1

	// Algorithm 3 of https://eprint.iacr.org/2022/198
	var v int
	// Lines 3 to 9 (negative set of a[i] = -g^{k} mod 2N)
	for i := Nhalf - 1; i > 0; i-- {
		if v, err = eval.evaluateFromDiscreteLogSets(GaloisElement, discreteLogSets, -i, v, acc, BRK); err != nil {
			return
		}
	}

	// Line 10 (0 in the negative set is 2N)
	if _, err = eval.evaluateFromDiscreteLogSets(GaloisElement, discreteLogSets, eval.paramsBR.N()<<1, 0, acc, BRK); err != nil {
		return
	}

	// Line 12
	// acc = acc(X^{-g})
	if err = eval.Automorphism(acc, eval.paramsBR.RingQ().NthRoot()-ring.GaloisGen, acc); err != nil {
		return
	}

	// Lines 13 - 19 (positive set of a[i] = g^{k} mod 2N)
	for i := Nhalf - 1; i > 0; i-- {
		if v, err = eval.evaluateFromDiscreteLogSets(GaloisElement, discreteLogSets, i, v, acc, BRK); err != nil {
			return
		}
	}

	// Lines 20 - 21 (0 in the positive set is 0)
	if _, err = eval.evaluateFromDiscreteLogSets(GaloisElement, discreteLogSets, 0, 0, acc, BRK); err != nil {
		return
	}

	return
}

// evaluateFromDiscreteLogSets loops of Algorithm 3 of https://eprint.iacr.org/2022/198
func (eval *Evaluator) evaluateFromDiscreteLogSets(GaloisElement func(k int) (galEl uint64), sets map[int][]int, k, v int, acc *rlwe.Ciphertext, BRK BlindRotationEvaluationKeySet) (int, error) {

	// Checks if k is in the discrete log sets
	if set, ok := sets[k]; ok {

		// First condition of line 7 or 17
		if v != 0 {

			if err := eval.Automorphism(acc, GaloisElement(v), acc); err != nil {
				return v, err
			}

			v = 0
		}

		for _, j := range set {

			brk, err := BRK.GetBlindRotationKey(j)
			if err != nil {
				return v, err
			}

			// acc = acc * RGSW(X^{s[j]})
			eval.ExternalProduct(acc, brk, acc)
		}
	}

	v++

	// Second and third conditions of line 7 or 17
	if v == windowSize || k == 1 {

		if err := eval.Automorphism(acc, GaloisElement(v), acc); err != nil {
			return v, err
		}

		v = 0
	}

	return v, nil
}

// getGaloisElementInverseMap generates a map [(+/-) g^{k} mod 2N] = +/- k
func getGaloisElementInverseMap(GaloisGen uint64, N int) (GaloisGenDiscreteLog map[uint64]int) {

	twoN := N << 1
	NHalf := N >> 1
	mask := uint64(twoN - 1)

	GaloisGenDiscreteLog = map[uint64]int{}

	var pow uint64 = 1
	for i := 0; i < NHalf; i++ {
		GaloisGenDiscreteLog[pow] = i
		GaloisGenDiscreteLog[uint64(twoN)-pow] = -i
		pow *= GaloisGen
		pow &= mask
	}

	return
}

// getDiscreteLogSets returns map[+/-k] = [i...] for a[0 <= i < N] = {(+/-) g^{k} mod 2N for +/- k}
func (eval *Evaluator) getDiscreteLogSets(a []uint64) (discreteLogSets map[int][]int) {

	GaloisGenDiscreteLog := eval.galoisGenDiscreteLog

	// Maps (2*N*a[i]/QLWE) to -N/2 < k <= N/2 for a[i] = (+/- 1) * g^{k}
	discreteLogSets = map[int][]int{}
	for i, ai := range a {

		if ai&1 != 1 && ai != 0 {
			panic("getDiscreteLogSets: a[i] is not odd and thus not an element of Z_{2N}^{*} -> a[i] = (+/- 1) * g^{k} does not exist.")
		}

		dlog := GaloisGenDiscreteLog[ai]

		if _, ok := discreteLogSets[dlog]; !ok {
			discreteLogSets[dlog] = []int{i}
		} else {
			discreteLogSets[dlog] = append(discreteLogSets[dlog], i)
		}
	}

	return
}

// modSwitchRLWETo2NLvl applies round(x * 2N / Q) to the coefficients of polQ and returns the result on pol2N.
// makeOdd ensures that output coefficients are odd by xoring with 1 (if not already zero).
func (eval *Evaluator) modSwitchRLWETo2NLvl(level int, polQ, pol2N ring.Poly, makeOdd bool) {
	coeffsBigint := make([]*big.Int, len(polQ.Coeffs[0]))

	ringQ := eval.paramsLWE.RingQ().AtLevel(level)

	ringQ.PolyToBigint(polQ, 1, coeffsBigint)

	QBig := ringQ.ModulusAtLevel[level]

	twoN := uint64(eval.paramsBR.N() << 1)
	twoNBig := bignum.NewInt(twoN)
	tmp := pol2N.Coeffs[0]
	N := ringQ.N()
	for i := 0; i < N; i++ {
		coeffsBigint[i].Mul(coeffsBigint[i], twoNBig)
		bignum.DivRound(coeffsBigint[i], QBig, coeffsBigint[i])
		tmp[i] = coeffsBigint[i].Uint64() & (twoN - 1)

		if makeOdd && tmp[i]&1 == 0 && tmp[i] != 0 {
			tmp[i] ^= 1
		}
	}
}
