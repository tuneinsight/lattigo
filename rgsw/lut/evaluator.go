package lut

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v4/rgsw"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// Evaluator is a struct that stores the necessary
// data to handle LWE <-> RLWE conversion and
// blind rotations.
type Evaluator struct {
	*rgsw.Evaluator
	paramsLUT rlwe.Parameters
	paramsLWE rlwe.Parameters

	poolMod2N [2]ring.Poly

	accumulator *rlwe.Ciphertext

	galoisGenDiscretLog map[uint64]int
}

// NewEvaluator instaniates a new Evaluator.
func NewEvaluator(paramsLUT, paramsLWE rlwe.Parameters) (eval *Evaluator) {
	eval = new(Evaluator)
	eval.Evaluator = rgsw.NewEvaluator(paramsLUT, nil)
	eval.paramsLUT = paramsLUT
	eval.paramsLWE = paramsLWE

	eval.poolMod2N = [2]ring.Poly{paramsLWE.RingQ().NewPoly(), paramsLWE.RingQ().NewPoly()}
	eval.accumulator = rlwe.NewCiphertext(paramsLUT, 1, paramsLUT.MaxLevel())
	eval.accumulator.IsNTT = true // This flag is always true

	// Generates a map for the discret log of (+/- 1) * GaloisGen^k for 0 <= k < N-1.
	// galoisGenDiscretLog: map[+/-G^{k} mod 2N] = k
	eval.galoisGenDiscretLog = getGaloisElementInverseMap(ring.GaloisGen, paramsLUT.N())

	return
}

// EvaluateAndRepack extracts on the fly LWE samples, evaluates the provided LUT on the LWE and repacks everything into a single rlwe.Ciphertext.
// lutPolyWithSlotIndex : a map with [slot_index] -> LUT
// repackIndex : a map with [slot_index_have] -> slot_index_want
func (eval *Evaluator) EvaluateAndRepack(ct *rlwe.Ciphertext, lutPolyWithSlotIndex map[int]*ring.Poly, repackIndex map[int]int, key BlindRotatationEvaluationKeySet, repackKey rlwe.EvaluationKeySet) (res *rlwe.Ciphertext, err error) {
	cts, err := eval.Evaluate(ct, lutPolyWithSlotIndex, key)

	if err != nil {
		return nil, err
	}

	ciphertexts := make(map[int]*rlwe.Ciphertext)

	for i := range cts {
		ciphertexts[repackIndex[i]] = cts[i]
	}

	eval.Evaluator = eval.Evaluator.WithKey(repackKey)

	return eval.Pack(ciphertexts, eval.paramsLUT.LogN(), true)
}

// Evaluate extracts on the fly LWE samples and evaluates the provided LUT on the LWE.
// lutPolyWithSlotIndex : a map with [slot_index] -> LUT
// Returns a map[slot_index] -> LUT(ct[slot_index])
func (eval *Evaluator) Evaluate(ct *rlwe.Ciphertext, lutPolyWithSlotIndex map[int]*ring.Poly, key BlindRotatationEvaluationKeySet) (res map[int]*rlwe.Ciphertext, err error) {

	evk, err := key.GetEvaluationKeySet()

	if err != nil {
		return nil, err
	}

	eval.Evaluator = eval.Evaluator.WithKey(evk)

	bRLWEMod2N := eval.poolMod2N[0]
	aRLWEMod2N := eval.poolMod2N[1]

	acc := eval.accumulator

	brk, err := key.GetBlingRotateKey(0)

	if err != nil {
		return nil, err
	}

	ringQLUT := eval.paramsLUT.RingQ().AtLevel(brk.LevelQ())
	ringQLWE := eval.paramsLWE.RingQ().AtLevel(ct.Level())

	if ct.IsNTT {
		ringQLWE.INTT(ct.Value[0], acc.Value[0])
		ringQLWE.INTT(ct.Value[1], acc.Value[1])
	} else {
		ring.CopyLvl(ct.Level(), ct.Value[0], acc.Value[0])
		ring.CopyLvl(ct.Level(), ct.Value[1], acc.Value[1])
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
	mask := uint64(ringQLUT.N()<<1) - 1
	for j := 1; j < NLWE; j++ {
		tmp0[j] = -tmp1[ringQLWE.N()-j] & mask
	}

	// Switch modulus from Q to 2N
	eval.modSwitchRLWETo2NLvl(ct.Level(), acc.Value[0], bRLWEMod2N, false)

	res = make(map[int]*rlwe.Ciphertext)

	var prevIndex int
	for index := 0; index < NLWE; index++ {

		if lutpoly, ok := lutPolyWithSlotIndex[index]; ok {

			mulBySmallMonomialMod2N(mask, aRLWEMod2N, index-prevIndex)
			prevIndex = index

			a := aRLWEMod2N.Coeffs[0]
			b := bRLWEMod2N.Coeffs[0][index]

			// Line 2 of Algorithm 7 of https://eprint.iacr.org/2022/198
			// Acc = (f(X^{-g}) * X^{-g * b}, 0)
			Xb := ringQLUT.NewMonomialXi(int(b))
			ringQLUT.NTT(Xb, Xb)
			ringQLUT.MForm(Xb, Xb)
			ringQLUT.MulCoeffsMontgomery(*lutpoly, Xb, acc.Value[1]) // use unused buffer because AutomorphismNTT is not in place
			ringQLUT.AutomorphismNTT(acc.Value[1], ringQLUT.NthRoot()-ring.GaloisGen, acc.Value[0])
			acc.Value[1].Zero()

			// Line 3 of Algorithm 7 https://eprint.iacr.org/2022/198 (Algorithm 3 of https://eprint.iacr.org/2022/198)
			if err = eval.BlindRotateCore(a, acc, key); err != nil {
				panic(err)
			}

			// f(X) * X^{b + <a, s>}
			res[index] = acc.CopyNew()

			if !eval.paramsLUT.NTTFlag() {
				ringQLUT.INTT(res[index].Value[0], res[index].Value[0])
				ringQLUT.INTT(res[index].Value[1], res[index].Value[1])
				res[index].IsNTT = false
			}
		}
	}

	return
}

// BlindRotateCore implements Algorithm 3 of https://eprint.iacr.org/2022/198
func (eval *Evaluator) BlindRotateCore(a []uint64, acc *rlwe.Ciphertext, evk BlindRotatationEvaluationKeySet) (err error) {

	// GaloisElement(k) = GaloisGen^{k} mod 2N
	GaloisElement := eval.paramsLUT.GaloisElement

	// Maps a[i] to (+/-) g^{k} mod 2N
	discretLogSets := eval.getDiscretLogSets(a)

	Nhalf := eval.paramsLUT.N() >> 1

	// Algorithm 3 of https://eprint.iacr.org/2022/198
	var v int
	// Lines 3 to 9 (negative set of a[i] = -g^{k} mod 2N)
	for i := Nhalf - 1; i > 0; i-- {
		if v, err = eval.evaluateFromDiscretLogSets(GaloisElement, discretLogSets, -i, v, acc, evk); err != nil {
			return
		}
	}

	// Line 10 (0 in the negative set is 2N)
	if _, err = eval.evaluateFromDiscretLogSets(GaloisElement, discretLogSets, eval.paramsLUT.N()<<1, 0, acc, evk); err != nil {
		return
	}

	// Line 12
	// acc = acc(X^{-g})
	if err = eval.Automorphism(acc, eval.paramsLUT.RingQ().NthRoot()-ring.GaloisGen, acc); err != nil {
		return
	}

	// Lines 13 - 19 (positive set of a[i] = g^{k} mod 2N)
	for i := Nhalf - 1; i > 0; i-- {
		if v, err = eval.evaluateFromDiscretLogSets(GaloisElement, discretLogSets, i, v, acc, evk); err != nil {
			return
		}
	}

	// Lines 20 - 21 (0 in the positive set is 0)
	if _, err = eval.evaluateFromDiscretLogSets(GaloisElement, discretLogSets, 0, 0, acc, evk); err != nil {
		return
	}

	return
}

// evaluateFromDiscretLogSets loops of Algorithm 3 of https://eprint.iacr.org/2022/198
func (eval *Evaluator) evaluateFromDiscretLogSets(GaloisElement func(k int) (galEl uint64), sets map[int][]int, k, v int, acc *rlwe.Ciphertext, evk BlindRotatationEvaluationKeySet) (int, error) {

	// Checks if k is in the discret log sets
	if set, ok := sets[k]; ok {

		// First condition of line 7 or 17
		if v != 0 {

			if err := eval.Automorphism(acc, GaloisElement(v), acc); err != nil {
				return v, err
			}

			v = 0
		}

		for _, j := range set {

			brk, err := evk.GetBlingRotateKey(j)
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
func getGaloisElementInverseMap(GaloisGen uint64, N int) (GaloisGenDiscretLog map[uint64]int) {

	twoN := N << 1
	NHalf := N >> 1
	mask := uint64(twoN - 1)

	GaloisGenDiscretLog = map[uint64]int{}

	var pow uint64 = 1
	for i := 0; i < NHalf; i++ {
		GaloisGenDiscretLog[pow] = i
		GaloisGenDiscretLog[uint64(twoN)-pow] = -i
		pow *= GaloisGen
		pow &= mask
	}

	return
}

// getDiscretLogSets returns map[+/-k] = [i...] for a[0 <= i < N] = {(+/-) g^{k} mod 2N for +/- k}
func (eval *Evaluator) getDiscretLogSets(a []uint64) (discretLogSets map[int][]int) {

	GaloisGenDiscretLog := eval.galoisGenDiscretLog

	// Maps (2*N*a[i]/QLWE) to -N/2 < k <= N/2 for a[i] = (+/- 1) * g^{k}
	discretLogSets = map[int][]int{}
	for i, ai := range a {

		dlog := GaloisGenDiscretLog[ai]

		if _, ok := discretLogSets[dlog]; !ok {
			discretLogSets[dlog] = []int{i}
		} else {
			discretLogSets[dlog] = append(discretLogSets[dlog], i)
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

	twoN := uint64(eval.paramsLUT.N() << 1)
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
