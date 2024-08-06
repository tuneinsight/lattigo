package rlwe

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils"
)

// RingPackingEvaluator is an evaluator for Ring-LWE packing operations.
// All fields of this struct are public, enabling custom instantiations.
type RingPackingEvaluator struct {
	*RingPackingEvaluationKey

	Evaluators map[int]*Evaluator

	//XPow2NTT: [1, x, x^2, x^4, ..., x^2^s] / (X^2^s +1)
	XPow2NTT map[int][]ring.Poly

	//XInvPow2NTT: [1, x^-1, x^-2, x^-4, ..., x^-2^s/2] / (X^2^s +1)
	XInvPow2NTT map[int][]ring.Poly
}

// NewRingPackingEvaluator instantiates a new RingPackingEvaluator from a RingPackingEvaluationKey.
func NewRingPackingEvaluator(evk *RingPackingEvaluationKey) *RingPackingEvaluator {

	Evaluators := map[int]*Evaluator{}
	XPow2NTT := map[int][]ring.Poly{}
	XInvPow2NTT := map[int][]ring.Poly{}

	minLogN := evk.MinLogN()
	maxLogN := evk.MaxLogN()

	levelQ := evk.Parameters[minLogN].GetRLWEParameters().MaxLevel()

	for i := minLogN; i < maxLogN+1; i++ {
		pi := evk.Parameters[i].GetRLWEParameters()
		Evaluators[i] = NewEvaluator(pi, nil)
		XPow2NTT[i] = GenXPow2NTT(pi.RingQ().AtLevel(levelQ), pi.LogN(), false)
		XInvPow2NTT[i] = GenXPow2NTT(pi.RingQ().AtLevel(levelQ), pi.LogN(), true)
	}

	return &RingPackingEvaluator{
		RingPackingEvaluationKey: evk,
		Evaluators:               Evaluators,
		XPow2NTT:                 XPow2NTT,
		XInvPow2NTT:              XInvPow2NTT,
	}
}

// ShallowCopy creates a shallow copy of this struct in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Evaluators can be used concurrently.
func (eval RingPackingEvaluator) ShallowCopy() *RingPackingEvaluator {

	Evaluators := map[int]*Evaluator{}
	for i := range eval.Evaluators {
		Evaluators[i] = eval.Evaluators[i].ShallowCopy()
	}

	return &RingPackingEvaluator{
		RingPackingEvaluationKey: eval.RingPackingEvaluationKey,
		Evaluators:               Evaluators,
		XPow2NTT:                 eval.XPow2NTT,
		XInvPow2NTT:              eval.XInvPow2NTT,
	}
}

// Extract takes as input a ciphertext encrypting P(X) = c[i] * X^i and returns a map of
// ciphertexts of degree eval.MinLogN(), each encrypting P(X) = c[i] * X^{0} for i in idx.
// All non-constant coefficients are zeroed and thus correctness is ensured if this method
// is composed with either Repack or RepackNaive.
func (eval RingPackingEvaluator) Extract(ct *Ciphertext, idx map[int]bool) (cts map[int]*Ciphertext, err error) {
	return eval.extract(ct, idx, false)
}

// ExtractNaive takes as input a ciphertext encrypting P(X) = c[i] * X^i and returns a map of
// ciphertexts of degree eval.MinLogN(), each encrypting P(X) = c[i] * X^{0} for i in idx.
// Non-constant coefficients are NOT zeroed thus correctness is only ensured if this method
// is composed with Repack.
//
// If eval.MinLogN() = eval.MaxLogN(), no evaluation keys are required for this method.
// If eval.MinLogN() < eval.MaxLogN(), only RingSwitchingKeys are required for this method.
func (eval RingPackingEvaluator) ExtractNaive(ct *Ciphertext, idx map[int]bool) (cts map[int]*Ciphertext, err error) {
	return eval.extract(ct, idx, true)
}

// If naive = false, then all non-constant coefficients are zeroed.
func (eval RingPackingEvaluator) extract(ct *Ciphertext, idx map[int]bool, naive bool) (cts map[int]*Ciphertext, err error) {

	logNMax := ct.LogN()
	logNMin := eval.MinLogN()
	level := ct.Level()

	logNFactor := logNMax - logNMin
	NFactor := 1 << logNFactor

	keys := utils.GetSortedKeys(idx)

	_, logGap, err := getMinimumGap(keys)

	if err != nil {
		return nil, fmt.Errorf("getMinimumGap: %w", err)
	}

	// First recursively splits the ciphertexts into smaller ciphertexts of half the ring
	// degree until the minimum ring degre is reached
	tmpCts := make(map[int]*Ciphertext)
	tmpCts[0] = ct.CopyNew()
	for i := 0; i < logNFactor; i++ {
		t := 1 << i

		// Each split of the ring divides the gap a factor of two
		logGap = utils.Max(0, logGap-1)

		for j := 0; j < t; j++ {

			if tmpCts[j] != nil {

				ctEvenNHalf := NewCiphertext(eval.Parameters[logNMax-i-1], 1, level)
				ctOddNHalf := NewCiphertext(eval.Parameters[logNMax-i-1], 1, level)

				if err = eval.Split(tmpCts[j], ctEvenNHalf, ctOddNHalf); err != nil {
					return nil, fmt.Errorf("eval.split(cts[%d]): %w", j, err)
				}

				tmpCts[j] = ctEvenNHalf
				tmpCts[j+t] = ctOddNHalf
			}
		}
	}

	gap := 1 << logGap

	// Applies the same split on the index map, but also update the
	// indexes to take into account the new ordering
	buckets := make(map[int][]int)
	for _, i := range keys {
		bucket := i & (NFactor - 1)
		buckets[bucket] = append(buckets[bucket], i/NFactor)
	}

	// For each small ciphertext, extracts the relevant values
	cts = make(map[int]*Ciphertext)
	for i := range buckets {
		var ciphertexts map[int]*Ciphertext

		if naive {

			ciphertexts = map[int]*Ciphertext{}
			for _, j := range buckets[i] {
				ciphertexts[j] = tmpCts[i].CopyNew()
			}

			XInvPow2NTT := eval.XInvPow2NTT[logNMin]

			ringQ := eval.Parameters[logNMin].GetRLWEParameters().RingQ().AtLevel(level)

			// Rotates ciphertexts to move c[i] * X^{i} -> c[i] * X^{0}
			// by sequentially multplying with the appropriate X^{-2^{i}}.
			for i := 0; i < logNMin; i++ {
				for j := range ciphertexts {
					if (j>>i)&1 == 1 {
						ct := ciphertexts[j]
						ringQ.MulCoeffsMontgomery(ct.Value[0], XInvPow2NTT[i], ct.Value[0])
						ringQ.MulCoeffsMontgomery(ct.Value[1], XInvPow2NTT[i], ct.Value[1])
					}
				}
			}

		} else {
			if ciphertexts, err = eval.Expand(tmpCts[i], logGap); err != nil {
				return nil, fmt.Errorf("evalN.expand(tmpCt[%d], %d): %w", i, logGap, err)
			}
		}

		for _, j := range buckets[i] {
			if ct, ok := ciphertexts[j]; ok {
				cts[i+j*NFactor] = ct
			} else {
				return nil, fmt.Errorf("invalid ciphertexts map: index i+j*(NFactor*gap)=%d is nil", i+j*(NFactor*gap))
			}

		}
	}

	return
}

// Split splits a ciphertext of degree N into two ciphertexts of degree N/2:
// ctN[X] = ctEvenNHalf[Y] + X * ctOddNHalf[Y] where Y = X^2.
func (eval RingPackingEvaluator) Split(ctN, ctEvenNHalf, ctOddNHalf *Ciphertext) (err error) {

	if eval.MinLogN() == eval.MaxLogN() {
		return fmt.Errorf("method is not supported when eval.MinLogN() == eval.MaxLogN()")
	}

	if ctN.LogN() <= eval.MinLogN() {
		return fmt.Errorf("ctN.Log() must be greater than eval.MinLogN()")
	}

	if ctEvenNHalf == nil {
		return fmt.Errorf("ctEvenNHalf cannot be nil")
	}

	if ctEvenNHalf.LogN() != ctN.LogN()-1 {
		return fmt.Errorf("ctEvenNHalf.LogN() must be equal to ctN.LogN()-1")
	}

	LogN := ctN.LogN()

	evalN := eval.Evaluators[LogN]
	evkNToNHalf := eval.RingSwitchingKeys[LogN][LogN-1]

	ctTmp := NewCiphertext(eval.Parameters[LogN], 1, ctN.Level())

	// SkN -> SkNHalf
	if err = evalN.ApplyEvaluationKey(ctN, evkNToNHalf, ctTmp); err != nil {
		return fmt.Errorf("ApplyEvaluationKey: %w", err)
	}

	r := eval.Parameters[LogN].GetRLWEParameters().RingQ().AtLevel(ctN.Level())

	// Maps to smaller ring degree X -> Y = X^{2}

	*ctEvenNHalf.MetaData = *ctN.MetaData
	SwitchCiphertextRingDegreeNTT(ctTmp.El(), r, ctEvenNHalf.El())
	ctEvenNHalf.LogDimensions.Cols--

	// Maps to smaller ring degree X -> Y = X^{2}
	if ctOddNHalf != nil {

		if ctOddNHalf.LogN() != ctN.LogN()-1 {
			return fmt.Errorf("ctOddNHalf.LogN() must be equal to ctN.LogN()-1")
		}

		*ctOddNHalf.MetaData = *ctN.MetaData
		r.MulCoeffsMontgomery(ctTmp.Value[0], eval.XInvPow2NTT[LogN][0], ctTmp.Value[0])
		r.MulCoeffsMontgomery(ctTmp.Value[1], eval.XInvPow2NTT[LogN][0], ctTmp.Value[1])
		SwitchCiphertextRingDegreeNTT(ctTmp.El(), r, ctOddNHalf.El())
		ctOddNHalf.LogDimensions.Cols--
	}

	return
}

// SplitNew splits a ciphertext of degree N into two ciphertexts of degree N/2:
// ctN[X] = ctEvenNHalf[Y] + X * ctOddNHalf[Y] where Y = X^2.
func (eval RingPackingEvaluator) SplitNew(ctN *Ciphertext) (ctEvenNHalf, ctOddNHalf *Ciphertext, err error) {
	if eval.MinLogN() == eval.MaxLogN() {
		return nil, nil, fmt.Errorf("method is not supported when eval.MinLogN() == eval.MaxLogN()")
	}
	LogN := ctN.LogN()
	ctEvenNHalf = NewCiphertext(eval.Parameters[LogN-1], 1, ctN.Level())
	ctOddNHalf = NewCiphertext(eval.Parameters[LogN-1], 1, ctN.Level())
	return ctEvenNHalf, ctOddNHalf, eval.Split(ctN, ctEvenNHalf, ctOddNHalf)
}

// Repack takes as input a map of ciphertext and repacks the constant coefficient each ciphertext
// into a single ciphertext of degree eval.MaxLogN() following the indexing of the map.
//
// For example, if cts = map[int]*Ciphertext{0:ct0, 1:ct1, 4:ct2}, then the method will return
// a ciphertext encrypting P(X) = ct0[0] + ct1[0] * X + ct2[0] * X^4.
//
// The method accepts ciphertexts of a ring degree between eval.MinLogN() and eval.MaxLogN().
//
// All non-constant coefficient are zeroed during the repacking, thus correctness is ensured if this
// method can be composed with either Extract or ExtractNaive.
func (eval RingPackingEvaluator) Repack(cts map[int]*Ciphertext) (ct *Ciphertext, err error) {
	return eval.repack(cts, false)
}

// RepackNaive takes as input a map of ciphertext and repacks the constant coefficient each ciphertext
// into a single ciphertext of degree eval.MaxLogN() following the indexing of the map.
//
// For example, if cts = map[int]*Ciphertext{0:ct0, 1:ct1, 4:ct2}, then the method will return
// a ciphertext encrypting P(X) = ct0[0] + ct1[0] * X + ct2[0] * X^4.
//
// The method accepts ciphertexts of a ring degree between eval.MinLogN() and eval.MaxLogN().
//
// If eval.MinLogN() = eval.MaxLogN(), no evaluation keys are required for this method.
// If eval.MinLogN() < eval.MaxLogN(), only RingSwitchingKeys are required for this method.
//
// Unlike Repack, non-constant coefficient are NOT zeroed during the repacking, thus correctness is only
// ensured if this method is composed with either Extract.
func (eval RingPackingEvaluator) RepackNaive(cts map[int]*Ciphertext) (ct *Ciphertext, err error) {
	return eval.repack(cts, true)
}

func (eval RingPackingEvaluator) repack(cts map[int]*Ciphertext, naive bool) (ct *Ciphertext, err error) {

	keys := utils.GetSortedKeys(cts)

	logNMin := cts[keys[0]].LogN()
	logNMax := eval.MaxLogN()
	level := cts[keys[0]].Level()

	logNFactor := logNMax - logNMin
	NFactor := 1 << logNFactor

	// List of map containing the repacking of cts
	ctsSmallN := make([]map[int]*Ciphertext, NFactor)
	for i := range ctsSmallN {
		ctsSmallN[i] = map[int]*Ciphertext{}
	}

	// Assigns to each map the corresponding ciphertext.
	// This takes into account the future merging, that merges
	// ciphertexts in a base-2 tree-like fashion by evaluating
	// ctN[X] = ctEvenNHalf[Y] + X * ctOddNHalf[Y] where Y = X^2.
	for _, i := range keys {
		ctsSmallN[i&(NFactor-1)][i/NFactor] = cts[i]
	}

	// Map of repacked ciphertext that will then be merged together.
	// Each merging takes two ciphertexts, doubles their ring degree
	// and adds them together.
	ctsLargeN := map[int]*Ciphertext{}
	for i := 0; i < NFactor; i++ {

		if naive {

			tmpCts := ctsSmallN[i]

			XPow2NTT := eval.XPow2NTT[logNMin]

			ringQ := eval.Parameters[logNMin].GetRLWEParameters().RingQ().AtLevel(level)

			for i := 0; i < logNMin; i++ {

				t := 1 << (logNMin - 1 - i)

				for jx, jy := 0, t; jx < t; jx, jy = jx+1, jy+1 {

					a := tmpCts[jx]
					b := tmpCts[jy]

					if b != nil {

						//X^(N/2^L)
						ringQ.MulCoeffsMontgomery(b.Value[0], XPow2NTT[len(XPow2NTT)-i-1], b.Value[0])
						ringQ.MulCoeffsMontgomery(b.Value[1], XPow2NTT[len(XPow2NTT)-i-1], b.Value[1])

						if a != nil {
							// a = a + b * X^{N/2^{i}}
							ringQ.Add(a.Value[0], b.Value[0], a.Value[0])
							ringQ.Add(a.Value[1], b.Value[1], a.Value[1])

						} else {
							// if ct[jx] == nil, then simply re-assigns
							tmpCts[jx] = tmpCts[jy]
						}

						tmpCts[jy] = nil
					}
				}
			}

			ctsLargeN[i] = tmpCts[0]

		} else {
			if len(ctsSmallN[i]) != 0 {
				if ctsLargeN[i], err = eval.Pack(ctsSmallN[i], logNMin, true); err != nil {
					return nil, fmt.Errorf("eval.pack(ctsSmallN[%d], logGap=%d, true): %w", i, logNMin, err)
				}
			}
		}
	}

	// Merges the cipehrtexts in a base-2 tree like fashion.
	for i := logNFactor - 1; i >= 0; i-- {
		t := 1 << i

		for j := 0; j < t; j++ {

			if ctsLargeN[j] != nil || ctsLargeN[j+1] != nil {

				ctN := NewCiphertext(eval.Parameters[logNMax-i], 1, level)

				if err = eval.Merge(ctsLargeN[j], ctsLargeN[j+t], ctN); err != nil {
					return nil, fmt.Errorf("eval.split(cts[%d]): %w", j, err)
				}

				ctsLargeN[j] = ctN
				ctsLargeN[j+t] = nil
			}
		}
	}

	return ctsLargeN[0], nil
}

// Merge merges two ciphertexts of degree N/2 into a ciphertext of degre N:
// ctN[X] = ctEvenNHalf[Y] + X * ctOddNHalf[Y] where Y = X^2.
func (eval RingPackingEvaluator) Merge(ctEvenNHalf, ctOddNHalf, ctN *Ciphertext) (err error) {

	if eval.MinLogN() == eval.MaxLogN() {
		return fmt.Errorf("method is not supported when eval.MinLogN() == eval.MaxLogN()")
	}

	if ctEvenNHalf == nil {
		return fmt.Errorf("ctEvenNHalf cannot be nil")
	}

	if ctEvenNHalf.LogN() >= eval.MaxLogN() {
		return fmt.Errorf("ctEvenNHalf.LogN() must be smaller than eval.MaxLogN()")
	}

	if ctN.LogN() != ctEvenNHalf.LogN()+1 {
		return fmt.Errorf("ctN.LogN() must be equal to ctEvenNHalf.LogN()+1")
	}

	if ctOddNHalf != nil {
		if ctEvenNHalf.LogN() != ctOddNHalf.LogN() {
			return fmt.Errorf("ctEvenNHalf.LogN() and ctOddNHalf.LogN() must be equal")
		}
	}

	LogN := ctN.LogN()

	evalN := eval.Evaluators[LogN]
	evkNHalfToN := eval.RingSwitchingKeys[LogN-1][LogN]
	r := eval.Parameters[LogN].GetRLWEParameters().RingQ().AtLevel(ctN.Level())

	ctTmp := NewCiphertext(eval.Parameters[LogN], 1, ctN.Level())

	*ctN.MetaData = *ctEvenNHalf.MetaData
	SwitchCiphertextRingDegreeNTT(ctEvenNHalf.El(), r, ctN.El())

	if ctOddNHalf != nil {
		SwitchCiphertextRingDegreeNTT(ctOddNHalf.El(), r, ctTmp.El())
		r.MulCoeffsMontgomeryThenAdd(ctTmp.Value[0], eval.XPow2NTT[LogN][0], ctN.Value[0])
		r.MulCoeffsMontgomeryThenAdd(ctTmp.Value[1], eval.XPow2NTT[LogN][0], ctN.Value[1])
	}

	// SkNHalf -> SkN
	if err = evalN.ApplyEvaluationKey(ctN, evkNHalfToN, ctN); err != nil {
		return fmt.Errorf("evalN.ApplyEvaluationKey(ctN, evkNToNHalf, ctN): %w", err)
	}

	ctN.LogDimensions.Cols++
	return
}

// MergeNew merges two ciphertexts of degree N/2 into a ciphertext of degre N:
// ctN[X] = ctEvenNHalf[Y] + X * ctOddNHalf[Y] where Y = X^2.
func (eval RingPackingEvaluator) MergeNew(ctEvenNHalf, ctOddNHalf *Ciphertext) (ctN *Ciphertext, err error) {

	if eval.MinLogN() == eval.MaxLogN() {
		return nil, fmt.Errorf("method is not supported when eval.MinLogN() == eval.MaxLogN()")
	}

	if ctEvenNHalf == nil {
		return nil, fmt.Errorf("ctEvenNHalf cannot be nil")
	}

	if ctEvenNHalf.LogN() >= eval.MaxLogN() {
		return nil, fmt.Errorf("ctEvenNHalf.LogN() must be smaller than eval.MaxLogN()")
	}

	ctN = NewCiphertext(eval.Parameters[ctEvenNHalf.LogN()+1], 1, ctEvenNHalf.Level())
	return ctN, eval.Merge(ctEvenNHalf, ctOddNHalf, ctN)
}

// Expand expands a RLWE Ciphertext encrypting P(X) = ci * X^i and returns a map of
// ciphertexts, each encrypting P(X) = ci * X^0, indexed by i, for 0<= i < 2^{logN}
// and i divisible by 2^{logGap}.
//
// This method is a used as a sub-routine of the Extract method.
//
// The method will return an error if:
//   - The input ciphertext degree is not one
//   - The ring type is not ring.Standard
func (eval RingPackingEvaluator) Expand(ct *Ciphertext, logGap int) (cts map[int]*Ciphertext, err error) {

	if ct.Degree() != 1 {
		return nil, fmt.Errorf("ct.Degree() != 1")
	}

	logN := ct.LogN()

	var params Parameters

	if p, ok := eval.Parameters[logN]; !ok {
		return nil, fmt.Errorf("eval.Parameters[%d] is nil", logN)
	} else {
		params = *p.GetRLWEParameters()
	}

	if eval.ExtractKeys == nil {
		return nil, fmt.Errorf("eval.ExtractKeys is nil")
	}

	var evk EvaluationKeySet
	if p, ok := eval.ExtractKeys[params.LogN()]; !ok {
		return nil, fmt.Errorf("eval.ExtractKeys[%d] is nil", params.LogN())
	} else {
		evk = p
	}

	evalN := eval.Evaluators[params.LogN()].WithKey(evk)
	xPow2 := eval.XInvPow2NTT[params.LogN()]

	level := ct.Level()

	ringQ := params.RingQ().AtLevel(level)

	if params.RingType() != ring.Standard {
		return nil, fmt.Errorf("method is only supported for ring.Type = ring.Standard (X^{-2^{i}} does not exist in the sub-ring Z[X + X^{-1}])")
	}

	cts = map[int]*Ciphertext{}
	cts[0] = ct.CopyNew()
	cts[0].LogDimensions = ring.Dimensions{Rows: 0, Cols: 0}

	if ct := cts[0]; !ct.IsNTT {
		ringQ.NTT(ct.Value[0], ct.Value[0])
		ringQ.NTT(ct.Value[1], ct.Value[1])
		ct.IsNTT = true
	}

	// Multiplies by 2^{-logN} mod Q
	NInv := new(big.Int).SetUint64(1 << logN)
	NInv.ModInverse(NInv, ringQ.ModulusAtLevel[level])

	ringQ.MulScalarBigint(cts[0].Value[0], NInv, cts[0].Value[0])
	ringQ.MulScalarBigint(cts[0].Value[1], NInv, cts[0].Value[1])

	gap := 1 << logGap

	tmp, err := NewCiphertextAtLevelFromPoly(level, []ring.Poly{evalN.BuffCt.Value[0], evalN.BuffCt.Value[1]})

	// Sanity check, this error should not happen unless the
	// evaluator's buffer thave been improperly tempered with.
	if err != nil {
		panic(err)
	}

	*tmp.MetaData = *ct.MetaData

	for i := 0; i < logN; i++ {

		n := 1 << i

		galEl := uint64(ringQ.N()/n + 1)

		for j := 0; j < n; j += gap {

			c0 := cts[j]

			// X -> X^{N/n + 1}
			//[a, b, c, d] -> [a, -b, c, -d]
			if err = evalN.Automorphism(c0, galEl, tmp); err != nil {
				return nil, fmt.Errorf("evalN.Automorphism(c0, galEl, tmp): %w", err)
			}

			if j+n/gap > 0 {

				c1 := cts[j].CopyNew()

				// Zeroes odd coeffs: [a, b, c, d] + [a, -b, c, -d] -> [2a, 0, 2b, 0]
				ringQ.Add(c0.Value[0], tmp.Value[0], c0.Value[0])
				ringQ.Add(c0.Value[1], tmp.Value[1], c0.Value[1])

				// Zeroes even coeffs: [a, b, c, d] - [a, -b, c, -d] -> [0, 2b, 0, 2d]
				ringQ.Sub(c1.Value[0], tmp.Value[0], c1.Value[0])
				ringQ.Sub(c1.Value[1], tmp.Value[1], c1.Value[1])

				// c1 * X^{-2^{i}}: [0, 2b, 0, 2d] * X^{-n} -> [2b, 0, 2d, 0]
				ringQ.MulCoeffsMontgomery(c1.Value[0], xPow2[i], c1.Value[0])
				ringQ.MulCoeffsMontgomery(c1.Value[1], xPow2[i], c1.Value[1])

				cts[j+n] = c1

			} else {

				// Zeroes odd coeffs: [a, b, c, d] + [a, -b, c, -d] -> [2a, 0, 2b, 0]
				ringQ.Add(c0.Value[0], tmp.Value[0], c0.Value[0])
				ringQ.Add(c0.Value[1], tmp.Value[1], c0.Value[1])
			}
		}
	}

	for _, ct := range cts {
		if ct != nil && !ct.IsNTT {
			ringQ.INTT(ct.Value[0], ct.Value[0])
			ringQ.INTT(ct.Value[1], ct.Value[1])
			ct.IsNTT = false
		}
	}
	return
}

// Pack packs a map of of ciphertexts, each encrypting Pi(X) = ci * X^{i} for 0 <= i * 2^{inputLogGap} < 2^{LogN}
// and indexed by j, for 0<= j < 2^{eval.MaxLogN()} and returns ciphertext encrypting P(X) = Pi(X) * X^i.
// zeroGarbageSlots: if set to true, slots which are not multiples of X^{2^{logGap}} will be zeroed during the procedure.
//
// The method will return an error if:
//   - The number of ciphertexts is 0
//   - Any input ciphertext degree is not one
//   - Gaps between ciphertexts is smaller than inputLogGap > N
//   - The ring type is not ring.Standard
//
// Example: we want to pack 4 ciphertexts into one, and keep only coefficients which are a multiple of X^{4}.
//
//	To do so, we must set logGap = 2.
//	Here the `X` slots are treated as garbage slots that we want to discard during the procedure.
//
//	input: map[int]{
//	   0: [x00, X, X, X, x01, X, X, X],   with logGap = 2
//	   1: [x10, X, X, X, x11, X, X, X],
//	   2: [x20, X, X, X, x21, X, X, X],
//	   3: [x30, X, X, X, x31, X, X, X],
//		}
//
//	 Step 1:
//	         map[0]: 2^{-1} * (map[0] + X^2 * map[2] + phi_{5^2}(map[0] - X^2 * map[2]) = [x00, X, x20, X, x01, X, x21, X]
//	         map[1]: 2^{-1} * (map[1] + X^2 * map[3] + phi_{5^2}(map[1] - X^2 * map[3]) = [x10, X, x30, X, x11, X, x31, X]
//	 Step 2:
//	         map[0]: 2^{-1} * (map[0] + X^1 * map[1] + phi_{5^4}(map[0] - X^1 * map[1]) = [x00, x10, x20, x30, x01, x11, x21, x22]
func (eval RingPackingEvaluator) Pack(cts map[int]*Ciphertext, inputLogGap int, zeroGarbageSlots bool) (ct *Ciphertext, err error) {

	if len(cts) == 0 {
		return nil, fmt.Errorf("len(cts) = 0")
	}

	keys := utils.GetSortedKeys(cts)

	logN := cts[keys[0]].LogN()

	var params Parameters

	if p, ok := eval.Parameters[logN]; !ok {
		return nil, fmt.Errorf("eval.Parameters[%d] is nil", logN)
	} else {
		params = *p.GetRLWEParameters()
	}

	if eval.RepackKeys == nil {
		return nil, fmt.Errorf("eval.RepackKeys is nil")
	}

	var evk EvaluationKeySet
	if p, ok := eval.RepackKeys[params.LogN()]; !ok {
		return nil, fmt.Errorf("eval.RepackKeys[%d] is nil", params.LogN())
	} else {
		evk = p
	}

	evalN := eval.Evaluators[params.LogN()].WithKey(evk)

	xPow2 := eval.XPow2NTT[params.LogN()]

	if params.RingType() != ring.Standard {
		return nil, fmt.Errorf("procedure is only supported for ring.Type = ring.Standard (X^{2^{i}} does not exist in the sub-ring Z[X + X^{-1}])")
	}

	level := cts[keys[0]].Level()

	var gap, logGap int

	if len(keys) > 1 {
		if gap, logGap, err = getMinimumGap(keys); err != nil {
			return nil, fmt.Errorf("getMinimumGap: %w", err)
		}
	} else {
		gap = params.N()
		logGap = params.LogN()
	}

	ringQ := params.RingQ().AtLevel(level)

	logStart := logN - inputLogGap
	logEnd := logN

	if !zeroGarbageSlots {
		if gap > 0 {
			logEnd -= logGap
		}
	}

	if logStart >= logEnd {
		return nil, fmt.Errorf("gaps between ciphertexts is smaller than inputLogGap > N")
	}

	NInv := new(big.Int).SetUint64(uint64(1 << (logEnd - logStart)))
	NInv.ModInverse(NInv, ringQ.ModulusAtLevel[level])

	for _, key := range keys {

		ct := cts[key]

		if ct.Degree() != 1 {
			return nil, fmt.Errorf("cts[%d].Degree() != 1", key)
		}

		if !ct.IsNTT {
			ringQ.NTT(ct.Value[0], ct.Value[0])
			ringQ.NTT(ct.Value[1], ct.Value[1])
			ct.IsNTT = true
		}

		ringQ.MulScalarBigint(ct.Value[0], NInv, ct.Value[0])
		ringQ.MulScalarBigint(ct.Value[1], NInv, ct.Value[1])
	}

	tmpa := &Ciphertext{}
	tmpa.Value = []ring.Poly{ringQ.NewPoly(), ringQ.NewPoly()}
	tmpa.MetaData = &MetaData{}
	tmpa.MetaData.IsNTT = true

	for i := logStart; i < logEnd; i++ {

		t := 1 << (logN - 1 - i)

		for jx, jy := 0, t; jx < t; jx, jy = jx+1, jy+1 {

			a := cts[jx]
			b := cts[jy]

			if b != nil {

				//X^(N/2^L)
				ringQ.MulCoeffsMontgomery(b.Value[0], xPow2[len(xPow2)-i-1], b.Value[0])
				ringQ.MulCoeffsMontgomery(b.Value[1], xPow2[len(xPow2)-i-1], b.Value[1])

				if a != nil {

					// tmpa = phi(a - b * X^{N/2^{i}}, 2^{i-1})
					ringQ.Sub(a.Value[0], b.Value[0], tmpa.Value[0])
					ringQ.Sub(a.Value[1], b.Value[1], tmpa.Value[1])

					// a = a + b * X^{N/2^{i}}
					ringQ.Add(a.Value[0], b.Value[0], a.Value[0])
					ringQ.Add(a.Value[1], b.Value[1], a.Value[1])

				} else {
					// if ct[jx] == nil, then simply re-assigns
					cts[jx] = cts[jy]
				}

				cts[jy] = nil
			}

			if a != nil {

				var galEl uint64

				if i == 0 {
					galEl = ringQ.NthRoot() - 1
				} else {
					galEl = params.GaloisElement(1 << (i - 1))
				}

				if b != nil {
					if err = evalN.Automorphism(tmpa, galEl, tmpa); err != nil {
						return nil, fmt.Errorf("evalN.Automorphism(tmpa, galEl, tmpa): %w", err)
					}
				} else {
					if err = evalN.Automorphism(a, galEl, tmpa); err != nil {
						return nil, fmt.Errorf("evalN.Automorphism(a, galEl, tmpa): %w", err)
					}
				}

				// a + b * X^{N/2^{i}} + phi(a - b * X^{N/2^{i}}, 2^{i-1})
				ringQ.Add(a.Value[0], tmpa.Value[0], a.Value[0])
				ringQ.Add(a.Value[1], tmpa.Value[1], a.Value[1])

			} else if b != nil {

				var galEl uint64

				if i == 0 {
					galEl = ringQ.NthRoot() - 1
				} else {
					galEl = params.GaloisElement(1 << (i - 1))
				}

				if err = evalN.Automorphism(b, galEl, tmpa); err != nil {
					return nil, fmt.Errorf("evalN.Automorphism(b, galEl, tmpa): %w", err)
				}

				// b * X^{N/2^{i}} - phi(b * X^{N/2^{i}}, 2^{i-1}))
				ringQ.Sub(b.Value[0], tmpa.Value[0], b.Value[0])
				ringQ.Sub(b.Value[1], tmpa.Value[1], b.Value[1])
			}
		}
	}

	return cts[0], nil
}

// GenXPow2NTT generates X^({-1 if div else 1} * {2^{0 <= i < LogN}}) in NTT.
func GenXPow2NTT(r *ring.Ring, logN int, div bool) (xPow []ring.Poly) {

	// Compute X^{-n} from 0 to LogN
	xPow = make([]ring.Poly, logN)

	moduli := r.ModuliChain()[:r.Level()+1]
	BRC := r.BRedConstants()

	var idx int
	for i := 0; i < logN; i++ {

		idx = 1 << i

		if div {
			idx = r.N() - idx
		}

		xPow[i] = r.NewPoly()

		if i == 0 {

			for j := range moduli {
				xPow[i].Coeffs[j][idx] = ring.MForm(1, moduli[j], BRC[j])
			}

			r.NTT(xPow[i], xPow[i])

		} else {
			r.MulCoeffsMontgomery(xPow[i-1], xPow[i-1], xPow[i]) // X^{n} = X^{1} * X^{n-1}
		}
	}

	if div {
		r.Neg(xPow[0], xPow[0])
	}

	return
}

func getMinimumGap(list []int) (gap, logGap int, err error) {

	// The loops over to find the smallest gap
	gap = 0x7fffffffffffffff // 2^{63}-1
	for i := 1; i < len(list); i++ {

		a, b := list[i-1], list[i]

		if a > b {
			return gap, logGap, fmt.Errorf("invalid index list: element must be sorted from smallest to largest")
		} else if a == b {
			return gap, logGap, fmt.Errorf("invalid index list: contains duplicated elements")
		}

		if tmp := b - a; tmp < gap {
			gap = tmp
		}

		if gap == 1 {
			break
		}
	}

	// Sets gap to the largest power-of-two that divides it.
	// We will then discart all coefficients that are not a
	// multiple of this gap (and thus possibly entire ciph-
	// ertexts).
	for gap&1 == 0 {
		logGap++
		gap >>= 1
	}

	return
}
