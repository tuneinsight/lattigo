package rlwe

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
	"math/bits"
)

// Evaluator is a struct that holds the necessary elements to execute general homomorphic
// operation on RLWE ciphertexts, such as automorphisms, key-switching and relinearization.
type Evaluator struct {
	*evaluatorBase
	*evaluatorBuffers

	Rlk             *RelinearizationKey
	Rtks            *RotationKeySet
	PermuteNTTIndex map[uint64][]uint64

	BasisExtender *ring.BasisExtender
	Decomposer    *ring.Decomposer
}

type evaluatorBase struct {
	params Parameters
}

type evaluatorBuffers struct {
	// BuffQ[0]/BuffP[0] : on the fly decomp(c2)
	// BuffQ[1-5]/BuffP[1-5] : available
	BuffQP        [6]ringqp.Poly
	BuffInvNTT    *ring.Poly
	BuffDecompQP  []ringqp.Poly // Memory Buff for the basis extension in hoisting
	BuffBitDecomp []uint64
}

func newEvaluatorBase(params Parameters) *evaluatorBase {
	ev := new(evaluatorBase)
	ev.params = params
	return ev
}

func newEvaluatorBuffers(params Parameters) *evaluatorBuffers {

	buff := new(evaluatorBuffers)
	decompRNS := params.DecompRNS(params.QCount()-1, params.PCount()-1)
	ringQP := params.RingQP()

	buff.BuffQP = [6]ringqp.Poly{ringQP.NewPoly(), ringQP.NewPoly(), ringQP.NewPoly(), ringQP.NewPoly(), ringQP.NewPoly(), ringQP.NewPoly()}

	buff.BuffInvNTT = params.RingQ().NewPoly()

	buff.BuffDecompQP = make([]ringqp.Poly, decompRNS)
	for i := 0; i < decompRNS; i++ {
		buff.BuffDecompQP[i] = ringQP.NewPoly()
	}

	buff.BuffBitDecomp = make([]uint64, params.RingQ().N)

	return buff
}

// NewEvaluator creates a new Evaluator.
func NewEvaluator(params Parameters, evaluationKey *EvaluationKey) (eval *Evaluator) {
	eval = new(Evaluator)
	eval.evaluatorBase = newEvaluatorBase(params)
	eval.evaluatorBuffers = newEvaluatorBuffers(params)

	if params.RingP() != nil {
		eval.BasisExtender = ring.NewBasisExtender(params.RingQ(), params.RingP())
		eval.Decomposer = ring.NewDecomposer(params.RingQ(), params.RingP())
	}

	if evaluationKey != nil {
		if evaluationKey.Rlk != nil {
			eval.Rlk = evaluationKey.Rlk
		}

		if evaluationKey.Rtks != nil {
			eval.Rtks = evaluationKey.Rtks
			eval.PermuteNTTIndex = *eval.permuteNTTIndexesForKey(eval.Rtks)
		}
	}

	return
}

// permuteNTTIndexesForKey generates pemutation indexes for automorphisms for ciphertext
// that are given in the NTT domain.
func (eval *Evaluator) permuteNTTIndexesForKey(rtks *RotationKeySet) *map[uint64][]uint64 {
	if rtks == nil {
		return &map[uint64][]uint64{}
	}
	permuteNTTIndex := make(map[uint64][]uint64, len(rtks.Keys))
	for galEl := range rtks.Keys {
		permuteNTTIndex[galEl] = eval.params.RingQ().PermuteNTTIndex(galEl)
	}
	return &permuteNTTIndex
}

// ShallowCopy creates a shallow copy of this Evaluator in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Evaluators can be used concurrently.
func (eval *Evaluator) ShallowCopy() *Evaluator {
	return &Evaluator{
		evaluatorBase:    eval.evaluatorBase,
		Decomposer:       eval.Decomposer,
		BasisExtender:    eval.BasisExtender.ShallowCopy(),
		evaluatorBuffers: newEvaluatorBuffers(eval.params),
		Rlk:              eval.Rlk,
		Rtks:             eval.Rtks,
		PermuteNTTIndex:  eval.PermuteNTTIndex,
	}
}

// WithKey creates a shallow copy of the receiver Evaluator for which the new EvaluationKey is evaluationKey
// and where the temporary buffers are shared. The receiver and the returned Evaluators cannot be used concurrently.
func (eval *Evaluator) WithKey(evaluationKey *EvaluationKey) *Evaluator {
	var indexes map[uint64][]uint64
	if evaluationKey.Rtks == eval.Rtks {
		indexes = eval.PermuteNTTIndex
	} else {
		indexes = *eval.permuteNTTIndexesForKey(evaluationKey.Rtks)
	}
	return &Evaluator{
		evaluatorBase:    eval.evaluatorBase,
		evaluatorBuffers: eval.evaluatorBuffers,
		Decomposer:       eval.Decomposer,
		BasisExtender:    eval.BasisExtender,
		Rlk:              evaluationKey.Rlk,
		Rtks:             evaluationKey.Rtks,
		PermuteNTTIndex:  indexes,
	}
}

// MergeRLWE merges a batch of RLWE, packing the first coefficient of each RLWE into a single RLWE.
// The operation will require N/gap + log(gap) key-switches, where gap is the minimum gap between
// two non-zero coefficients of the final ciphertext.
// The method takes as input a map of Ciphertexts, indexing in which coefficient, of the final
// ciphertext, the first coefficient of each ciphertext of the map must be packed.
func (eval *Evaluator) MergeRLWE(ciphertexts map[int]*Ciphertext) (ciphertext *Ciphertext) {

	params := eval.params
	ringQ := params.RingQ()

	// Compute X^{n} from 0 to LogN LUT
	xPow := make([]*ring.Poly, params.LogN())
	for i := 0; i < params.LogN(); i++ {
		xPow[i] = ringQ.NewPoly()
		if i == 0 {
			for j := 0; j < params.MaxLevel()+1; j++ {
				xPow[i].Coeffs[j][1<<i] = ring.MForm(1, ringQ.Modulus[j], ringQ.BredParams[j])
			}
			ringQ.NTT(xPow[i], xPow[i])
		} else {
			ringQ.MulCoeffsMontgomery(xPow[i-1], xPow[i-1], xPow[i]) // X^{n} = X^{1} * X^{n-1}
		}
	}

	var level int
	for i := range ciphertexts {
		level = ciphertexts[i].Level()
		break
	}

	NInv := ringQ.NttNInv
	Q := ringQ.Modulus
	mredParams := ringQ.MredParams

	// Multiplies by (Slots * N) ^-1 mod Q
	for i := range ciphertexts {
		if ciphertexts[i] != nil {
			v0, v1 := ciphertexts[i].Value[0], ciphertexts[i].Value[1]
			for j := 0; j < level+1; j++ {
				ring.MulScalarMontgomeryVec(v0.Coeffs[j], v0.Coeffs[j], NInv[j], Q[j], mredParams[j])
				ring.MulScalarMontgomeryVec(v1.Coeffs[j], v1.Coeffs[j], NInv[j], Q[j], mredParams[j])
			}
		}
	}

	ciphertextslist := make([]*Ciphertext, ringQ.N)

	for i := range ciphertexts {
		ciphertextslist[i] = ciphertexts[i]
	}

	if ciphertextslist[0] == nil {
		ciphertextslist[0] = NewCiphertextNTT(params, 1, level)
	}

	return eval.mergeRLWERecurse(ciphertextslist, xPow)
}

func (eval *Evaluator) mergeRLWERecurse(ciphertexts []*Ciphertext, xPow []*ring.Poly) *Ciphertext {

	ringQ := eval.params.RingQ()

	L := bits.Len64(uint64(len(ciphertexts))) - 1

	if L == 0 {
		return ciphertexts[0]
	}

	odd := make([]*Ciphertext, len(ciphertexts)>>1)
	even := make([]*Ciphertext, len(ciphertexts)>>1)

	for i := 0; i < len(ciphertexts)>>1; i++ {
		odd[i] = ciphertexts[2*i]
		even[i] = ciphertexts[2*i+1]
	}

	ctEven := eval.mergeRLWERecurse(odd, xPow)
	ctOdd := eval.mergeRLWERecurse(even, xPow)

	if ctEven == nil && ctOdd == nil {
		return nil
	}

	var tmpEven *Ciphertext
	if ctEven != nil {
		tmpEven = ctEven.CopyNew()
	}

	// ctOdd * X^(N/2^L)
	if ctOdd != nil {

		level := ctOdd.Level()

		//X^(N/2^L)
		ringQ.MulCoeffsMontgomeryLvl(level, ctOdd.Value[0], xPow[len(xPow)-L], ctOdd.Value[0])
		ringQ.MulCoeffsMontgomeryLvl(level, ctOdd.Value[1], xPow[len(xPow)-L], ctOdd.Value[1])

		if ctEven != nil {
			// ctEven + ctOdd * X^(N/2^L)
			ringQ.AddLvl(level, ctEven.Value[0], ctOdd.Value[0], ctEven.Value[0])
			ringQ.AddLvl(level, ctEven.Value[1], ctOdd.Value[1], ctEven.Value[1])

			// phi(ctEven - ctOdd * X^(N/2^L), 2^(L-2))
			ringQ.SubLvl(level, tmpEven.Value[0], ctOdd.Value[0], tmpEven.Value[0])
			ringQ.SubLvl(level, tmpEven.Value[1], ctOdd.Value[1], tmpEven.Value[1])
		}
	}

	if ctEven != nil {

		level := ctEven.Level()

		// if L-2 == -1, then gal = -1
		if L == 1 {
			eval.Automorphism(tmpEven, uint64(2*ringQ.N-1), tmpEven)
		} else {
			eval.Automorphism(tmpEven, eval.params.GaloisElementForColumnRotationBy(1<<(L-2)), tmpEven)
		}

		// ctEven + ctOdd * X^(N/2^L) + phi(ctEven - ctOdd * X^(N/2^L), 2^(L-2))
		ringQ.AddLvl(level, ctEven.Value[0], tmpEven.Value[0], ctEven.Value[0])
		ringQ.AddLvl(level, ctEven.Value[1], tmpEven.Value[1], ctEven.Value[1])
	}

	return ctEven
}
