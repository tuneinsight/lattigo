package rlwe

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// Operand is a common interface for Ciphertext and Plaintext types.
type Operand interface {
	El() *Ciphertext
	Degree() int
	Level() int
	GetMetaData() *MetaData
}

// Evaluator is a struct that holds the necessary elements to execute general homomorphic
// operation on RLWE ciphertexts, such as automorphisms, key-switching and relinearization.
type Evaluator struct {
	EvaluationKeySetInterface
	*evaluatorBase
	*evaluatorBuffers

	AutomorphismIndex map[uint64][]uint64

	BasisExtender *ring.BasisExtender
	Decomposer    *ring.Decomposer
}

type evaluatorBase struct {
	params Parameters
}

type evaluatorBuffers struct {
	BuffCt *Ciphertext
	// BuffQP[0-1]: Key-Switch output Key-Switch on the fly decomp(c2)
	// BuffQP[2-5]: Available
	BuffQP        [6]ringqp.Poly
	BuffInvNTT    *ring.Poly
	BuffDecompQP  []ringqp.Poly // Memory Buff for the basis extension in hoisting
	BuffBitDecomp []uint64
}

func newEvaluatorBase(params Parameters) *evaluatorBase {
	return &evaluatorBase{
		params: params,
	}
}

func newEvaluatorBuffers(params Parameters) *evaluatorBuffers {

	buff := new(evaluatorBuffers)
	decompRNS := params.DecompRNS(params.MaxLevelQ(), params.MaxLevelP())
	ringQP := params.RingQP()

	buff.BuffCt = NewCiphertext(params, 2, params.MaxLevel())

	buff.BuffQP = [6]ringqp.Poly{
		*ringQP.NewPoly(),
		*ringQP.NewPoly(),
		*ringQP.NewPoly(),
		*ringQP.NewPoly(),
		*ringQP.NewPoly(),
		*ringQP.NewPoly(),
	}

	buff.BuffInvNTT = params.RingQ().NewPoly()

	buff.BuffDecompQP = make([]ringqp.Poly, decompRNS)
	for i := 0; i < decompRNS; i++ {
		buff.BuffDecompQP[i] = *ringQP.NewPoly()
	}

	buff.BuffBitDecomp = make([]uint64, params.RingQ().N())

	return buff
}

// NewEvaluator creates a new Evaluator.
func NewEvaluator(params Parameters, evk EvaluationKeySetInterface) (eval *Evaluator) {
	eval = new(Evaluator)
	eval.evaluatorBase = newEvaluatorBase(params)
	eval.evaluatorBuffers = newEvaluatorBuffers(params)

	if params.RingP() != nil {
		eval.BasisExtender = ring.NewBasisExtender(params.RingQ(), params.RingP())
		eval.Decomposer = ring.NewDecomposer(params.RingQ(), params.RingP())
	}

	eval.EvaluationKeySetInterface = evk

	var AutomorphismIndex map[uint64][]uint64

	if !utils.IsNil(evk) {
		if galEls := evk.GetGaloisKeysList(); len(galEls) != 0 {
			AutomorphismIndex = make(map[uint64][]uint64)

			N := params.N()
			NthRoot := params.RingQ().NthRoot()

			for _, galEl := range galEls {
				AutomorphismIndex[galEl] = ring.AutomorphismNTTIndex(N, NthRoot, galEl)
			}
		}
	}

	eval.AutomorphismIndex = AutomorphismIndex

	return
}

// Parameters returns the parameters used to instantiate the target evaluator.
func (eval *Evaluator) Parameters() Parameters {
	return eval.params
}

// CheckAndGetGaloisKey returns an error if the GaloisKey for the given Galois element is missing or the EvaluationKey interface is nil.
func (eval *Evaluator) CheckAndGetGaloisKey(galEl uint64) (evk *GaloisKey, err error) {
	if eval.EvaluationKeySetInterface != nil {
		if evk, err = eval.GetGaloisKey(galEl); err != nil {
			return nil, fmt.Errorf("%w: key for galEl %d = 5^{%d} key is missing", err, galEl, eval.params.RotationFromGaloisElement(galEl))
		}
	} else {
		return nil, fmt.Errorf("evaluation key interface is nil")
	}

	return
}

// CheckAndGetRelinearizationKey returns an error if the RelinearizationKey is missing or the EvaluationKey interface is nil.
func (eval *Evaluator) CheckAndGetRelinearizationKey() (evk *RelinearizationKey, err error) {
	if eval.EvaluationKeySetInterface != nil {
		if evk, err = eval.GetRelinearizationKey(); err != nil {
			return nil, fmt.Errorf("%w: relineariztion key is missing", err)
		}
	} else {
		return nil, fmt.Errorf("evaluation key interface is nil")
	}

	return
}

// CheckBinary checks that:
//
//		Inputs are not nil
//		op0.Degree() + op1.Degree() != 0 (i.e at least one operand is a ciphertext)
//		opOut.Degree() >= opOutMinDegree
//		op0.IsNTT == op1.IsNTT == DefaultNTTFlag
//	 op0.EncodingDomain == op1.EncodingDomain
//
// The method will also update the MetaData of OpOut:
//
// IsNTT <- DefaultNTTFlag
// EncodingDomain <- op0.EncodingDomain
// LogSlots <- max(op0.LogSlots, op1.LogSlots)
//
// and returns max(op0.Degree(), op1.Degree(), opOut.Degree()) and min(op0.Level(), op1.Level(), opOut.Level())
func (eval *Evaluator) CheckBinary(op0, op1, opOut Operand, opOutMinDegree int) (degree, level int) {

	degree = utils.Max(op0.Degree(), op1.Degree())
	degree = utils.Max(degree, opOut.Degree())
	level = utils.Min(op0.Level(), op1.Level())
	level = utils.Min(level, opOut.Level())

	if op0 == nil || op1 == nil || opOut == nil {
		panic("op0, op1 and opOut cannot be nil")
	}

	if op0.Degree()+op1.Degree() == 0 {
		panic("op0 and op1 cannot be both plaintexts")
	}

	if op0.El().IsNTT != op1.El().IsNTT || op0.El().IsNTT != eval.params.DefaultNTTFlag() {
		panic(fmt.Sprintf("op0.El().IsNTT or op1.El().IsNTT != %t", eval.params.DefaultNTTFlag()))
	} else {
		opOut.El().IsNTT = op0.El().IsNTT
	}

	if op0.El().IsNTT != op1.El().IsNTT || op0.El().IsNTT != eval.params.DefaultNTTFlag() {
		panic(fmt.Sprintf("op0.El().IsNTT or op1.El().IsNTT != %t", eval.params.DefaultNTTFlag()))
	} else {
		opOut.El().IsNTT = op0.El().IsNTT
	}

	if op0.El().EncodingDomain != op1.El().EncodingDomain {
		panic("op1.El().EncodingDomain != op2.El().EncodingDomain")
	} else {
		opOut.El().EncodingDomain = op0.El().EncodingDomain
	}

	opOut.El().LogSlots = utils.MaxInt(op0.El().LogSlots, op1.El().LogSlots)

	return
}

// CheckUnary checks that op0 and opOut are not nil and that op0 respects the DefaultNTTFlag.
//
// The method will also update the metadata of opOut:
//
// IsNTT <- DefaultNTTFlag
// EncodingDomain <- op0.EncodingDomain
// LogSlots <- op0.LogSlots
//
// Also returns max(op0.Degree(), opOut.Degree()) and min(op0.Level(), opOut.Level()).
func (eval *Evaluator) CheckUnary(op0, opOut Operand) (degree, level int) {

	if op0 == nil || opOut == nil {
		panic("op0 and opOut cannot be nil")
	}

	if op0.El().IsNTT != eval.params.DefaultNTTFlag() {
		panic(fmt.Sprintf("op0.IsNTT() != %t", eval.params.DefaultNTTFlag()))
	} else {
		opOut.El().IsNTT = op0.El().IsNTT
	}

	opOut.El().EncodingDomain = op0.El().EncodingDomain

	opOut.El().LogSlots = op0.El().LogSlots

	return utils.MaxInt(op0.Degree(), opOut.Degree()), utils.MinInt(op0.Level(), opOut.Level())
}

// ShallowCopy creates a shallow copy of this Evaluator in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Evaluators can be used concurrently.
func (eval *Evaluator) ShallowCopy() *Evaluator {
	return &Evaluator{
		evaluatorBase:             eval.evaluatorBase,
		Decomposer:                eval.Decomposer,
		BasisExtender:             eval.BasisExtender.ShallowCopy(),
		evaluatorBuffers:          newEvaluatorBuffers(eval.params),
		EvaluationKeySetInterface: eval.EvaluationKeySetInterface,
		AutomorphismIndex:         eval.AutomorphismIndex,
	}
}

// WithKey creates a shallow copy of the receiver Evaluator for which the new EvaluationKey is evaluationKey
// and where the temporary buffers are shared. The receiver and the returned Evaluators cannot be used concurrently.
func (eval *Evaluator) WithKey(evk EvaluationKeySetInterface) *Evaluator {

	var AutomorphismIndex map[uint64][]uint64

	if galEls := evk.GetGaloisKeysList(); len(galEls) != 0 {
		AutomorphismIndex = make(map[uint64][]uint64)

		N := eval.params.N()
		NthRoot := eval.params.RingQ().NthRoot()

		for _, galEl := range galEls {
			AutomorphismIndex[galEl] = ring.AutomorphismNTTIndex(N, NthRoot, galEl)
		}
	}

	return &Evaluator{
		evaluatorBase:             eval.evaluatorBase,
		evaluatorBuffers:          eval.evaluatorBuffers,
		Decomposer:                eval.Decomposer,
		BasisExtender:             eval.BasisExtender,
		EvaluationKeySetInterface: evk,
		AutomorphismIndex:         AutomorphismIndex,
	}
}
