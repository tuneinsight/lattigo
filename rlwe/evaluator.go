package rlwe

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// Evaluator is a struct that holds the necessary elements to execute general homomorphic
// operation on RLWE ciphertexts, such as automorphisms, key-switching and relinearization.
type Evaluator struct {
	EvaluationKeySet
	*evaluatorBase
	*evaluatorBuffers

	AutomorphismIndex map[uint64][]uint64

	BasisExtender *ring.BasisExtender
	Decomposer    *ring.Decomposer
}

type evaluatorBase struct {
	params ParametersInterface
}

type evaluatorBuffers struct {
	BuffCt *Ciphertext
	// BuffQP[0-1]: Key-Switch output Key-Switch on the fly decomp(c2)
	// BuffQP[2-5]: Available
	BuffQP        [6]ringqp.Poly
	BuffInvNTT    ring.Poly
	BuffDecompQP  []ringqp.Poly // Memory Buff for the basis extension in hoisting
	BuffBitDecomp []uint64
}

func newEvaluatorBase(params ParametersInterface) *evaluatorBase {
	return &evaluatorBase{
		params: params,
	}
}

func newEvaluatorBuffers(params ParametersInterface) *evaluatorBuffers {

	buff := new(evaluatorBuffers)
	decompRNS := params.DecompRNS(params.MaxLevelQ(), 0)
	ringQP := params.RingQP()

	buff.BuffCt = NewCiphertext(params, 2, params.MaxLevel())

	buff.BuffQP = [6]ringqp.Poly{
		ringQP.NewPoly(),
		ringQP.NewPoly(),
		ringQP.NewPoly(),
		ringQP.NewPoly(),
		ringQP.NewPoly(),
		ringQP.NewPoly(),
	}

	buff.BuffInvNTT = params.RingQ().NewPoly()

	buff.BuffDecompQP = make([]ringqp.Poly, decompRNS)
	for i := 0; i < decompRNS; i++ {
		buff.BuffDecompQP[i] = ringQP.NewPoly()
	}

	buff.BuffBitDecomp = make([]uint64, params.RingQ().N())

	return buff
}

// NewEvaluator creates a new Evaluator.
func NewEvaluator(params ParametersInterface, evk EvaluationKeySet) (eval *Evaluator) {
	eval = new(Evaluator)
	eval.evaluatorBase = newEvaluatorBase(params)
	eval.evaluatorBuffers = newEvaluatorBuffers(params)

	if params.RingP() != nil {
		eval.BasisExtender = ring.NewBasisExtender(params.RingQ(), params.RingP())
		eval.Decomposer = ring.NewDecomposer(params.RingQ(), params.RingP())
	}

	eval.EvaluationKeySet = evk

	var AutomorphismIndex map[uint64][]uint64

	if !utils.IsNil(evk) {
		if galEls := evk.GetGaloisKeysList(); len(galEls) != 0 {
			AutomorphismIndex = make(map[uint64][]uint64)

			N := params.N()
			NthRoot := params.RingQ().NthRoot()

			var err error
			for _, galEl := range galEls {
				if AutomorphismIndex[galEl], err = ring.AutomorphismNTTIndex(N, NthRoot, galEl); err != nil {
					panic(err)
				}
			}
		}
	}

	eval.AutomorphismIndex = AutomorphismIndex

	return
}

// Parameters returns the parameters used to instantiate the target evaluator.
func (eval Evaluator) Parameters() ParametersInterface {
	return eval.params
}

// CheckAndGetGaloisKey returns an error if the GaloisKey for the given Galois element is missing or the EvaluationKey interface is nil.
func (eval Evaluator) CheckAndGetGaloisKey(galEl uint64) (evk *GaloisKey, err error) {
	if eval.EvaluationKeySet != nil {
		if evk, err = eval.GetGaloisKey(galEl); err != nil {
			return nil, fmt.Errorf("%w: key for galEl %d = 5^{%d} key is missing", err, galEl, eval.params.SolveDiscreteLogGaloisElement(galEl))
		}
	} else {
		return nil, fmt.Errorf("evaluation key interface is nil")
	}

	if eval.AutomorphismIndex == nil {
		eval.AutomorphismIndex = map[uint64][]uint64{}
	}

	if _, ok := eval.AutomorphismIndex[galEl]; !ok {
		if eval.AutomorphismIndex[galEl], err = ring.AutomorphismNTTIndex(eval.params.N(), eval.params.RingQ().NthRoot(), galEl); err != nil {
			panic(err)
		}
	}

	return
}

// CheckAndGetRelinearizationKey returns an error if the RelinearizationKey is missing or the EvaluationKey interface is nil.
func (eval Evaluator) CheckAndGetRelinearizationKey() (evk *RelinearizationKey, err error) {
	if eval.EvaluationKeySet != nil {
		if evk, err = eval.GetRelinearizationKey(); err != nil {
			return nil, fmt.Errorf("%w: relineariztion key is missing", err)
		}
	} else {
		return nil, fmt.Errorf("evaluation key interface is nil")
	}

	return
}

// InitOutputBinaryOp initializes the output Operand opOut for receiving the result of a binary operation over
// op0 and op1. The method also performs the following checks:
//
// 1. Inputs are not nil
// 2. MetaData are not nil
// 3. op0.Degree() + op1.Degree() != 0 (i.e at least one operand is a ciphertext)
// 4. op0.IsNTT == op1.IsNTT == DefaultNTTFlag
// 5. op0.IsBatched == op1.IsBatched
//
// The opOut metadata are initilized as:
// IsNTT <- DefaultNTTFlag
// IsBatched <- op0.IsBatched
// LogDimensions <- max(op0.LogDimensions, op1.LogDimensions)
//
// The method returns max(op0.Degree(), op1.Degree(), opOut.Degree()) and min(op0.Level(), op1.Level(), opOut.Level())
func (eval Evaluator) InitOutputBinaryOp(op0, op1 *Operand[ring.Poly], opInTotalMaxDegree int, opOut *Operand[ring.Poly]) (degree, level int, err error) {

	if op0 == nil || op1 == nil || opOut == nil {
		return 0, 0, fmt.Errorf("op0, op1 and opOut cannot be nil")
	}

	if op0.MetaData == nil || op1.MetaData == nil || opOut.MetaData == nil {
		return 0, 0, fmt.Errorf("op0, op1 and opOut MetaData cannot be nil")
	}

	degree = utils.Max(op0.Degree(), op1.Degree())
	degree = utils.Max(degree, opOut.Degree())
	level = utils.Min(op0.Level(), op1.Level())
	level = utils.Min(level, opOut.Level())

	totDegree := op0.Degree() + op1.Degree()

	if totDegree == 0 {
		return 0, 0, fmt.Errorf("op0 and op1 cannot be both plaintexts")
	}

	if totDegree > opInTotalMaxDegree {
		return 0, 0, fmt.Errorf("op0 and op1 total degree cannot exceed %d but is %d", opInTotalMaxDegree, totDegree)
	}

	if op0.El().IsNTT != op1.El().IsNTT || op0.El().IsNTT != eval.params.NTTFlag() {
		return 0, 0, fmt.Errorf("op0.El().IsNTT or op1.El().IsNTT != %t", eval.params.NTTFlag())
	} else {
		opOut.El().IsNTT = op0.El().IsNTT
	}

	if op0.El().IsBatched != op1.El().IsBatched {
		return 0, 0, fmt.Errorf("op1.El().IsBatched != opOut.El().IsBatched")
	} else {
		opOut.El().IsBatched = op0.El().IsBatched
	}

	opOut.El().LogDimensions.Rows = utils.Max(op0.El().LogDimensions.Rows, op1.El().LogDimensions.Rows)
	opOut.El().LogDimensions.Cols = utils.Max(op0.El().LogDimensions.Cols, op1.El().LogDimensions.Cols)

	return
}

// InitOutputUnaryOp initializes the output Operand opOut for receiving the result of a unary operation over
// op0. The method also performs the following checks:
//
// 1. Input and output are not nil
// 2. Inoutp and output Metadata are not nil
// 2. op0.IsNTT == DefaultNTTFlag
//
// The method will also update the metadata of opOut:
//
// IsNTT <- NTTFlag
// IsBatched <- op0.IsBatched
// LogDimensions <- op0.LogDimensions
//
// The method returns max(op0.Degree(), opOut.Degree()) and min(op0.Level(), opOut.Level()).
func (eval Evaluator) InitOutputUnaryOp(op0, opOut *Operand[ring.Poly]) (degree, level int, err error) {

	if op0 == nil || opOut == nil {
		return 0, 0, fmt.Errorf("op0 and opOut cannot be nil")
	}

	if op0.MetaData == nil || opOut.MetaData == nil {
		return 0, 0, fmt.Errorf("op0 and opOut MetaData cannot be nil")
	}

	if op0.El().IsNTT != eval.params.NTTFlag() {
		return 0, 0, fmt.Errorf("op0.IsNTT() != %t", eval.params.NTTFlag())
	} else {
		opOut.El().IsNTT = op0.El().IsNTT
	}

	opOut.El().IsBatched = op0.El().IsBatched
	opOut.El().LogDimensions = op0.El().LogDimensions

	return utils.Max(op0.Degree(), opOut.Degree()), utils.Min(op0.Level(), opOut.Level()), nil
}

// ShallowCopy creates a shallow copy of this Evaluator in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Evaluators can be used concurrently.
func (eval Evaluator) ShallowCopy() *Evaluator {
	return &Evaluator{
		evaluatorBase:     eval.evaluatorBase,
		Decomposer:        eval.Decomposer,
		BasisExtender:     eval.BasisExtender.ShallowCopy(),
		evaluatorBuffers:  newEvaluatorBuffers(eval.params),
		EvaluationKeySet:  eval.EvaluationKeySet,
		AutomorphismIndex: eval.AutomorphismIndex,
	}
}

// WithKey creates a shallow copy of the receiver Evaluator for which the new EvaluationKey is evaluationKey
// and where the temporary buffers are shared. The receiver and the returned Evaluators cannot be used concurrently.
func (eval Evaluator) WithKey(evk EvaluationKeySet) *Evaluator {

	var AutomorphismIndex map[uint64][]uint64

	if galEls := evk.GetGaloisKeysList(); len(galEls) != 0 {
		AutomorphismIndex = make(map[uint64][]uint64)

		N := eval.params.N()
		NthRoot := eval.params.RingQ().NthRoot()

		var err error
		for _, galEl := range galEls {
			if AutomorphismIndex[galEl], err = ring.AutomorphismNTTIndex(N, NthRoot, galEl); err != nil {
				panic(err)
			}
		}
	}

	return &Evaluator{
		evaluatorBase:     eval.evaluatorBase,
		evaluatorBuffers:  eval.evaluatorBuffers,
		Decomposer:        eval.Decomposer,
		BasisExtender:     eval.BasisExtender,
		EvaluationKeySet:  evk,
		AutomorphismIndex: AutomorphismIndex,
	}
}
