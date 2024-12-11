package rlwe

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/structs"
)

// Evaluator is a struct that holds the necessary elements to execute general homomorphic
// operation on RLWE ciphertexts, such as automorphisms, key-switching and relinearization.
type Evaluator struct {
	params Parameters
	EvaluationKeySet
	*EvaluatorBuffers

	automorphismIndex map[uint64][]uint64

	BasisExtender *ring.BasisExtender
	Decomposer    *ring.Decomposer
}

type EvaluatorBuffers struct {
	BuffQPPool  structs.BufferPool[*ringqp.Poly]
	BuffQPool   structs.BufferPool[*ring.Poly]
	BuffBitPool structs.BufferPool[*[]uint64]
	BuffCtPool  structs.BufferPool[*Ciphertext]
}

func newBuffer[T any](f func() T) structs.BufferPool[T] {
	// Uncomment to try with free lists instead of sync pool:
	// nbItemsInPool := 10
	// return structs.NewFreeList(nbItemsInPool, f)
	return structs.NewSyncPool(f)
}
func NewEvaluatorBuffersWithUintPool(params Parameters) *EvaluatorBuffers {
	buff := new(EvaluatorBuffers)
	ringQP := params.RingQP()

	buffUint := newBuffer(func() *[]uint64 {
		buff := make([]uint64, params.RingQ().N())
		return &buff
	})

	buff.BuffQPPool = structs.NewBuffFromUintPool(buffUint,
		func(bp structs.BufferPool[*[]uint64]) *ringqp.Poly {
			return ringQP.NewPolyQPFromUintPool(bp)
		},
		func(bp structs.BufferPool[*[]uint64], poly *ringqp.Poly) {
			ringqp.RecyclePolyQPFromUintPool(bp, poly)
		},
	)
	buff.BuffQPool = structs.NewBuffFromUintPool(buffUint,
		func(bp structs.BufferPool[*[]uint64]) *ring.Poly {
			return ring.NewPolyFromUintPool(bp, params.ringQ.N(), params.ringQ.Level())
		},
		func(bp structs.BufferPool[*[]uint64], poly *ring.Poly) {
			ring.RecyclePolyInUintPool(bp, poly)
		},
	)
	buff.BuffCtPool = structs.NewBuffFromUintPool(buffUint,
		func(bp structs.BufferPool[*[]uint64]) *Ciphertext {
			return NewCiphertextFromUintPool(bp, params, 2, params.MaxLevel())
		},
		func(bp structs.BufferPool[*[]uint64], ct *Ciphertext) {
			RecycleCiphertextInUintPool(bp, ct)
		},
	)
	buff.BuffBitPool = buffUint
	return buff
}

func NewEvaluatorBuffers(params Parameters) *EvaluatorBuffers {

	buff := new(EvaluatorBuffers)
	ringQP := params.RingQP()

	buff.BuffQPPool = newBuffer(func() *ringqp.Poly {
		poly := ringQP.NewPoly()
		return &poly
	})
	buff.BuffQPool = newBuffer(func() *ring.Poly {
		poly := params.RingQ().NewPoly()
		return &poly
	})
	buff.BuffCtPool = newBuffer(func() *Ciphertext {
		return NewCiphertext(params, 2, params.MaxLevel())
	})
	buff.BuffBitPool = newBuffer(func() *[]uint64 {
		buff := make([]uint64, params.RingQ().N())
		return &buff
	})
	return buff
}

// NewEvaluator creates a new [Evaluator].
func NewEvaluator(params ParameterProvider, evk EvaluationKeySet) (eval *Evaluator) {
	eval = new(Evaluator)
	p := params.GetRLWEParameters()
	eval.params = *p
	// All buffer use the same sync.Pool of *[]uint64
	eval.EvaluatorBuffers = NewEvaluatorBuffersWithUintPool(eval.params)
	// Uncomment following line to have one sync.Pool per buffer type
	// eval.EvaluatorBuffers = NewEvaluatorBuffers(eval.params)

	if p.RingP() != nil {
		eval.BasisExtender = ring.NewBasisExtender(p.RingQ(), p.RingP())
	}

	eval.Decomposer = ring.NewDecomposer(p.RingQ(), p.RingP())

	eval.EvaluationKeySet = evk

	var AutomorphismIndex map[uint64][]uint64

	if !utils.IsNil(evk) {
		if galEls := evk.GetGaloisKeysList(); len(galEls) != 0 {
			AutomorphismIndex = make(map[uint64][]uint64)

			N := p.N()
			NthRoot := p.RingQ().NthRoot()

			var err error
			for _, galEl := range galEls {
				if AutomorphismIndex[galEl], err = ring.AutomorphismNTTIndex(N, NthRoot, galEl); err != nil {
					// Sanity check, this error should not happen.
					panic(err)
				}
			}
		}
	}

	eval.automorphismIndex = AutomorphismIndex

	return
}

func (eval *Evaluator) GetRLWEParameters() *Parameters {
	return &eval.params
}

// CheckAndGetGaloisKey returns an error if the [GaloisKey] for the given Galois element is missing or the [EvaluationKey] interface is nil.
func (eval Evaluator) CheckAndGetGaloisKey(galEl uint64) (evk *GaloisKey, err error) {
	if eval.EvaluationKeySet != nil {
		if evk, err = eval.GetGaloisKey(galEl); err != nil {
			return nil, fmt.Errorf("%w: key for galEl %d = 5^{%d} key is missing", err, galEl, eval.params.SolveDiscreteLogGaloisElement(galEl))
		}
	} else {
		return nil, fmt.Errorf("evaluation key interface is nil")
	}

	if eval.automorphismIndex == nil {
		eval.automorphismIndex = map[uint64][]uint64{}
	}

	if _, ok := eval.automorphismIndex[galEl]; !ok {
		if eval.automorphismIndex[galEl], err = ring.AutomorphismNTTIndex(eval.params.N(), eval.params.RingQ().NthRoot(), galEl); err != nil {
			// Sanity check, this error should not happen.
			panic(err)
		}
	}

	return
}

// CheckAndGetRelinearizationKey returns an error if the [RelinearizationKey] is missing or the [EvaluationKey] interface is nil.
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

// InitOutputBinaryOp initializes the output [Element] opOut for receiving the result of a binary operation over
// op0 and op1. The method also performs the following checks:
//
//  1. Inputs are not nil
//  2. MetaData are not nil
//  3. op0.Degree() + op1.Degree() != 0 (i.e at least one [Element] is a ciphertext)
//  4. op0.IsNTT == op1.IsNTT == DefaultNTTFlag
//  5. op0.IsBatched == op1.IsBatched
//
// The opOut metadata are initilized as:
// IsNTT <- DefaultNTTFlag
// IsBatched <- op0.IsBatched
// LogDimensions <- max(op0.LogDimensions, op1.LogDimensions)
//
// The method returns max(op0.Degree(), op1.Degree(), opOut.Degree()) and min(op0.Level(), op1.Level(), opOut.Level())
func (eval Evaluator) InitOutputBinaryOp(op0, op1 *Element[ring.Poly], opInTotalMaxDegree int, opOut *Element[ring.Poly]) (degree, level int, err error) {

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

// InitOutputUnaryOp initializes the output [Element] opOut for receiving the result of a unary operation over
// op0. The method also performs the following checks:
//
//  1. Input and output are not nil
//  2. Inoutp and output Metadata are not nil
//  3. op0.IsNTT == DefaultNTTFlag
//
// The method will also update the metadata of opOut:
//
// IsNTT <- NTTFlag
// IsBatched <- op0.IsBatched
// LogDimensions <- op0.LogDimensions
//
// The method returns max(op0.Degree(), opOut.Degree()) and min(op0.Level(), opOut.Level()).
func (eval Evaluator) InitOutputUnaryOp(op0, opOut *Element[ring.Poly]) (degree, level int, err error) {

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

// ShallowCopy creates a shallow copy of this [Evaluator] in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// evaluators can be used concurrently.
func (eval Evaluator) ShallowCopy() *Evaluator {
	return &Evaluator{
		params:            eval.params,
		Decomposer:        eval.Decomposer,
		BasisExtender:     eval.BasisExtender,
		EvaluatorBuffers:  eval.EvaluatorBuffers,
		EvaluationKeySet:  eval.EvaluationKeySet,
		automorphismIndex: eval.automorphismIndex,
	}
}

// WithKey creates a shallow copy of the receiver [Evaluator] for which the new [EvaluationKey] is evaluationKey
// and where the temporary buffers are shared. The receiver and the returned evaluators cannot be used concurrently.
func (eval Evaluator) WithKey(evk EvaluationKeySet) *Evaluator {

	var AutomorphismIndex map[uint64][]uint64

	if galEls := evk.GetGaloisKeysList(); len(galEls) != 0 {
		AutomorphismIndex = make(map[uint64][]uint64)

		N := eval.params.N()
		NthRoot := eval.params.RingQ().NthRoot()

		var err error
		for _, galEl := range galEls {
			if AutomorphismIndex[galEl], err = ring.AutomorphismNTTIndex(N, NthRoot, galEl); err != nil {
				// Sanity check, this error should not happen.
				panic(err)
			}
		}
	}

	return &Evaluator{
		params:            eval.params,
		EvaluatorBuffers:  eval.EvaluatorBuffers,
		Decomposer:        eval.Decomposer,
		BasisExtender:     eval.BasisExtender,
		EvaluationKeySet:  evk,
		automorphismIndex: AutomorphismIndex,
	}
}

func (eval Evaluator) AutomorphismIndex(galEl uint64) []uint64 {
	return eval.automorphismIndex[galEl]
}

func (eval Evaluator) GetEvaluatorBuffer() *EvaluatorBuffers {
	return eval.EvaluatorBuffers
}

func (eval Evaluator) GetBuffQPPool() structs.BufferPool[*ringqp.Poly] {
	return eval.BuffQPPool
}

func (eval Evaluator) GetBuffCtPool() structs.BufferPool[*Ciphertext] {
	return eval.BuffCtPool
}

func (eval Evaluator) ModDownQPtoQNTT(levelQ, levelP int, p1Q, p1P, p2Q ring.Poly) {
	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, p1Q, p1P, p2Q)
}
