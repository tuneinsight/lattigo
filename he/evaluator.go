package he

import (
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// EvaluatorInterface defines a set of common and scheme agnostic homomorphic operations provided by an Evaluator struct.
type EvaluatorInterface interface {
	Add(op0 *rlwe.Ciphertext, op1 interface{}, opOut *rlwe.Ciphertext) (err error)
	Sub(op0 *rlwe.Ciphertext, op1 interface{}, opOut *rlwe.Ciphertext) (err error)
	Mul(op0 *rlwe.Ciphertext, op1 interface{}, opOut *rlwe.Ciphertext) (err error)
	MulNew(op0 *rlwe.Ciphertext, op1 interface{}) (opOut *rlwe.Ciphertext, err error)
	MulRelinNew(op0 *rlwe.Ciphertext, op1 interface{}) (opOut *rlwe.Ciphertext, err error)
	MulThenAdd(op0 *rlwe.Ciphertext, op1 interface{}, opOut *rlwe.Ciphertext) (err error)
	Relinearize(op0, op1 *rlwe.Ciphertext) (err error)
	Rescale(op0, op1 *rlwe.Ciphertext) (err error)
	Parameters() rlwe.ParametersInterface
}

type Evaluator struct {
	rlwe.Evaluator
}

func NewEvaluator(params rlwe.ParametersInterface, evk rlwe.EvaluationKeySet) (eval *Evaluator) {
	return &Evaluator{*rlwe.NewEvaluator(params, evk)}
}

func (eval Evaluator) WithKey(evk rlwe.EvaluationKeySet) *Evaluator {
	return &Evaluator{*eval.Evaluator.WithKey(evk)}
}

func (eval Evaluator) ShallowCopy() *Evaluator {
	return &Evaluator{*eval.Evaluator.ShallowCopy()}
}
