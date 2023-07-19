package he

import (
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

type Evaluator struct {
	rlwe.Evaluator
}

func NewEvaluator(params rlwe.ParametersInterface, evk rlwe.EvaluationKeySet) (eval *Evaluator) {
	return &Evaluator{*rlwe.NewEvaluator(params, evk)}
}

func (eval Evaluator) WithKey(evk rlwe.EvaluationKeySet) *Evaluator {
	return &Evaluator{*eval.Evaluator.WithKey(evk)}
}
