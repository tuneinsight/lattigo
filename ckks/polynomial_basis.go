package ckks

import (
	"github.com/tuneinsight/lattigo/v3/rlwe"
)

type EvaluatorInterface struct {
	Evaluator
	scale rlwe.Scale
}

func (eval *EvaluatorInterface) Add(op0, op1, op2 *rlwe.Ciphertext) {
	eval.Evaluator.Add(&Ciphertext{op0}, &Ciphertext{op1}, &Ciphertext{op2})
}

func (eval *EvaluatorInterface) Sub(op0, op1, op2 *rlwe.Ciphertext) {
	eval.Evaluator.Sub(&Ciphertext{op0}, &Ciphertext{op1}, &Ciphertext{op2})
}

func (eval *EvaluatorInterface) AddConst(op1 *rlwe.Ciphertext, constant interface{}, op2 *rlwe.Ciphertext) {
	eval.Evaluator.AddConst(&Ciphertext{op1}, constant, &Ciphertext{op2})
}

func (eval *EvaluatorInterface) Mul(op0, op1, op2 *rlwe.Ciphertext) {
	eval.Evaluator.Mul(&Ciphertext{op0}, &Ciphertext{op1}, &Ciphertext{op2})
}

func (eval *EvaluatorInterface) MulNew(op0, op1 *rlwe.Ciphertext) (op2 *rlwe.Ciphertext) {
	return eval.Evaluator.MulNew(&Ciphertext{op0}, &Ciphertext{op1}).Ciphertext
}

func (eval *EvaluatorInterface) MulAndAdd(op0, op1, op2 *rlwe.Ciphertext) {
	eval.Evaluator.MulAndAdd(&Ciphertext{op0}, &Ciphertext{op1}, &Ciphertext{op2})
}

func (eval *EvaluatorInterface) MulRelinNew(op0, op1 *rlwe.Ciphertext) (op2 *rlwe.Ciphertext) {
	return eval.Evaluator.MulRelinNew(&Ciphertext{op0}, &Ciphertext{op1}).Ciphertext
}

func (eval *EvaluatorInterface) Rescale(op0, op1 *rlwe.Ciphertext) (err error) {
	return eval.Evaluator.Rescale(&Ciphertext{op0}, eval.scale, &Ciphertext{op1})
}

func (eval *EvaluatorInterface) Relinearize(op0, op1 *rlwe.Ciphertext) {
	eval.Evaluator.Relinearize(&Ciphertext{op0}, &Ciphertext{op1})
}

func (eval *EvaluatorInterface) Parameters() rlwe.Parameters {
	return eval.Evaluator.Parameters().Parameters
}
