// Package bfv is a depreciated placeholder package wrapping the bgv package for backward compatibility. This package will be removed in the next major version.
package bfv

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/bgv"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

func NewPlaintext(params Parameters, level int) (pt *rlwe.Plaintext) {
	return rlwe.NewPlaintext(params.Parameters, level)
}

func NewCiphertext(params Parameters, degree, level int) (ct *rlwe.Ciphertext) {
	return rlwe.NewCiphertext(params.Parameters, degree, level)
}

func NewEncryptor(params Parameters, key interface{}) rlwe.Encryptor {
	return rlwe.NewEncryptor(params.Parameters, key)
}

func NewDecryptor(params Parameters, key *rlwe.SecretKey) rlwe.Decryptor {
	return rlwe.NewDecryptor(params.Parameters, key)
}

func NewKeyGenerator(params Parameters) *rlwe.KeyGenerator {
	return rlwe.NewKeyGenerator(params.Parameters)
}

func NewPRNGEncryptor(params Parameters, key *rlwe.SecretKey) rlwe.PRNGEncryptor {
	return rlwe.NewPRNGEncryptor(params.Parameters, key)
}

type Encoder struct {
	*bgv.Encoder
}

func NewEncoder(params Parameters) *Encoder {
	return &Encoder{bgv.NewEncoder(bgv.Parameters(params))}
}

type Evaluator struct {
	*bgv.Evaluator
}

func NewEvaluator(params Parameters, evk rlwe.EvaluationKeySetInterface) *Evaluator {
	return &Evaluator{bgv.NewEvaluator(bgv.Parameters(params), evk)}
}

func (eval *Evaluator) WithKey(evk rlwe.EvaluationKeySetInterface) *Evaluator {
	return &Evaluator{eval.Evaluator.WithKey(evk)}
}

func (eval *Evaluator) ShallowCopy() *Evaluator {
	return &Evaluator{eval.Evaluator.ShallowCopy()}
}

// Mul multiplies op0 with op1 without relinearization and returns the result in op2.
// The procedure will panic if either op0 or op1 are have a degree higher than 1.
// The procedure will panic if op2.Degree != op0.Degree + op1.Degree.
func (eval *Evaluator) Mul(op0 *rlwe.Ciphertext, op1 interface{}, op2 *rlwe.Ciphertext) {
	switch op1 := op1.(type) {
	case rlwe.Operand:
		eval.Evaluator.MulInvariant(op0, op1, op2)
	case uint64:
		eval.Evaluator.Mul(op0, op1, op0)
	default:
		panic(fmt.Sprintf("invalid op1.(Type), expected rlwe.Operand or uint64, but got %T", op1))
	}

}

func (eval *Evaluator) MulNew(op0 *rlwe.Ciphertext, op1 interface{}) (op2 *rlwe.Ciphertext) {
	switch op1 := op1.(type) {
	case rlwe.Operand:
		return eval.Evaluator.MulInvariantNew(op0, op1)
	case uint64:
		return eval.Evaluator.MulNew(op0, op1)
	default:
		panic(fmt.Sprintf("invalid op1.(Type), expected rlwe.Operand or uint64, but got %T", op1))
	}
}

func (eval *Evaluator) MulRelinNew(op0 *rlwe.Ciphertext, op1 interface{}) (op2 *rlwe.Ciphertext) {
	return eval.Evaluator.MulRelinInvariantNew(op0, op1)
}

func (eval *Evaluator) MulRelin(op0 *rlwe.Ciphertext, op1 interface{}, op2 *rlwe.Ciphertext) {
	eval.Evaluator.MulRelinInvariant(op0, op1, op2)
}

func (eval *Evaluator) Polynomial(input, pol interface{}) (opOut *rlwe.Ciphertext, err error) {
	return eval.Evaluator.Polynomial(input, pol, true, eval.Parameters().DefaultScale())
}

func NewLinearTransformEncoder(ecd *Encoder, diagonals map[int][]uint64) rlwe.LinearTransformEncoder {
	return bgv.NewLinearTransformEncoder(ecd.Encoder, diagonals)
}
