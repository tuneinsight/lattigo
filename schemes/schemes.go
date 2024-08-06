// Package schemes contains the implemented cryptosystems.
package schemes

import "github.com/tuneinsight/lattigo/v6/core/rlwe"

// Encoder is a scheme-agnostic encoding interface.
type Encoder interface {
	Encode(values interface{}, pt *rlwe.Plaintext) error
	Decode(pt *rlwe.Plaintext, values interface{}) error
	Embed(values interface{}, metadata *rlwe.MetaData, polyOut interface{}) error
}

// Evaluator is a scheme-agnostic evaluator interface.
type Evaluator interface {
	rlwe.ParameterProvider
	rlwe.EvaluatorProvider
	Add(op0 *rlwe.Ciphertext, op1 rlwe.Operand, opOut *rlwe.Ciphertext) (err error)
	AddNew(op0 *rlwe.Ciphertext, op1 rlwe.Operand) (opOut *rlwe.Ciphertext, err error)
	Sub(op0 *rlwe.Ciphertext, op1 rlwe.Operand, opOut *rlwe.Ciphertext) (err error)
	SubNew(op0 *rlwe.Ciphertext, op1 rlwe.Operand) (opOut *rlwe.Ciphertext, err error)
	Mul(op0 *rlwe.Ciphertext, op1 rlwe.Operand, opOut *rlwe.Ciphertext) (err error)
	MulNew(op0 *rlwe.Ciphertext, op1 rlwe.Operand) (opOut *rlwe.Ciphertext, err error)
	MulRelin(op0 *rlwe.Ciphertext, op1 rlwe.Operand, opOut *rlwe.Ciphertext) (err error)
	MulRelinNew(op0 *rlwe.Ciphertext, op1 rlwe.Operand) (opOut *rlwe.Ciphertext, err error)
	MulThenAdd(op0 *rlwe.Ciphertext, op1 rlwe.Operand, opOut *rlwe.Ciphertext) (err error)
	Relinearize(op0, op1 *rlwe.Ciphertext) (err error)
	Rescale(op0, op1 *rlwe.Ciphertext) (err error)
}
