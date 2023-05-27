// Package bfv is a depreciated placeholder package wrapping the bgv package for backward compatibility. This package will be removed in the next major version.
package bfv

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/bgv"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// NewPlaintext allocates a new rlwe.Plaintext.
//
// inputs:
// - params: bfv.Parameters
// - level: the level of the plaintext
//
// output: a newly allocated rlwe.Plaintext at the specified level.
func NewPlaintext(params Parameters, level int) (pt *rlwe.Plaintext) {
	return rlwe.NewPlaintext(params.Parameters, level)
}

// NewCiphertext allocates a new rlwe.Ciphertext.
//
// inputs:
// - params: bfv.Parameters
// - degree: the degree of the ciphertext
// - level: the level of the Ciphertext
//
// output: a newly allocated rlwe.Ciphertext of the specified degree and level.
func NewCiphertext(params Parameters, degree, level int) (ct *rlwe.Ciphertext) {
	return rlwe.NewCiphertext(params.Parameters, degree, level)
}

// NewEncryptor instantiates a new rlwe.Encryptor.
//
// inputs:
// - params: bfv.Parameters
// - key: *rlwe.SecretKey or *rlwe.PublicKey
//
// output: an rlwe.Encryptor instantiated with the provided key.
func NewEncryptor[T *rlwe.SecretKey | *rlwe.PublicKey](params Parameters, key T) rlwe.Encryptor {
	return rlwe.NewEncryptor(params.Parameters, key)
}

// NewPRNGEncryptor instantiates a new rlwe.PRNGEncryptor.
//
// inputs:
// - params: bfv.Parameters
// - key: *rlwe.SecretKey
//
// output: an rlwe.PRNGEncryptor instantiated with the provided key.
func NewPRNGEncryptor(params Parameters, key *rlwe.SecretKey) rlwe.PRNGEncryptor {
	return rlwe.NewPRNGEncryptor(params.Parameters, key)
}

// NewDecryptor instantiates a new rlwe.Decryptor.
//
// inputs:
// - params: bfv.Parameters
// - key: *rlwe.SecretKey
//
// output: an rlwe.Decryptor instantiated with the provided key.
func NewDecryptor(params Parameters, key *rlwe.SecretKey) rlwe.Decryptor {
	return rlwe.NewDecryptor(params.Parameters, key)
}

// NewKeyGenerator instantiates a new rlwe.KeyGenerator.
//
// inputs:
// - params: bfv.Parameters
//
// output: an rlwe.KeyGenerator.
func NewKeyGenerator(params Parameters) *rlwe.KeyGenerator {
	return rlwe.NewKeyGenerator(params.Parameters)
}

// Encoder is a structure that stores the parameters to encode values on a plaintext in a SIMD (Single-Instruction Multiple-Data) fashion.
type Encoder struct {
	*bgv.Encoder
}

// NewEncoder creates a new Encoder from the provided parameters.
func NewEncoder(params Parameters) *Encoder {
	return &Encoder{bgv.NewEncoder(bgv.Parameters(params))}
}

// Evaluator is a struct that holds the necessary elements to perform the homomorphic operations between ciphertexts and/or plaintexts.
// It also holds a memory buffer used to store intermediate computations.
type Evaluator struct {
	*bgv.Evaluator
}

// NewEvaluator creates a new Evaluator, that can be used to do homomorphic
// operations on ciphertexts and/or plaintexts. It stores a memory buffer
// and ciphertexts that will be used for intermediate values.
func NewEvaluator(params Parameters, evk rlwe.EvaluationKeySetInterface) *Evaluator {
	return &Evaluator{bgv.NewEvaluator(bgv.Parameters(params), evk)}
}

// WithKey creates a shallow copy of this Evaluator in which the read-only data-structures are
// shared with the receiver but the EvaluationKey is evaluationKey.
func (eval *Evaluator) WithKey(evk rlwe.EvaluationKeySetInterface) *Evaluator {
	return &Evaluator{eval.Evaluator.WithKey(evk)}
}

// ShallowCopy creates a shallow copy of this Evaluator in which the read-only data-structures are
// shared with the receiver.
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

// MulNew multiplies op0 with op1 without relinearization and returns the result in a new op2.
// The procedure will panic if either op0.Degree or op1.Degree > 1.
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

// MulRelinNew multiplies op0 with op1 with relinearization and returns the result in a new op2.
// The procedure will panic if either op0.Degree or op1.Degree > 1.
// The procedure will panic if the evaluator was not created with an relinearization key.
func (eval *Evaluator) MulRelinNew(op0 *rlwe.Ciphertext, op1 interface{}) (op2 *rlwe.Ciphertext) {
	return eval.Evaluator.MulRelinInvariantNew(op0, op1)
}

// MulRelin multiplies op0 with op1 with relinearization and returns the result in op2.
// The procedure will panic if either op0.Degree or op1.Degree > 1.
// The procedure will panic if op2.Degree != op0.Degree + op1.Degree.
// The procedure will panic if the evaluator was not created with an relinearization key.
func (eval *Evaluator) MulRelin(op0 *rlwe.Ciphertext, op1 interface{}, op2 *rlwe.Ciphertext) {
	eval.Evaluator.MulRelinInvariant(op0, op1, op2)
}

// Polynomial evaluates opOut = P(input).
//
// inputs:
// - input: *rlwe.Ciphertext or *rlwe.PoweBasis
// - pol: *polynomial.Polynomial, *rlwe.Polynomial or *rlwe.PolynomialVector
//
// output: an *rlwe.Ciphertext encrypting pol(input)
func (eval *Evaluator) Polynomial(input, pol interface{}) (opOut *rlwe.Ciphertext, err error) {
	return eval.Evaluator.Polynomial(input, pol, true, eval.Parameters().DefaultScale())
}

// NewLinearTransformEncoder returns new instance of an rlwe.LinearTransformEncoder.
// An rlwe.LinearTransformEncoder is given as input to rlwe
func NewLinearTransformEncoder(ecd *Encoder, diagonals map[int][]uint64) rlwe.LinearTransformEncoder {
	return bgv.NewLinearTransformEncoder(ecd.Encoder, diagonals)
}
