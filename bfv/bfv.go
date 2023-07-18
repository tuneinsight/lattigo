// Package bfv provides an RNS-accelerated implementation of the Fan-Vercauteren version of Brakerski's (BFV) scale-invariant homomorphic encryption scheme.
// The BFV scheme enables SIMD modular arithmetic over encrypted vectors or integers.
package bfv

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/bgv"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// NewPlaintext allocates a new rlwe.Plaintext.
//
// inputs:
//   - params: an rlwe.ParametersInterface interface
//   - level: the level of the plaintext
//
// output: a newly allocated rlwe.Plaintext at the specified level.
//
// Note: the user can update the field `MetaData` to set a specific scaling factor,
// plaintext dimensions (if applicable) or encoding domain, before encoding values
// on the created plaintext.
func NewPlaintext(params rlwe.ParametersInterface, level int) (pt *rlwe.Plaintext) {
	return rlwe.NewPlaintext(params, level)
}

// NewCiphertext allocates a new rlwe.Ciphertext.
//
// inputs:
//
//   - params: an rlwe.ParametersInterface interface
//
//   - degree: the degree of the ciphertext
//
//   - level: the level of the Ciphertext
//
// output: a newly allocated rlwe.Ciphertext of the specified degree and level.
func NewCiphertext(params rlwe.ParametersInterface, degree, level int) (ct *rlwe.Ciphertext) {
	return rlwe.NewCiphertext(params, degree, level)
}

// NewEncryptor instantiates a new rlwe.Encryptor.
//
// inputs:
//   - params: an rlwe.ParametersInterface interface
//   - key: *rlwe.SecretKey or *rlwe.PublicKey
//
// output: an rlwe.Encryptor instantiated with the provided key.
func NewEncryptor(params rlwe.ParametersInterface, key rlwe.EncryptionKey) (*rlwe.Encryptor, error) {
	return rlwe.NewEncryptor(params, key)
}

// NewDecryptor instantiates a new rlwe.Decryptor.
//
// inputs:
//   - params: an rlwe.ParametersInterface interface
//   - key: *rlwe.SecretKey
//
// output: an rlwe.Decryptor instantiated with the provided key.
func NewDecryptor(params rlwe.ParametersInterface, key *rlwe.SecretKey) (*rlwe.Decryptor, error) {
	return rlwe.NewDecryptor(params, key)
}

// NewKeyGenerator instantiates a new rlwe.KeyGenerator.
//
// inputs:
//   - params: an rlwe.ParametersInterface interface
//
// output: an rlwe.KeyGenerator.
func NewKeyGenerator(params rlwe.ParametersInterface) *rlwe.KeyGenerator {
	return rlwe.NewKeyGenerator(params)
}

// Encoder is a structure that stores the parameters to encode values on a plaintext in a SIMD (Single-Instruction Multiple-Data) fashion.
type Encoder struct {
	*bgv.Encoder
}

// NewEncoder creates a new Encoder from the provided parameters.
func NewEncoder(params Parameters) *Encoder {
	return &Encoder{bgv.NewEncoder(params.Parameters)}
}

// ShallowCopy creates a shallow copy of this Encoder in which the read-only data-structures are
// shared with the receiver.
func (e Encoder) ShallowCopy() *Encoder {
	return &Encoder{Encoder: e.Encoder.ShallowCopy()}
}

type encoder[T int64 | uint64, U *ring.Poly | ringqp.Poly | *rlwe.Plaintext] struct {
	*Encoder
}

func (e encoder[T, U]) Encode(values []T, metadata *rlwe.MetaData, output U) (err error) {
	return e.Encoder.Embed(values, false, metadata, output)
}

// Evaluator is a struct that holds the necessary elements to perform the homomorphic operations between ciphertexts and/or plaintexts.
// It also holds a memory buffer used to store intermediate computations.
type Evaluator struct {
	*bgv.Evaluator
}

// NewEvaluator creates a new Evaluator, that can be used to do homomorphic
// operations on ciphertexts and/or plaintexts. It stores a memory buffer
// and ciphertexts that will be used for intermediate values.
func NewEvaluator(params Parameters, evk rlwe.EvaluationKeySet) *Evaluator {
	return &Evaluator{bgv.NewEvaluator(params.Parameters, evk)}
}

// WithKey creates a shallow copy of this Evaluator in which the read-only data-structures are
// shared with the receiver but the EvaluationKey is evaluationKey.
func (eval Evaluator) WithKey(evk rlwe.EvaluationKeySet) *Evaluator {
	return &Evaluator{eval.Evaluator.WithKey(evk)}
}

// ShallowCopy creates a shallow copy of this Evaluator in which the read-only data-structures are
// shared with the receiver.
func (eval Evaluator) ShallowCopy() *Evaluator {
	return &Evaluator{eval.Evaluator.ShallowCopy()}
}

// Mul multiplies op0 with op1 without relinearization and returns the result in opOut.
// inputs:
//   - op0: an *rlwe.Ciphertext
//   - op1: an rlwe.OperandInterface[ring.Poly], an uint64 or an []uint64 slice (of size at most N where N is the smallest integer satisfying T = 1 mod 2N)
//   - opOut: an *rlwe.Ciphertext
//
// The procedure will return an error if either op0 or op1 are have a degree higher than 1.
// The procedure will return an error if opOut.Degree != op0.Degree + op1.Degree.
func (eval Evaluator) Mul(op0 *rlwe.Ciphertext, op1 interface{}, opOut *rlwe.Ciphertext) (err error) {
	switch op1 := op1.(type) {
	case rlwe.OperandInterface[ring.Poly], []uint64:
		return eval.Evaluator.MulScaleInvariant(op0, op1, opOut)
	case uint64, int64, int:
		return eval.Evaluator.Mul(op0, op1, op0)
	default:
		return fmt.Errorf("invalid op1.(Type), expected rlwe.OperandInterface[ring.Poly], []uint64 or uint64, int64, int, but got %T", op1)
	}

}

// MulNew multiplies op0 with op1 without relinearization and returns the result in a new opOut.
// inputs:
//   - op0: an *rlwe.Ciphertext
//   - op1: an rlwe.OperandInterface[ring.Poly], an uint64 or an []uint64 slice (of size at most N where N is the smallest integer satisfying T = 1 mod 2N)
//   - opOut: an *rlwe.Ciphertext
//
// The procedure will return an error if either op0.Degree or op1.Degree > 1.
func (eval Evaluator) MulNew(op0 *rlwe.Ciphertext, op1 interface{}) (opOut *rlwe.Ciphertext, err error) {
	switch op1 := op1.(type) {
	case rlwe.OperandInterface[ring.Poly], []uint64:
		return eval.Evaluator.MulScaleInvariantNew(op0, op1)
	case uint64, int64, int:
		return eval.Evaluator.MulNew(op0, op1)
	default:
		return nil, fmt.Errorf("invalid op1.(Type), expected rlwe.OperandInterface[ring.Poly], []uint64 or  uint64, int64, int, but got %T", op1)
	}
}

// MulRelinNew multiplies op0 with op1 with relinearization and returns the result in a new opOut.
// inputs:
//   - op0: an *rlwe.Ciphertext
//   - op1: an rlwe.OperandInterface[ring.Poly], an uint64 or an []uint64 slice (of size at most N where N is the smallest integer satisfying T = 1 mod 2N)
//   - opOut: an *rlwe.Ciphertext
//
// The procedure will return an error if either op0.Degree or op1.Degree > 1.
// The procedure will return an error if the evaluator was not created with an relinearization key.
func (eval Evaluator) MulRelinNew(op0 *rlwe.Ciphertext, op1 interface{}) (opOut *rlwe.Ciphertext, err error) {
	return eval.Evaluator.MulRelinScaleInvariantNew(op0, op1)
}

// MulRelin multiplies op0 with op1 with relinearization and returns the result in opOut.
// inputs:
//   - op0: an *rlwe.Ciphertext
//   - op1: an rlwe.OperandInterface[ring.Poly], an uint64 or an []uint64 slice (of size at most N where N is the smallest integer satisfying T = 1 mod 2N)
//   - opOut: an *rlwe.Ciphertext
//
// The procedure will return an error if either op0.Degree or op1.Degree > 1.
// The procedure will return an error if opOut.Degree != op0.Degree + op1.Degree.
// The procedure will return an error if the evaluator was not created with an relinearization key.
func (eval Evaluator) MulRelin(op0 *rlwe.Ciphertext, op1 interface{}, opOut *rlwe.Ciphertext) (err error) {
	return eval.Evaluator.MulRelinScaleInvariant(op0, op1, opOut)
}

// NewPowerBasis creates a new PowerBasis from the input ciphertext.
// The input ciphertext is treated as the base monomial X used to
// generate the other powers X^{n}.
func NewPowerBasis(ct *rlwe.Ciphertext) rlwe.PowerBasis {
	return rlwe.NewPowerBasis(ct, bignum.Monomial)
}

// Polynomial evaluates opOut = P(input).
//
// inputs:
//   - input: *rlwe.Ciphertext or *rlwe.PoweBasis
//   - pol: *bignum.Polynomial, *rlwe.Polynomial or *rlwe.PolynomialVector
//
// output: an *rlwe.Ciphertext encrypting pol(input)
func (eval Evaluator) Polynomial(input, pol interface{}) (opOut *rlwe.Ciphertext, err error) {
	return eval.Evaluator.Polynomial(input, pol, true, eval.Evaluator.Parameters().PlaintextScale())
}

type PolynomialEvaluator struct {
	bgv.PolynomialEvaluator
}

func NewPolynomialEvaluator(eval *Evaluator) *PolynomialEvaluator {
	return &PolynomialEvaluator{PolynomialEvaluator: *bgv.NewPolynomialEvaluator(eval.Evaluator, false)}
}

// NewLinearTransformation allocates a new LinearTransformation with zero values according to the parameters specified by the LinearTranfromationParameters.
func NewLinearTransformation[T int64 | uint64](params rlwe.ParametersInterface, lt rlwe.LinearTranfromationParameters[T]) rlwe.LinearTransformation {
	return rlwe.NewLinearTransformation(params, lt)
}

// EncodeLinearTransformation encodes a linear transformation on a pre-allocated linear transformation.
// The method will return an error if the non-zero diagonals between the pre-allocated linear transformation and the parameters of the linear transformation to encode do not match.
func EncodeLinearTransformation[T int64 | uint64](allocated rlwe.LinearTransformation, params rlwe.LinearTranfromationParameters[T], ecd *Encoder) (err error) {
	return rlwe.EncodeLinearTransformation[T](allocated, params, &encoder[T, ringqp.Poly]{ecd})
}
