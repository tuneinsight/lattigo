// Package bfv provides an RNS-accelerated implementation of the Fan-Vercauteren version of Brakerski's (BFV) scale-invariant homomorphic encryption scheme.
// The BFV scheme enables SIMD modular arithmetic over encrypted vectors or integers.
package bfv

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/schemes/bgv"
)

// NewPlaintext allocates a new rlwe.Plaintext from the BFV parameters, at the
// specified level. If the level argument is not provided, the plaintext is
// initialized at level params.MaxLevelQ().
//
// The plaintext is initialized with its metadata so that it can be passed to a,
// bfv.Encoder. Before doing so, the user can update the MetaData field to set
// a specific scaling factor,
// plaintext dimensions (if applicable) or encoding domain.
func NewPlaintext(params Parameters, level ...int) (pt *rlwe.Plaintext) {
	pt = rlwe.NewPlaintext(params, level...)
	pt.IsBatched = true
	pt.Scale = params.DefaultScale()
	pt.LogDimensions = params.LogMaxDimensions()
	return
}

// NewCiphertext allocates a new rlwe.Ciphertext from the BFV parameters,
// at the specified level and ciphertext degree. If the level argument is not
// provided, the ciphertext is initialized at level params.MaxLevelQ().
//
// To create a ciphertext for encrypting a new message, the ciphertext should be
// at degree 1.
func NewCiphertext(params Parameters, degree int, level ...int) (ct *rlwe.Ciphertext) {
	ct = rlwe.NewCiphertext(params, degree, level...)
	ct.IsBatched = true
	ct.Scale = params.DefaultScale()
	ct.LogDimensions = params.LogMaxDimensions()
	return
}

// NewEncryptor instantiates a new rlwe.Encryptor from the given BFV parameters and
// encryption key. This key can be either a *rlwe.SecretKey or a *rlwe.PublicKey.
func NewEncryptor(params Parameters, key rlwe.EncryptionKey) *rlwe.Encryptor {
	return rlwe.NewEncryptor(params, key)
}

// NewDecryptor instantiates a new rlwe.Decryptor from the given BFV parameters and
// secret decryption key.
func NewDecryptor(params Parameters, key *rlwe.SecretKey) *rlwe.Decryptor {
	return rlwe.NewDecryptor(params, key)
}

// NewKeyGenerator instantiates a new rlwe.KeyGenerator from the given
// BFV parameters.
func NewKeyGenerator(params Parameters) *rlwe.KeyGenerator {
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
//   - op1:
//   - rlwe.ElementInterface[ring.Poly]
//   - *big.Int, uint64, int64, int
//   - []uint64 or []int64 (of size at most N where N is the smallest integer satisfying PlaintextModulus = 1 mod 2N)
//   - opOut: an *rlwe.Ciphertext
//
// The procedure will return an error if either op0 or op1 are have a degree higher than 1.
// The procedure will return an error if opOut.Degree != op0.Degree + op1.Degree.
func (eval Evaluator) Mul(op0 *rlwe.Ciphertext, op1 rlwe.Operand, opOut *rlwe.Ciphertext) (err error) {
	switch op1 := op1.(type) {
	case rlwe.ElementInterface[ring.Poly], []uint64:
		return eval.Evaluator.MulScaleInvariant(op0, op1, opOut)
	case uint64, int64, int:
		return eval.Evaluator.Mul(op0, op1, op0)
	default:
		return fmt.Errorf("invalid op1.(Type), expected rlwe.ElementInterface[ring.Poly], []uint64 or uint64, int64, int, but got %T", op1)
	}

}

// MulNew multiplies op0 with op1 without relinearization and returns the result in a new opOut.
// inputs:
//   - op0: an *rlwe.Ciphertext
//   - op1:
//   - rlwe.ElementInterface[ring.Poly]
//   - *big.Int, uint64, int64, int
//   - []uint64 or []int64 (of size at most N where N is the smallest integer satisfying PlaintextModulus = 1 mod 2N)
//   - opOut: an *rlwe.Ciphertext
//
// The procedure will return an error if either op0.Degree or op1.Degree > 1.
func (eval Evaluator) MulNew(op0 *rlwe.Ciphertext, op1 rlwe.Operand) (opOut *rlwe.Ciphertext, err error) {
	switch op1 := op1.(type) {
	case rlwe.ElementInterface[ring.Poly], []uint64:
		return eval.Evaluator.MulScaleInvariantNew(op0, op1)
	case uint64, int64, int:
		return eval.Evaluator.MulNew(op0, op1)
	default:
		return nil, fmt.Errorf("invalid op1.(Type), expected rlwe.ElementInterface[ring.Poly], []uint64 or  uint64, int64, int, but got %T", op1)
	}
}

// MulRelinNew multiplies op0 with op1 with relinearization and returns the result in a new opOut.
// inputs:
//   - op0: an *rlwe.Ciphertext
//   - op1:
//   - rlwe.ElementInterface[ring.Poly]
//   - *big.Int, uint64, int64, int
//   - []uint64 or []int64 (of size at most N where N is the smallest integer satisfying PlaintextModulus = 1 mod 2N)
//   - opOut: an *rlwe.Ciphertext
//
// The procedure will return an error if either op0.Degree or op1.Degree > 1.
// The procedure will return an error if the evaluator was not created with an relinearization key.
func (eval Evaluator) MulRelinNew(op0 *rlwe.Ciphertext, op1 rlwe.Operand) (opOut *rlwe.Ciphertext, err error) {
	return eval.Evaluator.MulRelinScaleInvariantNew(op0, op1)
}

// MulRelin multiplies op0 with op1 with relinearization and returns the result in opOut.
// inputs:
//   - op0: an *rlwe.Ciphertext
//   - op1:
//   - rlwe.ElementInterface[ring.Poly]
//   - *big.Int, uint64, int64, int
//   - []uint64 or []int64 (of size at most N where N is the smallest integer satisfying PlaintextModulus = 1 mod 2N)
//   - opOut: an *rlwe.Ciphertext
//
// The procedure will return an error if either op0.Degree or op1.Degree > 1.
// The procedure will return an error if opOut.Degree != op0.Degree + op1.Degree.
// The procedure will return an error if the evaluator was not created with an relinearization key.
func (eval Evaluator) MulRelin(op0 *rlwe.Ciphertext, op1 rlwe.Operand, opOut *rlwe.Ciphertext) (err error) {
	return eval.Evaluator.MulRelinScaleInvariant(op0, op1, opOut)
}

// Rescale does nothing when instantiated with the BFV scheme.
func (eval Evaluator) Rescale(op0, op1 *rlwe.Ciphertext) (err error) {
	return nil
}
