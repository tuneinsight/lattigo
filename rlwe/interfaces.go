package rlwe

import (
	"github.com/tuneinsight/lattigo/v4/ring"

	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

// ParametersInterface defines a set of common and scheme agnostic methods provided by a Parameter struct.
type ParametersInterface interface {
	RingType() ring.Type
	N() int
	LogN() int
	PlaintextDimensions() [2]int
	PlaintextLogDimensions() [2]int
	PlaintextSlots() int
	PlaintextLogSlots() int
	PlaintextModulus() uint64
	PlaintextScale() Scale
	PlaintextPrecision() uint
	PlaintextScaleToModuliRatio() int
	MaxLevel() int
	MaxLevelQ() int
	MaxLevelP() int
	Q() []uint64
	P() []uint64
	QCount() int
	PCount() int
	QiOverflowMargin(levelQ int) int
	PiOverflowMargin(levelP int) int
	RingQ() *ring.Ring
	RingP() *ring.Ring
	RingQP() *ringqp.Ring
	DecompRNS(levelQ, levelP int) int
	DecompPw2(levelQ, levelP, Base2Decomposition int) int
	NTTFlag() bool
	Xe() ring.DistributionParameters
	Xs() ring.DistributionParameters
	XsHammingWeight() int
	GaloisElement(k int) (galEl uint64)
	GaloisElements(k []int) (galEls []uint64)
	GaloisElementsForLinearTransform(nonZeroDiagonals []int, LogSlots, LogBSGSRatio int) (galEls []uint64)
	SolveDiscreteLogGaloisElement(galEl uint64) (k int)
	ModInvGaloisElement(galEl uint64) (galElInv uint64)

	Equal(other ParametersInterface) bool
}

// DecryptorInterface is a generic RLWE decryption interface.
type DecryptorInterface interface {
	Decrypt(ct *Ciphertext, pt *Plaintext)
	DecryptNew(ct *Ciphertext) (pt *Plaintext)
	ShallowCopy() DecryptorInterface
	WithKey(sk *SecretKey) Decryptor
}

// EncryptorInterface a generic RLWE encryption interface.
type EncryptorInterface interface {
	Encrypt(pt *Plaintext, ct interface{})
	EncryptZero(ct interface{})

	EncryptZeroNew(level int) (ct *Ciphertext)
	EncryptNew(pt *Plaintext) (ct *Ciphertext)

	ShallowCopy() EncryptorInterface
	WithKey(key interface{}) EncryptorInterface
}

// PRNGEncryptorInterface is an interface for encrypting RLWE ciphertexts from a secret-key and
// a pre-determined PRNG. An Encryptor constructed from a secret-key complies to this
// interface.
type PRNGEncryptorInterface interface {
	EncryptorInterface
	WithPRNG(prng sampling.PRNG) PRNGEncryptorInterface
}

// EncoderInterface defines a set of common and scheme agnostic method provided by an Encoder struct.
type EncoderInterface[T any, U *ring.Poly | ringqp.Poly | *Plaintext] interface {
	Encode(values []T, metaData MetaData, output U) (err error)
	Parameters() ParametersInterface
}

// EvaluatorInterface defines a set of common and scheme agnostic homomorphic operations provided by an Evaluator struct.
type EvaluatorInterface interface {
	Add(op0 *Ciphertext, op1 interface{}, opOut *Ciphertext)
	Sub(op0 *Ciphertext, op1 interface{}, opOut *Ciphertext)
	Mul(op0 *Ciphertext, op1 interface{}, opOut *Ciphertext)
	MulNew(op0 *Ciphertext, op1 interface{}) (opOut *Ciphertext)
	MulRelinNew(op0 *Ciphertext, op1 interface{}) (opOut *Ciphertext)
	MulThenAdd(op0 *Ciphertext, op1 interface{}, opOut *Ciphertext)
	Relinearize(op0, op1 *Ciphertext)
	Rescale(op0, op1 *Ciphertext) (err error)
	Parameters() ParametersInterface
}

// PolynomialEvaluatorInterface defines the set of common and scheme agnostic homomorphic operations
// that are required for the encrypted evaluation of plaintext polynomial.
type PolynomialEvaluatorInterface interface {
	EvaluatorInterface
	EvaluatePolynomialVectorFromPowerBasis(targetLevel int, pol PolynomialVector, pb PowerBasis, targetScale Scale) (res *Ciphertext, err error)
}
