package rlwe

import (
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/ring/distribution"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
)

// ParametersInterface defines a set of common and scheme agnostic methods provided by a Parameter struct.
type ParametersInterface interface {
	RingType() ring.Type
	N() int
	LogN() int
	MaxSlots() [2]int
	MaxLogSlots() [2]int
	PlaintextModulus() uint64
	DefaultScale() Scale
	DefaultPrecision() uint
	DefaultScaleModuliRatio() int
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
	Pow2Base() int
	DecompPw2(levelQ, levelP int) int
	DefaultNTTFlag() bool
	Xe() distribution.Distribution
	Xs() distribution.Distribution
	XsHammingWeight() int
	GaloisElement(k int) (galEl uint64)
	GaloisElements(k []int) (galEls []uint64)
	SolveDiscretLogGaloisElement(galEl uint64) (k int)
	ModInvGaloisElement(galEl uint64) (galElInv uint64)

	Equal(other ParametersInterface) bool
}

// EncoderInterface defines a set of common and scheme agnostic method provided by an Encoder struct.
type EncoderInterface[T any, U *ring.Poly | ringqp.Poly | *Plaintext] interface {
	Encode(values []T, metaData MetaData, output U) (err error)
	Parameters() ParametersInterface
}

// EvaluatorInterface defines a set of common and scheme agnostic homomorphic operations provided by an Evaluator struct.
type EvaluatorInterface interface {
	Add(op0 *Ciphertext, op1 interface{}, op2 *Ciphertext)
	Sub(op0 *Ciphertext, op1 interface{}, op2 *Ciphertext)
	Mul(op0 *Ciphertext, op1 interface{}, op2 *Ciphertext)
	MulNew(op0 *Ciphertext, op1 interface{}) (op2 *Ciphertext)
	MulRelinNew(op0 *Ciphertext, op1 interface{}) (op2 *Ciphertext)
	MulThenAdd(op0 *Ciphertext, op1 interface{}, op2 *Ciphertext)
	Relinearize(op0, op1 *Ciphertext)
	Rescale(op0, op1 *Ciphertext) (err error)
	Parameters() ParametersInterface
}

// PolynomialEvaluatorInterface defines the set of common and scheme agnostic homomorphic operations
// that are required for the encrypted evaluation of plaintext polynomial.
type PolynomialEvaluatorInterface interface {
	EvaluatorInterface
	EvaluatePolynomialVectorFromPowerBasis(targetLevel int, pol *PolynomialVector, pb *PowerBasis, targetScale Scale) (res *Ciphertext, err error)
}
