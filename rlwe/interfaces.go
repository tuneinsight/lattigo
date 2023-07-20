package rlwe

import (
	"github.com/tuneinsight/lattigo/v4/ring"

	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
)

// ParametersInterface defines a set of common and scheme agnostic methods provided by a Parameter struct.
type ParametersInterface interface {
	RingType() ring.Type
	N() int
	LogN() int
	PlaintextDimensions() ring.Dimensions
	PlaintextLogDimensions() ring.Dimensions
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
	SolveDiscreteLogGaloisElement(galEl uint64) (k int)
	ModInvGaloisElement(galEl uint64) (galElInv uint64)

	Equal(other ParametersInterface) bool
}
