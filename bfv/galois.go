package bfv

import (
	"github.com/ldsec/lattigo/v2/ring"
)

type Rotation interface {
	GaloisElement() uint64
}

func NewRotationRight(params *Parameters, k uint64) Rotation {
	return &rotationRight{rotation{ring.ModExp(GaloisGen, k, 2*params.N())}, k}
}

func NewRotationLeft(params *Parameters, k uint64) Rotation {
	return &rotationLeft{rotation{ring.ModExp(GaloisGen, 2*params.N()-k, 2*params.N())}, k}
}

func NewRotationRow(params *Parameters) Rotation {
	return &rotationRow{rotation{2*params.N() - 1}}
}

type rotation struct {
	galEl uint64
}

type rotationLeft struct {
	rotation
	k uint64
}

type rotationRight struct {
	rotation
	k uint64
}

type rotationRow struct {
	rotation
}

func (r *rotation) GaloisElement() uint64 {
	return r.galEl
}

// RotationType is a type used to represent the rotations types.
type RotationType int

// Constants for rotation types
const (
	RotationRight = iota + 1
	RotationLeft
	RotationRow
)

func getGaloisElementForRotationRev(rotType RotationType, k uint64, n uint64) uint64 {
	switch rotType {
	case RotationLeft:
		return ring.ModExp(GaloisGen, 2*n-k, 2*n)
	case RotationRight:
		return ring.ModExp(GaloisGen, k, 2*n)
	case RotationRow:
		return 2*n - 1
	}
	return 0
}

func getGaloisElementForRotation(rotType RotationType, k uint64, n uint64) uint64 {
	switch rotType {
	case RotationLeft:
		return ring.ModExp(GaloisGen, k, 2*n)
	case RotationRight:
		return ring.ModExp(GaloisGen, 2*n-k, 2*n)
	case RotationRow:
		return 2*n - 1
	}
	return 0
}
