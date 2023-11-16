package utils

import (
	"unsafe"

	cs "golang.org/x/exp/constraints"
)

type Number interface {
	cs.Complex | cs.Float | cs.Integer
}

// Pointy creates a new T variable and returns its pointer.
func Pointy[T Number](x T) *T {
	return &x
}

// PointyIntToPointUint64 converts *int to *uint64.
func PointyIntToPointUint64(x *int) *uint64 {
	/* #nosec G103 -- behavior and consequences well understood, pointer type cast */
	return (*uint64)(unsafe.Pointer(uintptr(unsafe.Pointer(x))))
}
