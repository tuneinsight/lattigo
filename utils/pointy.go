package utils

import (
	"unsafe"
)

// PointyInt creates a new int variable and returns its pointer.
func PointyInt(x int) *int {
	return &x
}

// PointyUint64 creates a new uint64 variable and returns its pointer.
func PointyUint64(x uint64) *uint64 {
	return &x
}

// PointyIntToPointUint64 converts *int to *uint64.
func PointyIntToPointUint64(x *int) *uint64 {
	/* #nosec G103 -- behavior and consequences well understood */
	return (*uint64)(unsafe.Pointer(uintptr(unsafe.Pointer(x))))
}
