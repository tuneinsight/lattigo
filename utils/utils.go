package utils

import (
	"math/bits"
)

// EqualSliceUint64 checks the equality between two uint64 slices.
func EqualSliceUint64(a, b []uint64) (v bool) {
	v = true
	for i := range a {
		v = v && (a[i] == b[i])
	}
	return
}

// EqualSliceUint8 checks the equality between two uint8 slices.
func EqualSliceUint8(a, b []uint8) (v bool) {
	v = true
	for i := range a {
		v = v && (a[i] == b[i])
	}
	return
}

// IsInSliceUint64 checks if x is in slice.
func IsInSliceUint64(x uint64, slice []uint64) (v bool) {
	for i := range slice {
		v = v || (slice[i] == x)
	}
	return
}

// MinUint64 returns the minimum value of the input slice of uint64 values.
func MinUint64(a, b uint64) (r uint64) {
	if a <= b {
		return a
	}
	return b
}

// MaxUint64 returns the maximum value of the input slice of uint64 values.
func MaxUint64(a, b uint64) (r uint64) {
	if a >= b {
		return a
	}
	return b
}

// MaxFloat64 returns the maximum value of the input slice of uint64 values.
func MaxFloat64(a, b float64) (r float64) {
	if a >= b {
		return a
	}
	return b
}

// BitReverse64 returns the bit-reverse value of the input value, within a ring of 2^bitLen.
func BitReverse64(index, bitLen uint64) uint64 {
	return bits.Reverse64(index) >> (64 - bitLen)
}

// HammingWeight64 returns the hammingweight if the input value.
func HammingWeight64(x uint64) uint64 {
	x -= (x >> 1) & 0x5555555555555555
	x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333)
	x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0f
	return ((x * 0x0101010101010101) & 0xffffffffffffffff) >> 56
}
