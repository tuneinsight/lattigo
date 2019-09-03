package bfv

import (
	"bytes"
	"encoding/binary"
	"golang.org/x/crypto/blake2b"
)

// hash hashes a slice of uint64 values (data) and returns the digest in bytes.
func hash(data []uint64) (digest []byte, err error) {
	hash, err := blake2b.New512(nil)
	buff := make([]byte, 8)
	for _, x := range data {
		binary.BigEndian.PutUint64(buff, x)
		hash.Write(buff)
	}
	digest = hash.Sum(nil)
	return

}

// verifyHash compares to byte slices and return true if they are equal, else false.
func verifyHash(hash0, hash1 []byte) bool {
	if res := bytes.Compare(hash0, hash1); res != 0 {
		return false
	} else {
		return true
	}
}

// equalslice compares two slices of uint64 values, and return true if they are equal, else false.
func equalslice(a, b []uint64) bool {

	if len(a) != len(b) {
		return false
	}

	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// min returns the minimum value of the input slice of uint64 values.
func min(values []uint64) (r uint64) {
	r = values[0]
	for _, i := range values[1:] {
		if i < r {
			r = i
		}
	}
	return
}

// max returns the maximum value of the input slice of uint64 values.
func max(values []uint64) (r uint64) {
	r = values[0]
	for _, i := range values[1:] {
		if i > r {
			r = i
		}
	}
	return
}

// bitReverse64 returns the bit-reverse value of the input value, within a context of 2^bitLen.
func bitReverse64(index, bitLen uint64) uint64 {
	indexReverse := uint64(0)
	for i := uint64(0); i < bitLen; i++ {
		if (index>>i)&1 != 0 {
			indexReverse |= 1 << (bitLen - 1 - i)
		}
	}
	return indexReverse
}

// hammingWeight64 returns the hammingweight if the input value.
func hammingWeight64(x uint64) uint64 {
	x -= (x >> 1) & 0x5555555555555555
	x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333)
	x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0f
	return ((x * 0x0101010101010101) & 0xffffffffffffffff) >> 56
}
