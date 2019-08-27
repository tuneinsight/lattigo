package bfv

import (
	"bytes"
	"encoding/binary"
	"github.com/lca1/lattigo/ring"
	"golang.org/x/crypto/blake2b"
	"math/bits"
)

func Hash(data []uint64) (value []byte, err error) {
	hash, err := blake2b.New512(nil)
	buff := make([]byte, 8)
	for _, x := range data {
		binary.BigEndian.PutUint64(buff, x)
		hash.Write(buff)
	}
	value = hash.Sum(nil)
	return

}

func VerifyHash(hash0, hash1 []byte) bool {
	if res := bytes.Compare(hash0, hash1); res != 0 {
		return false
	} else {
		return true
	}
}

func EqualSlice(a, b []uint64) bool {

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

func min(values []uint64) (r uint64) {
	r = values[0]
	for _, i := range values[1:] {
		if i < r {
			r = i
		}
	}
	return
}

func max(values []uint64) (r uint64) {
	r = values[0]
	for _, i := range values[1:] {
		if i > r {
			r = i
		}
	}
	return
}

func bitReverse64(index, bitLen uint64) uint64 {
	indexReverse := uint64(0)
	for i := uint64(0); i < bitLen; i++ {
		if (index>>i)&1 != 0 {
			indexReverse |= 1 << (bitLen - 1 - i)
		}
	}
	return indexReverse
}

func hammingWeight64(x uint64) uint64 {
	x -= (x >> 1) & 0x5555555555555555
	x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333)
	x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0f
	return ((x * 0x0101010101010101) & 0xffffffffffffffff) >> 56
}

func modexp(x, e, p uint64) (result uint64) {
	params := ring.BRedParams(p)
	result = 1
	for i := e; i > 0; i >>= 1 {
		if i&1 == 1 {
			result = ring.BRed(result, x, p, params)
		}
		x = ring.BRed(x, x, p, params)
	}
	return result
}

// Returns (x*2^n)%q where x is in montgomery form
func PowerOf2(x, n, q, qInv uint64) (r uint64) {
	ahi, alo := x>>(64-n), x<<n
	R := alo * qInv
	H, _ := bits.Mul64(R, q)
	r = ahi - H + q
	if r >= q {
		r -= q
	}
	return
}
