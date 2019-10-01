package bfv

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
func min(a, b uint64) (r uint64) {
	if a <= b {
		return a
	}
	return b
}

// max returns the maximum value of the input slice of uint64 values.
func max(a, b uint64) (r uint64) {
	if a >= b {
		return a
	}
	return b
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
