package ckks

import (
	"errors"
	"github.com/ldsec/lattigo/ring"
	"math"
	"math/big"
	"math/bits"
)

// Multiplies x by 2^n and returns the result mod q
// Unaffected by overflows, arbitrary n allowed.
// Expects inputs in the range of uint64.
func scaleUp(x float64, n, q uint64) uint64 {

	var a uint64
	if n < 53 {
		a = n
	} else {
		a = 53
	}

	var is_negative bool

	if x < 0 {
		is_negative = true
		x *= -1
	}

	x_int := uint64(x)
	x_flo := x - float64(x_int)
	x_int %= q

	for i := uint64(0); i < a; i++ { // stops at 53 to avoid an overflow of the float
		x_int <<= 1
		x_flo *= 2
		if x_int >= q {
			x_int -= q
		}
	}

	x_int += uint64(math.Round(x_flo))
	x_int %= q

	if n > 53 { // continues with int and and float merged without a risk of overflow
		for i := a; i < n; i++ {
			x_int <<= 1
			if x_int >= q {
				x_int -= q
			}
		}
	}

	if is_negative {
		return q - x_int
	} else {
		return x_int
	}
}

// Divides x by n^2, returns a float
func scaleDown(coeff *ring.Int, n uint64) (x float64) {

	if n > 53 { // if n > float64 precision, then first reduce the integer to 53 bits, and the scales down
		coeff.Rsh(coeff, n-53)
		n -= (n - 53)
	}

	x, _ = new(big.Float).SetInt(&coeff.Value).Float64()

	x /= float64(uint64(1 << n))

	return
}

// Generates CKKS Primes given logQ = size of the primes, logN = size of N and level, the number
// of levels we require. Will return all the appropriate primes, up to the number of level, with the
// best avaliable precision for the given level.
func GenerateCKKSPrimes(logQ, logN, levels uint64) ([]uint64, error) {

	if logQ > 60 {
		return nil, errors.New("error : logQ must be between 1 and 62")
	}

	var x, y, Qpow2, _2N uint64

	primes := []uint64{}

	Qpow2 = 1 << logQ

	_2N = 2 << logN

	x = Qpow2 + 1
	y = Qpow2 + 1

	for true {

		if ring.IsPrime(y) {
			primes = append(primes, y)
			if uint64(len(primes)) == levels {
				return primes, nil
			}
		}

		y -= _2N

		if ring.IsPrime(x) {
			primes = append(primes, x)
			if uint64(len(primes)) == levels {
				return primes, nil
			}
		}

		x += _2N
	}

	return primes, nil
}

func equalslice64(a, b []uint64) bool {

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

func equalslice8(a, b []uint8) bool {

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
	return bits.Reverse64(index) >> (64 - bitLen)
}

func hammingWeight64(x uint64) uint64 {
	x -= (x >> 1) & 0x5555555555555555
	x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333)
	x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0f
	return ((x * 0x0101010101010101) & 0xffffffffffffffff) >> 56
}
