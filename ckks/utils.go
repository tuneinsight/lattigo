package ckks

import (
	"errors"
	"github.com/ldsec/lattigo/ring"
	"math/big"
	"math/rand"
)

func randomFloat(min, max float64) float64 {
	return min + rand.Float64()*(max-min)
}

func randomComplex(min, max float64) complex128 {
	return complex(randomFloat(min, max), randomFloat(min, max))
}

func scaleUpExact(value float64, n float64, q uint64) (res uint64) {

	var is_negative bool
	var xFlo *big.Float
	var xInt *big.Int

	is_negative = false
	if value < 0 {
		is_negative = true
		xFlo = big.NewFloat(-n * value)
	} else {
		xFlo = big.NewFloat(n * value)
	}

	xInt = ring.NewUint(0)
	xFlo.Int(xInt)
	xInt.Mod(xInt, ring.NewUint(q))

	res = xInt.Uint64()

	if is_negative {
		res = q - res
	}

	return
}

func scaleUpVecExact(values []float64, n float64, moduli []uint64, coeffs [][]uint64) {

	var is_negative bool
	var xFlo *big.Float
	var xInt *big.Int
	tmp := new(big.Int)

	for i := range values {

		if n*values[i] > 1.8446744073709552e+19 {

			is_negative = false
			if values[i] < 0 {
				is_negative = true
				xFlo = big.NewFloat(-n * values[i])
			} else {
				xFlo = big.NewFloat(n * values[i])
			}

			xInt = ring.NewUint(0)
			xFlo.Int(xInt)

			for j := range moduli {
				tmp.Mod(xInt, ring.NewUint(moduli[j]))
				if is_negative {
					coeffs[j][i] = moduli[j] - tmp.Uint64()
				} else {
					coeffs[j][i] = tmp.Uint64()
				}
			}
		} else {

			if values[i] < 0 {
				for j := range moduli {
					coeffs[j][i] = moduli[j] - (uint64(-n*values[i]) % moduli[j])
				}
			} else {
				for j := range moduli {
					coeffs[j][i] = uint64(n*values[i]) % moduli[j]
				}
			}
		}
	}

	return
}

func modVec(values []*big.Int, q uint64, coeffs []uint64) {
	tmp := ring.NewUint(0)
	for i := range values {
		coeffs[i] = tmp.Mod(values[i], ring.NewUint(q)).Uint64()
	}
}

// Divides x by n^2, returns a float
func scaleDown(coeff *big.Int, n float64) (x float64) {

	x, _ = new(big.Float).SetInt(coeff).Float64()
	x /= n

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

func sliceBitReverseInPlaceComplex128(slice []complex128, N uint64) {

	var bit uint64

	i, j := uint64(1), uint64(0)
	for i < N {

		bit = N >> 1

		for j >= bit {
			j -= bit
			bit >>= 1
		}

		j += bit

		if i < j {
			slice[i], slice[j] = slice[j], slice[i]
		}

		i++
	}
}
