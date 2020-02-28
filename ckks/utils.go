package ckks

import (
	"github.com/ldsec/lattigo/ring"
	"math/big"
	"math/cmplx"
	"math/rand"
)

func exp2pi(x complex128) complex128 {
	return cmplx.Exp(2 * 3.141592653589793 * complex(0, 1) * x)
}

func randomFloat(min, max float64) float64 {
	return min + rand.Float64()*(max-min)
}

func randomComplex(min, max float64) complex128 {
	return complex(randomFloat(min, max), randomFloat(min, max))
}

func scaleUpExact(value float64, n float64, q uint64) (res uint64) {

	var isNegative bool
	var xFlo *big.Float
	var xInt *big.Int

	isNegative = false
	if value < 0 {
		isNegative = true
		xFlo = big.NewFloat(-n * value)
	} else {
		xFlo = big.NewFloat(n * value)
	}

	xFlo.Add(xFlo, big.NewFloat(0.5))

	xInt = new(big.Int)
	xFlo.Int(xInt)
	xInt.Mod(xInt, ring.NewUint(q))

	res = xInt.Uint64()

	if isNegative {
		res = q - res
	}

	return
}

func scaleUpVecExact(values []float64, n float64, moduli []uint64, coeffs [][]uint64) {

	var isNegative bool
	var xFlo *big.Float
	var xInt *big.Int
	tmp := new(big.Int)

	for i := range values {

		if n*values[i] > 1.8446744073709552e+19 {

			isNegative = false
			if values[i] < 0 {
				isNegative = true
				xFlo = big.NewFloat(-n * values[i])
			} else {
				xFlo = big.NewFloat(n * values[i])
			}

			xFlo.Add(xFlo, big.NewFloat(0.5))

			xInt = new(big.Int)
			xFlo.Int(xInt)

			for j := range moduli {
				tmp.Mod(xInt, ring.NewUint(moduli[j]))
				if isNegative {
					coeffs[j][i] = moduli[j] - tmp.Uint64()
				} else {
					coeffs[j][i] = tmp.Uint64()
				}
			}
		} else {

			if values[i] < 0 {
				for j := range moduli {
					coeffs[j][i] = moduli[j] - (uint64(-n*values[i]+0.5) % moduli[j])
				}
			} else {
				for j := range moduli {
					coeffs[j][i] = uint64(n*values[i]+0.5) % moduli[j]
				}
			}
		}
	}

	return
}

func modVec(values []*big.Int, q uint64, coeffs []uint64) {
	tmp := new(big.Int)
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

func genBigIntChain(Q []uint64) (bigintChain []*big.Int) {

	bigintChain = make([]*big.Int, len(Q))
	bigintChain[0] = ring.NewUint(Q[0])
	for i := 1; i < len(Q); i++ {
		bigintChain[i] = ring.NewUint(Q[i])
		bigintChain[i].Mul(bigintChain[i], bigintChain[i-1])
	}
	return
}

// GenSwitchkeysRescalingParams generates the parameters for rescaling the switching keys
func GenSwitchkeysRescalingParams(Q, P []uint64) (params []uint64) {

	params = make([]uint64, len(Q))

	PBig := ring.NewUint(1)
	for _, pj := range P {
		PBig.Mul(PBig, ring.NewUint(pj))
	}

	tmp := ring.NewUint(0)

	for i := 0; i < len(Q); i++ {

		params[i] = tmp.Mod(PBig, ring.NewUint(Q[i])).Uint64()
		params[i] = ring.ModExp(params[i], Q[i]-2, Q[i])
		params[i] = ring.MForm(params[i], Q[i], ring.BRedParams(Q[i]))
	}

	return
}

// GenModuli generates the appropriate primes from the parameters using generateCKKSPrimes, such that all the primes are different.
func GenModuli(params *Parameters) (Q []uint64, P []uint64) {

	// Extracts all the different primes bit size and maps their number
	primesbitlen := make(map[uint64]uint64)
	for _, qi := range params.LogQi {

		primesbitlen[qi]++

		if qi > 60 {
			panic("cannot GenModuli: the provided LogQi must be smaller than 61")
		}
	}

	for _, pj := range params.LogPi {

		primesbitlen[pj]++

		if pj > 61 {
			panic("cannot GenModuli: the provided LogPi must be smaller than 62")
		}
	}

	// For each bit-size, finds that many primes
	primes := make(map[uint64][]uint64)
	for key, value := range primesbitlen {
		primes[key] = ring.GenerateNTTPrimes(key, params.LogN, value)
	}

	// Assigns the primes to the ckks moduli chain
	Q = make([]uint64, len(params.LogQi))
	for i, qi := range params.LogQi {
		Q[i] = primes[qi][0]
		primes[qi] = primes[qi][1:]
	}

	// Assigns the primes to the special primes list for the keys context
	P = make([]uint64, len(params.LogPi))
	for i, pj := range params.LogPi {
		P[i] = primes[pj][0]
		primes[pj] = primes[pj][1:]
	}

	return Q, P
}

func sliceBitReverseInPlaceComplex128(slice []complex128, N uint64) {

	var bit, j uint64

	for i := uint64(1); i < N; i++ {

		bit = N >> 1

		for j >= bit {
			j -= bit
			bit >>= 1
		}

		j += bit

		if i < j {
			slice[i], slice[j] = slice[j], slice[i]
		}
	}
}

func max(array []complex128) (m float64) {
	m = real(array[0])
	for _, i := range array[1:] {
		if real(i) > m {
			m = real(i)
		}
	}
	return
}

func min(array []complex128) (m float64) {
	m = real(array[0])
	for _, i := range array[1:] {
		if real(i) < m {
			m = real(i)
		}
	}
	return
}
