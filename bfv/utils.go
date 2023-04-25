package bfv

import (
	"math/big"

	"github.com/tuneinsight/lattigo/ring"
)

// GenLiftParams generates the lifting parameters.
func GenLiftParams(context *ring.Context, t uint64) (deltaMont []uint64) {

	delta := new(big.Int).Quo(context.ModulusBigint, ring.NewUint(t))

	deltaMont = make([]uint64, len(context.Modulus))

	tmp := new(big.Int)
	bredParams := context.GetBredParams()
	for i, Qi := range context.Modulus {
		deltaMont[i] = tmp.Mod(delta, ring.NewUint(Qi)).Uint64()
		deltaMont[i] = ring.MForm(deltaMont[i], Qi, bredParams[i])
	}

	return
}

// GenModuli generates the appropriate primes from the parameters using generateCKKSPrimes such that all primes are different.
func GenModuli(params *Parameters) (Q []uint64, P []uint64, QMul []uint64) {

	// Extracts all the different primes bit-size and maps their number
	primesbitlen := make(map[uint64]uint64)

	for _, qi := range params.LogQi {

		primesbitlen[qi]++

		if qi > 60 {
			panic("cannot GenModuli: the provided LogQi must be smaller than 61")
		}
	}

	for _, pj := range params.LogPi {

		primesbitlen[pj]++

		if pj > 60 {
			panic("cannot GenModuli: the provided LogPi must be smaller than 61")
		}
	}

	for _, qi := range params.LogQiMul {

		primesbitlen[qi]++

		if qi > 60 {
			panic("cannot GenModuli: the provided LogQiMul must be smaller than 61")
		}
	}

	// For each bit-size, it finds that many primes
	primes := make(map[uint64][]uint64)
	for key, value := range primesbitlen {
		primes[key] = ring.GenerateNTTPrimes(key, params.LogN, value)
	}

	// Assigns the primes to the CKKS moduli chain
	Q = make([]uint64, len(params.LogQi))
	for i, qi := range params.LogQi {
		Q[i] = primes[qi][0]
		primes[qi] = primes[qi][1:]
	}

	// Assigns the primes to the special primes list for the the keys context
	P = make([]uint64, len(params.LogPi))
	for i, pj := range params.LogPi {
		P[i] = primes[pj][0]
		primes[pj] = primes[pj][1:]
	}

	QMul = make([]uint64, len(params.LogQiMul))
	for i, qi := range params.LogQiMul {
		QMul[i] = primes[qi][0]
		primes[qi] = primes[qi][1:]
	}

	return Q, P, QMul
}
