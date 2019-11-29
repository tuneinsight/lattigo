package bfv

import (
	"github.com/ldsec/lattigo/ring"
	"math/big"
)

// GenLiftParams generates the lifting parameters
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
func GenModuli(params *Parameters) (Q1 []uint64, P []uint64, Q2 []uint64) {

	// Extracts all the different primes bit size and maps their number
	primesbitlen := make(map[uint64]uint64)
	for i, qi := range params.Q1 {

		primesbitlen[uint64(qi)]++

		if uint64(params.Q1[i]) > 60 {
			panic("provided moduli must be smaller than 61")
		}
	}

	for _, pj := range params.P {
		primesbitlen[uint64(pj)]++

		if uint64(pj) > 60 {
			panic("provided P must be smaller than 61")
		}
	}

	for i, qi := range params.Q2 {

		primesbitlen[uint64(qi)]++

		if uint64(params.Q2[i]) > 60 {
			panic("provided moduli must be smaller than 61")
		}
	}

	// For each bitsize, finds that many primes
	primes := make(map[uint64][]uint64)
	for key, value := range primesbitlen {
		primes[key] = ring.GenerateNTTPrimes(key, uint64(params.LogN), value)
	}

	// Assigns the primes to the ckks moduli chain
	Q1 = make([]uint64, len(params.Q1))
	for i, qi := range params.Q1 {
		Q1[i] = primes[uint64(params.Q1[i])][0]
		primes[uint64(qi)] = primes[uint64(qi)][1:]
	}

	// Assigns the primes to the special primes list for the the keyscontext
	P = make([]uint64, len(params.P))
	for i, pj := range params.P {
		P[i] = primes[uint64(pj)][0]
		primes[uint64(pj)] = primes[uint64(pj)][1:]
	}

	Q2 = make([]uint64, len(params.Q2))
	for i, qi := range params.Q2 {
		Q2[i] = primes[uint64(params.Q2[i])][0]
		primes[uint64(qi)] = primes[uint64(qi)][1:]
	}

	return Q1, P, Q2
}
