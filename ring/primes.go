package ring

import (
	"fmt"
	"math/bits"

	"github.com/ldsec/lattigo/v2/utils"
)

// IsPrime applies a Miller-Rabin test on the given uint64 variable, returning true if the input is probably prime, and false otherwise.
func IsPrime(num uint64) bool {

	if num < 2 {
		return false
	}

	for _, smallPrime := range smallPrimes {
		if num == smallPrime {
			return true
		}
	}

	for _, smallPrime := range smallPrimes {
		if num%smallPrime == 0 {
			return false
		}
	}

	s := num - 1
	k := 0
	for (s & 1) == 0 {
		s >>= 1
		k++
	}

	bredParams := BRedParams(num)
	var mask, b uint64
	mask = (1 << uint64(bits.Len64(num))) - 1

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	for trial := 0; trial < 50; trial++ {

		b = RandUniform(prng, num-1, mask)

		for b < 2 {
			b = RandUniform(prng, num-1, mask)
		}

		x := ModExp(b, s, num)

		if x != 1 {
			i := 0
			for x != num-1 {

				if i == k-1 {
					return false
				}

				i++
				x = BRed(x, x, num, bredParams)
			}
		}
	}
	return true
}

// GenerateNTTPrimes generates n NthRoot NTT friendly primes given logQ = size of the primes.
// It will return all the appropriate primes, up to the number of n, with the
// best available deviation from the base power of 2 for the given n.
func GenerateNTTPrimes(logQ, NthRoot, n uint64) (primes []uint64) {

	if logQ > 61 {
		panic("logQ must be between 1 and 61")
	}

	if logQ == 61 {
		return GenerateNTTPrimesP(logQ, NthRoot, n)
	}

	return GenerateNTTPrimesQ(logQ, NthRoot, n)
}

// NextNTTPrime returns the next NthRoot NTT prime after q.
// The input q must be itself an NTT prime for the given NthRoot.
func NextNTTPrime(q, NthRoot uint64) (qNext uint64, err error) {

	qNext = q + NthRoot

	for !IsPrime(qNext) {

		qNext += NthRoot

		if bits.Len64(qNext) > 61 {
			return 0, fmt.Errorf("Next NTT prime exceeds the maximum bit-size of 61 bits")
		}
	}

	return qNext, nil
}

// PreviousNTTPrime returns the previous NthRoot NTT prime after q.
// The input q must be itself an NTT prime for the given NthRoot.
func PreviousNTTPrime(q, NthRoot uint64) (qPrev uint64, err error) {

	if q < NthRoot {
		return 0, fmt.Errorf("Previous NTT prime is smaller than NthRoot")
	}

	qPrev = q - NthRoot

	for !IsPrime(qPrev) {

		if q < NthRoot {
			return 0, fmt.Errorf("Previous NTT prime is smaller than NthRoot")
		}

		qPrev -= NthRoot
	}

	return qPrev, nil
}

// GenerateNTTPrimesQ generates "levels" different NthRoot NTT-friendly
// primes starting from 2**LogQ and alternating between upward and downward.
func GenerateNTTPrimesQ(logQ, NthRoot, levels uint64) (primes []uint64) {

	var nextPrime, previousPrime, Qpow2 uint64
	var checkfornextprime, checkforpreviousprime bool

	primes = []uint64{}

	Qpow2 = 1 << logQ

	nextPrime = Qpow2 + 1
	previousPrime = Qpow2 + 1

	checkfornextprime = true
	checkforpreviousprime = true

	for true {

		if !(checkfornextprime || checkforpreviousprime) {
			panic("GenerateNTTPrimesQ error: cannot generate enough primes for the given parameters")
		}

		if checkfornextprime {

			if nextPrime > 0xffffffffffffffff-NthRoot {

				checkfornextprime = false

			} else {

				nextPrime += NthRoot

				if IsPrime(nextPrime) {

					primes = append(primes, nextPrime)

					if uint64(len(primes)) == levels {
						return
					}
				}
			}
		}

		if checkforpreviousprime {

			if previousPrime < NthRoot {

				checkforpreviousprime = false

			} else {

				previousPrime -= NthRoot

				if IsPrime(previousPrime) {

					primes = append(primes, previousPrime)

					if uint64(len(primes)) == levels {
						return
					}
				}
			}

		}
	}

	return
}

// GenerateNTTPrimesP generates "levels" different NthRoot NTT-friendly
// primes starting from 2**LogP and downward.
// Special case were primes close to 2^{LogP} but with a smaller bit-size than LogP are sought.
func GenerateNTTPrimesP(logP, NthRoot, n uint64) (primes []uint64) {

	var x, Ppow2 uint64

	primes = []uint64{}

	Ppow2 = 1 << logP

	x = Ppow2 + 1

	for true {

		// We start by subtracting 2N to ensure that the prime bit-length is smaller than LogP

		if x > NthRoot {

			x -= NthRoot

			if IsPrime(x) {

				primes = append(primes, x)

				if uint64(len(primes)) == n {
					return primes
				}
			}

		} else {
			panic("GenerateNTTPrimesP error: cannot generate enough primes for the given parameters")
		}
	}

	return
}
