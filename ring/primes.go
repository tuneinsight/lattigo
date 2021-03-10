package ring

import (
	"fmt"
	"math/bits"
)

// IsPrime applies the Baillie-PSW, which is 100% accurate for numbers bellow 2^64.
func IsPrime(x uint64) bool {
	return NewUint(x).ProbablyPrime(0)
}

// GenerateNTTPrimes generates n NthRoot NTT friendly primes given logQ = size of the primes.
// It will return all the appropriate primes, up to the number of n, with the
// best available deviation from the base power of 2 for the given n.
func GenerateNTTPrimes(logQ, NthRoot, n int) (primes []uint64) {

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
func NextNTTPrime(q uint64, NthRoot int) (qNext uint64, err error) {

	qNext = q + uint64(NthRoot)

	for !IsPrime(qNext) {

		qNext += uint64(NthRoot)

		if bits.Len64(qNext) > 61 {
			return 0, fmt.Errorf("next NTT prime exceeds the maximum bit-size of 61 bits")
		}
	}

	return qNext, nil
}

// PreviousNTTPrime returns the previous NthRoot NTT prime after q.
// The input q must be itself an NTT prime for the given NthRoot.
func PreviousNTTPrime(q uint64, NthRoot int) (qPrev uint64, err error) {

	if q < uint64(NthRoot) {
		return 0, fmt.Errorf("previous NTT prime is smaller than NthRoot")
	}

	qPrev = q - uint64(NthRoot)

	for !IsPrime(qPrev) {

		if q < uint64(NthRoot) {
			return 0, fmt.Errorf("previous NTT prime is smaller than NthRoot")
		}

		qPrev -= uint64(NthRoot)
	}

	return qPrev, nil
}

// GenerateNTTPrimesQ generates "levels" different NthRoot NTT-friendly
// primes starting from 2**LogQ and alternating between upward and downward.
func GenerateNTTPrimesQ(logQ, NthRoot, levels int) (primes []uint64) {

	var nextPrime, previousPrime, Qpow2 uint64
	var checkfornextprime, checkforpreviousprime bool

	primes = []uint64{}

	Qpow2 = uint64(1 << logQ)

	nextPrime = Qpow2 + 1
	previousPrime = Qpow2 + 1

	checkfornextprime = true
	checkforpreviousprime = true

	for {

		if !(checkfornextprime || checkforpreviousprime) {
			panic("generateNTTPrimesQ error: cannot generate enough primes for the given parameters")
		}

		if checkfornextprime {

			if nextPrime > 0xffffffffffffffff-uint64(NthRoot) {

				checkfornextprime = false

			} else {

				nextPrime += uint64(NthRoot)

				if IsPrime(nextPrime) {

					primes = append(primes, nextPrime)

					if len(primes) == levels {
						return
					}
				}
			}
		}

		if checkforpreviousprime {

			if previousPrime < uint64(NthRoot) {

				checkforpreviousprime = false

			} else {

				previousPrime -= uint64(NthRoot)

				if IsPrime(previousPrime) {

					primes = append(primes, previousPrime)

					if len(primes) == levels {
						return
					}
				}
			}

		}
	}
}

// GenerateNTTPrimesP generates "levels" different NthRoot NTT-friendly
// primes starting from 2**LogP and downward.
// Special case were primes close to 2^{LogP} but with a smaller bit-size than LogP are sought.
func GenerateNTTPrimesP(logP, NthRoot, n int) (primes []uint64) {

	var x, Ppow2 uint64

	primes = []uint64{}

	Ppow2 = uint64(1 << logP)

	x = Ppow2 + 1

	for {

		// We start by subtracting 2N to ensure that the prime bit-length is smaller than LogP

		if x > uint64(NthRoot) {

			x -= uint64(NthRoot)

			if IsPrime(x) {

				primes = append(primes, x)

				if len(primes) == n {
					return primes
				}
			}

		} else {
			panic("generateNTTPrimesP error: cannot generate enough primes for the given parameters")
		}
	}
}
