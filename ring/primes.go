package ring

import (
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// IsPrime applies the Baillie-PSW, which is 100% accurate for numbers bellow 2^64.
func IsPrime(x uint64) bool {
	return bignum.NewInt(x).ProbablyPrime(0)
}

// NTTFriendlyPrimesGenerator is a struct used to generate NTT friendly primes.
type NTTFriendlyPrimesGenerator struct {
	Size                           float64
	NextPrime, PrevPrime, NthRoot  uint64
	CheckNextPrime, CheckPrevPrime bool
}

// NewNTTFriendlyPrimesGenerator instantiates a new NTTFriendlyPrimesGenerator.
// Primes generated are of the form 2^{BitSize} +/- k * {NthRoot} + 1.
func NewNTTFriendlyPrimesGenerator(BitSize, NthRoot uint64) NTTFriendlyPrimesGenerator {

	CheckNextPrime := true
	CheckPrevPrime := true

	NextPrime := uint64(1<<BitSize) + 1
	PrevPrime := uint64(1<<BitSize) + 1

	if NextPrime > 0xffffffffffffffff-NthRoot {
		CheckNextPrime = false
	}

	if PrevPrime < NthRoot {
		CheckPrevPrime = false
	}

	PrevPrime -= NthRoot

	return NTTFriendlyPrimesGenerator{
		CheckNextPrime: CheckNextPrime,
		CheckPrevPrime: CheckPrevPrime,
		NthRoot:        NthRoot,
		NextPrime:      NextPrime,
		PrevPrime:      PrevPrime,
		Size:           float64(BitSize),
	}
}

// NextUpstreamPrimes returns the next k primes of the form 2^{BitSize} + k * {NthRoot} + 1.
func (n *NTTFriendlyPrimesGenerator) NextUpstreamPrimes(k int) (primes []uint64, err error) {
	primes = make([]uint64, k)

	for i := range primes {
		if primes[i], err = n.NextUpstreamPrime(); err != nil {
			return
		}
	}

	return
}

// NextDownstreamPrimes returns the next k primes of the form 2^{BitSize} - k * {NthRoot} + 1.
func (n *NTTFriendlyPrimesGenerator) NextDownstreamPrimes(k int) (primes []uint64, err error) {
	primes = make([]uint64, k)

	for i := range primes {
		if primes[i], err = n.NextDownstreamPrime(); err != nil {
			return
		}
	}

	return
}

// NextAlternatingPrimes returns the next k primes of the form 2^{BitSize} +/- k * {NthRoot} + 1.
func (n *NTTFriendlyPrimesGenerator) NextAlternatingPrimes(k int) (primes []uint64, err error) {
	primes = make([]uint64, k)

	for i := range primes {
		if primes[i], err = n.NextAlternatingPrime(); err != nil {
			return
		}
	}

	return
}

// NextUpstreamPrime returns the next prime of the form 2^{BitSize} + k * {NthRoot} + 1.
func (n *NTTFriendlyPrimesGenerator) NextUpstreamPrime() (uint64, error) {

	NextPrime := n.NextPrime
	NthRoot := n.NthRoot
	CheckNextPrime := n.CheckNextPrime
	Size := n.Size

	for {
		if CheckNextPrime {

			// Stops if the next prime would overlap with primes of the next bit-size or if an uint64 overflow would occur.
			if math.Log2(float64(NextPrime))-Size >= 0.5 {

				n.CheckNextPrime = false

				return 0, fmt.Errorf("cannot NextUpstreamPrime: prime list for upstream primes is exhausted (overlap with next bit-size or prime > 2^{64})")

			} else {

				if IsPrime(NextPrime) {

					n.NextPrime = NextPrime + NthRoot

					n.CheckNextPrime = CheckNextPrime

					return NextPrime, nil
				}

				NextPrime += NthRoot
			}
		}
	}
}

// NextDownstreamPrime returns the next prime of the form 2^{BitSize} - k * {NthRoot} + 1.
func (n *NTTFriendlyPrimesGenerator) NextDownstreamPrime() (uint64, error) {

	PrevPrime := n.PrevPrime
	NthRoot := n.NthRoot
	CheckPrevPrime := n.CheckPrevPrime
	Size := n.Size

	for {

		if CheckPrevPrime {

			// Stops if the next prime would overlap with the primes of the previous bit-size or if an uint64 overflow would occur.
			if Size-math.Log2(float64(PrevPrime)) >= 0.5 || PrevPrime < NthRoot {

				n.CheckPrevPrime = false

				return 0, fmt.Errorf("cannot NextDownstreamPrime: prime list for downstream primes is exhausted (overlap with previous bit-size or prime < NthRoot")

			} else {

				if IsPrime(PrevPrime) {

					n.PrevPrime = PrevPrime - NthRoot

					n.CheckPrevPrime = CheckPrevPrime

					return PrevPrime, nil
				}

				PrevPrime -= NthRoot
			}
		}
	}
}

// NextAlternatingPrime returns the next prime of the form 2^{BitSize} +/- k * {NthRoot} + 1.
func (n *NTTFriendlyPrimesGenerator) NextAlternatingPrime() (uint64, error) {

	NextPrime := n.NextPrime
	PrevPrime := n.PrevPrime

	NthRoot := n.NthRoot

	CheckNextPrime := n.CheckNextPrime
	CheckPrevPrime := n.CheckPrevPrime

	Size := n.Size

	for {

		if !(CheckNextPrime || CheckPrevPrime) {
			return 0, fmt.Errorf("cannot NextAlternatingPrime: prime list for both upstream and downstream primes is exhausted (overlap with previous/next bit-size or NthRoot > prime > 2^{64} ")
		}

		if CheckNextPrime {

			// Stops if the next prime would overlap with primes of the next bit-size or if an uint64 overflow would occure.
			if math.Log2(float64(NextPrime))-Size >= 0.5 || NextPrime > 0xffffffffffffffff-NthRoot {

				CheckNextPrime = false

			} else {

				if IsPrime(NextPrime) {

					n.NextPrime = NextPrime + NthRoot
					n.PrevPrime = PrevPrime

					n.CheckNextPrime = CheckNextPrime
					n.CheckPrevPrime = CheckPrevPrime

					return NextPrime, nil
				}

				NextPrime += NthRoot
			}
		}

		if CheckPrevPrime {

			// Stops if the next prime would overlap with the primes of the previous bit-size or if an uint64 overflow would occure.
			if Size-math.Log2(float64(PrevPrime)) >= 0.5 || PrevPrime < NthRoot {

				CheckPrevPrime = false

			} else {

				if IsPrime(PrevPrime) {

					n.NextPrime = NextPrime
					n.PrevPrime = PrevPrime - NthRoot

					n.CheckNextPrime = CheckNextPrime
					n.CheckPrevPrime = CheckPrevPrime

					return PrevPrime, nil
				}

				PrevPrime -= NthRoot
			}
		}
	}
}
