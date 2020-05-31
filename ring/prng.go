package ring

import (
	"encoding/binary"
	"github.com/ldsec/lattigo/utils"
	"math"
	"math/bits"
)

// CRPGenerator is the structure storing the parameters for deterministicaly securely
// generating random polynomials using the structure PRNG.
type CRPGenerator struct {
	prng    utils.PRNG
	context *Context
	sum     []byte
	masks   []uint64
}

// NewCRPGenerator creates a new CRPGenerator, that will deterministically and securely generate uniform polynomials
// in the domain of the input context using the hash function blake2b. The PRNG can be instantiated with a key, if no
// key is used, set key=nil.
func NewCRPGenerator(key []byte, context *Context) *CRPGenerator {
	var err error
	crpgenerator := new(CRPGenerator)

	if crpgenerator.prng, err = utils.NewKeyedPRNG(key); err != nil {
		panic(err)
	}

	crpgenerator.context = context
	crpgenerator.masks = make([]uint64, len(context.Modulus))

	for i, qi := range context.Modulus {
		crpgenerator.masks[i] = (1 << uint64(bits.Len64(qi))) - 1
	}

	crpgenerator.sum = make([]byte, context.N)

	return crpgenerator
}

// GetClock returns the current clock of the CRPGenerator.
func (crpgenerator *CRPGenerator) GetClock() uint64 {
	return crpgenerator.prng.GetClock()
}

// SetClock sets the clock of the CRPGenerator to the given input by clocking it until the
// clock cycle reaches the desired number. If the given input is smaller than the current clock,
// it will panic.
func (crpgenerator *CRPGenerator) SetClock(n uint64) {
	if err := crpgenerator.prng.SetClock(crpgenerator.sum, n); err != nil {
		panic(err)
	}
}

// ClockNew generates and returns a new uniform polynomial and increases the clock cycle by 1.
func (crpgenerator *CRPGenerator) ClockUniformNew() (crp *Poly) {
	crp = crpgenerator.context.NewPoly()
	crpgenerator.ClockUniform(crp)
	return
}

// Clock samples a random polynomial on the probided polynomial and increases the clock cycle by 1.
func (crpgenerator *CRPGenerator) ClockUniform(crp *Poly) {
	var maxBytes = crpgenerator.context.N
	var coeff, ptr uint64

	// Starts with 32 random bytes from the prng
	crpgenerator.prng.Clock(crpgenerator.sum)

	for i := uint64(0); i < crpgenerator.context.N; i++ {

		for j, qi := range crpgenerator.context.Modulus {

			for {
				// If not enough randomBytes, reads 64 more random bytes from the PRNG
				if ptr == maxBytes {
					crpgenerator.prng.Clock(crpgenerator.sum)
					ptr = 0
				}

				// Converts 4 randomBytes into an uint64 of at most 60 bits (maximum size of the modulus)
				coeff = binary.BigEndian.Uint64(crpgenerator.sum[ptr:ptr+8]) & crpgenerator.masks[j]
				ptr += 8

				// If coeff is smaller than qi, breaks
				if coeff < qi {
					break
				}
			}

			// Assigns the coeff to the polynomial
			crp.Coeffs[j][i] = coeff
		}
	}
}

func (crpgenerator *CRPGenerator) ClockGaussian(crp *Poly, sigma float64, bound uint64) {

	var coeffFlo float64
	var coeffInt uint64
	var sign uint64
	var ptr uint64

	context := crpgenerator.context

	crpgenerator.prng.Clock(crpgenerator.sum)

	for i := uint64(0); i < context.N; i++ {

		for {

			coeffFlo, sign, ptr = crpgenerator.normFloat64(ptr)

			if coeffInt = uint64(coeffFlo * sigma); coeffInt <= bound {
				break
			}
		}

		for j, qi := range context.Modulus {
			crp.Coeffs[j][i] = (coeffInt * sign) | (qi-coeffInt)*(sign^1)
		}
	}
}

func (crpgenerator *CRPGenerator) ClockGaussianAndAdd(crp *Poly, sigma float64, bound uint64) {

	var coeffFlo float64
	var coeffInt uint64
	var sign uint64
	var ptr uint64

	context := crpgenerator.context

	crpgenerator.prng.Clock(crpgenerator.sum)

	for i := uint64(0); i < context.N; i++ {

		for {

			coeffFlo, sign, ptr = crpgenerator.normFloat64(ptr)

			if coeffInt = uint64(coeffFlo * sigma); coeffInt <= bound {
				break
			}
		}

		for j, qi := range context.Modulus {
			crp.Coeffs[j][i] = CRed(crp.Coeffs[j][i]+((coeffInt*sign)|(qi-coeffInt)*(sign^1)), qi)
		}
	}
}

// NormFloat64 returns a normally distributed float64 in
// the range -math.MaxFloat64 through +math.MaxFloat64 inclusive,
// with standard normal distribution (mean = 0, stddev = 1).
// To produce a different normal distribution, callers can
// adjust the output using:
//
//  sample = NormFloat64() * desiredStdDev + desiredMean
// Algorithm adapted from https://golang.org/src/math/rand/normal.go
func (crpgenerator *CRPGenerator) normFloat64(ptr uint64) (float64, uint64, uint64) {

	var maxBytes = crpgenerator.context.N

	for {

		if ptr == maxBytes {
			crpgenerator.prng.Clock(crpgenerator.sum)
			ptr = 0
		}

		juint32 := binary.BigEndian.Uint32(crpgenerator.sum[ptr : ptr+4])
		ptr += 8

		j := int32(juint32 & 0x7fffffff)
		sign := uint64(juint32 >> 31)

		i := j & 0x7F

		x := float64(j) * float64(wn[i])

		// 1
		if uint32(j) < kn[i] {

			// This case should be hit better than 99% of the time.
			return x, sign, ptr
		}

		// 2
		if i == 0 {

			// This extra work is only required for the base strip.
			for {

				if ptr == maxBytes {
					crpgenerator.prng.Clock(crpgenerator.sum)
					ptr = 0
				}

				x = -math.Log(randFloat64(crpgenerator.sum[ptr:ptr+8])) * (1.0 / 3.442619855899)
				ptr += 8

				if ptr == maxBytes {
					crpgenerator.prng.Clock(crpgenerator.sum)
					ptr = 0
				}

				y := -math.Log(randFloat64(crpgenerator.sum[ptr : ptr+8]))
				ptr += 8

				if y+y >= x*x {
					break
				}
			}

			return x + 3.442619855899, sign, ptr
		}

		if ptr == maxBytes {
			crpgenerator.prng.Clock(crpgenerator.sum)
			ptr = 0
		}

		// 3
		if fn[i]+float32(randFloat64(crpgenerator.sum[ptr:ptr+8]))*(fn[i-1]-fn[i]) < float32(math.Exp(-0.5*x*x)) {
			return x, sign, ptr + 8
		}
		ptr += 8
	}
}
