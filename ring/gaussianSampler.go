package ring

import (
	"crypto/rand"
)

// SampleGaussian samples a truncated gaussian polynomial with variance
// sigma within the given bound using the Ziggurat algorithm.
func (context *Context) SampleGaussian(pol *Poly, sigma float64, bound uint64) {

	var coeffFlo float64
	var coeffInt uint64
	var sign uint64

	randomBytes := make([]byte, 1024)

	if _, err := rand.Read(randomBytes); err != nil {
		panic("crypto rand error")
	}

	for i := uint64(0); i < context.N; i++ {

		for {
			coeffFlo, sign, randomBytes = normFloat64(randomBytes)

			if coeffInt = uint64(coeffFlo * sigma); coeffInt <= bound {
				break
			}
		}

		for j := 0; j < len(pol.Coeffs); j++ {
			pol.Coeffs[j][i] = (coeffInt * sign) | (context.Modulus[j]-coeffInt)*(sign^1)
		}
	}
}

// SampleGaussianAndAdd adds on the input polynomial a truncated gaussian polynomial
// with variance sigma within the given bound using the Ziggurat algorithm.
func (context *Context) SampleGaussianAndAdd(pol *Poly, sigma float64, bound uint64) {

	var coeffFlo float64
	var coeffInt uint64
	var sign uint64

	randomBytes := make([]byte, 1024)

	if _, err := rand.Read(randomBytes); err != nil {
		panic("crypto rand error")
	}

	for i := uint64(0); i < context.N; i++ {

		for {
			coeffFlo, sign, randomBytes = normFloat64(randomBytes)

			if coeffInt = uint64(coeffFlo * sigma); coeffInt <= bound {
				break
			}
		}

		for j := 0; j < len(pol.Coeffs); j++ {
			pol.Coeffs[j][i] = CRed(pol.Coeffs[j][i]+((coeffInt*sign)|(context.Modulus[j]-coeffInt)*(sign^1)), context.Modulus[j])
		}
	}
}

// SampleGaussianNew samples a new truncated gaussian polynomial with
// variance sigma within the given bound using the Ziggurat algorithm.
func (context *Context) SampleGaussianNew(sigma float64, bound uint64) (pol *Poly) {
	pol = context.NewPoly()
	context.SampleGaussian(pol, sigma, bound)
	return
}

// SampleGaussianNTT samples a trucated gaussian polynomial in the NTT domain
// with variance sigma within the given bound using the Ziggurat algorithm.
func (context *Context) SampleGaussianNTT(pol *Poly, sigma float64, bound uint64) {
	context.SampleGaussian(pol, sigma, bound)
	context.NTT(pol, pol)
}

// SampleGaussianNTTNew samples a new trucated gaussian polynomial in the NTT domain
// with variance sigma within the given bound using the Ziggurat algorithm
func (context *Context) SampleGaussianNTTNew(sigma float64, bound uint64) (pol *Poly) {
	pol = context.SampleGaussianNew(sigma, bound)
	context.NTT(pol, pol)
	return
}

// Kept for backwards compatibility
type Sampler struct {
	context *Context
	sigma   float64
	bound   uint64
}

// Kept for backwards compatibility
func (context *Context) NewSampler(sigma float64, bound uint64) *Sampler {
	sampler := new(Sampler)
	sampler.context = context
	sampler.sigma = sigma
	sampler.bound = bound
	return sampler
}

// Kept for backwards compatibility
func (sampler *Sampler) SampleNew() *Poly {
	Pol := sampler.context.NewPoly()
	sampler.Sample(Pol)
	return Pol
}

// Kept for backwards compatibility
func (sampler *Sampler) Sample(Pol *Poly) {
	sampler.context.SampleGaussian(Pol, sampler.sigma, sampler.bound)
	return
}

// Kept for backwards compatibility
func (sampler *Sampler) SampleAndAdd(Pol *Poly) {
	sampler.context.SampleGaussianAndAdd(Pol, sampler.sigma, sampler.bound)
}

// Kept for backwards compatibility
func (sampler *Sampler) SampleNTTNew() *Poly {
	Pol := sampler.SampleNew()
	sampler.context.NTT(Pol, Pol)
	return Pol
}

// Kept for backwards compatibility
func (sampler *Sampler) SampleNTT(Pol *Poly) {
	sampler.Sample(Pol)
	sampler.context.NTT(Pol, Pol)
}

// kysampling use the binary expension and random bytes matrix to sample a discret gaussian value and its sign.
// Kept For ternarySampler functionality
func kysampling(M [][]uint8, randomBytes []byte, pointer uint8, bytePointer uint64, byteLength uint64) (uint64, uint64, []byte, uint8, uint64) {

	var sign uint8

	d := 0
	col := 0
	colLen := len(M)

	for {

		// Uses one random byte per cycle and cycle through the randombytes
		for i := pointer; i < 8; i++ {

			d = (d << 1) + 1 - int((uint8(randomBytes[bytePointer])>>i)&1)

			// There is small probability that it will get out of the bound, then
			// rerun until it gets a proper output
			if d > colLen-1 {
				return kysampling(M, randomBytes, i, bytePointer, byteLength)
			}

			for row := colLen - 1; row >= 0; row-- {

				d -= int(M[row][col])

				if d == -1 {

					// Sign
					if i == 7 {
						pointer = 0
						bytePointer++

						// If the last bit of the byte was read, samples a new byte for the sign
						if bytePointer >= byteLength {
							bytePointer = 0
							if _, err := rand.Read(randomBytes); err != nil {
								panic("crypto rand error")
							}
						}

						sign = uint8(randomBytes[bytePointer]) & 1

					} else {
						pointer = i
						// Else the sign is the next bit of the byte
						sign = uint8(randomBytes[bytePointer]>>(i+1)) & 1
					}

					return uint64(row), uint64(sign), randomBytes, pointer + 1, bytePointer
				}
			}

			col++
		}

		// Resets the bit pointer and discards the used byte
		pointer = 0
		bytePointer++

		// Sample new bytes if the last byte was discarded
		if bytePointer >= byteLength {
			bytePointer = 0
			if _, err := rand.Read(randomBytes); err != nil {
				panic("crypto rand error")
			}
		}

	}
}
