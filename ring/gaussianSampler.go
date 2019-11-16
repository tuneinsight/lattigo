package ring

import (
	"crypto/rand"
	"math"
)

type GaussiamSampler struct {
	context *Context
	sigma   float64
	bound   uint64
}

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

func (context *Context) SampleGaussianNew(sigma float64, bound uint64) (pol *Poly) {
	pol = context.NewPoly()
	context.SampleGaussian(pol, sigma, bound)
	return
}

func (context *Context) SampleGaussianNTT(pol *Poly, sigma float64, bound uint64) {
	context.SampleGaussian(pol, sigma, bound)
	context.NTT(pol, pol)
}

func (context *Context) SampleGaussianNTTNew(sigma float64, bound uint64) (pol *Poly) {
	pol = context.SampleGaussianNew(sigma, bound)
	context.NTT(pol, pol)
	return
}

// KYSampler is the structure holding the parameters for the gaussian sampling.
type KYSampler struct {
	context *Context
	sigma   float64
	bound   int
	Matrix  [][]uint8
}

// NewKYSampler creates a new KYSampler with sigma and bound that will be used to sample polynomial within the provided discret gaussian distribution.
func (context *Context) NewKYSampler(sigma float64, bound int) *KYSampler {
	kysampler := new(KYSampler)
	kysampler.context = context
	kysampler.sigma = sigma
	kysampler.bound = bound
	kysampler.Matrix = computeMatrix(sigma, bound)
	return kysampler
}

//gaussian computes (1/variange*sqrt(pi)) * exp((x^2) / (2*variance^2)),  2.50662827463100050241576528481104525300698674060993831662992357 = sqrt(2*pi)
func gaussian(x, sigma float64) float64 {
	return (1 / (sigma * 2.5066282746310007)) * math.Exp(-((math.Pow(x, 2)) / (2 * sigma * sigma)))
}

// computeMatrix computes the binary expension with precision x in bits of the normal distribution
// with sigma and bound. Returns a matrix of the form M = [[0,1,0,0,...],[0,0,1,,0,...]],
// where each row is the binary expension of the normal distribution of index(row) with sigma and bound (center=0).
func computeMatrix(sigma float64, bound int) [][]uint8 {
	var g float64
	var x uint64

	precision := uint64(56)

	M := make([][]uint8, bound)

	breakCounter := 0

	for i := 0; i < bound; i++ {

		g = gaussian(float64(i), sigma)

		if i == 0 {
			g *= math.Exp2(float64(precision) - 1)
		} else {
			g *= math.Exp2(float64(precision))
		}

		x = uint64(g)

		if x == 0 {
			break
		}

		M[i] = make([]uint8, precision-1)

		for j := uint64(0); j < precision-1; j++ {
			M[i][j] = uint8((x >> (precision - j - 2)) & 1)
		}

		breakCounter += 1
	}

	M = M[:breakCounter]

	return M
}

func kysampling(M [][]uint8, randomBytes []byte, pointer uint8) (uint64, uint64, []byte, uint8) {

	var sign uint8

	d := 0
	col := 0
	colLen := len(M)

	for {

		// Uses one random byte per cycle and cycle through the randombytes
		for i := pointer; i < 8; i++ {

			d = (d << 1) + 1 - int((uint8(randomBytes[0])>>i)&1)

			// There is small probability that it will get out of the bound, then
			// rerun until it gets a proper output
			if d > colLen-1 {
				return kysampling(M, randomBytes, i)
			}

			for row := colLen - 1; row >= 0; row-- {

				d -= int(M[row][col])

				if d == -1 {

					// Sign
					if i == 7 {
						pointer = 0
						// If the last bit of the byte was read, samples a new byte for the sign
						randomBytes = randomBytes[1:]

						if len(randomBytes) == 0 {
							randomBytes = make([]byte, 8)
							if _, err := rand.Read(randomBytes); err != nil {
								panic("crypto rand error")
							}
						}

						sign = uint8(randomBytes[0]) & 1

					} else {
						pointer = i
						// Else the sign is the next bit of the byte
						sign = uint8(randomBytes[0]>>(i+1)) & 1
					}

					return uint64(row), uint64(sign), randomBytes, pointer + 1
				}
			}

			col += 1
		}

		// Resets the bit pointer and discards the used byte
		pointer = 0
		randomBytes = randomBytes[1:]

		// Sample 8 new bytes if the last byte was discarded
		if len(randomBytes) == 0 {
			randomBytes = make([]byte, 8)
			if _, err := rand.Read(randomBytes); err != nil {
				panic("crypto rand error")
			}
		}

	}
}

// SampleNew samples a new polynomial with gaussian distribution given the target kys parameters.
func (kys *KYSampler) SampleNew() *Poly {
	Pol := kys.context.NewPoly()
	kys.Sample(Pol)
	return Pol
}

// SampleNew samples on the target polynomial coefficients with gaussian distribution given the target kys parameters.
func (kys *KYSampler) Sample(Pol *Poly) {

	var coeff uint64
	var sign uint64

	randomBytes := make([]byte, 8)
	pointer := uint8(0)

	if _, err := rand.Read(randomBytes); err != nil {
		panic("crypto rand error")
	}

	for i := uint64(0); i < kys.context.N; i++ {

		coeff, sign, randomBytes, pointer = kysampling(kys.Matrix, randomBytes, pointer)

		for j, qi := range kys.context.Modulus {
			Pol.Coeffs[j][i] = (coeff & (sign * 0xFFFFFFFFFFFFFFFF)) | ((qi - coeff) & ((sign ^ 1) * 0xFFFFFFFFFFFFFFFF))

		}
	}
}

// SampleNew samples on the target polynomial coefficients with gaussian distribution given the target kys parameters.
func (kys *KYSampler) SampleAndAddLvl(level uint64, Pol *Poly) {

	var coeff uint64
	var sign uint64

	randomBytes := make([]byte, 8)
	pointer := uint8(0)

	if _, err := rand.Read(randomBytes); err != nil {
		panic("crypto rand error")
	}

	for i := uint64(0); i < kys.context.N; i++ {

		coeff, sign, randomBytes, pointer = kysampling(kys.Matrix, randomBytes, pointer)

		for j := uint64(0); j < level+1; j++ {
			Pol.Coeffs[j][i] = CRed(Pol.Coeffs[j][i]+((coeff*sign)|(kys.context.Modulus[j]-coeff)*(sign^1)), kys.context.Modulus[j])
		}
	}
}

func (kys *KYSampler) SampleAndAdd(Pol *Poly) {
	kys.SampleAndAddLvl(uint64(len(kys.context.Modulus))-1, Pol)
}

// SampleNTTNew samples a polynomial with gaussian distribution given the target kys context and apply the NTT.
func (kys *KYSampler) SampleNTTNew() *Poly {
	Pol := kys.SampleNew()
	kys.context.NTT(Pol, Pol)
	return Pol
}

// SampleNew samples on the target polynomial coefficients with gaussian distribution given the target kys parameters,and applies the NTT.
func (kys *KYSampler) SampleNTT(Pol *Poly) {
	kys.Sample(Pol)
	kys.context.NTT(Pol, Pol)
}
