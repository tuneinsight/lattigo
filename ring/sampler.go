package ring

import (
	"crypto/rand"
	"encoding/binary"
	"math"
	"math/bits"
)

// This code is to sample a value from discrete gaussian distribution.
// All the algorithms originate from https://eprint.iacr.org/2013/383.pdf

// GaussSampling returns a value sampled from discrete gaussian distribution,
// originates from Algorithm 11 & 12 in the paper.
func GaussSampling(sigma float64) int64 {
	const sigma2 = 0.8493218                // sigma2 = sqrt(1/(2*ln2)) -- page 28
	k := uint64(math.Round(sigma / sigma2)) // sigma = k * sigma2 -- page 29
	var x, y, z uint64
	var b1, b2 bool
	for {
		x = binaryGauss()
		y = randUniform(k)
		b1 = bernoulliExp(y*(y+2*k*x), 2*sigma*sigma)
		if b1 {
			z = k*x + y
			b2 = bernoulli(0.5)
			if z != 0 || b2 {
				if b2 {
					return int64(z)
				} else {
					return -int64(z)
				}
			}
		}
	}
}

// binaryGauss returns a random uint value drawn from binary gaussian distribution,
// originates from Algorithm 10 in the paper.
func binaryGauss() uint64 {
	if bernoulli(0.5) == false {
		return 0
	}
	// i < 16 represents infinite
	for i := 1; i < 16; i++ {
		randomBits := randInt64(uint64(2*i - 1))
		if randomBits != 0 {
			return binaryGauss()
		}
		if randomBits == 0 {
			return uint64(i)
		}
	}
	return 0
}

// bernoulli returns a random bool value drawn from exponential bernoulli distribution
// originates from Algorithm 8 in the paper.
func bernoulliExp(x uint64, f float64) bool {
	xBinary, xBitlen := NewInt(int64(x)).Bits()
	if xBitlen == 0 {
		return true
	}
	for i := xBitlen; i > 0; i-- {
		if xBinary[i-1] == 1 {
			c := math.Exp(-math.Exp2(float64(i)) / f)
			if bernoulli(c) == false {
				return false
			}
		}
	}
	return true
}

// bernoulli returns a random bool value drawn from bernoulli distribution
func bernoulli(p float64) bool {
	pInt := uint64(p * (1 << 31))
	randomInt := randInt32(0xFFFFFFFF)
	if randomInt < pInt {
		return true
	}
	return false
}

type KYSampler struct {
	context *Context
	sigma   float64
	bound   int
	Matrix  [][]uint8
}

func (context *Context) NewKYSampler(sigma float64, bound int) *KYSampler {
	kysampler := new(KYSampler)
	kysampler.context = context
	kysampler.sigma = sigma
	kysampler.bound = bound
	kysampler.Matrix = computeMatrix(sigma, bound)
	return kysampler
}

//2.50662827463100050241576528481104525300698674060993831662992357 = sqrt(2*pi)
func gaussian(x, sigma float64) float64 {
	return (1 / (sigma * 2.5066282746310007)) * math.Exp(-((math.Pow(x, 2)) / (2 * sigma * sigma)))
}

// Computes the binary expension with precision x in bits of the normal distribution
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

func (kys *KYSampler) SampleNew() *Poly {
	Pol := kys.context.NewPoly()
	kys.Sample(Pol)
	return Pol
}

func (kys *KYSampler) Sample(Pol *Poly) {

	var coeff uint64
	var sign uint64

	randomBytes := make([]byte, 8)
	pointer := uint8(0)

	if _, err := rand.Read(randomBytes); err != nil {
		panic("crypto rand error")
	}

	for i := uint64(0); i < kys.context.N; i++ {

		coeff, sign, randomBytes, pointer = KYSampling(kys.Matrix, randomBytes, pointer)

		for j, qi := range kys.context.Modulus {

			Pol.Coeffs[j][i] = (coeff & (sign * 0xFFFFFFFFFFFFFFFF)) | ((qi - coeff) & ((sign ^ 1) * 0xFFFFFFFFFFFFFFFF))

		}
	}
}

func (kys *KYSampler) SampleNTTNew() *Poly {
	Pol := kys.SampleNew()
	kys.context.NTT(Pol, Pol)
	return Pol
}

func (kys *KYSampler) SampleNTT(Pol *Poly) {
	kys.Sample(Pol)
	kys.context.NTT(Pol, Pol)
}

func KYSampling(M [][]uint8, randomBytes []byte, pointer uint8) (uint64, uint64, []byte, uint8) {

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
				return KYSampling(M, randomBytes, i)
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

type TernarySampler struct {
	context          *Context
	Matrix           [][]uint64
	MatrixMontgomery [][]uint64
}

func (context *Context) NewTernarySampler() *TernarySampler {

	sampler := new(TernarySampler)
	sampler.context = context

	sampler.Matrix = make([][]uint64, len(context.Modulus))
	sampler.MatrixMontgomery = make([][]uint64, len(context.Modulus))

	for i, Qi := range context.Modulus {

		sampler.Matrix[i] = make([]uint64, 3)
		sampler.Matrix[i][0] = Qi - 1
		sampler.Matrix[i][1] = 0
		sampler.Matrix[i][2] = 1

		sampler.MatrixMontgomery[i] = make([]uint64, 3)
		sampler.MatrixMontgomery[i][0] = MForm(Qi-1, Qi, context.bredParams[i])
		sampler.MatrixMontgomery[i][1] = 0
		sampler.MatrixMontgomery[i][2] = MForm(1, Qi, context.bredParams[i])
	}

	return sampler
}

func (sampler *TernarySampler) SampleNew() (pol *Poly) {
	pol = sampler.context.NewPoly()
	sampler.Sample(pol)
	return pol
}

func (sampler *TernarySampler) Sample(pol *Poly) {
	var coeff uint64
	for j := uint64(0); j < sampler.context.N; j++ {
		coeff = randUint3()
		for i := range sampler.context.Modulus {
			pol.Coeffs[i][j] = sampler.Matrix[i][coeff]
		}
	}
}

func (sampler *TernarySampler) SampleMontgomeryNew() (pol *Poly) {
	pol = sampler.context.NewPoly()
	sampler.SampleMontgomery(pol)
	return pol
}

func (sampler *TernarySampler) SampleMontgomery(pol *Poly) {
	var coeff uint64
	for j := uint64(0); j < sampler.context.N; j++ {
		coeff = randUint3()
		for i := range sampler.context.Modulus {
			pol.Coeffs[i][j] = sampler.MatrixMontgomery[i][coeff]
		}
	}
}

func (sampler *TernarySampler) SampleNTTNew() (pol *Poly) {
	pol = sampler.context.NewPoly()
	sampler.Sample(pol)
	sampler.context.NTT(pol, pol)
	return pol
}

func (sampler *TernarySampler) SampleNTT(pol *Poly) {
	sampler.Sample(pol)
	sampler.context.NTT(pol, pol)
}

func (sampler *TernarySampler) SampleMontgomeryNTTNew() (pol *Poly) {
	pol = sampler.SampleMontgomeryNew()
	sampler.context.NTT(pol, pol)
	return pol
}

func (sampler *TernarySampler) SampleMontgomeryNTT(pol *Poly) {
	sampler.SampleMontgomery(pol)
	sampler.context.NTT(pol, pol)
}

// TODO : precompute the mask in the ring context
// randUniform returns a uniformly distributed value in {0,v)
func randUniform(v uint64) uint64 {
	var length, mask, randomInt uint64

	length = uint64(bits.Len64(v))

	mask = (1 << length) - 1

	for {
		randomInt = randInt64(mask)
		if randomInt < v {
			return randomInt
		}
	}
}

func RandUniform(v uint64) uint64 {
	var length, mask, randomInt uint64

	length = uint64(bits.Len64(v))

	mask = (1 << length) - 1

	for {
		randomInt = randInt64(mask)
		if randomInt < v {
			return randomInt
		}
	}
}

// Samples an bit and a sign with rejection sampling (25% chance of failure)
// with probabilities :
// Pr[int = 0 : 1/3 ; int = 1 : 2/3]
// Pr[sign = 1 : 1/2; sign = 0 : 1/2]
func randInt3() (uint64, uint64) {
	var randomInt uint64

	for {
		randomInt = randInt8()
		if (randomInt & 3) < 3 {
			// (a|b) is 1/3 = 0 and 2/3 = 1 if (a<<1 | b) in [0,1,2]
			return ((randomInt >> 1) | randomInt) & 1, (randomInt >> 2) & 1
		}
	}
}

func randUint3() uint64 {
	var randomInt uint64
	for {
		randomInt = randInt8() & 3
		if randomInt < 3 {
			return randomInt
		}
	}
}

// TODO : precompute the mask in the ring context
// randInt generates a random uint32 value of given length
func randInt8() uint64 {

	// generate random 4 bytes
	randomBytes := make([]byte, 1)
	if _, err := rand.Read(randomBytes); err != nil {
		panic("crypto rand error")
	}
	// return required bits
	return uint64(randomBytes[0])
}

// TODO : precompute the mask in the ring context
// randInt generates a random uint32 value of given length
func randInt32(mask uint64) uint64 {

	// generate random 4 bytes
	randomBytes := make([]byte, 4)
	_, err := rand.Read(randomBytes)
	if err != nil {
		panic("crypto rand error")
	}

	// convert 4 bytes to a uint32
	randomUint32 := uint64(binary.BigEndian.Uint32(randomBytes))

	// return required bits
	return mask & randomUint32
}

// TODO : precompute the mask in the ring context
// randInt generates a random uint64 value of given length
func randInt64(mask uint64) uint64 {

	// generate random 8 bytes
	randomBytes := make([]byte, 8)
	_, err := rand.Read(randomBytes)
	if err != nil {
		panic("crypto rand error")
	}

	// convert 8 bytes to a uint64
	randomUint64 := binary.BigEndian.Uint64(randomBytes)

	// return required bits
	return mask & randomUint64
}
