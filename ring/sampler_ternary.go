package ring

import (
	"math"
	"math/bits"

	"github.com/ldsec/lattigo/utils"
)

const precision = uint64(56)

type baseSampler struct {
	prng     utils.PRNG
	baseRing *Ring
}

// TernarySampler is the state of a polynomial sampler in the ternary distribution.
type TernarySampler struct {
	baseSampler
	matrixProba  [2][precision - 1]uint8
	matrixValues [][3]uint64
	p            float64
	hw           uint64
	sample       func(poly *Poly)
}

// NewTernarySampler creates a new instance of TernarySampler from a PRNG, the ring definition and the distribution
// parameters: p is the probability of a coefficient being 0, (1-p)/2 is the probability of 1 and -1. If montgomery
// is set to true, polynomials read from this sampler are in Montgomery form.
func NewTernarySampler(prng utils.PRNG, baseRing *Ring, p float64, montgomery bool) *TernarySampler {
	ternarySampler := new(TernarySampler)
	ternarySampler.baseRing = baseRing
	ternarySampler.prng = prng
	ternarySampler.p = p
	ternarySampler.sample = ternarySampler.sampleProba

	ternarySampler.initialiseMatrix(montgomery)

	if p != 0.5 {
		ternarySampler.computeMatrixTernary(p)
	}

	return ternarySampler
}

// NewTernarySampler creates a new instance of TernarySampler from a PRNG, the ring definition and the desired
// hamming weight for the output polynomials. If montgomery is set to true, polynomials read from this sampler
// are in Montgomery form.
func NewTernarySamplerSparse(prng utils.PRNG, baseRing *Ring, hw uint64, montgomery bool) *TernarySampler {
	ternarySampler := new(TernarySampler)
	ternarySampler.baseRing = baseRing
	ternarySampler.prng = prng
	ternarySampler.hw = hw
	ternarySampler.sample = ternarySampler.sampleSparse

	ternarySampler.initialiseMatrix(montgomery)

	return ternarySampler
}

// Read samples a polynomial into pol.
func (ts *TernarySampler) Read(pol *Poly) {
	ts.sample(pol)
}

// ReadNew allocates and samples a polynomial.
func (ts *TernarySampler) ReadNew() (pol *Poly) {
	pol = ts.baseRing.NewPoly()
	ts.sample(pol)
	return pol
}

func (ternarySampler *TernarySampler) initialiseMatrix(montgomery bool) {
	ternarySampler.matrixValues = make([][3]uint64, len(ternarySampler.baseRing.Modulus))

	for i, Qi := range ternarySampler.baseRing.Modulus {

		ternarySampler.matrixValues[i][0] = 0

		if montgomery {
			ternarySampler.matrixValues[i][1] = MForm(1, Qi, ternarySampler.baseRing.BredParams[i])
			ternarySampler.matrixValues[i][2] = MForm(Qi-1, Qi, ternarySampler.baseRing.BredParams[i])
		} else {
			ternarySampler.matrixValues[i][1] = 1
			ternarySampler.matrixValues[i][2] = Qi - 1
		}
	}
}

func (ts *TernarySampler) computeMatrixTernary(p float64) {
	var g float64
	var x uint64

	g = p
	g *= math.Exp2(float64(precision))
	x = uint64(g)

	for j := uint64(0); j < precision-1; j++ {
		ts.matrixProba[0][j] = uint8((x >> (precision - j - 1)) & 1)
	}

	g = 1 - p
	g *= math.Exp2(float64(precision))
	x = uint64(g)

	for j := uint64(0); j < precision-1; j++ {
		ts.matrixProba[1][j] = uint8((x >> (precision - j - 1)) & 1)
	}

}

func (ternarySampler *TernarySampler) sampleProba(pol *Poly) {

	if ternarySampler.p == 0 {
		panic("cannot sample -> p = 0")
	}

	var coeff uint64
	var sign uint64
	var index uint64

	if ternarySampler.p == 0.5 {

		randomBytesCoeffs := make([]byte, ternarySampler.baseRing.N>>3)
		randomBytesSign := make([]byte, ternarySampler.baseRing.N>>3)

		ternarySampler.prng.Clock(randomBytesCoeffs)

		ternarySampler.prng.Clock(randomBytesSign)

		for i := uint64(0); i < ternarySampler.baseRing.N; i++ {
			coeff = uint64(uint8(randomBytesCoeffs[i>>3])>>(i&7)) & 1
			sign = uint64(uint8(randomBytesSign[i>>3])>>(i&7)) & 1

			index = (coeff & (sign ^ 1)) | ((sign & coeff) << 1)

			for j := range ternarySampler.baseRing.Modulus {
				pol.Coeffs[j][i] = ternarySampler.matrixValues[j][index] //(coeff & (sign^1)) | (qi - 1) * (sign & coeff)
			}
		}

	} else {

		randomBytes := make([]byte, ternarySampler.baseRing.N)

		pointer := uint8(0)
		bytePointer := uint64(0)

		ternarySampler.prng.Clock(randomBytes)

		for i := uint64(0); i < ternarySampler.baseRing.N; i++ {

			coeff, sign, randomBytes, pointer, bytePointer = ternarySampler.kysampling(ternarySampler.prng, randomBytes, pointer, bytePointer, ternarySampler.baseRing.N)

			index = (coeff & (sign ^ 1)) | ((sign & coeff) << 1)

			for j := range ternarySampler.baseRing.Modulus {
				pol.Coeffs[j][i] = ternarySampler.matrixValues[j][index] //(coeff & (sign^1)) | (qi - 1) * (sign & coeff)
			}
		}
	}
}

func (ternarySampler *TernarySampler) sampleSparse(pol *Poly) {

	if ternarySampler.hw > ternarySampler.baseRing.N {
		ternarySampler.hw = ternarySampler.baseRing.N
	}

	var mask, j uint64
	var coeff uint8

	index := make([]uint64, ternarySampler.baseRing.N)
	for i := uint64(0); i < ternarySampler.baseRing.N; i++ {
		index[i] = i
	}

	randomBytes := make([]byte, (uint64(math.Ceil(float64(ternarySampler.hw) / 8.0)))) // We sample ceil(hw/8) bytes
	pointer := uint8(0)

	ternarySampler.prng.Clock(randomBytes)

	for i := uint64(0); i < ternarySampler.hw; i++ {
		mask = (1 << uint64(bits.Len64(ternarySampler.baseRing.N-i))) - 1 // rejection sampling of a random variable between [0, len(index)]

		j = randInt32(ternarySampler.prng, mask)
		for j >= ternarySampler.baseRing.N-i {
			j = randInt32(ternarySampler.prng, mask)
		}

		coeff = (uint8(randomBytes[0]) >> (i & 7)) & 1 // random binary digit [0, 1] from the random bytes
		for i := range ternarySampler.baseRing.Modulus {
			pol.Coeffs[i][index[j]] = ternarySampler.matrixValues[i][coeff]
		}

		// Removes the element in position j of the slice (order not preserved)
		index[j] = index[len(index)-1]
		index = index[:len(index)-1]

		pointer++

		if pointer == 8 {
			randomBytes = randomBytes[1:]
			pointer = 0
		}
	}
}

// kysampling use the binary expension and random bytes matrix to sample a discret gaussian value and its sign.
func (ts *TernarySampler) kysampling(prng utils.PRNG, randomBytes []byte, pointer uint8, bytePointer uint64, byteLength uint64) (uint64, uint64, []byte, uint8, uint64) {

	var sign uint8

	d := 0
	col := 0
	colLen := len(ts.matrixProba)

	for {

		// Uses one random byte per cycle and cycle through the randombytes
		for i := pointer; i < 8; i++ {

			d = (d << 1) + 1 - int((uint8(randomBytes[bytePointer])>>i)&1)

			// There is small probability that it will get out of the bound, then
			// rerun until it gets a proper output
			if d > colLen-1 {
				return ts.kysampling(prng, randomBytes, i, bytePointer, byteLength)
			}

			for row := colLen - 1; row >= 0; row-- {

				d -= int(ts.matrixProba[row][col])

				if d == -1 {

					// Sign
					if i == 7 {
						pointer = 0
						// If the last bit of the array was read, samples a new one
						bytePointer++

						if bytePointer >= byteLength {
							bytePointer = 0
							prng.Clock(randomBytes)
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
		// If the last bit of the array was read, samples a new one
		bytePointer++

		if bytePointer >= byteLength {
			bytePointer = 0
			prng.Clock(randomBytes)
		}

	}
}
