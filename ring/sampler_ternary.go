package ring

import (
	"math"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

// TernarySampler keeps the state of a polynomial sampler in the ternary distribution.
type TernarySampler struct {
	baseSampler
	matrixProba  [2][precision - 1]uint8
	matrixValues [][3]uint64
	p            float64
	hw           int
	sample       func(poly *Poly)
}

// NewTernarySampler creates a new instance of TernarySampler from a PRNG, the ring definition and the distribution
// parameters: p is the probability of a coefficient being 0, (1-p)/2 is the probability of 1 and -1. If "montgomery"
// is set to true, polynomials read from this sampler are in Montgomery form.
func NewTernarySampler(prng utils.PRNG, baseRing *Ring, X *UniformTernary, montgomery bool) (ts *TernarySampler) {
	ts = new(TernarySampler)
	ts.baseRing = baseRing
	ts.prng = prng
	ts.p = X.P
	ts.sampleLvl = ts.sampleProbaLvl
	ts.sampleLvlAndAddLvl = ts.sampleProbaAndAddLvl

	ts.initializeMatrix(montgomery)

	if ts.p != 0.5 {
		ts.computeMatrixTernary(ts.p)
	}

	return
}

// AtLevel returns an instance of the target TernarySampler that operates at the target level.
// This instance is not thread safe and cannot be used concurrently to the base instance.
func (ts *TernarySampler) AtLevel(level int) Sampler {
	return &TernarySampler{
		baseSampler:  ts.baseSampler.AtLevel(level),
		matrixProba:  ts.matrixProba,
		matrixValues: ts.matrixValues,
		p:            ts.p,
		hw:           ts.hw,
		sample:       ts.sample,
	}
}

// NewTernarySamplerWithHammingWeight creates a new instance of a fixed-hamming-weight TernarySampler from a PRNG, the ring definition and the desired
// hamming weight for the output polynomials. If "montgomery" is set to true, polynomials read from this sampler
// are in Montgomery form.
func NewTernarySamplerWithHammingWeight(prng utils.PRNG, baseRing *Ring, X *SparseTernary, montgomery bool) (ts *TernarySampler) {
	ts = new(TernarySampler)
	ts.baseRing = baseRing
	ts.prng = prng
	ts.hw = X.HammingWeight
	ts.sampleLvl = ts.sampleSparseLvl
	ts.sampleLvlAndAddLvl = ts.sampleSparseAndAddLvl
	ts.initializeMatrix(montgomery)
	return
}

// Read samples a polynomial into pol.
func (ts *TernarySampler) Read(pol *Poly) {
	ts.sample(pol)
}

// ReadNew allocates and samples a polynomial at the max level.
func (ts *TernarySampler) ReadNew() (pol *Poly) {
	pol = ts.baseRing.NewPoly()
	ts.sample(pol)
	return pol
}

func (ts *TernarySampler) initializeMatrix(montgomery bool) {
	ts.matrixValues = make([][3]uint64, ts.baseRing.ModuliChainLength())

	// [0] = 0
	// [1] = 1 * 2^64 mod qi
	// [2] = (qi - 1) * 2^64 mod qi

	for i, s := range ts.baseRing.SubRings {

		modulus := s.Modulus
		brc := s.BRedConstant

		ts.matrixValues[i][0] = 0

		if montgomery {
			ts.matrixValues[i][1] = MForm(1, modulus, brc)
			ts.matrixValues[i][2] = MForm(modulus-1, modulus, brc)
		} else {
			ts.matrixValues[i][1] = 1
			ts.matrixValues[i][2] = Table.Modulus - 1
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

func (ts *TernarySampler) sampleProba(pol *Poly) {

	if ts.p == 0 {
		panic("cannot sample -> p = 0")
	}

	var coeff uint64
	var sign uint64
	var index uint64

	level := ts.baseRing.level

	N := ts.baseRing.N()

	m := ts.matrixValues

	if ts.p == 0.5 {

		randomBytesCoeffs := make([]byte, N>>3)
		randomBytesSign := make([]byte, N>>3)

		if _, err := ts.prng.Read(randomBytesCoeffs); err != nil {
			panic(err)
		}

		if _, err := ts.prng.Read(randomBytesSign); err != nil {
			panic(err)
		}

		for i := 0; i < N; i++ {
			coeff = uint64(uint8(randomBytesCoeffs[i>>3])>>(i&7)) & 1
			sign = uint64(uint8(randomBytesSign[i>>3])>>(i&7)) & 1

			index = (coeff & (sign ^ 1)) | ((sign & coeff) << 1)

			for j := 0; j < level+1; j++ {
				pol.Coeffs[j][i] = lut[j][index]
			}
		}

	} else {

		randomBytes := make([]byte, N)

		pointer := uint8(0)
		var bytePointer int

		if _, err := ts.prng.Read(randomBytes); err != nil {
			panic(err)
		}

		for i := 0; i < N; i++ {

			coeff, sign, randomBytes, pointer, bytePointer = ts.kysampling(ts.prng, randomBytes, pointer, bytePointer, N)

			index = (coeff & (sign ^ 1)) | ((sign & coeff) << 1)

			for j := 0; j < level+1; j++ {
				pol.Coeffs[j][i] = lut[j][index]
			}
		}
	}
}

func (ts *TernarySampler) sampleSparse(pol *Poly) {

	if ts.p == 0 {
		panic("cannot sample -> p = 0")
	}

	var coeff uint64
	var sign uint64
	var index uint64

	coeffs := pol.Coeffs

	moduli := ts.baseRing.Moduli()[:level+1]

	N := ts.baseRing.N()

	m := ts.matrixValues

	if ts.p == 0.5 {

		randomBytesCoeffs := make([]byte, N>>3)
		randomBytesSign := make([]byte, N>>3)

		ts.prng.Read(randomBytesCoeffs)

		ts.prng.Read(randomBytesSign)

		for i := 0; i < N; i++ {
			coeff = uint64(uint8(randomBytesCoeffs[i>>3])>>(i&7)) & 1
			sign = uint64(uint8(randomBytesSign[i>>3])>>(i&7)) & 1

			index = (coeff & (sign ^ 1)) | ((sign & coeff) << 1)

			for j, qi := range moduli {
				coeffs[j][i] = CRed(coeffs[j][i]+m[j][index], qi)
			}
		}

	} else {

		randomBytes := make([]byte, N)

		pointer := uint8(0)
		var bytePointer int

		ts.prng.Read(randomBytes)

		for i := 0; i < N; i++ {

			coeff, sign, randomBytes, pointer, bytePointer = ts.kysampling(ts.prng, randomBytes, pointer, bytePointer, N)

			index = (coeff & (sign ^ 1)) | ((sign & coeff) << 1)

			for j, qi := range moduli {
				coeffs[j][i] = CRed(coeffs[j][i]+m[j][index], qi)
			}
		}
	}
}

func (ts *TernarySampler) sampleSparseLvl(level int, pol *Poly) {

	N := ts.baseRing.N()

	if ts.hw > N {
		ts.hw = N
	}

	var mask, j uint64
	var coeff uint8

	index := make([]int, N)
	for i := 0; i < N; i++ {
		index[i] = i
	}

	randomBytes := make([]byte, (uint64(math.Ceil(float64(ts.hw) / 8.0)))) // We sample ceil(hw/8) bytes
	pointer := uint8(0)

	if _, err := ts.prng.Read(randomBytes); err != nil {
		panic(err)
	}

	level := ts.baseRing.level

	for i := 0; i < hw; i++ {
		mask = (1 << uint64(bits.Len64(uint64(N-i)))) - 1 // rejection sampling of a random variable between [0, len(index)]

		j = randInt32(ts.prng, mask)
		for j >= uint64(N-i) {
			j = randInt32(ts.prng, mask)
		}

		coeff = (uint8(randomBytes[0]) >> (i & 7)) & 1 // random binary digit [0, 1] from the random bytes (0 = 1, 1 = -1)
		for k := 0; k < level+1; k++ {
			pol.Coeffs[k][index[j]] = ts.matrixValues[k][coeff+1]
		}

		// Remove the element in position j of the slice (order not preserved)
		index[j] = index[len(index)-1]
		index = index[:len(index)-1]

		pointer++

		if pointer == 8 {
			randomBytes = randomBytes[1:]
			pointer = 0
		}
	}

	for _, i := range index {
		for k := 0; k < level+1; k++ {
			pol.Coeffs[k][i] = 0
		}
	}
}

func (ts *TernarySampler) sampleSparseAndAddLvl(level int, pol *Poly) {

	N := ts.baseRing.N()

	if ts.hw > N {
		ts.hw = N
	}

	var mask, j uint64
	var coeff uint8

	index := make([]int, N)
	for i := 0; i < N; i++ {
		index[i] = i
	}

	randomBytes := make([]byte, (uint64(math.Ceil(float64(ts.hw) / 8.0)))) // We sample ceil(hw/8) bytes
	pointer := uint8(0)

	ts.prng.Read(randomBytes)

	coeffs := pol.Coeffs

	moduli := ts.baseRing.Moduli()[:level+1]

	m := ts.matrixValues

	for i := 0; i < ts.hw; i++ {
		mask = (1 << uint64(bits.Len64(uint64(N-i)))) - 1 // rejection sampling of a random variable between [0, len(index)]

		j = randInt32(ts.prng, mask)
		for j >= uint64(N-i) {
			j = randInt32(ts.prng, mask)
		}

		coeff = (uint8(randomBytes[0]) >> (i & 7)) & 1 // random binary digit [0, 1] from the random bytes (0 = 1, 1 = -1)
		for k, qi := range moduli {
			coeffs[k][index[j]] = CRed(coeffs[k][index[j]]+m[k][coeff+1], qi)
		}

		// Remove the element in position j of the slice (order not preserved)
		index[j] = index[len(index)-1]
		index = index[:len(index)-1]

		pointer++

		if pointer == 8 {
			randomBytes = randomBytes[1:]
			pointer = 0
		}
	}
}

// kysampling uses the binary expansion and random bytes matrix to sample a discrete Gaussian value and its sign.
func (ts *TernarySampler) kysampling(prng sampling.PRNG, randomBytes []byte, pointer uint8, bytePointer, byteLength int) (uint64, uint64, []byte, uint8, int) {

	var sign uint8

	d := 0
	col := 0
	colLen := len(ts.matrixProba)

	for {

		// Use one random byte per cycle and cycle through the randomBytes
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
						// If the last bit of the array was read, sample a new one
						bytePointer++

						if bytePointer >= byteLength {
							bytePointer = 0
							if _, err := prng.Read(randomBytes); err != nil {
								panic(err)
							}
						}

						sign = uint8(randomBytes[bytePointer]) & 1

					} else {
						pointer = i
						// Otherwise, the sign is the next bit of the byte
						sign = uint8(randomBytes[bytePointer]>>(i+1)) & 1
					}

					return uint64(row), uint64(sign), randomBytes, pointer + 1, bytePointer
				}
			}

			col++
		}

		// Reset the bit pointer and discard the used byte
		pointer = 0
		// If the last bit of the array was read, sample a new one
		bytePointer++

		if bytePointer >= byteLength {
			bytePointer = 0
			if _, err := prng.Read(randomBytes); err != nil {
				panic(err)
			}
		}
	}
}
