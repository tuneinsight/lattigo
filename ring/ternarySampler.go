package ring

import (
	"crypto/rand"
	"math"
	"math/bits"
)

func computeMatrixTernary(p float64) (M [][]uint8) {
	var g float64
	var x uint64

	precision := uint64(56)

	M = make([][]uint8, 2)

	g = p
	g *= math.Exp2(float64(precision))
	x = uint64(g)

	M[0] = make([]uint8, precision-1)
	for j := uint64(0); j < precision-1; j++ {
		M[0][j] = uint8((x >> (precision - j - 1)) & 1)
	}

	g = 1 - p
	g *= math.Exp2(float64(precision))
	x = uint64(g)

	M[1] = make([]uint8, precision-1)
	for j := uint64(0); j < precision-1; j++ {
		M[1][j] = uint8((x >> (precision - j - 1)) & 1)
	}

	return M
}

// SampleMontgomeryNew samples coefficients with ternary distribution in montgomery form on the target polynomial.
func (context *Context) sampleTernary(samplerMatrix [][]uint64, p float64, pol *Poly) {

	if p == 0 {
		panic("cannot sample -> p = 0")
	}

	var coeff uint64
	var sign uint64
	var index uint64

	if p == 0.5 {

		randomBytesCoeffs := make([]byte, context.N>>3)
		randomBytesSign := make([]byte, context.N>>3)

		if _, err := rand.Read(randomBytesCoeffs); err != nil {
			panic("crypto rand error")
		}

		if _, err := rand.Read(randomBytesSign); err != nil {
			panic("crypto rand error")
		}

		for i := uint64(0); i < context.N; i++ {
			coeff = uint64(uint8(randomBytesCoeffs[i>>3])>>(i&7)) & 1
			sign = uint64(uint8(randomBytesSign[i>>3])>>(i&7)) & 1

			index = (coeff & (sign ^ 1)) | ((sign & coeff) << 1)

			for j := range context.Modulus {
				pol.Coeffs[j][i] = samplerMatrix[j][index] //(coeff & (sign^1)) | (qi - 1) * (sign & coeff)
			}
		}

	} else {

		matrix := computeMatrixTernary(p)

		randomBytes := make([]byte, 8)

		pointer := uint8(0)

		if _, err := rand.Read(randomBytes); err != nil {
			panic("crypto rand error")
		}

		for i := uint64(0); i < context.N; i++ {

			coeff, sign, randomBytes, pointer = kysampling(matrix, randomBytes, pointer)

			index = (coeff & (sign ^ 1)) | ((sign & coeff) << 1)

			for j := range context.Modulus {
				pol.Coeffs[j][i] = samplerMatrix[j][index] //(coeff & (sign^1)) | (qi - 1) * (sign & coeff)
			}
		}
	}
}

func (context *Context) SampleTernaryUniform(pol *Poly) {
	context.sampleTernary(context.matrixTernary, 1.0/3.0, pol)
}

func (context *Context) SampleTernary(pol *Poly, p float64) {
	context.sampleTernary(context.matrixTernary, p, pol)
}

func (context *Context) SampleTernaryMontgomery(pol *Poly, p float64) {
	context.sampleTernary(context.matrixTernaryMontgomery, p, pol)
}

// SampleNew samples a new polynomial with ternary distribution.
func (context *Context) SampleTernaryNew(p float64) (pol *Poly) {
	pol = context.NewPoly()
	context.SampleTernary(pol, p)
	return pol
}

// SampleMontgomeryNew samples a new polynomial with ternary distribution in montgomery form.
func (context *Context) SampleTernaryMontgomeryNew(p float64) (pol *Poly) {
	pol = context.NewPoly()
	context.SampleTernaryMontgomery(pol, p)
	return
}

// SampleNTTNew samples a new polynomial with ternary distribution in the NTT domain.
func (context *Context) SampleTernaryNTTNew(p float64) (pol *Poly) {
	pol = context.NewPoly()
	context.SampleTernary(pol, p)
	context.NTT(pol, pol)
	return
}

// SampleNTT samples coefficients with ternary distribution in the NTT domain on the target polynomial.
func (context *Context) SampleTernaryNTT(pol *Poly, p float64) {
	context.SampleTernary(pol, p)
	context.NTT(pol, pol)
}

// SampleNTTNew samples a new polynomial with ternary distribution in the NTT domain and in montgomery form.
func (context *Context) SampleTernaryMontgomeryNTTNew(p float64) (pol *Poly) {
	pol = context.SampleTernaryMontgomeryNew(p)
	context.NTT(pol, pol)
	return
}

// SampleNTT samples coefficients with ternary distribution in the NTT domain and in montgomery form on the target polynomial.
func (context *Context) SampleTernaryMontgomeryNTT(pol *Poly, p float64) {
	context.SampleTernaryMontgomery(pol, p)
	context.NTT(pol, pol)
}

// Samples a keys with distribution [-1, 1] = [1/2, 1/2] with exactly hw non zero coefficients
func (context *Context) sampleTernarySparse(samplerMatrix [][]uint64, pol *Poly, hw uint64) {

	if hw > context.N {
		hw = context.N
	}

	var mask, j uint64
	var coeff uint8

	index := make([]uint64, context.N)
	for i := uint64(0); i < context.N; i++ {
		index[i] = i
	}

	randomBytes := make([]byte, (uint64(math.Ceil(float64(hw) / 8.0)))) // We sample ceil(hw/8) bytes
	pointer := uint8(0)

	if _, err := rand.Read(randomBytes); err != nil {
		panic("crypto rand error")
	}

	for i := uint64(0); i < hw; i++ {
		mask = (1 << uint64(bits.Len64(context.N-i))) - 1 // rejection sampling of a random variable between [0, len(index)]

		j = randInt32(mask)
		for j >= context.N-i {
			j = randInt32(mask)
		}

		coeff = (uint8(randomBytes[0]) >> (i & 7)) & 1 // random binary digit [0, 1] from the random bytes
		for i := range context.Modulus {
			pol.Coeffs[i][index[j]] = samplerMatrix[i][coeff]
		}

		// Removes the element in position j of the slice (order not preserved)
		index[j] = index[len(index)-1]
		index = index[:len(index)-1]

		pointer += 1

		if pointer == 8 {
			randomBytes = randomBytes[1:]
			pointer = 0
		}
	}
}

func (context *Context) SampleTernarySparse(pol *Poly, hw uint64) {
	context.sampleTernarySparse(context.matrixTernary, pol, hw)
}

func (context *Context) SampleTernarySparseNew(hw uint64) (pol *Poly) {
	pol = context.NewPoly()
	context.SampleTernarySparse(pol, hw)
	return pol
}

func (context *Context) SampleTernarySparseNTT(pol *Poly, hw uint64) {
	context.SampleTernarySparse(pol, hw)
	context.NTT(pol, pol)
}

func (context *Context) SampleTernarySparseNTTNew(hw uint64) (pol *Poly) {
	pol = context.NewPoly()
	context.sampleTernarySparse(context.matrixTernaryMontgomery, pol, hw)
	context.NTT(pol, pol)
	return pol
}

func (context *Context) SampleTernarySparseMontgomery(pol *Poly, hw uint64) {
	context.sampleTernarySparse(context.matrixTernaryMontgomery, pol, hw)
}

func (context *Context) SampleSparseMontgomeryNew(hw uint64) (pol *Poly) {
	pol = context.NewPoly()
	context.sampleTernarySparse(context.matrixTernaryMontgomery, pol, hw)
	return pol
}

func (context *Context) SampleTernarySparseMontgomeryNTTNew(hw uint64) (pol *Poly) {
	pol = context.SampleSparseMontgomeryNew(hw)
	context.NTT(pol, pol)
	return pol
}

func (context *Context) SampleTernarySparseMontgomeryNTT(pol *Poly, hw uint64) {
	context.SampleTernarySparseMontgomery(pol, hw)
	context.NTT(pol, pol)
}
