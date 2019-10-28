package ring

import (
	"crypto/rand"
	"errors"
	"math"
	"math/bits"
)

// TernarySampler is the structure holding the parameters for sampling polynomials of the form [-1, 0, 1].
type TernarySampler struct {
	context          *Context
	Matrix           [][]uint64
	MatrixMontgomery [][]uint64

	KYMatrix [][]uint8
}

// NewTernarySampler creates a new TernarySampler from the target context.
func (context *Context) NewTernarySampler() *TernarySampler {

	sampler := new(TernarySampler)
	sampler.context = context

	sampler.Matrix = make([][]uint64, len(context.Modulus))
	sampler.MatrixMontgomery = make([][]uint64, len(context.Modulus))

	for i, Qi := range context.Modulus {

		sampler.Matrix[i] = make([]uint64, 3)
		sampler.Matrix[i][0] = 0
		sampler.Matrix[i][1] = 1
		sampler.Matrix[i][2] = Qi - 1

		sampler.MatrixMontgomery[i] = make([]uint64, 3)
		sampler.MatrixMontgomery[i][0] = 0
		sampler.MatrixMontgomery[i][1] = MForm(1, Qi, context.bredParams[i])
		sampler.MatrixMontgomery[i][2] = MForm(Qi-1, Qi, context.bredParams[i])
	}

	return sampler
}

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
func (sampler *TernarySampler) sample(samplerMatrix [][]uint64, p float64, pol *Poly) (err error) {

	if p == 0 {
		return errors.New("cannot sample -> p = 0")
	}

	var coeff uint64
	var sign uint64
	var index uint64

	if p == 0.5 {

		randomBytesCoeffs := make([]byte, sampler.context.N>>3)
		randomBytesSign := make([]byte, sampler.context.N>>3)

		if _, err := rand.Read(randomBytesCoeffs); err != nil {
			panic("crypto rand error")
		}

		if _, err := rand.Read(randomBytesSign); err != nil {
			panic("crypto rand error")
		}

		for i := uint64(0); i < sampler.context.N; i++ {
			coeff = uint64(uint8(randomBytesCoeffs[i>>3])>>(i&7)) & 1
			sign = uint64(uint8(randomBytesSign[i>>3])>>(i&7)) & 1

			index = (coeff & (sign ^ 1)) | ((sign & coeff) << 1)

			for j := range sampler.context.Modulus {
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

		for i := uint64(0); i < sampler.context.N; i++ {

			coeff, sign, randomBytes, pointer = kysampling(matrix, randomBytes, pointer)

			index = (coeff & (sign ^ 1)) | ((sign & coeff) << 1)

			for j := range sampler.context.Modulus {
				pol.Coeffs[j][i] = samplerMatrix[j][index] //(coeff & (sign^1)) | (qi - 1) * (sign & coeff)
			}
		}
	}

	return nil
}

func (sampler *TernarySampler) SampleUniform(pol *Poly) {
	_ = sampler.sample(sampler.Matrix, 1.0/3.0, pol)
}

func (sampler *TernarySampler) Sample(p float64, pol *Poly) (err error) {
	if err = sampler.sample(sampler.Matrix, p, pol); err != nil {
		return err
	}
	return nil
}

func (sampler *TernarySampler) SampleMontgomery(p float64, pol *Poly) (err error) {
	if err = sampler.sample(sampler.MatrixMontgomery, p, pol); err != nil {
		return err
	}
	return nil
}

// SampleNew samples a new polynomial with ternary distribution.
func (sampler *TernarySampler) SampleNew(p float64) (pol *Poly, err error) {
	pol = sampler.context.NewPoly()
	if err = sampler.Sample(p, pol); err != nil {
		return nil, err
	}
	return pol, nil
}

// SampleMontgomeryNew samples a new polynomial with ternary distribution in montgomery form.
func (sampler *TernarySampler) SampleMontgomeryNew(p float64) (pol *Poly, err error) {
	pol = sampler.context.NewPoly()
	if err = sampler.SampleMontgomery(p, pol); err != nil {
		return nil, err
	}
	return pol, nil
}

// SampleNTTNew samples a new polynomial with ternary distribution in the NTT domain.
func (sampler *TernarySampler) SampleNTTNew(p float64) (pol *Poly, err error) {
	pol = sampler.context.NewPoly()
	if err = sampler.Sample(p, pol); err != nil {
		return nil, err
	}
	sampler.context.NTT(pol, pol)
	return pol, nil
}

// SampleNTT samples coefficients with ternary distribution in the NTT domain on the target polynomial.
func (sampler *TernarySampler) SampleNTT(p float64, pol *Poly) (err error) {
	if err = sampler.Sample(p, pol); err != nil {
		return err
	}
	sampler.context.NTT(pol, pol)

	return nil
}

// SampleNTTNew samples a new polynomial with ternary distribution in the NTT domain and in montgomery form.
func (sampler *TernarySampler) SampleMontgomeryNTTNew(p float64) (pol *Poly, err error) {
	if pol, err = sampler.SampleMontgomeryNew(p); err != nil {
		return nil, err
	}
	sampler.context.NTT(pol, pol)
	return pol, nil
}

// SampleNTT samples coefficients with ternary distribution in the NTT domain and in montgomery form on the target polynomial.
func (sampler *TernarySampler) SampleMontgomeryNTT(p float64, pol *Poly) (err error) {
	if err = sampler.SampleMontgomery(p, pol); err != nil {
		return err
	}
	sampler.context.NTT(pol, pol)
	return nil
}

// Samples a keys with distribution [-1, 1] = [1/2, 1/2] with exactly hw non zero coefficients
func (sampler *TernarySampler) SampleSparse(pol *Poly, hw uint64) {

	if hw > sampler.context.N {
		hw = sampler.context.N
	}

	var mask, j uint64
	var coeff uint8

	index := make([]uint64, sampler.context.N)
	for i := uint64(0); i < sampler.context.N; i++ {
		index[i] = i
	}

	randomBytes := make([]byte, (uint64(math.Ceil(float64(hw) / 8.0)))) // We sample ceil(hw/8) bytes
	pointer := uint8(0)

	if _, err := rand.Read(randomBytes); err != nil {
		panic("crypto rand error")
	}

	for i := uint64(0); i < hw; i++ {
		mask = (1 << uint64(bits.Len64(sampler.context.N-i))) - 1 // rejection sampling of a random variable between [0, len(index)]

		j = randInt32(mask)
		for j >= sampler.context.N-i {
			j = randInt32(mask)
		}

		coeff = (uint8(randomBytes[0]) >> (i & 7)) & 1 // random binary digit [0, 1] from the random bytes
		for i := range sampler.context.Modulus {
			pol.Coeffs[i][index[j]] = sampler.Matrix[i][coeff]
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

func (sampler *TernarySampler) SampleSparseMontgomeryNew(hw uint64) (pol *Poly) {
	pol = sampler.context.NewPoly()
	sampler.SampleSparseMontgomery(pol, hw)
	return pol
}

// Samples a keys with distribution [-1, 1] = [1/2, 1/2] with exactly hw non zero coefficients
func (sampler *TernarySampler) SampleSparseMontgomery(pol *Poly, hw uint64) {

	if hw > sampler.context.N {
		hw = sampler.context.N
	}

	var mask, j uint64
	var coeff uint8

	index := make([]uint64, sampler.context.N)
	for i := uint64(0); i < sampler.context.N; i++ {
		index[i] = i
	}

	randomBytes := make([]byte, (uint64(math.Ceil(float64(hw) / 8.0)))) // We sample ceil(hw/8) bytes
	pointer := uint8(0)

	if _, err := rand.Read(randomBytes); err != nil {
		panic("crypto rand error")
	}

	for i := uint64(0); i < hw; i++ {

		mask = (1 << uint64(bits.Len64(sampler.context.N-i))) - 1 // rejection sampling of a random variable between [0, len(index)]

		j = randInt32(mask)
		for j >= sampler.context.N-i {
			j = randInt32(mask)
		}

		coeff = (uint8(randomBytes[0]) >> (i & 7)) & 1 // random binary digit [0, 1] from the random bytes
		for i := range sampler.context.Modulus {
			pol.Coeffs[i][index[j]] = sampler.MatrixMontgomery[i][coeff]
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

func (sampler *TernarySampler) SampleSparseNTTNew(hw uint64) (pol *Poly) {
	pol = sampler.context.NewPoly()
	sampler.SampleSparse(pol, hw)
	sampler.context.NTT(pol, pol)
	return pol
}

func (sampler *TernarySampler) SampleSparseNTT(hw uint64, pol *Poly) {
	sampler.SampleSparse(pol, hw)
	sampler.context.NTT(pol, pol)
}

func (sampler *TernarySampler) SampleSparseMontgomeryNTTNew(hw uint64) (pol *Poly) {
	pol = sampler.SampleSparseMontgomeryNew(hw)
	sampler.context.NTT(pol, pol)
	return pol
}

func (sampler *TernarySampler) SampleSparseMontgomeryNTT(hw uint64, pol *Poly) {
	sampler.SampleSparseMontgomery(pol, hw)
	sampler.context.NTT(pol, pol)
}

func (sampler *TernarySampler) SampleSparseNew(hw uint64) (pol *Poly) {
	pol = sampler.context.NewPoly()
	sampler.SampleSparse(pol, hw)
	return pol
}
