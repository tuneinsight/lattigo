package ring

import (
	"crypto/rand"
)

// SampleGaussianLvl samples a truncated gaussian polynomial with variance
// sigma of moduli 0 to level within the given bound using the Ziggurat algorithm.
func (context *Context) SampleGaussianLvl(level uint64, pol *Poly, sigma float64, bound uint64) {

	var coeffFlo float64
	var coeffInt uint64
	var sign uint64
	var ptr uint64

	randomBytes := make([]byte, context.N)

	if _, err := rand.Read(randomBytes); err != nil {
		panic("crypto rand error")
	}

	for i := uint64(0); i < context.N; i++ {

		for {
			coeffFlo, sign, randomBytes, ptr = normFloat64(ptr, randomBytes)

			if coeffInt = uint64(coeffFlo * sigma); coeffInt <= bound {
				break
			}
		}

		for j, qi := range context.Modulus[:level+1] {
			pol.Coeffs[j][i] = (coeffInt * sign) | (qi-coeffInt)*(sign^1)
		}
	}
}

// SampleGaussianAndAddLvl adds on the input polynomial a truncated gaussian polynomial of moduli 0 to level
// with variance sigma within the given bound using the Ziggurat algorithm.
func (context *Context) SampleGaussianAndAddLvl(level uint64, pol *Poly, sigma float64, bound uint64) {

	var coeffFlo float64
	var coeffInt uint64
	var sign uint64
	var ptr uint64

	randomBytes := make([]byte, context.N)

	if _, err := rand.Read(randomBytes); err != nil {
		panic("crypto rand error")
	}

	for i := uint64(0); i < context.N; i++ {

		for {
			coeffFlo, sign, randomBytes, ptr = normFloat64(ptr, randomBytes)

			if coeffInt = uint64(coeffFlo * sigma); coeffInt <= bound {
				break
			}
		}

		for j, qi := range context.Modulus[:level+1] {
			pol.Coeffs[j][i] = CRed(pol.Coeffs[j][i]+((coeffInt*sign)|(qi-coeffInt)*(sign^1)), qi)
		}
	}
}

// SampleGaussianNew samples a new truncated gaussian polynomial with
// variance sigma within the given bound using the Ziggurat algorithm.
func (context *Context) SampleGaussianNew(sigma float64, bound uint64) (pol *Poly) {
	pol = context.NewPoly()
	context.SampleGaussianLvl(uint64(len(context.Modulus)-1), pol, sigma, bound)
	return
}

// SampleGaussianNTTLvl samples a trucated gaussian polynomial in the NTT domain of moduli 0 to level
// with variance sigma within the given bound using the Ziggurat algorithm.
func (context *Context) SampleGaussianNTTLvl(level uint64, pol *Poly, sigma float64, bound uint64) {
	context.SampleGaussianLvl(level, pol, sigma, bound)
	context.NTT(pol, pol)
}

// SampleGaussianNTTNew samples a new trucated gaussian polynomial in the NTT domain
// with variance sigma within the given bound using the Ziggurat algorithm
func (context *Context) SampleGaussianNTTNew(sigma float64, bound uint64) (pol *Poly) {
	pol = context.SampleGaussianNew(sigma, bound)
	context.NTT(pol, pol)
	return
}
