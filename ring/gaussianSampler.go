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

// kysampling use the binary expension and random bytes matrix to sample a discret gaussian value and its sign.
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
						// If the last bit of the array was read, samples a new one
						bytePointer++

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
		// If the last bit of the array was read, samples a new one
		bytePointer++

		if bytePointer >= byteLength {
			bytePointer = 0
			if _, err := rand.Read(randomBytes); err != nil {
				panic("crypto rand error")
			}
		}

	}
}
