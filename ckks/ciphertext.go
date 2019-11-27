package ckks

import (
	"github.com/ldsec/lattigo/ring"
)

// Ciphertext is a BigPoly of degree > 0.
type Ciphertext struct {
	*ckksElement
}

// NewCiphertextStruct returns a new Ciphertext element.
func NewCiphertextStruct() (ciphertext *Ciphertext) {
	return &Ciphertext{&ckksElement{}}
}

func NewRingContext(params *Parameters) *ring.Context {
	scalechain := make([]float64, len(params.Modulichain))

	// Extracts all the different primes bit size and maps their number
	primesbitlen := make(map[uint64]uint64)
	for i, qi := range params.Modulichain {

		primesbitlen[uint64(qi)]++

		if uint64(params.Modulichain[i]) > 60 {
			panic("provided moduli must be smaller than 61")
		}
	}

	for _, pj := range params.P {
		primesbitlen[uint64(pj)]++

		if uint64(pj) > 60 {
			panic("provided P must be smaller than 61")
		}
	}

	// For each bitsize, finds that many primes
	primes := make(map[uint64][]uint64)
	for key, value := range primesbitlen {
		primes[key] = GenerateCKKSPrimes(key, uint64(params.LogN), value)
	}

	// Assigns the primes to the ckks moduli chain
	moduli := make([]uint64, len(params.Modulichain))
	for i, qi := range params.Modulichain {
		moduli[i] = primes[uint64(params.Modulichain[i])][0]
		primes[uint64(qi)] = primes[uint64(qi)][1:]

		scalechain[i] = float64(moduli[i])
	}

	ringCtx := ring.NewContext()
	ringCtx.SetParameters(1<<params.LogN, moduli)

	err := ringCtx.GenNTTParams()
	if err != nil {
		panic(err)
	}

	return ringCtx
}

// NewCiphertext creates a new ciphertext parameterized by degree, level and scale.
func NewCiphertext(degree uint64, level uint64, scale float64, ringCtx *ring.Context) *Ciphertext {
	ciphertext := &Ciphertext{&ckksElement{}}

	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = ringCtx.NewPolyLvl(level)
	}

	ciphertext.scale = scale
	ciphertext.isNTT = true

	return ciphertext
}

// NewRandomCiphertext generates a new uniformely distributed ciphertext of degree, level and scale.
func NewRandomCiphertext(degree, level uint64, scale float64, ringCtx *ring.Context) (ciphertext *Ciphertext) {
	ciphertext = &Ciphertext{&ckksElement{}}

	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = ringCtx.NewUniformPolyLvl(level)
	}

	ciphertext.scale = scale
	ciphertext.isNTT = true

	return ciphertext
}
