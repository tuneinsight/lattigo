package rlwe

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// Plaintext is a common base type for RLWE plaintexts.
type Plaintext struct {
	Value *ring.Poly
}

// Ciphertext is a generic type for RLWE ciphertext.
type Ciphertext struct {
	Value []*ring.Poly
}

// AdditiveShare is a type for storing additively shared values in Z_Q[X] (RNS domain)
type AdditiveShare struct {
	Value ring.Poly
}

// AdditiveShareBigint is a type for storing additively shared values
// in Z (positional domain)
type AdditiveShareBigint struct {
	Value []*big.Int
}

// NewAdditiveShare instantiate a new additive share struct for the ring defined
// by the given parameters at maximum level.
func NewAdditiveShare(params Parameters) *AdditiveShare {
	return &AdditiveShare{Value: *ring.NewPoly(params.N(), 1)}
}

// NewAdditiveShareAtLevel instantiate a new additive share struct for the ring defined
// by the given parameters at level `level`.
func NewAdditiveShareAtLevel(params Parameters, level int) *AdditiveShare {
	return &AdditiveShare{Value: *ring.NewPoly(params.N(), level+1)}
}

// NewAdditiveShareBigint instantiate a new additive share struct composed of "n" big.Int elements
func NewAdditiveShareBigint(params Parameters, n int) *AdditiveShareBigint {
	v := make([]*big.Int, n)
	for i := range v {
		v[i] = new(big.Int)
	}
	return &AdditiveShareBigint{Value: v}
}

// NewPlaintext creates a new Plaintext at level `level` from the parameters.
func NewPlaintext(params Parameters, level int) *Plaintext {
	return &Plaintext{Value: ring.NewPoly(params.N(), level+1)}
}

// Degree returns the degree of the target element.
func (pt Plaintext) Degree() int {
	return 0
}

// Level returns the level of the target element.
func (pt Plaintext) Level() int {
	return len(pt.Value.Coeffs) - 1
}

// El returns the plaintext as a new `Element` for which the value points
// to the receiver `Value` field.
func (pt Plaintext) El() *Ciphertext {
	return &Ciphertext{Value: []*ring.Poly{pt.Value}}
}

// Copy copies the `other` plaintext value into the reciever plaintext.
func (pt *Plaintext) Copy(other *Plaintext) {
	if other != nil && other.Value != nil {
		pt.Value.Copy(other.Value)
	}
}

// NewCiphertext returns a new Element with zero values.
func NewCiphertext(params Parameters, degree, level int) *Ciphertext {
	el := new(Ciphertext)
	el.Value = make([]*ring.Poly, degree+1)
	for i := 0; i < degree+1; i++ {
		el.Value[i] = ring.NewPoly(params.N(), level+1)
	}
	return el
}

// NewCiphertextNTT returns a new Element with zero values and the NTT flags set.
func NewCiphertextNTT(params Parameters, degree, level int) *Ciphertext {
	el := new(Ciphertext)
	el.Value = make([]*ring.Poly, degree+1)
	for i := 0; i < degree+1; i++ {
		el.Value[i] = ring.NewPoly(params.N(), level+1)
		el.Value[i].IsNTT = true
	}
	return el
}

// NewCiphertextRandom generates a new uniformly distributed Ciphertext of degree, level and scale.
func NewCiphertextRandom(prng utils.PRNG, params Parameters, degree, level int) (ciphertext *Ciphertext) {
	ciphertext = NewCiphertext(params, degree, level)
	PopulateElementRandom(prng, params, ciphertext)
	return
}

// SetValue sets the input slice of polynomials as the value of the target element.
func (el *Ciphertext) SetValue(value []*ring.Poly) {
	el.Value = value
}

// Degree returns the degree of the target element.
func (el *Ciphertext) Degree() int {
	return len(el.Value) - 1
}

// Level returns the level of the target element.
func (el *Ciphertext) Level() int {
	return len(el.Value[0].Coeffs) - 1
}

// Resize resizes the degree of the target element.
// Sets the NTT flag of the added poly equal to the NTT flag
// to the poly at degree zero.
func (el *Ciphertext) Resize(params Parameters, degree int) {
	if el.Degree() > degree {
		el.Value = el.Value[:degree+1]
	} else if el.Degree() < degree {
		for el.Degree() < degree {
			el.Value = append(el.Value, []*ring.Poly{new(ring.Poly)}...)
			el.Value[el.Degree()].Coeffs = make([][]uint64, el.Level()+1)
			for i := 0; i < el.Level()+1; i++ {
				el.Value[el.Degree()].Coeffs[i] = make([]uint64, params.N())
				el.Value[el.Degree()].IsNTT = el.Value[0].IsNTT
			}
		}
	}
}

// SwitchCiphertextRingDegreeNTT changes the ring degree of ctIn to the one of ctOut.
// Maps Y^{N/n} -> X^{N} or X^{N} -> Y^{N/n}.
// If the ring degree of ctOut is larger than the one of ctIn, then the ringQ of ctIn
// must be provided (else a nil pointer).
// The ctIn must be in the NTT domain and ctOut will be in the NTT domain.
func SwitchCiphertextRingDegreeNTT(ctIn *Ciphertext, ringQSmallDim, ringQLargeDim *ring.Ring, ctOut *Ciphertext) {

	NIn, NOut := len(ctIn.Value[0].Coeffs[0]), len(ctOut.Value[0].Coeffs[0])

	if NIn > NOut {
		gap := NIn / NOut
		buff := make([]uint64, NIn)
		for i := range ctOut.Value {
			for j := range ctOut.Value[i].Coeffs {
				tmpIn, tmpOut := ctIn.Value[i].Coeffs[j], ctIn.Value[i].Coeffs[j]
				ringQLargeDim.InvNTTSingle(j, tmpIn, buff)
				for w0, w1 := 0, 0; w0 < NOut; w0, w1 = w0+1, w1+gap {
					tmpOut[w0] = buff[w1]
				}
				ringQSmallDim.NTTSingle(j, tmpOut, tmpOut)
			}
		}
	} else {
		for i := range ctOut.Value {
			ring.MapSmallDimensionToLargerDimensionNTT(ctIn.Value[i], ctOut.Value[i])
		}
	}
}

// SwitchCiphertextRingDegree changes the ring degree of ctIn to the one of ctOut.
// Maps Y^{N/n} -> X^{N} or X^{N} -> Y^{N/n}.
// If the ring degree of ctOut is larger than the one of ctIn, then the ringQ of ctIn
// must be provided (else a nil pointer).
func SwitchCiphertextRingDegree(ctIn *Ciphertext, ctOut *Ciphertext) {

	NIn, NOut := len(ctIn.Value[0].Coeffs[0]), len(ctOut.Value[0].Coeffs[0])

	gapIn, gapOut := NOut/NIn, 1
	if NIn > NOut {
		gapIn, gapOut = 1, NIn/NOut
	}

	for i := range ctOut.Value {
		for j := range ctOut.Value[i].Coeffs {
			tmp0, tmp1 := ctOut.Value[i].Coeffs[j], ctIn.Value[i].Coeffs[j]
			for w0, w1 := 0, 0; w0 < NOut; w0, w1 = w0+gapIn, w1+gapOut {
				tmp0[w0] = tmp1[w1]
			}
		}
	}
}

// CopyNew creates a new element as a copy of the target element.
func (el *Ciphertext) CopyNew() *Ciphertext {

	ctxCopy := new(Ciphertext)

	ctxCopy.Value = make([]*ring.Poly, el.Degree()+1)
	for i := range el.Value {
		ctxCopy.Value[i] = el.Value[i].CopyNew()
	}

	return ctxCopy
}

// Copy copies the input element and its parameters on the target element.
func (el *Ciphertext) Copy(ctxCopy *Ciphertext) {

	if el != ctxCopy {
		for i := range ctxCopy.Value {
			el.Value[i].Copy(ctxCopy.Value[i])
		}
	}
}

// El returns a pointer to this Element
func (el *Ciphertext) El() *Ciphertext {
	return el
}

// RLWEElement returns a pointer to this Element
func (el *Ciphertext) RLWEElement() *Ciphertext {
	return el
}

// GetSmallestLargest returns the provided element that has the smallest degree as a first
// returned value and the largest degree as second return value. If the degree match, the
// order is the same as for the input.
func GetSmallestLargest(el0, el1 *Ciphertext) (smallest, largest *Ciphertext, sameDegree bool) {
	switch {
	case el0.Degree() > el1.Degree():
		return el1, el0, false
	case el0.Degree() < el1.Degree():
		return el0, el1, false
	}
	return el0, el1, true
}

// PopulateElementRandom creates a new rlwe.Element with random coefficients
func PopulateElementRandom(prng utils.PRNG, params Parameters, el *Ciphertext) {
	sampler := ring.NewUniformSampler(prng, params.RingQ())
	for i := range el.Value {
		sampler.Read(el.Value[i])
	}
}
