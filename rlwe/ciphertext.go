package rlwe

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// Ciphertext is a generic type for RLWE ciphertexts.
type Ciphertext struct {
	Scale
	Value []*ring.Poly
}

// NewCiphertext returns a new Element with zero values.
func NewCiphertext(params Parameters, degree, level int) *Ciphertext {
	el := new(Ciphertext)
	el.Value = make([]*ring.Poly, degree+1)
	for i := 0; i < degree+1; i++ {
		el.Value[i] = ring.NewPoly(params.N(), level)
	}
	return el
}

// NewCiphertextNTT returns a new Element with zero values and the NTT flags set.
func NewCiphertextNTT(params Parameters, degree, level int) *Ciphertext {
	el := new(Ciphertext)
	el.Value = make([]*ring.Poly, degree+1)
	for i := 0; i < degree+1; i++ {
		el.Value[i] = ring.NewPoly(params.N(), level)
		el.Value[i].IsNTT = true
	}
	return el
}

// NewCiphertextAtLevelFromPoly constructs a new Ciphetext at a specific level
// where the message is set to the passed poly. No checks are performed on poly and
// the returned Ciphertext will share its backing array of coefficient.
func NewCiphertextAtLevelFromPoly(level int, poly [2]*ring.Poly) *Ciphertext {
	v0, v1 := new(ring.Poly), new(ring.Poly)
	v0.Coeffs, v1.Coeffs = poly[0].Coeffs[:level+1], poly[1].Coeffs[:level+1]
	v0.Buff, v1.Buff = poly[0].Buff[:poly[0].N()*(level+1)], poly[1].Buff[:poly[1].N()*(level+1)]
	return &Ciphertext{Value: []*ring.Poly{v0, v1}}
}

// NewCiphertextRandom generates a new uniformly distributed Ciphertext of degree, level.
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
func (el *Ciphertext) Resize(degree, level int) {

	if el.Level() != level {
		for i := range el.Value {
			el.Value[i].Resize(level)
		}
	}

	if el.Degree() > degree {
		el.Value = el.Value[:degree+1]
	} else if el.Degree() < degree {
		for el.Degree() < degree {
			el.Value = append(el.Value, []*ring.Poly{ring.NewPoly(el.Value[0].N(), level)}...)
			el.Value[el.Degree()].IsNTT = el.Value[0].IsNTT
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
				tmpIn, tmpOut := ctIn.Value[i].Coeffs[j], ctOut.Value[i].Coeffs[j]
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

	ctxCopy.Scale = el.Scale.CopyNew()

	return ctxCopy
}

// Copy copies the input element and its parameters on the target element.
func (el *Ciphertext) Copy(ctxCopy *Ciphertext) {
	if el != ctxCopy {
		for i := range ctxCopy.Value {
			el.Value[i].Copy(ctxCopy.Value[i])
		}

		el.Scale = ctxCopy.Scale.CopyNew()
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
