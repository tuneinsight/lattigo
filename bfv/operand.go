package bfv

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// Operand is a common interface for Ciphertext and Plaintext.
type Operand interface {
	El() *Element
	Degree() int
}

// Element is a common struct for Plaintexts and Ciphertexts. It stores a value
// as a slice of polynomials, and an isNTT flag that indicates if the element is in the NTT domain.
type Element struct {
	value []*ring.Poly
}

func getSmallestLargest(el0, el1 *Element) (smallest, largest *Element, sameDegree bool) {
	switch {
	case el0.Degree() > el1.Degree():
		return el1, el0, false
	case el0.Degree() < el1.Degree():
		return el0, el1, false
	}
	return el0, el1, true
}

func newCiphertextElement(params *Parameters, degree int) *Element {
	el := new(Element)
	el.value = make([]*ring.Poly, degree+1)
	for i := 0; i < degree+1; i++ {
		el.value[i] = ring.NewPoly(params.N(), params.QiCount())
	}
	return el
}

func newPlaintextElement(params *Parameters) *Element {
	el := new(Element)
	el.value = []*ring.Poly{ring.NewPoly(params.N(), params.QiCount())}
	return el
}

func newPlaintextRingTElement(params *Parameters) *Element {
	el := new(Element)
	el.value = []*ring.Poly{ring.NewPoly(params.N(), 1)}
	return el
}

func newPlaintextMulElement(params *Parameters) *Element {
	el := new(Element)
	el.value = []*ring.Poly{ring.NewPoly(params.N(), params.QiCount())}
	return el
}

// NewElementRandom creates a new Element with random coefficients
func populateElementRandom(prng utils.PRNG, params *Parameters, el *Element) {

	ringQ, err := ring.NewRing(params.N(), params.qi)
	if err != nil {
		panic(err)
	}
	sampler := ring.NewUniformSampler(prng, ringQ)
	for i := range el.value {
		sampler.Read(el.value[i])
	}
}

// Value returns the value of the target Element (as a slice of polynomials in CRT form).
func (el *Element) Value() []*ring.Poly {
	return el.value
}

// SetValue assigns the input slice of polynomials to the target Element value.
func (el *Element) SetValue(value []*ring.Poly) {
	el.value = value
}

// Degree returns the degree of the target Element.
func (el *Element) Degree() int {
	return len(el.value) - 1
}

// Level returns the level of the target element.
func (el *Element) Level() int {
	return len(el.value[0].Coeffs) - 1
}

// Resize resizes the degree of the target element.
func (el *Element) Resize(params *Parameters, degree int) {
	if el.Degree() > degree {
		el.value = el.value[:degree+1]
	} else if el.Degree() < degree {
		for el.Degree() < degree {
			el.value = append(el.value, []*ring.Poly{new(ring.Poly)}...)
			el.value[el.Degree()].Coeffs = make([][]uint64, el.Level()+1)
			for i := 0; i < el.Level()+1; i++ {
				el.value[el.Degree()].Coeffs[i] = make([]uint64, params.N())
			}
		}
	}
}

// CopyNew creates a new Element which is a copy of the target Element, and returns the value as
// a Element.
func (el *Element) CopyNew() *Element {

	ctxCopy := new(Element)

	ctxCopy.value = make([]*ring.Poly, el.Degree()+1)
	for i := range el.value {
		ctxCopy.value[i] = el.value[i].CopyNew()
	}

	return ctxCopy
}

// Copy copies the value and parameters of the input on the target Element.
func (el *Element) Copy(ctxCopy *Element) {
	if el != ctxCopy {
		for i := range ctxCopy.Value() {
			el.Value()[i].Copy(ctxCopy.Value()[i])
		}
	}
}

// El sets the target Element type to Element.
func (el *Element) El() *Element {
	return el
}

// Ciphertext sets the target Element type to Ciphertext.
func (el *Element) Ciphertext() *Ciphertext {
	return &Ciphertext{el}
}

// Plaintext sets the target Element type to Plaintext.
func (el *Element) Plaintext() *Plaintext {
	return &Plaintext{el, el.value[0]}
}
