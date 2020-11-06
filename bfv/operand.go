package bfv

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// Operand is a common interface for Ciphertext and Plaintext.
type Operand interface {
	El() *Element
	Degree() uint64
}

// Element is a common struct for Plaintexts and Ciphertexts. It stores a value
// as a slice of polynomials, and an isNTT flag that indicates if the element is in the NTT domain.
type Element struct {
	value []*ring.Poly
	isNTT bool // If true, then in the NTT domain, else not
	inZQ  bool // If true, then in Z_{Q}, else in Z_{t}
}

// NewElement creates a new Element of the target degree with zero values.
func NewElement(params *Parameters, degree uint64, inZQ bool) *Element {

	el := new(Element)
	el.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		if inZQ {
			el.value[i] = ring.NewPoly(params.N(), params.QiCount())
		} else {
			el.value[i] = ring.NewPoly(params.N(), 1)
		}
	}
	el.isNTT = true
	el.inZQ = inZQ
	return el
}

// NewElementRandom creates a new Element with random coefficients
func NewElementRandom(prng utils.PRNG, params *Parameters, degree uint64, inZQ bool) *Element {

	ringQ, err := ring.NewRing(params.N(), params.qi)
	if err != nil {
		panic(err)
	}
	sampler := ring.NewUniformSampler(prng, ringQ)
	el := NewElement(params, degree, inZQ)
	for i := uint64(0); i < degree+1; i++ {
		sampler.Read(el.value[i])
	}
	return el
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
func (el *Element) Degree() uint64 {
	return uint64(len(el.value) - 1)
}

// Resize resizes the target Element degree to the degree given as input. If the input degree is bigger, then
// it will append new empty polynomials; if the degree is smaller, it will delete polynomials until the degree matches
// the input degree.
func (el *Element) Resize(params *Parameters, degree uint64) {
	if el.Degree() > degree {
		el.value = el.value[:degree+1]
	} else if el.Degree() < degree {
		for el.Degree() < degree {
			el.value = append(el.value, []*ring.Poly{new(ring.Poly)}...)
			el.value[el.Degree()].Coeffs = make([][]uint64, len(params.qi))
			for i := 0; i < len(params.qi); i++ {
				el.value[el.Degree()].Coeffs[i] = make([]uint64, uint64(1<<params.logN))
			}
		}
	}
}

// IsNTT returns true if the target Element is in the NTT domain, and false otherwise.
func (el *Element) IsNTT() bool {
	return el.isNTT
}

// SetIsNTT assigns the input Boolean value to the isNTT flag of the target Element.
func (el *Element) SetIsNTT(value bool) {
	el.isNTT = value
}

// CopyNew creates a new Element which is a copy of the target Element, and returns the value as
// a Element.
func (el *Element) CopyNew() *Element {

	ctxCopy := new(Element)

	ctxCopy.value = make([]*ring.Poly, el.Degree()+1)
	for i := range el.value {
		ctxCopy.value[i] = el.value[i].CopyNew()
	}
	ctxCopy.isNTT = el.isNTT

	return ctxCopy
}

// Copy copies the value and parameters of the input on the target Element.
func (el *Element) Copy(ctxCopy *Element) {
	if el != ctxCopy {
		for i := range ctxCopy.Value() {
			el.Value()[i].Copy(ctxCopy.Value()[i])
		}
		el.isNTT = ctxCopy.isNTT
	}
}

// NTT puts the target Element in the NTT domain and sets its isNTT flag to true. If it is already in the NTT domain, does nothing.
func (el *Element) NTT(ringQ *ring.Ring, c *Element) {
	if el.Degree() != c.Degree() {
		panic("cannot NTT: receiver element invalid degree (degrees do not match)")
	}
	if el.IsNTT() != true {
		for i := range el.value {
			ringQ.NTT(el.Value()[i], c.Value()[i])
		}
		c.SetIsNTT(true)
	}
}

// InvNTT puts the target Element outside of the NTT domain, and sets its isNTT flag to false. If it is not in the NTT domain, it does nothing.
func (el *Element) InvNTT(ringQ *ring.Ring, c *Element) {
	if el.Degree() != c.Degree() {
		panic("cannot InvNTT: receiver element invalid degree (degrees do not match)")
	}
	if el.IsNTT() != false {
		for i := range el.value {
			ringQ.InvNTT(el.Value()[i], c.Value()[i])
		}
		c.SetIsNTT(false)
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
