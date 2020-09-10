package ckks

import (
	"errors"

	"github.com/ldsec/lattigo/ring"
)

// Operand is a common interface for Ciphertext and Plaintext types.
type Operand interface {
	Element() *CkksElement
	Degree() uint64
	Level() uint64
	Scale() float64
}

// CkksElement is a generic type for ciphertext and plaintexts
type CkksElement struct {
	value []*ring.Poly
	scale float64
	isNTT bool
}

// NewCkksElement returns a new CkksElement with zero values.
func NewCkksElement() *CkksElement {
	return &CkksElement{}
}

// Value returns the slice of polynomials of the target element.
func (el *CkksElement) Value() []*ring.Poly {
	return el.value
}

// SetValue sets the input slice of polynomials as the value of the target element.
func (el *CkksElement) SetValue(value []*ring.Poly) {
	el.value = value
}

// Degree returns the degree of the target element.
func (el *CkksElement) Degree() uint64 {
	return uint64(len(el.value) - 1)
}

// Level returns the level of the target element.
func (el *CkksElement) Level() uint64 {
	return uint64(len(el.value[0].Coeffs) - 1)
}

// Scale returns the scale of the target element.
func (el *CkksElement) Scale() float64 {
	return el.scale
}

// SetScale sets the scale of the the target element to the input scale.
func (el *CkksElement) SetScale(scale float64) {
	el.scale = scale
}

// MulScale multiplies the scale of the target element with the input scale.
func (el *CkksElement) MulScale(scale float64) {
	el.scale *= scale
}

// DivScale divides the scale of the target element by the input scale.
func (el *CkksElement) DivScale(scale float64) {
	el.scale /= scale
}

// Resize resizes the degree of the target element.
func (el *CkksElement) Resize(params *Parameters, degree uint64) {
	if el.Degree() > degree {
		el.value = el.value[:degree+1]
	} else if el.Degree() < degree {
		for el.Degree() < degree {
			el.value = append(el.value, []*ring.Poly{new(ring.Poly)}...)
			el.value[el.Degree()].Coeffs = make([][]uint64, el.Level()+1)
			for i := uint64(0); i < el.Level()+1; i++ {
				el.value[el.Degree()].Coeffs[i] = make([]uint64, params.N())
			}
		}
	}
}

// IsNTT returns the value of the NTT flag of the target element.
func (el *CkksElement) IsNTT() bool {
	return el.isNTT
}

// SetIsNTT sets the value of the NTT flag of the target element with the input value.
func (el *CkksElement) SetIsNTT(value bool) {
	el.isNTT = value
}

// NTT puts the target element in the NTT domain and sets its isNTT flag to true. If it is already in the NTT domain, it does nothing.
func (el *CkksElement) NTT(ringQ *ring.Ring, c *CkksElement) error {
	if el.Degree() != c.Degree() {
		return errors.New("error: receiver element has invalid degree (it does not match)")
	}
	if el.IsNTT() != true {
		for i := range el.value {
			ringQ.NTTLvl(el.Level(), el.Value()[i], c.Value()[i])
		}
		c.SetIsNTT(true)
	}
	return nil
}

// InvNTT puts the target element outside of the NTT domain, and sets its isNTT flag to false. If it is not in the NTT domain, it does nothing.
func (el *CkksElement) InvNTT(ringQ *ring.Ring, c *CkksElement) error {
	if el.Degree() != c.Degree() {
		return errors.New("error: receiver element invalid degree (it does not match)")
	}
	if el.IsNTT() != false {
		for i := range el.value {
			ringQ.InvNTTLvl(el.Level(), el.Value()[i], c.Value()[i])
		}
		c.SetIsNTT(false)
	}
	return nil
}

// CopyNew creates a new element as a copy of the target element.
func (el *CkksElement) CopyNew() *CkksElement {

	ctxCopy := new(CkksElement)

	ctxCopy.value = make([]*ring.Poly, el.Degree()+1)
	for i := range el.value {
		ctxCopy.value[i] = el.value[i].CopyNew()
	}

	ctxCopy.CopyParams(el)

	return ctxCopy
}

// Copy copies the input element and its parameters on the target element.
func (el *CkksElement) Copy(ctxCopy *CkksElement) (err error) {

	if el != ctxCopy {
		for i := range ctxCopy.Value() {
			el.value[i].Copy(ctxCopy.Value()[i])
		}

		el.CopyParams(ctxCopy)
	}
	return nil
}

// CopyParams copies the input element parameters on the target element
func (el *CkksElement) CopyParams(CkksElement *CkksElement) {
	el.SetScale(CkksElement.Scale())
	el.SetIsNTT(CkksElement.IsNTT())
}

// Element sets the target element type to Element.
func (el *CkksElement) Element() *CkksElement {
	return el
}

// Ciphertext sets the target element type to Ciphertext.
func (el *CkksElement) Ciphertext() *Ciphertext {
	return &Ciphertext{el}
}

// Plaintext sets the target element type to Plaintext.
func (el *CkksElement) Plaintext() *Plaintext {
	return &Plaintext{el, el.value[0]}
}
