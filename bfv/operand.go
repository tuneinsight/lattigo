package bfv

import (
	"errors"
	"github.com/ldsec/lattigo/ring"
)

type Operand interface {
	Element() *bfvElement
	Degree() uint64
}

// bfvElement is a common struct between plaintexts and ciphertexts. It stores a value
// as a slice of polynomials, and an isNTT flag indicatig if the element is in the NTT domain.
type bfvElement struct {
	value []*ring.Poly
	isNTT bool
}

// NewCiphertext creates a new empty ciphertext of degree degree.
func (bfvcontext *BfvContext) NewBfvElement(degree uint64) *bfvElement {
	el := &bfvElement{}
	el.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		el.value[i] = bfvcontext.contextQ.NewPoly()
	}
	el.isNTT = false

	return el
}

// Value returns the value of the target ciphertext (as a slice of polynomials in CRT form).
func (el *bfvElement) Value() []*ring.Poly {
	return el.value
}

// SetValue assigns the input slice of polynomials to the target ciphertext value.
func (el *bfvElement) SetValue(value []*ring.Poly) {
	el.value = value
}

// Degree returns the degree of the target ciphertext.
func (el *bfvElement) Degree() uint64 {
	return uint64(len(el.value) - 1)
}

// Resize resizes the target ciphertext degree to the degree given as input. If the input degree is bigger then
// it will append new empty polynomials, if the degree is smaller, it will delete polynomials until the degree matches
// the input degree.
func (el *bfvElement) Resize(bfvcontext *BfvContext, degree uint64) {
	if el.Degree() > degree {
		el.value = el.value[:degree]
	} else if el.Degree() < degree {
		for el.Degree() < degree {
			el.value = append(el.value, []*ring.Poly{bfvcontext.contextQ.NewPoly()}...)
		}
	}
}

// IsNTT returns true if the target ciphertext is in the NTT domain, else false.
func (el *bfvElement) IsNTT() bool {
	return el.isNTT
}

// SetIsNTT assigns the input bolean value to the isNTT flag of the target ciphertext.
func (el *bfvElement) SetIsNTT(value bool) {
	el.isNTT = value
}

// CopyNew creates a new ciphertext which is a copy of the target ciphertext. Returns the value as
// a Element.
func (el *bfvElement) CopyNew() *bfvElement {

	ctxCopy := new(bfvElement)

	ctxCopy.value = make([]*ring.Poly, el.Degree()+1)
	for i := range el.value {
		ctxCopy.value[i] = el.value[i].CopyNew()
	}
	ctxCopy.isNTT = el.isNTT

	return ctxCopy
}

// Copy copies the value and parameters of the input on the target ciphertext.
func (el *bfvElement) Copy(ctxCopy *bfvElement) error {

	for i := range ctxCopy.Value() {
		ctxCopy.Value()[i].Copy(el.Value()[i])
	}
	ctxCopy.SetIsNTT(el.IsNTT())

	return nil
}

// NTT puts the target ciphertext in the NTT domain and sets its isNTT flag to true. If it is already in the NTT domain, does nothing.
func (el *bfvElement) NTT(bfvcontext *BfvContext, c *bfvElement) error {
	if el.Degree() != c.Degree() {
		return errors.New("error : receiver element invalide degree (does not match)")
	}
	if el.IsNTT() != true {
		for i := range el.value {
			bfvcontext.contextQ.NTT(el.Value()[i], c.Value()[i])
		}
		c.SetIsNTT(true)
	}
	return nil
}

// InvNTT puts the target ciphertext outside of the NTT domain, and sets its isNTT flag to false. If it is not in the NTT domain, does nothing.
func (el *bfvElement) InvNTT(bfvcontext *BfvContext, c *bfvElement) error {
	if el.Degree() != c.Degree() {
		return errors.New("error : receiver element invalide degree (does not match)")
	}
	if el.IsNTT() != false {
		for i := range el.value {
			bfvcontext.contextQ.InvNTT(el.Value()[i], c.Value()[i])
		}
		c.SetIsNTT(false)
	}
	return nil
}

func (el *bfvElement) Element() *bfvElement {
	return el
}

func (el *bfvElement) Ciphertext() *Ciphertext {
	return &Ciphertext{el}
}

func (el *bfvElement) Plaintext() *Plaintext {
	if len(el.value) != 1 {
		panic("not a plaintext element")
	}
	return &Plaintext{el, el.value[0]}
}
