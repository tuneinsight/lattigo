package bfv

import (
	"github.com/tuneinsight/lattigo/ring"
)

// Operand is a common interface for Ciphertext and Plaintext.
type Operand interface {
	Element() *bfvElement
	Degree() uint64
}

// bfvElement is a common struct for Plaintexts and Ciphertexts. It stores a value
// as a slice of polynomials, and an isNTT flag that indicates if the element is in the NTT domain.
type bfvElement struct {
	value []*ring.Poly
	isNTT bool
}

// newBfvElement creates a new bfvElement of the target degree with zero values.
func newBfvElement(params *Parameters, degree uint64) *bfvElement {

	if !params.isValid {
		panic("cannot newBfvElement: params not valid (check if they were generated properly)")
	}

	el := new(bfvElement)
	el.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		el.value[i] = ring.NewPoly(1<<params.LogN, uint64(len(params.LogQi)))
	}
	el.isNTT = true
	return el
}

func newBfvElementRandom(params *Parameters, degree uint64) *bfvElement {

	if !params.isValid {
		panic("cannot newBfvElementRandom: params not valid (check if they were generated properly)")
	}

	el := new(bfvElement)
	el.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		el.value[i] = ring.NewPolyUniform(1<<params.LogN, uint64(len(params.LogQi)))
	}
	el.isNTT = true
	return el
}

// Value returns the value of the target bfvElement (as a slice of polynomials in CRT form).
func (el *bfvElement) Value() []*ring.Poly {
	return el.value
}

// SetValue assigns the input slice of polynomials to the target bfvElement value.
func (el *bfvElement) SetValue(value []*ring.Poly) {
	el.value = value
}

// Degree returns the degree of the target bfvElement.
func (el *bfvElement) Degree() uint64 {
	return uint64(len(el.value) - 1)
}

// Resize resizes the target bfvElement degree to the degree given as input. If the input degree is bigger, then
// it will append new empty polynomials; if the degree is smaller, it will delete polynomials until the degree matches
// the input degree.
func (el *bfvElement) Resize(params *Parameters, degree uint64) {
	if el.Degree() > degree {
		el.value = el.value[:degree+1]
	} else if el.Degree() < degree {
		for el.Degree() < degree {
			el.value = append(el.value, []*ring.Poly{new(ring.Poly)}...)
			el.value[el.Degree()].Coeffs = make([][]uint64, len(params.LogQi))
			for i := 0; i < len(params.LogQi); i++ {
				el.value[el.Degree()].Coeffs[i] = make([]uint64, uint64(1<<params.LogN))
			}
		}
	}
}

// IsNTT returns true if the target bfvElement is in the NTT domain, and false otherwise.
func (el *bfvElement) IsNTT() bool {
	return el.isNTT
}

// SetIsNTT assigns the input Boolean value to the isNTT flag of the target bfvElement.
func (el *bfvElement) SetIsNTT(value bool) {
	el.isNTT = value
}

// CopyNew creates a new bfvElement which is a copy of the target bfvElement, and returns the value as
// a bfvElement.
func (el *bfvElement) CopyNew() *bfvElement {

	ctxCopy := new(bfvElement)

	ctxCopy.value = make([]*ring.Poly, el.Degree()+1)
	for i := range el.value {
		ctxCopy.value[i] = el.value[i].CopyNew()
	}
	ctxCopy.isNTT = el.isNTT

	return ctxCopy
}

// Copy copies the value and parameters of the input on the target bfvElement.
func (el *bfvElement) Copy(ctxCopy *bfvElement) {
	if el != ctxCopy {
		for i := range ctxCopy.Value() {
			el.Value()[i].Copy(ctxCopy.Value()[i])
		}
		el.isNTT = ctxCopy.isNTT
	}
}

// NTT puts the target bfvElement in the NTT domain and sets its isNTT flag to true. If it is already in the NTT domain, does nothing.
func (el *bfvElement) NTT(context *ring.Context, c *bfvElement) {
	if el.Degree() != c.Degree() {
		panic("cannot NTT: receiver element invalid degree (degrees do not match)")
	}
	if el.IsNTT() != true {
		for i := range el.value {
			context.NTT(el.Value()[i], c.Value()[i])
		}
		c.SetIsNTT(true)
	}
}

// InvNTT puts the target bfvElement outside of the NTT domain, and sets its isNTT flag to false. If it is not in the NTT domain, it does nothing.
func (el *bfvElement) InvNTT(context *ring.Context, c *bfvElement) {
	if el.Degree() != c.Degree() {
		panic("cannot InvNTT: receiver element invalid degree (degrees do not match)")
	}
	if el.IsNTT() != false {
		for i := range el.value {
			context.InvNTT(el.Value()[i], c.Value()[i])
		}
		c.SetIsNTT(false)
	}
}

func (el *bfvElement) Element() *bfvElement {
	return el
}

func (el *bfvElement) Ciphertext() *Ciphertext {
	return &Ciphertext{el}
}

func (el *bfvElement) Plaintext() *Plaintext {
	return &Plaintext{el, el.value[0]}
}
