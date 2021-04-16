package rlwe

import (
	"fmt"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// Element is a generic type for ciphertext and plaintexts
type Element struct {
	Value []*ring.Poly
	isNTT bool
}

// NewElement returns a new Element with zero values.
func NewElement(params *Parameters, degree uint64) *Element {
	el := new(Element)
	el.Value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		el.Value[i] = ring.NewPoly(params.N(), params.QiCount())
	}
	return el
}

// SetValue sets the input slice of polynomials as the value of the target element.
func (el *Element) SetValue(value []*ring.Poly) {
	el.Value = value
}

// Degree returns the degree of the target element.
func (el *Element) Degree() uint64 {
	return uint64(len(el.Value) - 1)
}

// Level returns the level of the target element.
func (el *Element) Level() uint64 {
	return uint64(len(el.Value[0].Coeffs) - 1)
}

// Resize resizes the degree of the target element.
func (el *Element) Resize(params *Parameters, degree uint64) {
	if el.Degree() > degree {
		el.Value = el.Value[:degree+1]
	} else if el.Degree() < degree {
		for el.Degree() < degree {
			el.Value = append(el.Value, []*ring.Poly{new(ring.Poly)}...)
			el.Value[el.Degree()].Coeffs = make([][]uint64, el.Level()+1)
			for i := uint64(0); i < el.Level()+1; i++ {
				el.Value[el.Degree()].Coeffs[i] = make([]uint64, params.N())
			}
		}
	}
}

// IsNTT returns the value of the NTT flag of the target element.
func (el *Element) IsNTT() bool {
	return el.isNTT
}

// SetIsNTT sets the value of the NTT flag of the target element with the input value.
func (el *Element) SetIsNTT(value bool) {
	el.isNTT = value
}

// NTT puts the target element in the NTT domain and sets its isNTT flag to true. If it is already in the NTT domain, it does nothing.
func (el *Element) NTT(ringQ *ring.Ring, c *Element) {
	if el.Degree() != c.Degree() {
		panic(fmt.Errorf("error: receiver element has invalid degree (it does not match)"))
	}
	if !el.IsNTT() {
		for i := range el.Value {
			ringQ.NTTLvl(el.Level(), el.Value[i], c.Value[i])
		}
		c.SetIsNTT(true)
	}
}

// InvNTT puts the target element outside of the NTT domain, and sets its isNTT flag to false. If it is not in the NTT domain, it does nothing.
func (el *Element) InvNTT(ringQ *ring.Ring, c *Element) {
	if el.Degree() != c.Degree() {
		panic(fmt.Errorf("error: receiver element invalid degree (it does not match)"))
	}
	if el.IsNTT() {
		for i := range el.Value {
			ringQ.InvNTTLvl(el.Level(), el.Value[i], c.Value[i])
		}
		c.SetIsNTT(false)
	}
}

// CopyNew creates a new element as a copy of the target element.
func (el *Element) CopyNew() *Element {

	ctxCopy := new(Element)

	ctxCopy.Value = make([]*ring.Poly, el.Degree()+1)
	for i := range el.Value {
		ctxCopy.Value[i] = el.Value[i].CopyNew()
	}

	ctxCopy.CopyParams(el)

	return ctxCopy
}

// Copy copies the input element and its parameters on the target element.
func (el *Element) Copy(ctxCopy *Element) {

	if el != ctxCopy {
		for i := range ctxCopy.Value {
			el.Value[i].Copy(ctxCopy.Value[i])
		}

		el.CopyParams(ctxCopy)
	}
}

// CopyParams copies the input element parameters on the target element
func (el *Element) CopyParams(Element *Element) {
	el.SetIsNTT(Element.IsNTT())
}

// El sets the target element type to Element.
func (el *Element) El() *Element {
	return el
}

func GetSmallestLargest(el0, el1 *Element) (smallest, largest *Element, sameDegree bool) {
	switch {
	case el0.Degree() > el1.Degree():
		return el1, el0, false
	case el0.Degree() < el1.Degree():
		return el0, el1, false
	}
	return el0, el1, true
}

// NewElementRandom creates a new rlwe.Element with random coefficients
func PopulateElementRandom(prng utils.PRNG, params *Parameters, el *Element) {

	ringQ, err := ring.NewRing(params.N(), params.qi)
	if err != nil {
		panic(err)
	}
	sampler := ring.NewUniformSampler(prng, ringQ)
	for i := range el.Value {
		sampler.Read(el.Value[i])
	}
}
