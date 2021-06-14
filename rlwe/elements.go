package rlwe

import (
	"fmt"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// Ciphertext is a common interface for RLWE ciphertexts.
type Ciphertext interface {
	RLWEElement() *Element
}

// Plaintext is a common base type for RLWE plaintexts.
type Plaintext struct {
	Value *ring.Poly
	IsNTT bool
}

// AdditiveShare is a type for storing additively shared values in the
// polynomial ring.
type AdditiveShare struct {
	Value ring.Poly
}

// NewAdditiveShare instantiate a new additive share struct for the ring defined
// by the given parameters at maximum level.
func NewAdditiveShare(params Parameters) AdditiveShare {
	return AdditiveShare{Value: *ring.NewPoly(params.N(), 1)}
}

// NewAdditiveShareAtLevel instantiate a new additive share struct for the ring defined
// by the given parameters at level `level`.
func NewAdditiveShareAtLevel(params Parameters, level int) AdditiveShare {
	return AdditiveShare{Value: *ring.NewPoly(params.N(), level+1)}
}

// NewPlaintext creates a new Plaintext at maximum level from the parameters.
func NewPlaintext(params Parameters) *Plaintext {
	return &Plaintext{Value: ring.NewPoly(params.N(), params.QCount())}
}

// NewPlaintextAtLevel creates a new Plaintext at level `level` from the parameters.
func NewPlaintextAtLevel(params Parameters, level int) *Plaintext {
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
func (pt Plaintext) El() *Element {
	return &Element{Value: []*ring.Poly{pt.Value}, IsNTT: true}
}

// Copy copies the `other` plaintext value into the reciever plaintext.
func (pt *Plaintext) Copy(other *Plaintext) {
	if other != nil && other.Value != nil {
		pt.Value.Copy(other.Value)
	}
}

// Element is a generic type for ciphertext and plaintexts
type Element struct {
	Value []*ring.Poly
	IsNTT bool
}

// NewElement returns a new Element with zero values.
func NewElement(params Parameters, degree int) *Element {
	return NewElementAtLevel(params, degree, params.QCount()-1)
}

// NewElementAtLevel returns a new Element with zero values.
func NewElementAtLevel(params Parameters, degree, level int) *Element {
	el := new(Element)
	el.Value = make([]*ring.Poly, degree+1)
	for i := 0; i < degree+1; i++ {
		el.Value[i] = ring.NewPoly(params.N(), level+1)
	}
	return el
}

// SetValue sets the input slice of polynomials as the value of the target element.
func (el *Element) SetValue(value []*ring.Poly) {
	el.Value = value
}

// Degree returns the degree of the target element.
func (el *Element) Degree() int {
	return len(el.Value) - 1
}

// Level returns the level of the target element.
func (el *Element) Level() int {
	return len(el.Value[0].Coeffs) - 1
}

// Resize resizes the degree of the target element.
func (el *Element) Resize(params Parameters, degree int) {
	if el.Degree() > degree {
		el.Value = el.Value[:degree+1]
	} else if el.Degree() < degree {
		for el.Degree() < degree {
			el.Value = append(el.Value, []*ring.Poly{new(ring.Poly)}...)
			el.Value[el.Degree()].Coeffs = make([][]uint64, el.Level()+1)
			for i := 0; i < el.Level()+1; i++ {
				el.Value[el.Degree()].Coeffs[i] = make([]uint64, params.N())
			}
		}
	}
}

// NTT puts the target element in the NTT domain and sets its isNTT flag to true. If it is already in the NTT domain, it does nothing.
func (el *Element) NTT(ringQ *ring.Ring, c *Element) {
	if el.Degree() != c.Degree() {
		panic(fmt.Errorf("error: receiver element has invalid degree (it does not match)"))
	}
	if !el.IsNTT {
		for i := range el.Value {
			ringQ.NTTLvl(el.Level(), el.Value[i], c.Value[i])
		}
		c.IsNTT = true
	}
}

// InvNTT puts the target element outside of the NTT domain, and sets its isNTT flag to false. If it is not in the NTT domain, it does nothing.
func (el *Element) InvNTT(ringQ *ring.Ring, c *Element) {
	if el.Degree() != c.Degree() {
		panic(fmt.Errorf("error: receiver element invalid degree (it does not match)"))
	}
	if el.IsNTT {
		for i := range el.Value {
			ringQ.InvNTTLvl(el.Level(), el.Value[i], c.Value[i])
		}
		c.IsNTT = false
	}
}

// CopyNew creates a new element as a copy of the target element.
func (el *Element) CopyNew() *Element {

	ctxCopy := new(Element)

	ctxCopy.Value = make([]*ring.Poly, el.Degree()+1)
	for i := range el.Value {
		ctxCopy.Value[i] = el.Value[i].CopyNew()
	}

	ctxCopy.IsNTT = el.IsNTT

	return ctxCopy
}

// Copy copies the input element and its parameters on the target element.
func (el *Element) Copy(ctxCopy *Element) {

	if el != ctxCopy {
		for i := range ctxCopy.Value {
			el.Value[i].Copy(ctxCopy.Value[i])
		}

		el.IsNTT = ctxCopy.IsNTT
	}
}

// El returns a pointer to this Element
func (el *Element) El() *Element {
	return el
}

// RLWEElement returns a pointer to this Element
func (el *Element) RLWEElement() *Element {
	return el
}

// GetSmallestLargest returns the provided element that has the smallest degree as a first
// returned value and the largest degree as second return value. If the degree match, the
// order is the same as for the input.
func GetSmallestLargest(el0, el1 *Element) (smallest, largest *Element, sameDegree bool) {
	switch {
	case el0.Degree() > el1.Degree():
		return el1, el0, false
	case el0.Degree() < el1.Degree():
		return el0, el1, false
	}
	return el0, el1, true
}

// PopulateElementRandom creates a new rlwe.Element with random coefficients
func PopulateElementRandom(prng utils.PRNG, params Parameters, el *Element) {

	ringQ, err := ring.NewRing(params.N(), params.Q())
	if err != nil {
		panic(err)
	}
	sampler := ring.NewUniformSampler(prng, ringQ)
	for i := range el.Value {
		sampler.Read(el.Value[i])
	}
}
