package ckks

import (
	"github.com/ldsec/lattigo/v2/rlwe"
)

// Element is a generic type for ciphertext and plaintexts
type Element struct {
	rlwe.Element
	scale float64
}

func NewElement(params *Parameters, degree, level int, scale float64) *Element {
	return &Element{*rlwe.NewElementAtLevel(params.RLWEParameters(), degree, level), scale}
}

func (el *Element) El() *Element {
	return el
}

// Level returns the level of the target element.
func (el *Element) Level() int {
	return len(el.Value[0].Coeffs) - 1
}

// Scale returns the scale of the target element.
func (el *Element) Scale() float64 {
	return el.scale
}

// IsNTT returns true if the underlying rlwe.Element is in the NTT domain.
func (el *Element) IsNTT() bool {
	return el.Element.IsNTT
}

// SetScale sets the scale of the the target element to the input scale.
func (el *Element) SetScale(scale float64) {
	el.scale = scale
}

// MulScale multiplies the scale of the target element with the input scale.
func (el *Element) MulScale(scale float64) {
	el.scale *= scale
}

// DivScale divides the scale of the target element by the input scale.
func (el *Element) DivScale(scale float64) {
	el.scale /= scale
}

// Resize resizes the degree of the target element.
func (el *Element) Resize(params *Parameters, degree int) {
	el.Element.Resize(params.RLWEParameters(), degree)
}

func (el *Element) Copy(other *Element) {
	el.Element.Copy(&other.Element)
	el.scale = el.scale
}

func (el *Element) CopyNew() *Element {
	return &Element{*el.Element.CopyNew(), el.scale}
}
