package ckks

import (
	"github.com/ldsec/lattigo/v2/rlwe"
)

// Element is a generic type for ciphertext and plaintexts
type Element struct {
	rlwe.Element
	scale float64
}

func NewElement(params *Parameters, degree, level uint64, scale float64) *Element {
	return &Element{*rlwe.NewElementAtLevel(params.RLWEParameters(), degree, level), scale}
}

func (el *Element) El() *Element {
	return el
}

// Scale returns the scale of the target element.
func (el *Element) Scale() float64 {
	return el.scale
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
func (el *Element) Resize(params *Parameters, degree uint64) {
	el.Element.Resize(params.RLWEParameters(), degree)
}

func (el *Element) Copy(other *Element) {
	el.Element.Copy(&other.Element)
	el.scale = other.scale
}

func (el *Element) CopyNew() *Element {
	return &Element{*el.Element.CopyNew(), el.scale}
}
