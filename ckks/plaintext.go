package ckks

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
)

// Plaintext is is a Element with only one Poly.
type Plaintext struct {
	*rlwe.Plaintext
	Scale float64
}

// NewPlaintext creates a new Plaintext of level level and scale scale.
func NewPlaintext(params Parameters, level int, scale float64) *Plaintext {
	pt := &Plaintext{Plaintext: rlwe.NewPlaintext(params.Parameters, level), Scale: scale}
	pt.Value.IsNTT = true
	return pt
}

// ScalingFactor returns the scaling factor of the plaintext
func (p *Plaintext) ScalingFactor() float64 {
	return p.Scale
}

// SetScalingFactor sets the scaling factor of the target plaintext
func (p *Plaintext) SetScalingFactor(scale float64) {
	p.Scale = scale
}

// NewPlaintextAtLevelFromPoly construct a new Plaintext at a specific level
// where the message is set to the passed poly. No checks are performed on poly and
// the returned Plaintext will share its backing array of coefficient.
func NewPlaintextAtLevelFromPoly(level int, poly *ring.Poly) *Plaintext {
	v0 := new(ring.Poly)
	v0.IsNTT = true
	v0.Coeffs = poly.Coeffs[:level+1]
	return &Plaintext{Plaintext: &rlwe.Plaintext{Value: v0}, Scale: 0}
}
