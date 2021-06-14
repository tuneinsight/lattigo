package ckks

import (
	"github.com/ldsec/lattigo/v2/rlwe"
)

// Plaintext is is a Element with only one Poly.
type Plaintext struct {
	*rlwe.Plaintext
	Scale float64
}

// NewPlaintext creates a new Plaintext of level level and scale scale.
func NewPlaintext(params Parameters, level int, scale float64) *Plaintext {
	pt := &Plaintext{Plaintext: rlwe.NewPlaintextAtLevel(params.Parameters, level), Scale: scale}
	pt.IsNTT = true
	return pt
}

func (p *Plaintext) ScalingFactor() float64 {
	return p.Scale
}

func (p *Plaintext) SetScalingFactor(scale float64) {
	p.Scale = scale
}
