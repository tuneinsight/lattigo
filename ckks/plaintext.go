package ckks

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)


// Plaintext is is a Element with only one Poly.
type Plaintext struct {
	Real *rlwe.Plaintext
	Imag *rlwe.Plaintext
	IsComplex bool
	Scale float64
}

// NewPlaintext creates a new Plaintext of level level and scale scale.
func NewPlaintext(params Parameters, level int, scale float64) *Plaintext {
	pt := &Plaintext{Real: rlwe.NewPlaintext(params.Parameters, level), Imag:nil, IsComplex:true, Scale: scale}
	pt.Real.Value.IsNTT = true
	return pt
}

func (p *Plaintext) Degree() int {
	return 0
}

func (p *Plaintext) El() *rlwe.Ciphertext{
	return &rlwe.Ciphertext{Value: []*ring.Poly{p.Real.Value}}
}

func (p *Plaintext) Level() int {
	if p.Imag == nil{
		return p.Real.Level()
	}else{
		return utils.MinInt(p.Real.Level(), p.Imag.Level())
	}
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
	return &Plaintext{Real: &rlwe.Plaintext{Value: v0}, Scale: 0}
}
