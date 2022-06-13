package bgv

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
)

// Plaintext is is a Element with only one Poly.
type Plaintext struct {
	*rlwe.Plaintext
	Scale uint64
}

// NewPlaintext creates a new Plaintext of level level and scale scale.
func NewPlaintext(params Parameters, level int, scale uint64) *Plaintext {
	pt := &Plaintext{Plaintext: rlwe.NewPlaintext(params.Parameters, level), Scale: scale}
	pt.Value.IsNTT = true
	return pt
}

// NewPlaintextAtLevelFromPoly construct a new Plaintext at a specific level
// where the message is set to the passed poly. No checks are performed on poly and
// the returned Plaintext will share its backing array of coefficient.
func NewPlaintextAtLevelFromPoly(level int, poly *ring.Poly) *Plaintext {
	pt := rlwe.NewPlaintextAtLevelFromPoly(level, poly)
	pt.Value.IsNTT = true
	return &Plaintext{Plaintext: pt, Scale: 1}
}

// ScalingFactor returns the scaling factor of the target plaintext.
func (p *Plaintext) ScalingFactor() uint64 {
	return p.Scale
}
