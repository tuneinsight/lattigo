package bgv

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
)

// Plaintext is is a Element with only one Poly.
type Plaintext struct {
	*rlwe.Plaintext
}

// NewPlaintext creates a new Plaintext of level level and scale scale.
func NewPlaintext(params Parameters, level int, scale rlwe.Scale) *Plaintext {
	pt := &Plaintext{Plaintext: rlwe.NewPlaintext(params.Parameters, level)}
	pt.Plaintext.Scale = CheckScaleType(scale).CopyNew()
	pt.Value.IsNTT = true
	return pt
}

// NewPlaintextAtLevelFromPoly construct a new Plaintext at a specific level
// where the message is set to the passed poly. No checks are performed on poly and
// the returned Plaintext will share its backing array of coefficient.
func NewPlaintextAtLevelFromPoly(level int, poly *ring.Poly) *Plaintext {
	pt := rlwe.NewPlaintextAtLevelFromPoly(level, poly)
	pt.Value.IsNTT = true
	pt.Scale = NewScale(0)
	return &Plaintext{Plaintext: pt}
}

func (pt *Plaintext) Scale() rlwe.Scale {
	return pt.Plaintext.Scale
}
