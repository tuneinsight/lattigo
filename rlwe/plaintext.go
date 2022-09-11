package rlwe

import (
	"github.com/tuneinsight/lattigo/v3/ring"
)

// Plaintext is a common base type for RLWE plaintexts.
type Plaintext struct {
	Scale
	Value *ring.Poly
}

// NewPlaintext creates a new Plaintext at level `level` from the parameters.
func NewPlaintext(params Parameters, level int) *Plaintext {
	return &Plaintext{Value: ring.NewPoly(params.N(), level)}
}

// NewPlaintextAtLevelFromPoly constructs a new Plaintext at a specific level
// where the message is set to the passed poly. No checks are performed on poly and
// the returned Plaintext will share its backing array of coefficients.
func NewPlaintextAtLevelFromPoly(level int, poly *ring.Poly) *Plaintext {
	if len(poly.Coeffs) < level+1 {
		panic("cannot NewPlaintextAtLevelFromPoly: provided ring.Poly level is too small")
	}
	v0 := new(ring.Poly)
	v0.Coeffs = poly.Coeffs[:level+1]
	v0.Buff = poly.Buff[:poly.N()*(level+1)]
	return &Plaintext{Value: v0}
}

// Degree returns the degree of the target element.
func (pt *Plaintext) Degree() int {
	return 0
}

// Level returns the level of the target element.
func (pt *Plaintext) Level() int {
	return len(pt.Value.Coeffs) - 1
}

// El returns the plaintext as a new `Element` for which the value points
// to the receiver `Value` field.
func (pt *Plaintext) El() *Ciphertext {
	return &Ciphertext{Value: []*ring.Poly{pt.Value}, Scale: pt.Scale}
}

// Copy copies the `other` plaintext value into the reciever plaintext.
func (pt *Plaintext) Copy(other *Plaintext) {
	if other != nil && other.Value != nil {
		pt.Value.Copy(other.Value)

		if pt.Scale != nil {
			other.Scale = pt.Scale.CopyNew()
		}
	}
}
