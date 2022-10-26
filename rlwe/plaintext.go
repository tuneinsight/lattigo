package rlwe

import (
	"github.com/tuneinsight/lattigo/v4/ring"
)

// Plaintext is a common base type for RLWE plaintexts.
type Plaintext struct {
	MetaData
	Value *ring.Poly
}

// NewPlaintext creates a new Plaintext at level `level` from the parameters.
func NewPlaintext(params Parameters, level int) (pt *Plaintext) {
	return &Plaintext{Value: ring.NewPoly(params.N(), level), MetaData: params.DefaultMetaData()}
}

// NewPlaintextAtLevelFromPoly constructs a new Plaintext at a specific level
// where the message is set to the passed poly. No checks are performed on poly and
// the returned Plaintext will share its backing array of coefficients.
// Returned plaintext's MetaData is empty.
func NewPlaintextAtLevelFromPoly(level int, poly *ring.Poly) (pt *Plaintext) {
	if len(poly.Coeffs) < level+1 {
		panic("cannot NewPlaintextAtLevelFromPoly: provided ring.Poly level is too small")
	}
	v0 := new(ring.Poly)
	v0.Coeffs = poly.Coeffs[:level+1]
	v0.Buff = poly.Buff[:poly.N()*(level+1)]
	return &Plaintext{Value: v0}
}

// Degree returns the degree of the target Plaintext.
func (pt *Plaintext) Degree() int {
	return 0
}

// Level returns the level of the target Plaintext.
func (pt *Plaintext) Level() int {
	return len(pt.Value.Coeffs) - 1
}

// GetScale gets the scale of the target Plaintext.
func (pt *Plaintext) GetScale() Scale {
	return pt.Scale
}

// SetScale sets the scale of the target Plaintext.
func (pt *Plaintext) SetScale(scale Scale) {
	pt.Scale = scale
}

// El returns the plaintext as a new `Element` for which the value points
// to the receiver `Value` field.
func (pt *Plaintext) El() *Ciphertext {
	return &Ciphertext{Value: []*ring.Poly{pt.Value}, MetaData: pt.MetaData}
}

// Copy copies the `other` plaintext value into the reciever plaintext.
func (pt *Plaintext) Copy(other *Plaintext) {
	if other != nil && other.Value != nil {
		pt.Value.Copy(other.Value)
		pt.MetaData = other.MetaData
	}
}
