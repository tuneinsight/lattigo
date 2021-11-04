package ckks

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// Ciphertext is *ring.Poly array representing a polynomial of degree > 0 with coefficients in R_Q.
type Ciphertext struct {
	*rlwe.Ciphertext
	Scale float64
}

// NewCiphertext creates a new Ciphertext parameterized by degree, level and scale.
func NewCiphertext(params Parameters, degree, level int, scale float64) (ciphertext *Ciphertext) {
	ciphertext = &Ciphertext{Ciphertext: rlwe.NewCiphertext(params.Parameters, degree, level)}
	for _, pol := range ciphertext.Value {
		pol.IsNTT = true
	}
	ciphertext.Scale = scale
	return ciphertext
}

// NewCiphertextRandom generates a new uniformly distributed Ciphertext of degree, level and scale.
func NewCiphertextRandom(prng utils.PRNG, params Parameters, degree, level int, scale float64) (ciphertext *Ciphertext) {
	return &Ciphertext{rlwe.NewCiphertextRandom(prng, params.Parameters, degree, level), scale}
}

// NewCiphertextAtLevelFromPoly construct a new Ciphetext at a specific level
// where the message is set to the passed poly. No checks are performed on poly and
// the returned Ciphertext will share its backing array of coefficient.
func NewCiphertextAtLevelFromPoly(level int, poly [2]*ring.Poly) *Ciphertext {
	v0, v1 := new(ring.Poly), new(ring.Poly)
	v0.IsNTT, v1.IsNTT = true, true
	v0.Coeffs, v1.Coeffs = poly[0].Coeffs[:level+1], poly[1].Coeffs[:level+1]
	return &Ciphertext{Ciphertext: &rlwe.Ciphertext{Value: []*ring.Poly{v0, v1}}, Scale: 0}
}

// ScalingFactor returns the scaling factor of the ciphertext
func (ct *Ciphertext) ScalingFactor() float64 {
	return ct.Scale
}

// SetScalingFactor sets the scaling factor of the ciphertext
func (ct *Ciphertext) SetScalingFactor(scale float64) {
	ct.Scale = scale
}

// Copy copies the given ciphertext ctp into the receiver ciphertext.
func (ct *Ciphertext) Copy(ctp *Ciphertext) {
	ct.Ciphertext.Copy(ctp.Ciphertext)
	ct.Scale = ctp.Scale
}

// CopyNew makes a deep copy of the receiver ciphertext and returns it.
func (ct *Ciphertext) CopyNew() (ctc *Ciphertext) {
	ctc = &Ciphertext{Ciphertext: ct.Ciphertext.CopyNew(), Scale: ct.Scale}
	return
}
