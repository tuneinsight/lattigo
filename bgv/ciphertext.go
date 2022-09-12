package bgv

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// Ciphertext is *ring.Poly array representing a polynomial of degree > 0 with coefficients in R_Q.
type Ciphertext struct {
	*rlwe.Ciphertext
}

// NewCiphertext creates a new Ciphertext parameterized by degree, level and scale.
func NewCiphertext(params Parameters, degree, level int) (ciphertext *Ciphertext) {
	ciphertext = &Ciphertext{Ciphertext: rlwe.NewCiphertextNTT(params.Parameters, degree, level)}
	ciphertext.Ciphertext.Scale = NewScale(params, 1)
	return
}

// NewCiphertextRandom generates a new uniformly distributed Ciphertext of degree, level and scale.
func NewCiphertextRandom(prng utils.PRNG, params Parameters, degree, level int) (ciphertext *Ciphertext) {
	ciphertext = &Ciphertext{rlwe.NewCiphertextRandom(prng, params.Parameters, degree, level)}
	ciphertext.Ciphertext.Scale = NewScale(params, 1)
	for i := range ciphertext.Value {
		ciphertext.Value[i].IsNTT = true
	}
	return
}

func (ct *Ciphertext) Scale() rlwe.Scale {
	return ct.Ciphertext.Scale
}

// CopyNew creates a fresh copy of the target ciphertext.
func (ct *Ciphertext) CopyNew() *Ciphertext {
	return &Ciphertext{ct.Ciphertext.CopyNew()}
}

// Copy copies the input ciphertexton the target ciphertext.
func (ct *Ciphertext) Copy(input *Ciphertext) {
	ct.Ciphertext.Copy(input.Ciphertext)
}

func (ct *Ciphertext) UnmarshalBinary(data []byte) (err error) {
	ct.Ciphertext = &rlwe.Ciphertext{Scale: &Scale{}}
	return ct.Ciphertext.UnmarshalBinary(data)
}

// NewCiphertextAtLevelFromPoly construct a new Ciphetext at a specific level
// where the message is set to the passed poly. No checks are performed on poly and
// the returned Ciphertext will share its backing array of coefficient.
func NewCiphertextAtLevelFromPoly(level int, poly [2]*ring.Poly) *Ciphertext {
	ct := rlwe.NewCiphertextAtLevelFromPoly(level, poly)
	ct.Value[0].IsNTT, ct.Value[1].IsNTT = true, true
	return &Ciphertext{Ciphertext: ct}
}
