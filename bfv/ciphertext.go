package bfv

import (
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// Ciphertext is a *ring.Poly array representing a polynomial of degree > 0 with coefficients in R_Q.
type Ciphertext struct {
	*rlwe.Ciphertext
}

// NewCiphertext creates a new ciphertext parameterized by degree and at the max level.
func NewCiphertext(params Parameters, degree int) (ciphertext *Ciphertext) {
	return &Ciphertext{rlwe.NewCiphertext(params.Parameters, degree, params.MaxLevel())}
}

// NewCiphertextLvl creates a new ciphertext parameterized by degree and level.
func NewCiphertextLvl(params Parameters, degree, level int) (ciphertext *Ciphertext) {
	return &Ciphertext{rlwe.NewCiphertext(params.Parameters, degree, level)}
}

// NewCiphertextRandom generates a new uniformly distributed ciphertext of given degree at maximum level.
func NewCiphertextRandom(prng utils.PRNG, params Parameters, degree int) (ciphertext *Ciphertext) {
	return &Ciphertext{rlwe.NewCiphertextRandom(prng, params.Parameters, degree, params.MaxLevel())}
}

// NewCiphertextRandomLvl generates a new uniformly distributed ciphertext of given degree and level.
func NewCiphertextRandomLvl(prng utils.PRNG, params Parameters, degree, level int) (ciphertext *Ciphertext) {
	return &Ciphertext{rlwe.NewCiphertextRandom(prng, params.Parameters, degree, level)}
}

// CopyNew creates a deep copy of the receiver ciphertext and returns it.
func (ct *Ciphertext) CopyNew() *Ciphertext {
	return &Ciphertext{ct.Ciphertext.CopyNew()}
}

// MarshalBinary encodes a Ciphertext in a byte slice.
func (ct *Ciphertext) MarshalBinary() (data []byte, err error) {
	return ct.Ciphertext.MarshalBinary()
}

// UnmarshalBinary decodes a previously marshaled Ciphertext in the target Ciphertext.
func (ct *Ciphertext) UnmarshalBinary(data []byte) (err error) {
	ct.Ciphertext = new(rlwe.Ciphertext)
	return ct.Ciphertext.UnmarshalBinary(data)
}

// GetDataLen returns the length in bytes of the target Ciphertext.
func (ct *Ciphertext) GetDataLen(WithMetaData bool) (dataLen int) {
	return ct.Ciphertext.GetDataLen(WithMetaData)
}
