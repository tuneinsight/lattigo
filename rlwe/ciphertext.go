package rlwe

import (
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

// Ciphertext is a generic type for RLWE ciphertexts.
type Ciphertext struct {
	OperandQ
}

// NewCiphertext returns a new Ciphertext with zero values and an associated
// MetaData set to the Parameters default value.
func NewCiphertext(params ParametersInterface, degree, level int) (ct *Ciphertext) {
	op := *NewOperandQ(params, degree, level)
	op.Scale = params.DefaultScale()
	op.LogSlots = params.MaxLogSlots()
	return &Ciphertext{op}
}

// NewCiphertextAtLevelFromPoly constructs a new Ciphertext at a specific level
// where the message is set to the passed poly. No checks are performed on poly and
// the returned Ciphertext will share its backing array of coefficients.
// Returned Ciphertext's MetaData is empty.
func NewCiphertextAtLevelFromPoly(level int, poly []*ring.Poly) *Ciphertext {
	return &Ciphertext{*NewOperandQAtLevelFromPoly(level, poly)}
}

// NewCiphertextRandom generates a new uniformly distributed Ciphertext of degree, level.
func NewCiphertextRandom(prng sampling.PRNG, params ParametersInterface, degree, level int) (ciphertext *Ciphertext) {
	ciphertext = NewCiphertext(params, degree, level)
	PopulateElementRandom(prng, params, ciphertext.El())
	return
}

// CopyNew creates a new element as a copy of the target element.
func (ct *Ciphertext) CopyNew() *Ciphertext {
	return &Ciphertext{OperandQ: *ct.OperandQ.CopyNew()}
}

// Copy copies the input element and its parameters on the target element.
func (ct *Ciphertext) Copy(ctxCopy *Ciphertext) {
	ct.OperandQ.Copy(&ctxCopy.OperandQ)
}

// Equal performs a deep equal.
func (ct *Ciphertext) Equal(other *Ciphertext) bool {
	return ct.OperandQ.Equal(&other.OperandQ)
}
