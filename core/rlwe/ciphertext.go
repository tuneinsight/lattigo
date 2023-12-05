package rlwe

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
)

// Ciphertext is a generic type for RLWE ciphertexts.
type Ciphertext struct {
	Element[ring.Poly]
}

// NewCiphertext returns a new Ciphertext with zero values and an associated
// MetaData set to the Parameters default value.
func NewCiphertext(params ParameterProvider, degree int, level ...int) (ct *Ciphertext) {
	op := *NewElement(params, degree, level...)
	return &Ciphertext{op}
}

// NewCiphertextAtLevelFromPoly constructs a new Ciphertext at a specific level
// where the message is set to the passed poly. No checks are performed on poly and
// the returned Ciphertext will share its backing array of coefficients.
// Returned Ciphertext's MetaData is allocated but empty	.
func NewCiphertextAtLevelFromPoly(level int, poly []ring.Poly) (*Ciphertext, error) {

	operand, err := NewElementAtLevelFromPoly(level, poly)

	if err != nil {
		return nil, fmt.Errorf("cannot NewCiphertextAtLevelFromPoly: %w", err)
	}

	operand.MetaData = &MetaData{}

	return &Ciphertext{*operand}, nil
}

// NewCiphertextRandom generates a new uniformly distributed Ciphertext of degree, level.
func NewCiphertextRandom(prng sampling.PRNG, params ParameterProvider, degree, level int) (ciphertext *Ciphertext) {
	ciphertext = NewCiphertext(params, degree, level)
	PopulateElementRandom(prng, params, ciphertext.El())
	return
}

// CopyNew creates a new element as a copy of the target element.
func (ct Ciphertext) CopyNew() *Ciphertext {
	return &Ciphertext{Element: *ct.Element.CopyNew()}
}

// Copy copies the input element and its parameters on the target element.
func (ct Ciphertext) Copy(ctxCopy *Ciphertext) {
	ct.Element.Copy(&ctxCopy.Element)
}

// Equal performs a deep equal.
func (ct Ciphertext) Equal(other *Ciphertext) bool {
	return ct.Element.Equal(&other.Element)
}
