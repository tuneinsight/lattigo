package rlwe

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

// Ciphertext is a generic type for RLWE ciphertexts.
type Ciphertext struct {
	Element[ring.Poly]
}

// NewCiphertext returns a new [Ciphertext] with zero values and an associated
// MetaData set to the Parameters default value.
func NewCiphertext(params ParameterProvider, degree int, level ...int) (ct *Ciphertext) {
	op := *NewElement(params, degree, level...)
	return &Ciphertext{op}
}

// NewCiphertextAtLevelFromPoly constructs a new [Ciphertext] at a specific level
// where the message is set to the passed poly. No checks are performed on poly and
// the returned [Ciphertext] will share its backing array of coefficients.
// Returned [Ciphertext]'s MetaData is allocated but empty.
func NewCiphertextAtLevelFromPoly(level int, poly []ring.Poly) (*Ciphertext, error) {

	operand, err := NewElementAtLevelFromPoly(level, poly)

	if err != nil {
		return nil, fmt.Errorf("cannot NewCiphertextAtLevelFromPoly: %w", err)
	}

	operand.MetaData = &MetaData{}

	return &Ciphertext{*operand}, nil
}

// NewCiphertextRandom generates a new uniformly distributed [Ciphertext] of degree, level.
func NewCiphertextRandom(prng sampling.PRNG, params ParameterProvider, degree, level int) (ciphertext *Ciphertext) {
	ciphertext = NewCiphertext(params, degree, level)
	PopulateElementRandom(prng, params, ciphertext.El())
	return
}

// Plaintext casts the target ciphertext into a plaintext type.
// This method is allocation free.
func (ct Ciphertext) Plaintext() *Plaintext {
	return &Plaintext{
		Element: Element[ring.Poly]{
			Value:    ct.Element.Value[:1],
			MetaData: ct.MetaData,
		},
		Value: ct.Element.Value[0],
	}
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

// NewCiphertextFromUintPool returns a new [*Ciphertext], built from backing []uint64 arrays obtained from the pool stored in the ring passed with the params.
// After use, the [Ciphertext] should be recycled using the [RecycleCiphertextInUintPool] method.
func NewCiphertextFromUintPool(params ParameterProvider, degree int, levelQ int) *Ciphertext {
	p := params.GetRLWEParameters()

	ringQ := p.RingQ().AtLevel(levelQ)

	Value := make([]ring.Poly, degree+1)
	for i := range Value {
		Value[i] = *ringQ.NewPolyFromUintPool()
	}

	el := Element[ring.Poly]{
		Value: Value,
		MetaData: &MetaData{
			CiphertextMetaData: CiphertextMetaData{
				IsNTT: p.NTTFlag(),
			},
		},
	}
	return &Ciphertext{el}
}

// RecycleCiphertextInUintPool takes a reference to a [Ciphertext] and recycles its backing []uint64 arrays
// (i.e. they are returned to a pool). The input [Ciphertext] must not be used after calling this method.
func RecycleCiphertextInUintPool(params ParameterProvider, ct *Ciphertext) {
	ringQ := params.GetRLWEParameters().ringQ
	for i := range ct.Value {
		ringQ.RecyclePolyInUintPool(&ct.Value[i])
	}
	ct = nil
}
