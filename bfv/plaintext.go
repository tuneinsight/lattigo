package bfv

import (
	"github.com/ldsec/lattigo/ring"
)

// Plaintext is a BigPoly of degree 0.
type Plaintext struct {
	*bfvElement
	value *ring.Poly
}

// NewPlaintext creates a new plaintext from the target context.
func NewPlaintextFromParams(params *Parameters) *Plaintext {
	plaintext := &Plaintext{&bfvElement{}, nil}

	plaintext.bfvElement.value = []*ring.Poly{ring.NewPoly(params.N, uint64(len(params.Q1)))}
	plaintext.value = plaintext.bfvElement.value[0]
	plaintext.isNTT = false

	return plaintext
}

// Lift scales the coefficient of the plaintext by Q/t (ciphertext modulus / plaintext modulus) and switches
// its modulus from t to Q.
func (P *Plaintext) Lift(bfvContext *Context) {
	ringContext := bfvContext.contextQ1
	for j := uint64(0); j < ringContext.N; j++ {
		for i := len(ringContext.Modulus) - 1; i >= 0; i-- {
			P.value.Coeffs[i][j] = ring.MRed(P.value.Coeffs[0][j], bfvContext.deltaMont[i], ringContext.Modulus[i], ringContext.GetMredParams()[i])
		}
	}
}

// NTTPlainModulus applies the NTT on a plaintext within the plaintext modulus.
func (P *Plaintext) NTTPlainModulus(context *ring.Context) {

	if context.AllowsNTT() != true {
		panic("plaintext context doesn't allow a valid NTT")
	}

	context.NTT(P.value, P.value)
}

// InvNTTPlainModulus applies the NTT on a plaintext within the plaintext modulus.
func (P *Plaintext) InvNTTPlainModulus(context *ring.Context) {

	if context.AllowsNTT() != true {
		panic("plaintext context doesn't allow a valid InvNTT")
	}

	context.InvNTT(P.value, P.value)
}
