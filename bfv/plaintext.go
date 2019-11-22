package bfv

import (
	"github.com/ldsec/lattigo/ring"
	"math/bits"
)

// Plaintext is a BigPoly of degree 0.
type Plaintext struct {
	*bfvElement
	value *ring.Poly
}

// NewPlaintext creates a new plaintext from the target context.
func (context *Context) NewPlaintext() *Plaintext {

	plaintext := &Plaintext{&bfvElement{}, nil}
	plaintext.bfvElement.value = []*ring.Poly{context.contextQ.NewPoly()}
	plaintext.value = plaintext.bfvElement.value[0]
	plaintext.isNTT = false

	return plaintext
}

// NewRandomPlaintextCoeffs generates a slice of random coefficient sampled within the plaintext modulus.
func (context *Context) NewRandomPlaintextCoeffs() (coeffs []uint64) {
	coeffs = make([]uint64, context.n)
	mask := uint64(1<<uint64(bits.Len64(context.t))) - 1
	for i := uint64(0); i < context.n; i++ {
		coeffs[i] = ring.RandUniform(context.t, mask)
	}
	return
}

// SetCoefficientsInt64 sets the coefficients of a plaintext with the provided slice of int64.
func (P *Plaintext) setCoefficientsInt64(context *Context, coeffs []int64) {
	for i, coeff := range coeffs {
		for j := range context.contextQ.Modulus {
			P.value.Coeffs[j][i] = uint64((coeff%int64(context.t))+int64(context.t)) % context.t
		}
	}
}

// SetCoefficientsInt64 sets the coefficients of a plaintext with the provided slice of uint64.
func (P *Plaintext) setCoefficientsUint64(context *Context, coeffs []uint64) {

	for i, coeff := range coeffs {
		for j := range context.contextQ.Modulus {
			P.value.Coeffs[j][i] = coeff % context.t
		}
	}
}

// GetCoefficients returns the coefficients of the plaintext in their CRT representation (double-slice).
func (P *Plaintext) GetCoefficients() [][]uint64 {
	return P.value.GetCoefficients()
}

// Lift scales the coefficient of the plaintext by Q/t (ciphertext modulus / plaintext modulus) and switches
// its modulus from t to Q.
func (P *Plaintext) Lift(context *encoderContext) {
	ringContext := context.contextQ
	for j := uint64(0); j < ringContext.N; j++ {
		for i := len(ringContext.Modulus) - 1; i >= 0; i-- {
			P.value.Coeffs[i][j] = ring.MRed(P.value.Coeffs[0][j], context.deltaMont[i], ringContext.Modulus[i], ringContext.GetMredParams()[i])
		}
	}
}

// NTTPlainModulus applies the NTT on a plaintext within the plaintext modulus.
func (P *Plaintext) NTTPlainModulus(context *Context) {

	if context.contextT.AllowsNTT() != true {
		panic("plaintext context doesn't allow a valid NTT")
	}

	context.contextT.NTT(P.value, P.value)
}

// InvNTTPlainModulus applies the NTT on a plaintext within the plaintext modulus.
func (P *Plaintext) InvNTTPlainModulus(context *encoderContext) {

	if context.contextT.AllowsNTT() != true {
		panic("plaintext context doesn't allow a valid InvNTT")
	}

	context.contextT.InvNTT(P.value, P.value)
}
