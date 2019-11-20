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

// NewPlaintext creates a new plaintext from the target bfvcontext.
func (bfvcontext *BfvContext) NewPlaintext() *Plaintext {

	plaintext := &Plaintext{&bfvElement{}, nil}
	plaintext.bfvElement.value = []*ring.Poly{bfvcontext.contextQ.NewPoly()}
	plaintext.value = plaintext.bfvElement.value[0]
	plaintext.isNTT = false

	return plaintext
}

// NewRandomPlaintextCoeffs generates a slice of random coefficient sampled within the plaintext modulus.
func (bfvcontext *BfvContext) NewRandomPlaintextCoeffs() (coeffs []uint64) {
	coeffs = make([]uint64, bfvcontext.n)
	mask := uint64(1<<uint64(bits.Len64(bfvcontext.t))) - 1
	for i := uint64(0); i < bfvcontext.n; i++ {
		coeffs[i] = ring.RandUniform(bfvcontext.t, mask)
	}
	return
}

// SetCoefficientsInt64 sets the coefficients of a plaintext with the provided slice of int64.
func (P *Plaintext) setCoefficientsInt64(bfvcontext *BfvContext, coeffs []int64) {
	for i, coeff := range coeffs {
		for j := range bfvcontext.contextQ.Modulus {
			P.value.Coeffs[j][i] = uint64((coeff%int64(bfvcontext.t))+int64(bfvcontext.t)) % bfvcontext.t
		}
	}
}

// SetCoefficientsInt64 sets the coefficients of a plaintext with the provided slice of uint64.
func (P *Plaintext) setCoefficientsUint64(bfvcontext *BfvContext, coeffs []uint64) {

	for i, coeff := range coeffs {
		for j := range bfvcontext.contextQ.Modulus {
			P.value.Coeffs[j][i] = coeff % bfvcontext.t
		}
	}
}

// GetCoefficients returns the coefficients of the plaintext in their CRT representation (double-slice).
func (P *Plaintext) GetCoefficients() [][]uint64 {
	return P.value.GetCoefficients()
}

// Lift scales the coefficient of the plaintext by Q/t (ciphertext modulus / plaintext modulus) and switches
// its modulus from t to Q.
func (P *Plaintext) Lift(bfvcontext *BfvContext) {
	context := bfvcontext.contextQ
	for j := uint64(0); j < bfvcontext.n; j++ {
		for i := len(context.Modulus) - 1; i >= 0; i-- {
			P.value.Coeffs[i][j] = ring.MRed(P.value.Coeffs[0][j], bfvcontext.deltaMont[i], context.Modulus[i], context.GetMredParams()[i])
		}
	}
}

// NTTPlainModulus applies the NTT on a plaintext within the plaintext modulus.
func (P *Plaintext) NTTPlainModulus(bfvcontext *BfvContext) {

	if bfvcontext.contextT.AllowsNTT() != true {
		panic("plaintext context doesn't allow a valid NTT")
	}

	bfvcontext.contextT.NTT(P.value, P.value)
}

// InvNTTPlainModulus applies the NTT on a plaintext within the plaintext modulus.
func (P *Plaintext) InvNTTPlainModulus(bfvcontext *BfvContext) {

	if bfvcontext.contextT.AllowsNTT() != true {
		panic("plaintext context doesn't allow a valid InvNTT")
	}

	bfvcontext.contextT.InvNTT(P.value, P.value)
}
