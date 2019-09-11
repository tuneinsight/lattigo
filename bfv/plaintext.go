package bfv

import (
	"errors"
	"github.com/lca1/lattigo/ring"
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

// Resize does nothing on a plaintext since its always of degree 0.
func (P *Plaintext) Resize(bfvcontext *BfvContext, degree uint64) {

}

// Add adds p0 to p1 within the plaintext modulus and returns the result on the target plaintext.
func (P *Plaintext) Add(bfvcontext *BfvContext, p0, p1 *Plaintext) {
	bfvcontext.contextT.Add(p0.value, p1.value, P.value)
}

// Sub subtracts p1 to p0 within the plaintext modulus and returns the result on the target plaintext.
func (P *Plaintext) Sub(bfvcontext *BfvContext, p0, p1 *Plaintext) {
	bfvcontext.contextT.Sub(p0.value, p1.value, P.value)
}

// Mul multiplies p0 by p1 within the plaintext modulus and returns the result on the target plaintext.
func (P *Plaintext) Mul(bfvcontext *BfvContext, p0, p1 *Plaintext) {

	// Checks if the plaintext contexts has been validated (allowing NTT)
	// Else performe the multiplication with a naive convolution

	if bfvcontext.contextT.IsValidated() {
		bfvcontext.contextT.MulPoly(p0.value, p1.value, P.value)
	} else {
		bfvcontext.contextT.MulPolyNaive(p0.value, p1.value, P.value)
	}
}

// EMBInv applies the InvNTT on a plaintext within the plaintext modulus.
func (P *Plaintext) EMBInv(bfvcontext *BfvContext) error {

	if bfvcontext.contextT.IsValidated() != true {
		return errors.New("plaintext context doesn't allow a valid NTT")
	}

	bfvcontext.contextT.NTT(P.value, P.value)

	return nil
}

// EMB applies the NTT on a plaintext within the plaintext modulus.
func (P *Plaintext) EMB(bfvcontext *BfvContext) error {
	if bfvcontext.contextT.IsValidated() != true {
		return errors.New("plaintext context doesn't allow a valid InvNTT")
	}

	bfvcontext.contextT.InvNTT(P.value, P.value)

	return nil
}
