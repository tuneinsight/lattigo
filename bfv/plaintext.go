package bfv

import (
	"errors"
	"github.com/lca1/lattigo/ring"
	"math/bits"
)

// Plaintext is a BigPoly of degree 0.
type Plaintext BigPoly

// NewPlaintext creates a new plaintext from the target bfvcontext.
func (bfvcontext *BfvContext) NewPlaintext() *Plaintext {

	plaintext := new(Plaintext)
	plaintext.value = []*ring.Poly{bfvcontext.contextQ.NewPoly()}
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
			P.value[0].Coeffs[j][i] = uint64((coeff%int64(bfvcontext.t))+int64(bfvcontext.t)) % bfvcontext.t
		}
	}
}

// SetCoefficientsInt64 sets the coefficients of a plaintext with the provided slice of uint64.
func (P *Plaintext) setCoefficientsUint64(bfvcontext *BfvContext, coeffs []uint64) {

	for i, coeff := range coeffs {
		for j := range bfvcontext.contextQ.Modulus {
			P.value[0].Coeffs[j][i] = coeff % bfvcontext.t
		}
	}
}

// GetCoefficients returns the coefficients of the plaintext in their CRT representation (double-slice).
func (P *Plaintext) GetCoefficients() [][]uint64 {
	return P.value[0].GetCoefficients()
}

// Value returns the coefficients of the plaintext a slice of polynomial (which is always of length 1).
func (P *Plaintext) Value() []*ring.Poly {
	return P.value
}

// SetValue sets the coefficents of the plaintext with the input slice of polynomials.
func (P *Plaintext) SetValue(value []*ring.Poly) {
	P.value = value
}

// IsNTT returns true if the plaintext is in the NTT domain, else false.
func (P *Plaintext) IsNTT() bool {
	return P.isNTT
}

// SetisNTT sets the plaintext isNTT flag to the input value.
func (P *Plaintext) SetIsNTT(value bool) {
	P.isNTT = value
}

// Degree returns the degree of the plaintext (which is always zero).
func (P *Plaintext) Degree() uint64 {
	return uint64(len(P.value) - 1)
}

// Lift scales the coefficient of the plaintext by Q/t (ciphertext modulus / plaintext modulus) and switches
// its modulus from t to Q.
func (P *Plaintext) Lift(bfvcontext *BfvContext) {
	context := bfvcontext.contextQ
	for j := uint64(0); j < bfvcontext.n; j++ {
		for i := len(context.Modulus) - 1; i >= 0; i-- {
			P.value[0].Coeffs[i][j] = ring.MRed(P.value[0].Coeffs[0][j], bfvcontext.DeltaMont[i], context.Modulus[i], context.GetMredParams()[i])
		}
	}
}

// Resize does nothing on a plaintext since its always of degree 0.
func (P *Plaintext) Resize(bfvcontext *BfvContext, degree uint64) {

}

// Add adds p0 to p1 within the plaintext modulus and returns the result on the target plaintext.
func (P *Plaintext) Add(bfvcontext *BfvContext, p0, p1 *Plaintext) {
	bfvcontext.contextT.Add(p0.value[0], p1.value[0], P.value[0])
}

// Sub subtracts p1 to p0 within the plaintext modulus and returns the result on the target plaintext.
func (P *Plaintext) Sub(bfvcontext *BfvContext, p0, p1 *Plaintext) {
	bfvcontext.contextT.Sub(p0.value[0], p1.value[0], P.value[0])
}

// Mul multiplies p0 by p1 within the plaintext modulus and returns the result on the target plaintext.
func (P *Plaintext) Mul(bfvcontext *BfvContext, p0, p1 *Plaintext) {

	// Checks if the plaintext contexts has been validated (allowing NTT)
	// Else performe the multiplication with a naive convolution

	if bfvcontext.contextT.IsValidated() {
		bfvcontext.contextT.MulPoly(p0.value[0], p1.value[0], P.value[0])
	} else {
		bfvcontext.contextT.MulPolyNaive(p0.value[0], p1.value[0], P.value[0])
	}
}

// NTT puts a lifted plaintext in the NTT domain (the NTT is applied within the ciphertet modulus), sets its isNTT flag to true, and returns the result on p. If the isNTT flag is true does nothing.
func (P *Plaintext) NTT(bfvcontext *BfvContext, p BfvElement) error {
	if P.Degree() != p.Degree() {
		return errors.New("error : receiver element invalide degree (does not match)")
	}
	if P.IsNTT() != true {
		for i := range P.value {
			bfvcontext.contextQ.NTT(P.Value()[i], p.Value()[i])
		}
		p.SetIsNTT(true)
	}
	return nil
}

// InvNTT puts a lifted plaintext outside of the NTT domain (the InvNTT is applied within the ciphertext modulus), and sets its isNTT flag to false, and returns the result on p. If the isNTT flag is flase, does nothing.
func (P *Plaintext) InvNTT(bfvcontext *BfvContext, p BfvElement) error {
	if P.Degree() != p.Degree() {
		return errors.New("error : receiver element invalide degree (does not match)")
	}
	if P.IsNTT() != false {
		for i := range P.value {
			bfvcontext.contextQ.InvNTT(P.Value()[i], p.Value()[i])
		}
		p.SetIsNTT(false)
	}
	return nil
}

// EMBInv applies the InvNTT on a plaintext within the plaintext modulus.
func (P *Plaintext) EMBInv(bfvcontext *BfvContext) error {

	if bfvcontext.contextT.IsValidated() != true {
		return errors.New("plaintext context doesn't allow a valid NTT")
	}

	bfvcontext.contextT.NTT(P.value[0], P.value[0])

	return nil
}

// EMB applies the NTT on a plaintext within the plaintext modulus.
func (P *Plaintext) EMB(bfvcontext *BfvContext) error {
	if bfvcontext.contextT.IsValidated() != true {
		return errors.New("plaintext context doesn't allow a valid InvNTT")
	}

	bfvcontext.contextT.InvNTT(P.value[0], P.value[0])

	return nil
}

// CopyNew creates a new element which is a copy of the target plaintext.
func (P *Plaintext) CopyNew() BfvElement {

	PCopy := new(Plaintext)

	PCopy.value = make([]*ring.Poly, P.Degree()+1)
	for i := range P.value {
		PCopy.value[i] = P.value[i].CopyNew()
	}

	PCopy.isNTT = P.isNTT

	return PCopy
}

// Copy copies the target plaintext on the input element.
func (P *Plaintext) Copy(PCopy BfvElement) error {

	for i := range P.value {
		PCopy.Value()[i].Copy(P.Value()[i])
	}
	P.SetIsNTT(P.IsNTT())

	return nil
}
