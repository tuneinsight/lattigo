package bfv

import (
	"errors"
	"github.com/ldsec/lattigo/ring"
)

// The plaintext is a ring of N coefficients with two contexts.
// The first context is defined by the BFV parameters. The second
// context defines a NTT around its modulus it it permits it.

type Plaintext BigPoly

// NewCiphertext creates a new ciphertext
func (bfvcontext *BfvContext) NewPlaintext() *Plaintext {

	plaintext := new(Plaintext)

	plaintext.bfvcontext = bfvcontext
	plaintext.value = []*ring.Poly{bfvcontext.contextQ.NewPoly()}
	plaintext.isNTT = false

	return plaintext
}

func (bfvcontext *BfvContext) NewRandomPlaintextCoeffs() (coeffs []uint64) {
	coeffs = make([]uint64, bfvcontext.n)
	for i := uint64(0); i < bfvcontext.n; i++ {
		coeffs[i] = ring.RandUniform(bfvcontext.t)
	}
	return
}

func (P *Plaintext) SetCoefficientsInt64(coeffs []int64) {
	for i, coeff := range coeffs {
		for j := range P.bfvcontext.contextQ.Modulus {
			P.value[0].Coeffs[j][i] = uint64((coeff%int64(P.bfvcontext.t))+int64(P.bfvcontext.t)) % P.bfvcontext.t
		}
	}
}

func (P *Plaintext) SetCoefficientsUint64(coeffs []uint64) {

	for i, coeff := range coeffs {
		for j := range P.bfvcontext.contextQ.Modulus {
			P.value[0].Coeffs[j][i] = coeff % P.bfvcontext.t
		}
	}
}

func (P *Plaintext) GetCoefficients() [][]uint64 {
	return P.value[0].GetCoefficients()
}

func (P *Plaintext) BfvContext() *BfvContext {
	return P.bfvcontext
}

func (P *Plaintext) Value() []*ring.Poly {
	return P.value
}

func (P *Plaintext) SetValue(value []*ring.Poly) {
	P.value = value
}

func (P *Plaintext) IsNTT() bool {
	return P.isNTT
}

func (P *Plaintext) SetIsNTT(value bool) {
	P.isNTT = value
}

func (P *Plaintext) Degree() uint64 {
	return uint64(len(P.value) - 1)
}

func (P *Plaintext) Lift() {
	context := P.bfvcontext.contextQ
	for j := uint64(0); j < P.bfvcontext.n; j++ {
		for i := len(context.Modulus) - 1; i >= 0; i-- {
			P.value[0].Coeffs[i][j] = ring.MRed(P.value[0].Coeffs[0][j], P.bfvcontext.DeltaMont[i], context.Modulus[i], context.GetMredParams()[i])
		}
	}
}

func (P *Plaintext) Resize(degree uint64) {
	if P.Degree() > degree {
		P.value = P.value[:degree]
	} else if P.Degree() < degree {
		for P.Degree() < degree {
			P.value = append(P.value, []*ring.Poly{P.bfvcontext.contextQ.NewPoly()}...)
		}
	}
}

func (P *Plaintext) Add(p0, p1 *Plaintext) {
	P.bfvcontext.contextT.Add(p0.value[0], p1.value[0], P.value[0])
}

func (P *Plaintext) Sub(p0, p1 *Plaintext) {
	P.bfvcontext.contextT.Sub(p0.value[0], p1.value[0], P.value[0])
}

func (P *Plaintext) Mul(p0, p1 *Plaintext) {

	// Checks if the plaintext contexts has been validated (allowing NTT)
	// Else performe the multiplication with a naive convolution

	if P.bfvcontext.contextT.IsValidated() {
		P.bfvcontext.contextT.MulPoly(p0.value[0], p1.value[0], P.value[0])
	} else {
		P.bfvcontext.contextT.MulPolyNaive(p0.value[0], p1.value[0], P.value[0])
	}
}

func (P *Plaintext) NTT(p BfvElement) error {
	if P.Degree() != p.Degree() {
		return errors.New("error : receiver element invalide degree (does not match)")
	}
	if P.IsNTT() != true {
		for i := range P.value {
			P.bfvcontext.contextQ.NTT(P.Value()[i], p.Value()[i])
		}
		p.SetIsNTT(true)
	}
	return nil
}

func (P *Plaintext) InvNTT(p BfvElement) error {
	if P.Degree() != p.Degree() {
		return errors.New("error : receiver element invalide degree (does not match)")
	}
	if P.IsNTT() != false {
		for i := range P.value {
			P.bfvcontext.contextQ.InvNTT(P.Value()[i], p.Value()[i])
		}
		p.SetIsNTT(false)
	}
	return nil
}

func (P *Plaintext) EMBInv() error {

	if P.bfvcontext.contextT.IsValidated() != true {
		return errors.New("plaintext context doesn't allow a valid NTT")
	}

	P.bfvcontext.contextT.NTT(P.value[0], P.value[0])

	return nil
}

func (P *Plaintext) EMB() error {
	if P.bfvcontext.contextT.IsValidated() != true {
		return errors.New("plaintext context doesn't allow a valid InvNTT")
	}

	P.bfvcontext.contextT.InvNTT(P.value[0], P.value[0])

	return nil
}

// Creates a new ciphertext of the same format and coefficients.
func (P *Plaintext) CopyNew() BfvElement {

	PCopy := new(Plaintext)

	PCopy.value = make([]*ring.Poly, P.Degree()+1)
	for i := range P.value {
		PCopy.value[i] = P.value[i].CopyNew()
	}
	PCopy.bfvcontext = P.bfvcontext
	PCopy.isNTT = P.isNTT

	return PCopy
}

// Copies the value of the ciphertext on a reciever ciphertext of the same format
func (P *Plaintext) Copy(PCopy BfvElement) error {

	if !checkContext([]BfvElement{P, PCopy}) {
		return errors.New("input ciphertext are not using the same bfvcontext")
	}

	for i := range P.value {
		PCopy.Value()[i].Copy(P.Value()[i])
	}
	P.SetIsNTT(P.IsNTT())

	return nil
}
