package ckks

import (
	"github.com/lca1/lattigo/ring"
)

// Plaintext is BigPoly of degree 0.
type Plaintext BigPoly

// NewPlaintext creates a new plaintext of level level and scale scale.
func (ckkscontext *CkksContext) NewPlaintext(level uint64, scale uint64) *Plaintext {
	plaintext := new(Plaintext)
	plaintext.value = []*ring.Poly{ckkscontext.contextLevel[level].NewPoly()}
	plaintext.scale = scale
	plaintext.currentModulus = ring.Copy(ckkscontext.contextLevel[level].ModulusBigint)
	plaintext.isNTT = true
	return plaintext
}

// Value returns the value (polynomial) of the plaintext.
func (P *Plaintext) Value() []*ring.Poly {
	return P.value
}

// SetValue sets the value (polynomial) of the plaintext to the provided value.
func (P *Plaintext) SetValue(value []*ring.Poly) {
	P.value = value
}

// Resize does nothing on a plaintext since it is always of degree 0.
func (P *Plaintext) Resize(ckkscontext *CkksContext, degree uint64) {

}

// CurrentModulus returns the current modulus of the plaintext.
// This variable is only used during the decoding.
func (P *Plaintext) CurrentModulus() *ring.Int {
	return P.currentModulus
}

// SetCurrentModulus sets the current modulus to the provided values.
// This variable is only used during the decoding.
func (P *Plaintext) SetCurrentModulus(modulus *ring.Int) {
	P.currentModulus = ring.Copy(modulus)
}

// Degree returns the degree of the plaintext,
// this value should always be zero.
func (P *Plaintext) Degree() uint64 {
	return uint64(len(P.value) - 1)
}

// Level returns the current level of the plaintext.
func (P *Plaintext) Level() uint64 {
	return uint64(len(P.value[0].Coeffs) - 1)
}

// Scale returns the current scale of the plaintext (in log2).
func (P *Plaintext) Scale() uint64 {
	return P.scale
}

// SetScale sets the scale of the plaintext to the provided value (in log2).
func (P *Plaintext) SetScale(scale uint64) {
	P.scale = scale
}

// IsNTT returns true or false depending on if the plaintext is in the NTT domain or not.
func (P *Plaintext) IsNTT() bool {
	return P.isNTT
}

// SetIsNTT sets the isNTT value of the plaintext to the provided value.
func (P *Plaintext) SetIsNTT(isNTT bool) {
	P.isNTT = isNTT
}

// NTT applies the NTT transform to a plaintext and returns the result on the receiver element.
// Can only be used if the plaintext is not already in the NTT domain.
func (P *Plaintext) NTT(ckkscontext *CkksContext, ct0 CkksElement) {

	if P.isNTT != true {
		for i := range ct0.Value() {
			ckkscontext.contextLevel[P.Level()].NTT(P.value[i], ct0.Value()[i])
		}
		ct0.SetIsNTT(true)
	}
}

// InvNTT applies the inverse NTT transform to a plaintext and returns the result on the receiver element.
// Can only be used it the plaintext is in the NTT domain
func (P *Plaintext) InvNTT(ckkscontext *CkksContext, ct0 CkksElement) {

	if P.isNTT != false {
		for i := range ct0.Value() {
			ckkscontext.contextLevel[P.Level()].InvNTT(P.value[i], ct0.Value()[i])
		}
		ct0.SetIsNTT(false)
	}
}

// CopyNew creates a new plaintext with the same value and same parameters.
func (P *Plaintext) CopyNew() CkksElement {
	PCopy := new(Plaintext)
	PCopy.value = make([]*ring.Poly, 1)
	PCopy.value[0] = P.value[0].CopyNew()
	P.CopyParams(PCopy)
	return PCopy
}

// Copy copies the value and parameters of the target plaintext ot the receiver plaintext.
func (P *Plaintext) Copy(PCopy CkksElement) error {
	P.value[0].Copy(PCopy.Value()[0])
	P.CopyParams(PCopy)
	return nil
}

// CopyParams copies the parameters of the target plaintext to the receiver plaintext.
func (P *Plaintext) CopyParams(ckkselement CkksElement) {
	ckkselement.SetCurrentModulus(P.CurrentModulus())
	ckkselement.SetScale(P.Scale())
	ckkselement.SetIsNTT(P.IsNTT())
}
