package bfv

import (
	"errors"
	"github.com/lca1/lattigo/ring"
)

// Ciphertext is a BigPoly of degree > 0
type Ciphertext BigPoly

// NewCiphertext creates a new empty ciphertext of degree degree.
func (bfvcontext *BfvContext) NewCiphertext(degree uint64) *Ciphertext {
	ciphertext := new(Ciphertext)
	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = bfvcontext.contextQ.NewPoly()
	}
	ciphertext.isNTT = false

	return ciphertext
}

// NewCiphertext creates a new empty ciphertext of degree degree in the extended ciphertext context (Q + P).
func (bfvcontext *BfvContext) NewCiphertextBig(degree uint64) *Ciphertext {
	ciphertext := new(Ciphertext)
	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = bfvcontext.contextQP.NewPoly()
	}
	ciphertext.isNTT = false

	return ciphertext
}

// NewRandomCiphertext creates a new ciphertext with uniform coefficients.
func (bfvcontext *BfvContext) NewRandomCiphertext(degree uint64) *Ciphertext {
	ciphertext := new(Ciphertext)
	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = bfvcontext.contextQ.NewUniformPoly()
	}
	ciphertext.isNTT = false

	return ciphertext
}

// Value returns the value of the target ciphertext (as a slice of polynomials in CRT form).
func (ctx *Ciphertext) Value() []*ring.Poly {
	return ctx.value
}

// SetValue assigns the input slice of polynomial to the target ciphertext value.
func (ctx *Ciphertext) SetValue(value []*ring.Poly) {
	ctx.value = value
}

// Degree returns the target ciphertext degree.
func (ctx *Ciphertext) Degree() uint64 {
	return uint64(len(ctx.value) - 1)
}

// Resize resizes the target ciphertext degree to the degree given as input. If the input degree is bigger then
// it will append new empty polynomials, if the degree is smaller, it will delete polynomials until the degree matches
// the input degree.
func (ctx *Ciphertext) Resize(bfvcontext *BfvContext, degree uint64) {
	if ctx.Degree() > degree {
		ctx.value = ctx.value[:degree]
	} else if ctx.Degree() < degree {
		for ctx.Degree() < degree {
			ctx.value = append(ctx.value, []*ring.Poly{bfvcontext.contextQ.NewPoly()}...)
		}
	}
}

// IsNTT returns true if the target ciphertext is in the NTT domain, else false.
func (ctx *Ciphertext) IsNTT() bool {
	return ctx.isNTT
}

// SetIsNTT assigns the input bolean value to the isNTT flag of the target ciphertext.
func (ctx *Ciphertext) SetIsNTT(value bool) {
	ctx.isNTT = value
}

// CopyNew creates a new ciphertext which is a copy of the target ciphertext. Returns the value as
// a BfvElement.
func (ctx *Ciphertext) CopyNew() BfvElement {

	ctxCopy := new(Ciphertext)

	ctxCopy.value = make([]*ring.Poly, ctx.Degree()+1)
	for i := range ctx.value {
		ctxCopy.value[i] = ctx.value[i].CopyNew()
	}
	ctxCopy.isNTT = ctx.isNTT

	return ctxCopy
}

// Copy copies the value of the target ciphertext on the reciever ciphertext.
func (ctx *Ciphertext) Copy(ctxCopy BfvElement) error {

	for i := range ctxCopy.Value() {
		ctxCopy.Value()[i].Copy(ctx.Value()[i])
	}
	ctxCopy.SetIsNTT(ctx.IsNTT())

	return nil
}

// NTT puts the target ciphertext in the NTT domain and sets its isNTT flag to true. If it is already in the NTT domain, does nothing.
func (ctx *Ciphertext) NTT(bfvcontext *BfvContext, c BfvElement) error {
	if ctx.Degree() != c.Degree() {
		return errors.New("error : receiver element invalide degree (does not match)")
	}
	if ctx.IsNTT() != true {
		for i := range ctx.value {
			bfvcontext.contextQ.NTT(ctx.Value()[i], c.Value()[i])
		}
		c.SetIsNTT(true)
	}
	return nil
}

// InvNTT puts the target ciphertext outside of the NTT domain, and sets its isNTT flag to false. If it is not in the NTT domain, does nothing.
func (ctx *Ciphertext) InvNTT(bfvcontext *BfvContext, c BfvElement) error {
	if ctx.Degree() != c.Degree() {
		return errors.New("error : receiver element invalide degree (does not match)")
	}
	if ctx.IsNTT() != false {
		for i := range ctx.value {
			bfvcontext.contextQ.InvNTT(ctx.Value()[i], c.Value()[i])
		}
		c.SetIsNTT(false)
	}
	return nil
}
