package bfv

import (
	"errors"
	"github.com/lca1/lattigo/ring"
)

type Ciphertext BigPoly

// NewCiphertext creates a new ciphertext of degree n
func (bfvcontext *BfvContext) NewCiphertext(degree uint64) *Ciphertext {
	ciphertext := new(Ciphertext)
	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = bfvcontext.contextQ.NewPoly()
	}
	ciphertext.isNTT = false

	return ciphertext
}

// NewCiphertext creates a new ciphertext of degree n
func (bfvcontext *BfvContext) NewCiphertextBig(degree uint64) *Ciphertext {
	ciphertext := new(Ciphertext)
	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = bfvcontext.contextQP.NewPoly()
	}
	ciphertext.isNTT = false

	return ciphertext
}

func (bfvcontext *BfvContext) NewRandomCiphertext(degree uint64) *Ciphertext {
	ciphertext := new(Ciphertext)
	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = bfvcontext.contextQ.NewUniformPoly()
	}
	ciphertext.isNTT = false

	return ciphertext
}

func (ctx *Ciphertext) Value() []*ring.Poly {
	return ctx.value
}

func (ctx *Ciphertext) SetValue(value []*ring.Poly) {
	ctx.value = value
}

func (ctx *Ciphertext) Degree() uint64 {
	return uint64(len(ctx.value) - 1)
}

func (ctx *Ciphertext) Resize(bfvcontext *BfvContext, degree uint64) {
	if ctx.Degree() > degree {
		ctx.value = ctx.value[:degree]
	} else if ctx.Degree() < degree {
		for ctx.Degree() < degree {
			ctx.value = append(ctx.value, []*ring.Poly{bfvcontext.contextQ.NewPoly()}...)
		}
	}
}

func (ctx *Ciphertext) IsNTT() bool {
	return ctx.isNTT
}

func (ctx *Ciphertext) SetIsNTT(value bool) {
	ctx.isNTT = value
}

// Creates a new ciphertext of the same format and coefficients.
func (ctx *Ciphertext) CopyNew() BfvElement {

	ctxCopy := new(Ciphertext)

	ctxCopy.value = make([]*ring.Poly, ctx.Degree()+1)
	for i := range ctx.value {
		ctxCopy.value[i] = ctx.value[i].CopyNew()
	}
	ctxCopy.isNTT = ctx.isNTT

	return ctxCopy
}

// Copies the value of the ciphertext on a reciever ciphertext of the same format
func (ctx *Ciphertext) Copy(ctxCopy BfvElement) error {

	for i := range ctxCopy.Value() {
		ctxCopy.Value()[i].Copy(ctx.Value()[i])
	}
	ctxCopy.SetIsNTT(ctx.IsNTT())

	return nil
}

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
