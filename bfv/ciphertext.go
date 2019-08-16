package bfv

import (
	"errors"
	"github.com/lca1/lattigo-private/ring"
)

type Ciphertext BigPoly

// NewCiphertext creates a new ciphertext of degree n
func (bfvcontext *BfvContext) NewCiphertext(degree uint64) *Ciphertext {
	ciphertext := new(Ciphertext)
	ciphertext.bfvcontext = bfvcontext
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
	ciphertext.bfvcontext = bfvcontext

	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = bfvcontext.contextQP.NewPoly()
	}
	ciphertext.isNTT = false

	return ciphertext
}

func (bfvcontext *BfvContext) NewRandomCiphertext(degree uint64) *Ciphertext {
	ciphertext := new(Ciphertext)
	ciphertext.bfvcontext = bfvcontext

	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = bfvcontext.contextQ.NewUniformPoly()
	}
	ciphertext.isNTT = false

	return ciphertext
}

func (ctx *Ciphertext) BfvContext() *BfvContext {
	return ctx.bfvcontext
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

func (ctx *Ciphertext) Resize(degree uint64) {
	if ctx.Degree() > degree {
		ctx.value = ctx.value[:degree]
	} else if ctx.Degree() < degree {
		for ctx.Degree() < degree {
			ctx.value = append(ctx.value, []*ring.Poly{ctx.bfvcontext.contextQ.NewPoly()}...)
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
	ctxCopy.bfvcontext = ctx.bfvcontext
	ctxCopy.isNTT = ctx.isNTT

	return ctxCopy
}

// Copies the value of the ciphertext on a reciever ciphertext of the same format
func (ctx *Ciphertext) Copy(ctxCopy BfvElement) error {

	if !checkContext([]BfvElement{ctx, ctxCopy}) {
		return errors.New("input ciphertext are not using the same bfvcontext")
	}

	for i := range ctxCopy.Value() {
		ctxCopy.Value()[i].Copy(ctx.Value()[i])
	}
	ctxCopy.SetIsNTT(ctx.IsNTT())

	return nil
}

func (ctx *Ciphertext) NTT(c BfvElement) error {
	if ctx.Degree() != c.Degree() {
		return errors.New("error : receiver element invalide degree (does not match)")
	}
	if ctx.IsNTT() != true {
		for i := range ctx.value {
			ctx.bfvcontext.contextQ.NTT(ctx.Value()[i], c.Value()[i])
		}
		c.SetIsNTT(true)
	}
	return nil
}

func (ctx *Ciphertext) InvNTT(c BfvElement) error {
	if ctx.Degree() != c.Degree() {
		return errors.New("error : receiver element invalide degree (does not match)")
	}
	if ctx.IsNTT() != false {
		for i := range ctx.value {
			ctx.bfvcontext.contextQ.InvNTT(ctx.Value()[i], c.Value()[i])
		}
		c.SetIsNTT(false)
	}
	return nil
}
