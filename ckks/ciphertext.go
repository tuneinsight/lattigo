package ckks

import (
	"github.com/lca1/lattigo/ring"
)

type Ciphertext BigPoly

// NewCiphertext creates a new ciphertext parametised by degree, level and scale.
func (ckkscontext *CkksContext) NewCiphertext(degree uint64, level uint64, scale uint64) *Ciphertext {
	ciphertext := new(Ciphertext)
	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = ckkscontext.contextLevel[level].NewPoly()
	}
	ciphertext.scale = scale
	ciphertext.currentModulus = ring.Copy(ckkscontext.contextLevel[level].ModulusBigint)
	ciphertext.isNTT = true

	return ciphertext
}

// Value returns the ciphertext polynomial.
func (ctx *Ciphertext) Value() []*ring.Poly {
	return ctx.value
}

// SetValue assigns a polynomial to a ciphertext value.
func (ctx *Ciphertext) SetValue(value []*ring.Poly) {
	ctx.value = value
}

// Resize resize the degree of a ciphertext, by allocating or removing polynomial.
// To be used to format a receiver of the wrong degree to a receiver of the correct degree.
func (ctx *Ciphertext) Resize(ckkscontext *CkksContext, degree uint64) {
	if ctx.Degree() > degree {
		ctx.value = ctx.value[:degree+1]
	} else if ctx.Degree() < degree {
		for ctx.Degree() < degree {
			ctx.value = append(ctx.value, []*ring.Poly{ckkscontext.contextLevel[ctx.Level()].NewPoly()}...)
		}
	}
}

// CurrentModulus returns the current modulus of the ciphertext.
// This value is only used at the decryption process.
func (ctx *Ciphertext) CurrentModulus() *ring.Int {
	return ctx.currentModulus
}

// SetCurrentModulus assigns a new modulus to the ciphertext. This
// value is only used at the decryption process.
func (ctx *Ciphertext) SetCurrentModulus(modulus *ring.Int) {
	ctx.currentModulus = ring.Copy(modulus)
}

// Degree returns the current degree of the ciphertext.
func (ctx *Ciphertext) Degree() uint64 {
	return uint64(len(ctx.value) - 1)
}

// Level returns the current level of the ciphertext.
func (ctx *Ciphertext) Level() uint64 {
	return uint64(len(ctx.value[0].Coeffs) - 1)
}

// Scale returns the current scale of the ciphertext.
func (ctx *Ciphertext) Scale() uint64 {
	return ctx.scale
}

// SetScale assigns a new scale to the ciphertext.
// Assagning a new scale to a ciphertext can be used
// to obtain free division or multiplication by 2^n.
func (ctx *Ciphertext) SetScale(scale uint64) {
	ctx.scale = scale
}

// IsNTT returns true if the ciphertext is in NTT form,
// else returns false.
func (ctx *Ciphertext) IsNTT() bool {
	return ctx.isNTT
}

// SetIsNTT assigns a new value to the isNTT variable of the
// ciphertext.
func (ctx *Ciphertext) SetIsNTT(isNTT bool) {
	ctx.isNTT = isNTT
}

// NTT computes the NTT of the ciphertext and copies it on the receiver element.
// Can only be used if the reference ciphertext is not already in the NTT domain.
func (ctx *Ciphertext) NTT(ckkscontext *CkksContext, ct0 CkksElement) {

	if ctx.isNTT != true {
		for i := range ct0.Value() {
			ckkscontext.contextLevel[ctx.Level()].NTT(ctx.value[i], ct0.Value()[i])
		}
		ct0.SetIsNTT(true)
	}
}

// InvNTT computes the inverse NTT of the ciphertext and copies it on the receiver element.
// Can only be used if the reference ciphertext is in the NTT domain.
func (ctx *Ciphertext) InvNTT(ckkscontext *CkksContext, ct0 CkksElement) {

	if ctx.isNTT != false {
		for i := range ct0.Value() {
			ckkscontext.contextLevel[ctx.Level()].InvNTT(ctx.value[i], ct0.Value()[i])
		}
		ct0.SetIsNTT(false)
	}
}

// NewRandoMCiphertext generates a new uniformely distributed ciphertext of degree, level and scale.
func (ckkscontext *CkksContext) NewRandomCiphertext(degree, level, scale uint64) *Ciphertext {
	ciphertext := new(Ciphertext)

	ciphertext.value = make([]*ring.Poly, degree+1)
	for i := uint64(0); i < degree+1; i++ {
		ciphertext.value[i] = ckkscontext.contextLevel[level].NewUniformPoly()
	}

	ciphertext.scale = scale
	ciphertext.currentModulus = ring.Copy(ckkscontext.contextLevel[level].ModulusBigint)
	ciphertext.isNTT = true

	return ciphertext
}

// CopyNew generates a new copy of the reference ciphertext, with the
// same value and same parameters.
func (ctx *Ciphertext) CopyNew() CkksElement {

	ctxCopy := new(Ciphertext)

	ctxCopy.value = make([]*ring.Poly, ctx.Degree()+1)
	for i := range ctx.value {
		ctxCopy.value[i] = ctx.value[i].CopyNew()
	}

	ctx.CopyParams(ctxCopy)

	return ctxCopy
}

// Copy copies the reference ciphertext and its parameters on the receiver ciphertext.
func (ctx *Ciphertext) Copy(ctxCopy CkksElement) error {

	for i := range ctxCopy.Value() {
		ctx.value[i].Copy(ctxCopy.Value()[i])
	}

	ctx.CopyParams(ctxCopy)

	return nil
}

// CopyParams copies the reference ciphertext parameters on the
// receiver ciphertext.
func (ctx *Ciphertext) CopyParams(ckkselement CkksElement) {
	ckkselement.SetCurrentModulus(ctx.CurrentModulus())
	ckkselement.SetScale(ctx.Scale())
	ckkselement.SetIsNTT(ctx.IsNTT())
}

