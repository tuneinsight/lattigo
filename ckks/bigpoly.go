package ckks

import (
	"github.com/ldsec/lattigo/ring"
)

// BigPoly is a common struct between plaintexts and ciphertexts. It stores a value
// as a slice of polynomials, and an isNTT flag indicatig if the element is in the NTT domain.
type BigPoly struct {
	value          []*ring.Poly
	scale          uint64
	currentModulus *ring.Int
	isNTT          bool
}

// CkksElement is an interface implementing common methodes on plaintext and ciphertexts.
type CkksElement interface {
	Value() []*ring.Poly
	SetValue([]*ring.Poly)
	Scale() uint64
	SetScale(uint64)
	CurrentModulus() *ring.Int
	SetCurrentModulus(*ring.Int)
	CopyParams(CkksElement)
	Resize(*CkksContext, uint64)
	CopyNew() CkksElement
	Copy(CkksElement) error
	Degree() uint64
	Level() uint64
	NTT(*CkksContext, CkksElement)
	InvNTT(*CkksContext, CkksElement)
	IsNTT() bool
	SetIsNTT(bool)
}
