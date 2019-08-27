package ckks

import (
	"github.com/lca1/lattigo/ring"
)

type BigPoly struct {
	value          []*ring.Poly
	scale          uint64
	currentModulus *ring.Int
	isNTT          bool
}

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
