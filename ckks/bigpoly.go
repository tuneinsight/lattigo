package ckks

import (
	"github.com/lca1/lattigo/ring"
)

type BigPoly struct {
	value          []*ring.Poly
	ckkscontext    *CkksContext
	scale          uint64
	currentModulus *ring.Int
	isNTT          bool
	isComplex      bool
}

type CkksElement interface {
	Value() []*ring.Poly
	SetValue([]*ring.Poly)
	CkksContext() *CkksContext
	Scale() uint64
	SetScale(uint64)
	CurrentModulus() *ring.Int
	SetCurrentModulus(*ring.Int)
	CopyParams(CkksElement)
	Resize(uint64)
	CopyNew() CkksElement
	Copy(CkksElement) error
	Degree() uint64
	Level() uint64
	NTT(CkksElement)
	InvNTT(CkksElement)
	IsNTT() bool
	SetIsNTT(bool)
}
