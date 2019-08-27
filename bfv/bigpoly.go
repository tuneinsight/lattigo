package bfv

import (
	"github.com/lca1/lattigo/ring"
)

type BigPoly struct {
	value []*ring.Poly
	isNTT bool
}

type BfvElement interface {
	Value() []*ring.Poly
	SetValue([]*ring.Poly)
	Resize(*BfvContext, uint64)
	CopyNew() BfvElement
	Copy(BfvElement) error
	Degree() uint64
	NTT(*BfvContext, BfvElement) error
	InvNTT(*BfvContext, BfvElement) error
	IsNTT() bool
	SetIsNTT(bool)
}
