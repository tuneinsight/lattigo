package bfv

import (
	"github.com/lca1/lattigo/ring"
)

// BigPoly is a common struct between plaintexts and ciphertexts. It stores a value
// as a slice of polynomials, and an isNTT flag indicatig if the element is in the NTT domain.
type BigPoly struct {
	value []*ring.Poly
	isNTT bool
}

// BfvElement is an interface implementing common methodes on plaintext and ciphertexts.
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
