// Package rlwe implements the generic cryptographic functionalities and operations that are common to R-LWE schemes.
// The other implemented schemes extend this package with their specific operations and structures.
package rlwe

import (
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
)

type EvaluatorProvider interface {
	GetBuffQP() [6]ringqp.Poly
	GetBuffCt() *Ciphertext
	GetBuffDecompQP() []ringqp.Poly
	DecomposeNTT(level, levelP, pCount int, c1 ring.Poly, isNTT bool, BuffDecompQP []ringqp.Poly)
	CheckAndGetGaloisKey(galEl uint64) (evk *GaloisKey, err error)
	GadgetProductLazy(levelQ int, cx ring.Poly, gadgetCt *GadgetCiphertext, ct *Element[ringqp.Poly]) (err error)
	GadgetProductHoistedLazy(levelQ int, BuffQPDecompQP []ringqp.Poly, gadgetCt *GadgetCiphertext, ct *Element[ringqp.Poly]) (err error)
	AutomorphismHoistedLazy(levelQ int, ctIn *Ciphertext, c1DecompQP []ringqp.Poly, galEl uint64, ctQP *Element[ringqp.Poly]) (err error)
	ModDownQPtoQNTT(levelQ, levelP int, p1Q, p1P, p2Q ring.Poly)
	AutomorphismIndex(uint64) []uint64
}
