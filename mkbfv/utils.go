package mkbfv

import "github.com/ldsec/lattigo/v2/ring"

// Dot computes the dot product of two decomposed polynomials in R_Q^d
func Dot(decpoly1 *MKDecomposedPoly, decpoly2 *MKDecomposedPoly, accPoly *ring.Poly, ringQP *ring.Ring, beta uint64) {
	for l := uint64(0); l < beta; l++ {
		ringQP.MulCoeffsMontgomeryAndAdd(decpoly1.poly[l], decpoly2.poly[l], accPoly)
	}
}
