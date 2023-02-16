package ring

import (
	"github.com/tuneinsight/lattigo/v4/utils"
)

// UnfoldConjugateInvariantToStandard maps the compressed representation (N/2 coefficients)
// of Z_Q[X+X^-1]/(X^2N + 1) to full representation in Z_Q[X]/(X^2N+1).
// Requires degree(polyConjugateInvariant) = 2*degree(polyStd).
// Requires that polyStd and polyConjugateInvariant share the same moduli.
func (r *Ring) UnfoldConjugateInvariantToStandard(polyConjugateInvariant, polyStd *Poly) {

	if 2*len(polyConjugateInvariant.Coeffs[0]) != len(polyStd.Coeffs[0]) {
		panic("cannot UnfoldConjugateInvariantToStandard: Ring degree of polyConjugateInvariant must be twice the ring degree of polyStd")
	}

	N := len(polyConjugateInvariant.Coeffs[0])

	for i := 0; i < r.level+1; i++ {
		tmp2, tmp1 := polyStd.Coeffs[i], polyConjugateInvariant.Coeffs[i]
		copy(tmp2, tmp1)
		for idx, jdx := N-1, N; jdx < 2*N; idx, jdx = idx-1, jdx+1 {
			tmp2[jdx] = tmp1[idx]
		}
	}
}

// FoldStandardToConjugateInvariant folds [X]/(X^N+1) to [X+X^-1]/(X^N+1) in compressed form (N/2 coefficients).
// Requires degree(polyConjugateInvariant) = 2*degree(polyStd).
// Requires that polyStd and polyConjugateInvariant share the same moduli.
func (r *Ring) FoldStandardToConjugateInvariant(polyStandard *Poly, permuteNTTIndexInv []uint64, polyConjugateInvariant *Poly) {

	if len(polyStandard.Coeffs[0]) != 2*len(polyConjugateInvariant.Coeffs[0]) {
		panic("cannot FoldStandardToConjugateInvariant: Ring degree of p2 must be 2N and ring degree of p1 must be N")
	}

	N := r.N()

	level := r.level

	r.PermuteNTTWithIndex(polyStandard, permuteNTTIndexInv, polyConjugateInvariant)

	for i, s := range r.SubRings[:level+1] {
		s.Add(polyConjugateInvariant.Coeffs[i][:N], polyStandard.Coeffs[i][:N], polyConjugateInvariant.Coeffs[i][:N])
	}
}

// PadDefaultRingToConjugateInvariant converts a polynomial in Z[X]/(X^N +1) to a polynomial in Z[X+X^-1]/(X^2N+1).
func PadDefaultRingToConjugateInvariant(p1 *Poly, ringQ *Ring, IsNTT bool, p2 *Poly) {

	if p1 == p2 {
		panic("cannot PadDefaultRingToConjugateInvariant: p1 == p2 but method cannot be used in place")
	}

	level := utils.MinInt(p1.Level(), p2.Level())
	n := len(p1.Coeffs[0])

	for i := 0; i < level+1; i++ {
		qi := ringQ.SubRings[i].Modulus

		if len(p2.Coeffs[i]) != 2*len(p1.Coeffs[i]) {
			panic("cannot PadDefaultRingToConjugateInvariant: p2 degree must be twice the one of p1")
		}

		copy(p2.Coeffs[i], p1.Coeffs[i])

		tmp := p2.Coeffs[i]
		if IsNTT {
			for j := 0; j < n; j++ {
				tmp[n-j-1] = tmp[j]
			}
		} else {
			tmp[0] = 0
			for j := 1; j < n; j++ {
				tmp[n-j] = qi - tmp[j]
			}
		}
	}
}
