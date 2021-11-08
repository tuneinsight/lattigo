package ring

import (
	"github.com/ldsec/lattigo/v2/utils"
)

// UnfoldConjugateInvariantToStandard maps the compressed representation of Z_Q[X+X^-1]/(X^2N + 1) to full representation in Z_Q[X]/(X^2N+1).
// Requires degree(polyConjugateInvariant) = 2*degree(polyStd).
// Requires that polyStd and polyConjugateInvariant share the same moduli.
func (r *Ring) UnfoldConjugateInvariantToStandard(level int, polyConjugateInvariant, polyStd *Poly) {

	if 2*len(polyConjugateInvariant.Coeffs[0]) != len(polyStd.Coeffs[0]) {
		panic("Ring degree of polyConjugateInvariant must be twice the ring degree of polyStd")
	}

	N := len(polyConjugateInvariant.Coeffs[0])

	for i := 0; i < level+1; i++ {
		tmp2, tmp1 := polyStd.Coeffs[i], polyConjugateInvariant.Coeffs[i]
		copy(tmp2, tmp1)
		for idx, jdx := N-1, N; jdx < 2*N; idx, jdx = idx-1, jdx+1 {
			tmp2[jdx] = tmp1[idx]
		}
	}
}

// FoldStandardToConjugateInvariant folds [X] to [X+X^-1] in compressed form.
// Requires degree(polyConjugateInvariant) = 2*degree(polyStd).
// Requires that polyStd and polyConjugateInvariant share the same moduli.
func (r *Ring) FoldStandardToConjugateInvariant(level int, polyStandard *Poly, permuteNTTIndexInv []uint64, polyConjugateInvariant *Poly) {

	if len(polyStandard.Coeffs[0]) != 2*len(polyConjugateInvariant.Coeffs[0]) {
		panic("Ring degree of p2 must be 2N and ring degree of p1 must be N")
	}

	r.PermuteNTTWithIndexLvl(level, polyStandard, permuteNTTIndexInv, polyConjugateInvariant)
	for i := 0; i < level+1; i++ {
		AddVec(polyConjugateInvariant.Coeffs[i][:r.N], polyStandard.Coeffs[i][:r.N], polyConjugateInvariant.Coeffs[i][:r.N], r.Modulus[i])
	}
}

// PadDefaultRingToConjuateInvariant converts a polynomial in Z[X]/(X^N +1) to a polynomial in Z[X+X^-1]/(X^2N+1).
// Conversion assumes polynomials are outside of the NTT domain.
// Default ring Z[X]/(X^N +1) and ConjugateInvariant ring Z[X+X^-1]/(X^2N+1) must share the same moduli.
func PadDefaultRingToConjuateInvariant(p1 *Poly, ringStd, ringConjInv *Ring, p2 *Poly) {

	level := utils.MinInt(p1.Level(), p2.Level())
	n := len(p1.Coeffs[0])

	for i := 0; i < level+1; i++ {
		qi := ringStd.Modulus[i]
		if qi != ringConjInv.Modulus[i] {
			panic("p1 and p2 rings must share the same moduli")
		}

		if len(p2.Coeffs[i]) != 2*len(p1.Coeffs[i]) {
			panic("p2 degree must be twice the one of p1")
		}

		tmp := p2.Coeffs[i]
		tmp[0] = 0
		for j := 1; j < n; j++ {
			tmp[n-j] = qi - tmp[j]
		}
	}
}
