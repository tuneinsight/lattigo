package ring

// UnfoldConjugateInvariantToStandard maps the compressed representation (N/2 coefficients)
// of Z_Q[X+X^-1]/(X^2N + 1) to full representation in Z_Q[X]/(X^2N+1).
// Requires degree(polyConjugateInvariant) = 2*degree(polyStandard).
// Requires that polyStandard and polyConjugateInvariant share the same moduli.
func (r Ring) UnfoldConjugateInvariantToStandard(polyConjugateInvariant, polyStandard Poly) {

	// Sanity check
	if 2*polyConjugateInvariant.N() != polyStandard.N() {
		panic("cannot UnfoldConjugateInvariantToStandard: Ring degree of polyConjugateInvariant must be twice the ring degree of polyStandard")
	}

	N := polyConjugateInvariant.N()

	for i := 0; i < r.level+1; i++ {
		tmp2, tmp1 := polyStandard.Coeffs[i], polyConjugateInvariant.Coeffs[i]
		copy(tmp2, tmp1)
		for idx, jdx := N-1, N; jdx < 2*N; idx, jdx = idx-1, jdx+1 {
			tmp2[jdx] = tmp1[idx]
		}
	}
}

// FoldStandardToConjugateInvariant folds [X]/(X^N+1) to [X+X^-1]/(X^N+1) in compressed form (N/2 coefficients).
// Requires degree(polyConjugateInvariant) = 2*degree(polyStandard).
// Requires that polyStandard and polyConjugateInvariant share the same moduli.
func (r Ring) FoldStandardToConjugateInvariant(polyStandard Poly, permuteNTTIndexInv []uint64, polyConjugateInvariant Poly) {

	// Sanity check
	if polyStandard.N() != 2*polyConjugateInvariant.N() {
		panic("cannot FoldStandardToConjugateInvariant: Ring degree of polyStandard must be 2N and ring degree of polyConjugateInvariant must be N")
	}

	N := r.N()

	level := r.level

	r.AutomorphismNTTWithIndex(polyStandard, permuteNTTIndexInv, polyConjugateInvariant)

	for i, s := range r.SubRings[:level+1] {
		s.Add(polyConjugateInvariant.Coeffs[i][:N], polyStandard.Coeffs[i][:N], polyConjugateInvariant.Coeffs[i][:N])
	}
}

// PadDefaultRingToConjugateInvariant converts a polynomial in Z[X]/(X^N +1) to a polynomial in Z[X+X^-1]/(X^2N+1).
func (r Ring) PadDefaultRingToConjugateInvariant(polyStandard Poly, IsNTT bool, polyConjugateInvariant Poly) {

	// Sanity check
	if polyConjugateInvariant.N() != 2*polyStandard.N() {
		panic("cannot PadDefaultRingToConjugateInvariant: polyConjugateInvariant degree must be twice the one of polyStandard")
	}

	N := polyStandard.N()

	for i := 0; i < r.level+1; i++ {
		qi := r.SubRings[i].Modulus

		copy(polyConjugateInvariant.Coeffs[i], polyStandard.Coeffs[i])

		tmp := polyConjugateInvariant.Coeffs[i]
		if IsNTT {
			for j := 0; j < N; j++ {
				tmp[N-j-1] = tmp[j]
			}
		} else {
			tmp[0] = 0
			for j := 1; j < N; j++ {
				tmp[N-j] = qi - tmp[j]
			}
		}
	}
}
