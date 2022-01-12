package ring

// NumberTheoreticTransformer is an interface to provide
// flexibility on what type of NTT is used by the struct Ring.
type NumberTheoreticTransformer interface {
	Forward(r *Ring, p1, p2 *Poly)
	ForwardLvl(r *Ring, level int, p1, p2 *Poly)
	ForwardLazy(r *Ring, p1, p2 *Poly)
	ForwardLazyLvl(r *Ring, level int, p1, p2 *Poly)
	Backward(r *Ring, p1, p2 *Poly)
	BackwardLvl(r *Ring, level int, p1, p2 *Poly)
	BackwardLazy(r *Ring, p1, p2 *Poly)
	BackwardLazyLvl(r *Ring, level int, p1, p2 *Poly)
	ForwardVec(r *Ring, level int, p1, p2 []uint64)
	ForwardLazyVec(r *Ring, level int, p1, p2 []uint64)
	BackwardVec(r *Ring, level int, p1, p2 []uint64)
	BackwardLazyVec(r *Ring, level int, p1, p2 []uint64)
}

// NumberTheoreticTransformerStandard computes the standard nega-cyclic NTT in the ring Z[X]/(X^N+1).
type NumberTheoreticTransformerStandard struct {
}

// Forward writes the forward NTT in Z[X]/(X^N+1) of p1 on p2.
func (rntt NumberTheoreticTransformerStandard) Forward(r *Ring, p1, p2 *Poly) {
	for x := range r.Modulus {
		NTT(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsi[x], r.Modulus[x], r.MredParams[x], r.BredParams[x])
	}
}

// ForwardLvl writes the forward NTT in Z[X]/(X^N+1) of p1 on p2.
// Only computes the NTT for the first level+1 moduli.
func (rntt NumberTheoreticTransformerStandard) ForwardLvl(r *Ring, level int, p1, p2 *Poly) {
	for x := 0; x < level+1; x++ {
		NTT(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsi[x], r.Modulus[x], r.MredParams[x], r.BredParams[x])
	}
}

// ForwardLazy writes the forward NTT in Z[X]/(X^N+1) of p1 on p2.
// Returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerStandard) ForwardLazy(r *Ring, p1, p2 *Poly) {
	for x := range r.Modulus {
		NTTLazy(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsi[x], r.Modulus[x], r.MredParams[x], r.BredParams[x])
	}
}

// ForwardLazyLvl writes the forward NTT in Z[X]/(X^N+1) of p1 on p2.
// Only computes the NTT for the first level+1 moduli and returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerStandard) ForwardLazyLvl(r *Ring, level int, p1, p2 *Poly) {
	for x := 0; x < level+1; x++ {
		NTTLazy(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsi[x], r.Modulus[x], r.MredParams[x], r.BredParams[x])
	}
}

// Backward writes the backward NTT in Z[X]/(X^N+1) on p2.
func (rntt NumberTheoreticTransformerStandard) Backward(r *Ring, p1, p2 *Poly) {
	for x := range r.Modulus {
		InvNTT(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsiInv[x], r.NttNInv[x], r.Modulus[x], r.MredParams[x])
	}
}

// BackwardLvl writes the backward NTT in Z[X]/(X^N+1) on p2.
// Only computes the NTT for the first level+1 moduli.
func (rntt NumberTheoreticTransformerStandard) BackwardLvl(r *Ring, level int, p1, p2 *Poly) {
	for x := 0; x < level+1; x++ {
		InvNTT(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsiInv[x], r.NttNInv[x], r.Modulus[x], r.MredParams[x])
	}
}

// BackwardLazy writes the backward NTT in Z[X]/(X^N+1) on p2.
// Returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerStandard) BackwardLazy(r *Ring, p1, p2 *Poly) {
	for x := range r.Modulus {
		InvNTTLazy(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsiInv[x], r.NttNInv[x], r.Modulus[x], r.MredParams[x])
	}
}

// BackwardLazyLvl writes the backward NTT in Z[X]/(X^N+1) on p2.
// Only computes the NTT for the first level+1 moduli and returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerStandard) BackwardLazyLvl(r *Ring, level int, p1, p2 *Poly) {
	for x := 0; x < level+1; x++ {
		InvNTTLazy(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsiInv[x], r.NttNInv[x], r.Modulus[x], r.MredParams[x])
	}
}

// ForwardVec writes the forward NTT in Z[X]/(X^N+1) of the i-th level of p1 on the i-th level of p2.
func (rntt NumberTheoreticTransformerStandard) ForwardVec(r *Ring, level int, p1, p2 []uint64) {
	NTT(p1, p2, r.N, r.NttPsi[level], r.Modulus[level], r.MredParams[level], r.BredParams[level])
}

// ForwardLazyVec writes the forward NTT in Z[X]/(X^N+1) of the i-th level of p1 on the i-th level of p2.
// Returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerStandard) ForwardLazyVec(r *Ring, level int, p1, p2 []uint64) {
	NTTLazy(p1, p2, r.N, r.NttPsi[level], r.Modulus[level], r.MredParams[level], r.BredParams[level])
}

// BackwardVec writes the backward NTT in Z[X]/(X^N+1) of the i-th level of p1 on the i-th level of p2.
func (rntt NumberTheoreticTransformerStandard) BackwardVec(r *Ring, level int, p1, p2 []uint64) {
	InvNTT(p1, p2, r.N, r.NttPsiInv[level], r.NttNInv[level], r.Modulus[level], r.MredParams[level])
}

// BackwardLazyVec writes the backward NTT in Z[X]/(X^N+1) of the i-th level of p1 on the i-th level of p2.
// Returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerStandard) BackwardLazyVec(r *Ring, level int, p1, p2 []uint64) {
	InvNTTLazy(p1, p2, r.N, r.NttPsiInv[level], r.NttNInv[level], r.Modulus[level], r.MredParams[level])
}

// NumberTheoreticTransformerConjugateInvariant computes the NTT in the ring Z[X+X^-1]/(X^2N+1).
// Z[X+X^-1]/(X^2N+1) is a closed sub-ring of Z[X]/(X^2N+1). Note that the input polynomial only needs to be size N
// since the right half does not provide any additional information.
// See "Approximate Homomorphic Encryption over the Conjugate-invariant Ring", https://eprint.iacr.org/2018/952.
// The implemented approach is more efficient than the one proposed in the referenced work.
// It avoids the linear map Z[X + X^-1]/(X^2N + 1) <-> Z[X]/(X^N - 1) by instead directly computing the left
// half of the NTT of Z[X + X^-1]/(X^2N + 1) since the right half provides no additional information, which
// allows to (re)use nega-cyclic NTT.
type NumberTheoreticTransformerConjugateInvariant struct {
}

// Forward writes the forward NTT in Z[X+X^-1]/(X^2N+1) on p2.
func (rntt NumberTheoreticTransformerConjugateInvariant) Forward(r *Ring, p1, p2 *Poly) {
	for x := range r.Modulus {
		NTTConjugateInvariant(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsi[x], r.Modulus[x], r.MredParams[x], r.BredParams[x])
	}
}

// ForwardLvl writes the forward NTT in Z[X+X^-1]/(X^2N+1) on p2.
// Only computes the NTT for the first level+1 moduli.
func (rntt NumberTheoreticTransformerConjugateInvariant) ForwardLvl(r *Ring, level int, p1, p2 *Poly) {
	for x := 0; x < level+1; x++ {
		NTTConjugateInvariant(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsi[x], r.Modulus[x], r.MredParams[x], r.BredParams[x])
	}
}

// ForwardLazy writes the forward NTT in Z[X+X^-1]/(X^2N+1) on p2.
// Returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerConjugateInvariant) ForwardLazy(r *Ring, p1, p2 *Poly) {
	for x := range r.Modulus {
		NTTConjugateInvariantLazy(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsi[x], r.Modulus[x], r.MredParams[x], r.BredParams[x])
	}
}

// ForwardLazyLvl writes the forward NTT in Z[X+X^-1]/(X^2N+1) on p2.
// Only computes the NTT for the first level+1 moduli and returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerConjugateInvariant) ForwardLazyLvl(r *Ring, level int, p1, p2 *Poly) {
	for x := 0; x < level+1; x++ {
		NTTConjugateInvariantLazy(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsi[x], r.Modulus[x], r.MredParams[x], r.BredParams[x])
	}
}

// Backward writes the backward NTT in Z[X+X^-1]/(X^2N+1) on p2.
func (rntt NumberTheoreticTransformerConjugateInvariant) Backward(r *Ring, p1, p2 *Poly) {
	for x := range r.Modulus {
		InvNTTConjugateInvariant(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsiInv[x], r.NttNInv[x], r.Modulus[x], r.MredParams[x])
	}
}

// BackwardLvl writes the backward NTT in Z[X+X^-1]/(X^2N+1) on p2.
// Only computes the NTT for the first level+1 moduli.
func (rntt NumberTheoreticTransformerConjugateInvariant) BackwardLvl(r *Ring, level int, p1, p2 *Poly) {
	for x := 0; x < level+1; x++ {
		InvNTTConjugateInvariant(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsiInv[x], r.NttNInv[x], r.Modulus[x], r.MredParams[x])
	}
}

// BackwardLazy writes the backward NTT in Z[X+X^-1]/(X^2N+1) on p2.
// Returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerConjugateInvariant) BackwardLazy(r *Ring, p1, p2 *Poly) {
	for x := range r.Modulus {
		InvNTTConjugateInvariantLazy(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsiInv[x], r.NttNInv[x], r.Modulus[x], r.MredParams[x])
	}
}

// BackwardLazyLvl writes the backward NTT in Z[X+X^-1]/(X^2N+1) on p2.
// Only computes the NTT for the first level+1 moduli and returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerConjugateInvariant) BackwardLazyLvl(r *Ring, level int, p1, p2 *Poly) {
	for x := 0; x < level+1; x++ {
		InvNTTConjugateInvariantLazy(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsiInv[x], r.NttNInv[x], r.Modulus[x], r.MredParams[x])
	}
}

// ForwardVec writes the forward NTT in Z[X+X^-1]/(X^2N+1) of the i-th level of p1 on the i-th level of p2.
func (rntt NumberTheoreticTransformerConjugateInvariant) ForwardVec(r *Ring, level int, p1, p2 []uint64) {
	NTTConjugateInvariant(p1, p2, r.N, r.NttPsi[level], r.Modulus[level], r.MredParams[level], r.BredParams[level])
}

// ForwardLazyVec writes the forward NTT in Z[X+X^-1]/(X^2N+1) of the i-th level of p1 on the i-th level of p2.
// Returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerConjugateInvariant) ForwardLazyVec(r *Ring, level int, p1, p2 []uint64) {
	NTTConjugateInvariantLazy(p1, p2, r.N, r.NttPsi[level], r.Modulus[level], r.MredParams[level], r.BredParams[level])
}

// BackwardVec writes the backward NTT in Z[X+X^-1]/(X^2N+1) of the i-th level of p1 on the i-th level of p2.
func (rntt NumberTheoreticTransformerConjugateInvariant) BackwardVec(r *Ring, level int, p1, p2 []uint64) {
	InvNTTConjugateInvariant(p1, p2, r.N, r.NttPsiInv[level], r.NttNInv[level], r.Modulus[level], r.MredParams[level])
}

// BackwardLazyVec writes the backward NTT in Z[X+X^-1]/(X^2N+1) of the i-th level of p1 on the i-th level of p2.
// Returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerConjugateInvariant) BackwardLazyVec(r *Ring, level int, p1, p2 []uint64) {
	InvNTTConjugateInvariantLazy(p1, p2, r.N, r.NttPsiInv[level], r.NttNInv[level], r.Modulus[level], r.MredParams[level])
}
