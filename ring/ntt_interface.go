package ring

// NumberTheoreticTransformer is an interface to provide
// flexibility on what type of NTT is used by the struct Ring.
type NumberTheoreticTransformer interface {
	Forward(r *Ring, p1, p2 *Poly)
	ForwardLazy(r *Ring, p1, p2 *Poly)
	Backward(r *Ring, p1, p2 *Poly)
	BackwardLazy(r *Ring, p1, p2 *Poly)
	ForwardVec(table *Table, p1, p2 []uint64)
	ForwardLazyVec(table *Table, p1, p2 []uint64)
	BackwardVec(table *Table, p1, p2 []uint64)
	BackwardLazyVec(table *Table, p1, p2 []uint64)
}

// NumberTheoreticTransformerStandard computes the standard nega-cyclic NTT in the ring Z[X]/(X^N+1).
type NumberTheoreticTransformerStandard struct {
}

// Forward writes the forward NTT in Z[X]/(X^N+1) of p1 on p2.
func (rntt NumberTheoreticTransformerStandard) Forward(r *Ring, p1, p2 *Poly) {
	for x := range r.Tables[:r.level+1] {
		NTT(r.Tables[x], p1.Coeffs[x], p2.Coeffs[x])
	}
}

// ForwardLazy writes the forward NTT in Z[X]/(X^N+1) of p1 on p2.
// Returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerStandard) ForwardLazy(r *Ring, p1, p2 *Poly) {
	for x := range r.Tables[:r.level+1] {
		NTTLazy(r.Tables[x], p1.Coeffs[x], p2.Coeffs[x])
	}
}

// Backward writes the backward NTT in Z[X]/(X^N+1) on p2.
func (rntt NumberTheoreticTransformerStandard) Backward(r *Ring, p1, p2 *Poly) {
	for x := range r.Tables[:r.level+1] {
		InvNTT(r.Tables[x], p1.Coeffs[x], p2.Coeffs[x])
	}
}

// BackwardLazy writes the backward NTT in Z[X]/(X^N+1) on p2.
// Returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerStandard) BackwardLazy(r *Ring, p1, p2 *Poly) {
	for x := range r.Tables[:r.level+1] {
		InvNTTLazy(r.Tables[x], p1.Coeffs[x], p2.Coeffs[x])
	}
}

// ForwardVec writes the forward NTT in Z[X]/(X^N+1) of p1 on p2.
func (rntt NumberTheoreticTransformerStandard) ForwardVec(table *Table, p1, p2 []uint64) {
	NTT(table, p1, p2)
}

// ForwardLazyVec writes the forward NTT in Z[X]/(X^N+1) of p1 on p2.
// Returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerStandard) ForwardLazyVec(table *Table, p1, p2 []uint64) {
	NTTLazy(table, p1, p2)
}

// BackwardVec writes the backward NTT in Z[X]/(X^N+1) of p1 on p2.
func (rntt NumberTheoreticTransformerStandard) BackwardVec(table *Table, p1, p2 []uint64) {
	InvNTT(table, p1, p2)
}

// BackwardLazyVec writes the backward NTT in Z[X]/(X^N+1) p1 on p2.
// Returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerStandard) BackwardLazyVec(table *Table, p1, p2 []uint64) {
	InvNTTLazy(table, p1, p2)
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
	for x := range r.Tables[:r.level+1] {
		NTTConjugateInvariant(r.Tables[x], p1.Coeffs[x], p2.Coeffs[x])
	}
}

// ForwardLazy writes the forward NTT in Z[X+X^-1]/(X^2N+1) on p2.
// Returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerConjugateInvariant) ForwardLazy(r *Ring, p1, p2 *Poly) {
	for x := range r.Tables[:r.level+1] {
		NTTConjugateInvariantLazy(r.Tables[x], p1.Coeffs[x], p2.Coeffs[x])
	}
}

// Backward writes the backward NTT in Z[X+X^-1]/(X^2N+1) on p2.
func (rntt NumberTheoreticTransformerConjugateInvariant) Backward(r *Ring, p1, p2 *Poly) {
	for x := range r.Tables[:r.level+1] {
		InvNTTConjugateInvariant(r.Tables[x], p1.Coeffs[x], p2.Coeffs[x])
	}
}

// BackwardLazy writes the backward NTT in Z[X+X^-1]/(X^2N+1) on p2.
// Returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerConjugateInvariant) BackwardLazy(r *Ring, p1, p2 *Poly) {
	for x := range r.Tables[:r.level+1] {
		InvNTTConjugateInvariantLazy(r.Tables[x], p1.Coeffs[x], p2.Coeffs[x])
	}
}

// ForwardVec writes the forward NTT in Z[X+X^-1]/(X^2N+1) of p1 on p2.
func (rntt NumberTheoreticTransformerConjugateInvariant) ForwardVec(table *Table, p1, p2 []uint64) {
	NTTConjugateInvariant(table, p1, p2)
}

// ForwardLazyVec writes the forward NTT in Z[X+X^-1]/(X^2N+1) of p1 on p2.
// Returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerConjugateInvariant) ForwardLazyVec(table *Table, p1, p2 []uint64) {
	NTTConjugateInvariantLazy(table, p1, p2)
}

// BackwardVec writes the backward NTT in Z[X+X^-1]/(X^2N+1) of p1 on p2.
func (rntt NumberTheoreticTransformerConjugateInvariant) BackwardVec(table *Table, p1, p2 []uint64) {
	InvNTTConjugateInvariant(table, p1, p2)
}

// BackwardLazyVec writes the backward NTT in Z[X+X^-1]/(X^2N+1) of p1 on p2.
// Returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerConjugateInvariant) BackwardLazyVec(table *Table, p1, p2 []uint64) {
	InvNTTConjugateInvariantLazy(table, p1, p2)
}
