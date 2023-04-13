package rlwe

// RelinearizationKey is type of evaluation key used for ciphertext multiplication compactness.
// The Relinearization key encrypts s^{2} under s and is used to homomorphically re-encrypt the
// degree 2 term of a ciphertext (the term that decrypt with s^{2}) into a degree 1 term
// (a term that decrypts with s).
type RelinearizationKey struct {
	EvaluationKey
}

// NewRelinearizationKey allocates a new RelinearizationKey with zero coefficients.
func NewRelinearizationKey(params Parameters) *RelinearizationKey {
	return &RelinearizationKey{EvaluationKey: *NewEvaluationKey(params, params.MaxLevelQ(), params.MaxLevelP())}
}

// CopyNew creates a deep copy of the object and returns it.
func (rlk *RelinearizationKey) CopyNew() *RelinearizationKey {
	return &RelinearizationKey{EvaluationKey: *rlk.EvaluationKey.CopyNew()}
}

// Equal performs a deep equal.
func (rlk *RelinearizationKey) Equal(other *RelinearizationKey) bool {
	return rlk.EvaluationKey.Equal(&other.EvaluationKey)
}
