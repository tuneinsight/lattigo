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

// Equals returs true if the to objects are equal.
func (rlk *RelinearizationKey) Equals(other *RelinearizationKey) bool {
	return rlk.EvaluationKey.Equals(&other.EvaluationKey)
}

// CopyNew creates a deep copy of the object and returns it.
func (rlk *RelinearizationKey) CopyNew() *RelinearizationKey {
	return &RelinearizationKey{EvaluationKey: *rlk.EvaluationKey.CopyNew()}
}

// MarshalBinarySize returns the length in bytes that the object requires to be marshaled.
func (rlk *RelinearizationKey) MarshalBinarySize() (dataLen int) {
	return rlk.EvaluationKey.MarshalBinarySize()
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (rlk *RelinearizationKey) MarshalBinary() (data []byte, err error) {
	return rlk.EvaluationKey.MarshalBinary()
}

// Read encodes the object into a binary form on a preallocated slice of bytes
// and returns the number of bytes written.
func (rlk *RelinearizationKey) Read(data []byte) (ptr int, err error) {
	return rlk.EvaluationKey.Read(data)
}

// UnmarshalBinary decodes a slice of bytes generated by MarshalBinary
// or Read on the object.
func (rlk *RelinearizationKey) UnmarshalBinary(data []byte) (err error) {
	_, err = rlk.Write(data)
	return
}

// Write decodes a slice of bytes generated by MarshalBinary or
// Read on the object and returns the number of bytes read.
func (rlk *RelinearizationKey) Write(data []byte) (ptr int, err error) {
	return rlk.EvaluationKey.Write(data)
}