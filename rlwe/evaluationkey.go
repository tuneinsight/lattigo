package rlwe

// EvaluationKey is a public key indended to be used during the evaluation phase of a homomorphic circuit.
// It provides a one way public and non-interactive re-encryption from a ciphertext encrypted under `skIn`
// to a ciphertext encrypted under `skOut`.
//
// Such re-encryption is for example used for:
//
// - Homomorphic relinearization: re-encryption of a quadratic ciphertext (that requires (1, sk sk^2) to be decrypted)
// to a linear ciphertext (that required (1, sk) to be decrypted). In this case skIn = sk^2 an skOut = sk.
//
// - Homomorphic automorphisms: an automorphism in the ring Z[X]/(X^{N}+1) is defined as pi_k: X^{i} -> X^{i^k} with
// k coprime to 2N. Pi_sk is for exampled used during homomorphic slot rotations. Applying pi_k to a ciphertext encrypted
// under sk generates a new ciphertext encrypted under pi_k(sk), and an Evaluationkey skIn = pi_k(sk) to skOut = sk
// is used to bring it back to its original key.
type EvaluationKey struct {
	GadgetCiphertext
}

// NewEvaluationKey returns a new EvaluationKey with pre-allocated zero-value
func NewEvaluationKey(params Parameters, levelQ, levelP int) *EvaluationKey {
	return &EvaluationKey{GadgetCiphertext: *NewGadgetCiphertext(
		params,
		levelQ,
		levelP,
		params.DecompRNS(levelQ, levelP),
		params.DecompPw2(levelQ, levelP),
	)}
}

// CopyNew creates a deep copy of the target EvaluationKey and returns it.
func (evk *EvaluationKey) CopyNew() *EvaluationKey {
	return &EvaluationKey{GadgetCiphertext: *evk.GadgetCiphertext.CopyNew()}
}
