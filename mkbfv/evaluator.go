package mkbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
)

// MKEvaluator is a wrapper for the bfv evaluator
type MKEvaluator interface {
	Add(c1 *MKCiphertext, c2 *MKCiphertext, params *bfv.Parameters) *MKCiphertext
	MultSharedRelinKey(c1 *MKCiphertext, c2 *MKCiphertext, relinKey *MKRelinearizationKey, params *bfv.Parameters) *MKCiphertext
	MultRelinDynamic(c1 *MKCiphertext, c2 *MKCiphertext, evalKeys []*MKEvaluationKey, publicKeys []*MKPublicKey, params *bfv.Parameters) *MKCiphertext
}

type mkEvaluator struct {
	bfvEval bfv.Evaluator
	ringQ   *ring.Ring
}

// NewMKEvaluator returns an evaluator for the multi key bfv scheme.
func NewMKEvaluator(params *bfv.Parameters) MKEvaluator {
	r := GetRingQ(params)
	return &mkEvaluator{bfvEval: bfv.NewEvaluator(params, bfv.EvaluationKey{}), ringQ: r}
}

// Add adds the cyphertexts component wise and expend their list of involved peers
func (eval *mkEvaluator) Add(c1 *MKCiphertext, c2 *MKCiphertext, params *bfv.Parameters) *MKCiphertext {

	out := NewMKCiphertext(c1.peerIDs, eval.ringQ)

	eval.bfvEval.Add(c1.ciphertexts.Element, c2.ciphertexts.Element, out.ciphertexts)

	return out
}

// MultSharedRelinKey will compute the homomorphic multiplication and relinearize the resulting cyphertext using pre computed Relin key
func (eval *mkEvaluator) MultSharedRelinKey(c1 *MKCiphertext, c2 *MKCiphertext, relinKey *MKRelinearizationKey, params *bfv.Parameters) *MKCiphertext {

	out := NewMKCiphertext(c1.peerIDs, eval.ringQ)

	eval.bfvEval.Mul(c1.ciphertexts.Element.Ciphertext(), c2.ciphertexts.Element, out.ciphertexts)

	// Call Relin on the resulting ciphertext
	RelinearizationWithSharedRelinKey(relinKey, out)

	return out
}

// MultRelinDynamic will compute the homomorphic multiplication and relinearize the resulting cyphertext using dynamic relin
func (eval *mkEvaluator) MultRelinDynamic(c1 *MKCiphertext, c2 *MKCiphertext, evalKeys []*MKEvaluationKey, publicKeys []*MKPublicKey, params *bfv.Parameters) *MKCiphertext {

	out := NewMKCiphertext(c1.peerIDs, eval.ringQ)

	eval.bfvEval.Mul(c1.ciphertexts.Element.Ciphertext(), c2.ciphertexts.Element, out.ciphertexts)

	// Call Relin alg 2
	RelinearizationOnTheFly(evalKeys, publicKeys, out, params)

	return out
}
