package mkbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
)

// MKEvaluator is a wrapper for the bfv evaluator
type MKEvaluator interface {
	Add(c1 *MKCiphertext, c2 *MKCiphertext, out *MKCiphertext, ringQ *ring.Ring, params *bfv.Parameters)
	MultSharedRelinKey(c1 *MKCiphertext, c2 *MKCiphertext, out *MKCiphertext, relinKey *MKRelinearizationKey, params *bfv.Parameters)
	MultRelinDynamic(c1 *MKCiphertext, c2 *MKCiphertext, out *MKCiphertext, evalKeys []*MKEvaluationKey, publicKeys []*MKPublicKey, params *bfv.Parameters)
}

type mkEvaluator struct {
	bfvEval bfv.Evaluator
}

// NewMKEvaluator returns an evaluator for the multi key bfv scheme.
func NewMKEvaluator(params *bfv.Parameters) MKEvaluator {

	return &mkEvaluator{bfvEval: bfv.NewEvaluator(params, bfv.EvaluationKey{})}
}

// Add adds the cyphertexts component wise and expend their list of involved peers
func (eval *mkEvaluator) Add(c1 *MKCiphertext, c2 *MKCiphertext, out *MKCiphertext, ringQ *ring.Ring, params *bfv.Parameters) {

	eval.bfvEval.Add(c1.ciphertexts.Element, c2.ciphertexts.Element, out.ciphertexts)

}

// MultSharedRelinKey will compute the homomorphic multiplication and relinearize the resulting cyphertext using pre computed Relin key
func (eval *mkEvaluator) MultSharedRelinKey(c1 *MKCiphertext, c2 *MKCiphertext, out *MKCiphertext, relinKey *MKRelinearizationKey, params *bfv.Parameters) {

	eval.bfvEval.Mul(c1.ciphertexts.Element.Ciphertext(), c2.ciphertexts.Element, out.ciphertexts)

	// Call Relin on the resulting ciphertext
	RelinearizationWithSharedRelinKey(relinKey, out)

}

// MultRelinDynamic will compute the homomorphic multiplication and relinearize the resulting cyphertext using dynamic relin
func (eval *mkEvaluator) MultRelinDynamic(c1 *MKCiphertext, c2 *MKCiphertext, out *MKCiphertext, evalKeys []*MKEvaluationKey, publicKeys []*MKPublicKey, params *bfv.Parameters) {

	eval.bfvEval.Mul(c1.ciphertexts.Element.Ciphertext(), c2.ciphertexts.Element, out.ciphertexts)

	// Call Relin alg 2
	RelinearizationOnTheFly(evalKeys, publicKeys, out, params)
}
