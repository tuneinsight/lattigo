package mkbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
)

// Add adds the cyphertexts component wise and expend their list of involved peers
func Add(c1 *MKCiphertext, c2 *MKCiphertext, out *MKCiphertext, ringQ *ring.Ring, params *bfv.Parameters) {

	eval := bfv.NewEvaluator(params, bfv.EvaluationKey{}) // declare evaluator to get access to homomorphic operations

	eval.Add(c1.ciphertexts.Element, c2.ciphertexts.Element, out.ciphertexts)

}

// MultSharedRelinKey will compute the homomorphic multiplication and relinearize the resulting cyphertext using pre computed Relin key
func MultSharedRelinKey(c1 *MKCiphertext, c2 *MKCiphertext, out *MKCiphertext, relinKey *MKRelinearizationKey, params *bfv.Parameters) {

	eval := bfv.NewEvaluator(params, bfv.EvaluationKey{})

	eval.Mul(c1.ciphertexts.Element.Ciphertext(), c2.ciphertexts.Element, out.ciphertexts)

	// Call Relin on the resulting ciphertext
	RelinearizationWithSharedRelinKey(relinKey, out)

}

// MultRelinDynamic will compute the homomorphic multiplication and relinearize the resulting cyphertext using dynamic relin
func MultRelinDynamic() {
	// TODO: implement multiplication

	// Call Relin alg 2
}
