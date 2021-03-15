package mkbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
)

// RelinearizationWithSharedRelinKey implements the algorithm 1 in section 3.3.1 of the Chen paper
// It does relin using a precomputed shared key
func RelinearizationWithSharedRelinKey(relinKey *MKRelinearizationKey, ciphertext *MKCiphertext) {
	// TODO: implement Algorithm 1
}

// RelinearizationOnTheFly implements the algorithm 2 in section 3.3.1 of the Chen paper
// It does relin directly by linearizing each entry of the extended ciphertext and stores it in cPrime (of size k+1)
// There are (k+1)**2 ciphertexts, and k pairs of (evaluation keys Di,bi)
func RelinearizationOnTheFly(evaluationKeys []*MKEvaluationKey, publicKeys []*MKPublicKey, ciphertexts [][]*MKCiphertext, cPrime []*MKCiphertext, params *bfv.Parameters, evaluator *mkEvaluator) {

	ringQP := GetRingQP(params)
	ringQ := GetRingQ(params)
	k := uint64(len(evaluationKeys))
	c0Prime := ciphertexts[0][0]
	cPrime2D := make([][]*MKCiphertext, k+1)
	cjPrime := make([]*MKCiphertext, k+1) // why is there a j_prime?
	cPrime[0] = c0Prime                   // is equality like this okay??
	for i := uint64(1); i < k; i++ {
		evaluator.Add(ciphertexts[0][i], ciphertexts[i][0], cPrime[i], ringQ, params) // first for loop in algo2
	}
	for i := uint64(1); i < k; i++ {
		di0 := evaluationKeys[i].key[0]
		di1 := evaluationKeys[i].key[1]
		di2 := evaluationKeys[i].key[2]
		for j := uint64(1); j < k; j++ {
			Dot(GInverse(ciphertexts[i][j].ciphertexts.Element.Value(), params), publicKeys[j].key[0], cPrime2D[i][j], ringQ) // need to exact the polynomial from ci,j and from cPrime[i][j]
			evaluator.Add(c0Prime, Dot(GInverse(cPrime2D[i][j].ciphertexts.Element.Value(), params), di0, c0Prime, ringQ))
			evaluator.Add(cPrime[i], Dot(GInverse(cPrime2D[i][j].ciphertexts.Element.Value(), params), di1, cPrime[i], ringQ))
			evaluator.Add(cjPrime[j], Dot(GInverse(ciphertexts[i][j].ciphertexts.Element.Value()), di2, cjPrime[j], ringQ))
		}
	}

}
