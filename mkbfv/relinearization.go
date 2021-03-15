package mkbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
)

// RelinearizationWithSharedRelinKey implements the algorithm 1 in section 3.3.1 of the Chen paper
// It does relin using a precomputed shared key
func RelinearizationWithSharedRelinKey(relinKey *MKRelinearizationKey, ciphertext *MKCiphertext) {
	// TODO: implement Algorithm 1
}

// RelinearizationOnTheFly implements the algorithm 2 in section 3.3.1 of the Chen paper
// It does relin directly by linearizing each entry of the extended ciphertext and stores it in cPrime (of size k+1)
// There are (k+1)**2 ciphertexts, and k pairs of (evaluation keys Di,bi)
func RelinearizationOnTheFly(evaluationKeys []*MKEvaluationKey, publicKeys []*MKPublicKey, ciphertexts *MKCiphertext, params *bfv.Parameters) {

	ringQ := GetRingQ(params)
	k := uint64(len(evaluationKeys))
	res := make([]*ring.Poly, k+1)

	for i := uint64(0); i < k+1; i++ {
		res[i] = ringQ.NewPoly()
	}

	cipherParts := ciphertexts.ciphertexts.Value()

	c0Prime := cipherParts[0]
	tmpIJ := ringQ.NewPoly()
	tmpC0Prime := ringQ.NewPoly()
	tmpCIPrime := ringQ.NewPoly()

	for i := uint64(1); i <= k; i++ {
		ringQ.Add(cipherParts[i], cipherParts[(k+1)*i], res[i]) // first for loop in algo2
	}

	for i := uint64(1); i <= k; i++ { // second loop in alg2

		di0 := evaluationKeys[i].key[0]
		di1 := evaluationKeys[i].key[1]
		di2 := evaluationKeys[i].key[2]

		for j := uint64(1); j <= k; j++ {

			index := i*(k+1) + j

			Dot(GInverse(cipherParts[index], params), publicKeys[j].key[0], tmpIJ, ringQ) // line 6

			decomposedTmp := GInverse(tmpIJ, params) // inverse and matric mult (line 7)

			Dot(decomposedTmp, di0, tmpC0Prime, ringQ)
			Dot(decomposedTmp, di1, tmpCIPrime, ringQ)

			ringQ.Add(c0Prime, tmpC0Prime, c0Prime)
			ringQ.Add(res[i], tmpCIPrime, res[i])

			Dot(GInverse(cipherParts[index], params), di2, tmpIJ, ringQ) // inverse and dot product (line8)
			ringQ.Add(res[j], tmpIJ, res[j])

		}
	}

	res[0] = c0Prime
	ciphertexts.ciphertexts.SetValue(res)
}
