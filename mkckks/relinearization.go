package mkckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
)

// RelinearizationWithSharedRelinKey implements the algorithm 1 in section 3.3.1 of the Chen paper
// It does relin using a precomputed shared key
func RelinearizationWithSharedRelinKey(relinKey *MKRelinearizationKey, ciphertext *MKCiphertext) {
	// TODO: implement Algorithm 1
}

/*
// RelinearizationOnTheFly implements the algorithm 2 in section 3.3.1 of the Chen paper
// It does relin directly by linearizing each entry of the extended ciphertext and stores it in cPrime (of size k+1)
// There are (k+1)**2 ciphertexts, and k pairs of (evaluation keys Di,bi)
func RelinearizationOnTheFly(evaluationKeys []*MKEvaluationKey, publicKeys []*MKPublicKey, ciphertexts *MKCiphertext, params *ckks.Parameters) {

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

			decomposedIJ := GInverse(cipherParts[i*(k+1)+j], params)

			Dot(decomposedIJ, publicKeys[j].key[0], tmpIJ, ringQ) // line 6

			decomposedTmp := GInverse(tmpIJ, params) // inverse and matric mult (line 7)

			Dot(decomposedTmp, di0, tmpC0Prime, ringQ)
			Dot(decomposedTmp, di1, tmpCIPrime, ringQ)

			ringQ.Add(c0Prime, tmpC0Prime, c0Prime)
			ringQ.Add(res[i], tmpCIPrime, res[i])

			Dot(decomposedIJ, di2, tmpIJ, ringQ) // inverse and dot product (line8)
			ringQ.Add(res[j], tmpIJ, res[j])

		}
	}

	res[0] = c0Prime
	ciphertexts.ciphertexts.SetValue(res)
}
*/

// RelinearizationOnTheFly implements the algorithm 3 in the appendix of the Chen paper
// It does relin directly by linearizing each entry of the extended ciphertext and stores it in cPrime (of size k+1)
// There are (k+1)**2 ciphertexts, and k pairs of (evaluation keys Di,bi)
func RelinearizationOnTheFly(evaluationKeys []*MKEvaluationKey, publicKeys []*MKPublicKey, ciphertexts *MKCiphertext, params *ckks.Parameters) {
	ringQP := GetRingQP(params)
	ringQ := GetRingQ(params)
	ringP := GetRingP(params)

	var baseconverter *ring.FastBasisExtender
	baseconverter = ring.NewFastBasisExtender(ringQ, ringP)
	levelQ := uint64(len(ringQ.Modulus) - 1)

	k := uint64(len(evaluationKeys))
	restmpQ := make([]*ring.Poly, k+1)
	res := make([]*ring.Poly, k+1)
	restmpP := make([]*ring.Poly, k+1)

	for i := uint64(0); i < k+1; i++ {
		restmpQ[i] = ringQ.NewPoly()
		restmpP[i] = ringP.NewPoly()

		res[i] = ringQ.NewPoly()
	}

	cipherParts := ciphertexts.ciphertexts.Value()
	tmpC0Q := ringQP.NewPoly()
	tmpC0P := ringQP.NewPoly()

	tmpCiQ := ringQP.NewPoly()
	tmpCiP := ringQP.NewPoly()

	tmpIJQ := ringQP.NewPoly()
	tmpIJP := ringQP.NewPoly()

	for i := uint64(1); i < k; i++ {

		di0 := evaluationKeys[i-1].key[0]
		di1 := evaluationKeys[i-1].key[1]
		di2 := evaluationKeys[i-1].key[2]

		for j := uint64(1); j < k; j++ {

			cIJDoublePrimeQ := ringQP.NewPoly()
			cIJDoublePrimeP := ringQP.NewPoly()
			decomposedIJQ, decomposedIJP := GInverse(cipherParts[i*(k+1)+j], params)

			Dot(decomposedIJQ, publicKeys[j].key[0], cIJDoublePrimeQ, ringQ) // line 3
			Dot(decomposedIJP, publicKeys[j].key[0], cIJDoublePrimeP, ringP)
			cIJPrime := ringQ.NewPoly()

			baseconverter.ModDownSplitNTTPQ(levelQ, cIJDoublePrimeQ, cIJDoublePrimeP, cIJPrime) // line 4

			decomposedTmpQ, decomposedTmpP := GInverse(cIJPrime, params) // inverse and matrix mult (line 5)

			Dot(decomposedTmpQ, di0, tmpC0Q, ringQ)
			Dot(decomposedTmpP, di0, tmpC0P, ringP)

			Dot(decomposedTmpQ, di1, tmpCiQ, ringQ)
			Dot(decomposedTmpP, di1, tmpCiP, ringP)

			ringQ.Add(restmpQ[0], tmpC0Q, restmpQ[0])
			ringQ.Add(restmpQ[i], tmpCiQ, restmpQ[i])

			ringP.Add(restmpP[0], tmpC0P, restmpP[0])
			ringP.Add(restmpP[i], tmpCiP, restmpP[i])

			Dot(decomposedIJQ, di2, tmpIJQ, ringQ) // line 6 of algorithm
			Dot(decomposedIJP, di2, tmpIJP, ringP)

			ringQ.Add(restmpQ[j], tmpIJQ, restmpQ[j])
			ringP.Add(restmpP[j], tmpIJP, restmpP[j])
		}
	}

	tmpModDown := ringQ.NewPoly()

	baseconverter.ModDownSplitNTTPQ(levelQ, restmpQ[0], restmpP[0], tmpModDown)
	ringQ.Add(cipherParts[0], tmpModDown, res[0])

	for i := uint64(1); i <= k; i++ {

		ringQ.Add(cipherParts[i], cipherParts[(k+1)*i], res[i])
		tmpModDown := ringQ.NewPoly()
		baseconverter.ModDownSplitNTTPQ(levelQ, restmpQ[i], restmpP[i], tmpModDown)
		ringQ.Add(res[i], tmpModDown, res[i])
	}

	ciphertexts.ciphertexts.SetValue(res)
}
