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

	ringQ := GetRingQ(params)
	ringP := GetRingP(params)
	var baseconverter *ring.FastBasisExtender
	baseconverter = ring.NewFastBasisExtender(ringQ, ringP)
	level := ciphertexts.ciphertexts.Level()

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

	for i := uint64(1); i <= k; i++ {

		di0 := evaluationKeys[i-1].key[0]
		di1 := evaluationKeys[i-1].key[1]
		di2 := evaluationKeys[i-1].key[2]

		for j := uint64(1); j <= k; j++ {

			decomposedIJQ, decomposedIJP := GInverse(cipherParts[i*(k+1)+j], params, level) // line 3

			cIJtmpQ := DotLvl(level, decomposedIJQ, publicKeys[j-1].key[0], ringQ)
			cIJtmpP := Dot(decomposedIJP, publicKeys[j-1].key[0], ringP)

			cIJPrime := ringQ.NewPoly()

			baseconverter.ModDownSplitNTTPQ(level, cIJtmpQ, cIJtmpP, cIJPrime) // line 4

			decomposedTmpQ, decomposedTmpP := GInverse(cIJPrime, params, level) // inverse and matrix mult (line 5)

			tmpC0Q := DotLvl(level, decomposedTmpQ, di0, ringQ)
			tmpC0P := Dot(decomposedTmpP, di0, ringP)

			tmpCiQ := DotLvl(level, decomposedTmpQ, di1, ringQ)
			tmpCiP := Dot(decomposedTmpP, di1, ringP)

			ringQ.AddLvl(level, restmpQ[0], tmpC0Q, restmpQ[0])
			ringQ.AddLvl(level, restmpQ[i], tmpCiQ, restmpQ[i])
			ringQ.ReduceLvl(level, restmpQ[0], restmpQ[0])
			ringQ.ReduceLvl(level, restmpQ[i], restmpQ[i])

			ringP.Add(restmpP[0], tmpC0P, restmpP[0])
			ringP.Add(restmpP[i], tmpCiP, restmpP[i])
			ringP.Reduce(restmpP[0], restmpP[0])
			ringP.Reduce(restmpP[i], restmpP[i])

			tmpIJQ := DotLvl(level, decomposedIJQ, di2, ringQ) // line 6 of algorithm
			tmpIJP := Dot(decomposedIJP, di2, ringP)

			ringQ.AddLvl(level, restmpQ[j], tmpIJQ, restmpQ[j])
			ringP.Add(restmpP[j], tmpIJP, restmpP[j])
			ringP.Reduce(restmpP[j], restmpP[j])
			ringQ.ReduceLvl(level, restmpQ[j], restmpQ[j])

		}
	}

	tmpModDown := ringQ.NewPoly()

	baseconverter.ModDownSplitNTTPQ(level, restmpQ[0], restmpP[0], tmpModDown)
	ringQ.AddLvl(level, cipherParts[0], tmpModDown, res[0])
	ringQ.ReduceLvl(level, res[0], res[0])

	for i := uint64(1); i <= k; i++ {

		ringQ.AddLvl(level, cipherParts[i], cipherParts[(k+1)*i], res[i])

		baseconverter.ModDownSplitNTTPQ(level, restmpQ[i], restmpP[i], tmpModDown)
		ringQ.AddLvl(level, res[i], tmpModDown, res[i])
		ringQ.ReduceLvl(level, res[i], res[i])

	}

	ciphertexts.ciphertexts.SetValue(res)
}
