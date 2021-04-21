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

	baseconverter := ring.NewFastBasisExtender(ringQ, ringP)
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

		d0Q, d1Q, d2Q, d0P, d1P, d2P := prepareEvalKey(i, level, uint64(len(ringQ.Modulus)), params.Beta(), evaluationKeys)

		for j := uint64(1); j <= k; j++ {

			pkQ, pkP := preparePublicKey(j, level, uint64(len(ringQ.Modulus)), params.Beta(), publicKeys)

			decomposedIJQ, decomposedIJP := GInverse(cipherParts[i*(k+1)+j], params, level) // line 3

			cIJtmpQ := DotLvl(level, decomposedIJQ, pkQ, ringQ)
			cIJtmpP := Dot(decomposedIJP, pkP, ringP)

			//uncomment to test initialization
			/*restmpQ[i] = cIJtmpQ
			restmpP[i] = cIJtmpP
			restmpQ[j] = cIJtmpQ
			restmpP[j] = cIJtmpP*/

			cIJPrime := ringQ.NewPoly()

			baseconverter.ModDownSplitNTTPQ(level, cIJtmpQ, cIJtmpP, cIJPrime) // line 4

			//uncomment to test initialization
			/*res[i] = cIJPrime
			res[j] = cIJPrime*/

			decomposedTmpQ, decomposedTmpP := GInverse(cIJPrime, params, level) // inverse and matrix mult (line 5)

			tmpC0Q := DotLvl(level, decomposedTmpQ, d0Q, ringQ)
			tmpC0P := Dot(decomposedTmpP, d0P, ringP)

			tmpCiQ := DotLvl(level, decomposedTmpQ, d1Q, ringQ)
			tmpCiP := Dot(decomposedTmpP, d1P, ringP)

			ringQ.AddLvl(level, restmpQ[0], tmpC0Q, restmpQ[0])
			ringQ.AddLvl(level, restmpQ[i], tmpCiQ, restmpQ[i])

			ringP.Add(restmpP[0], tmpC0P, restmpP[0])
			ringP.Add(restmpP[i], tmpCiP, restmpP[i])

			tmpIJQ := DotLvl(level, decomposedIJQ, d2Q, ringQ) // line 6 of algorithm
			tmpIJP := Dot(decomposedIJP, d2P, ringP)

			ringQ.AddLvl(level, restmpQ[j], tmpIJQ, restmpQ[j])
			ringP.Add(restmpP[j], tmpIJP, restmpP[j])

		}
	}

	tmpModDown := ringQ.NewPoly()

	baseconverter.ModDownSplitNTTPQ(level, restmpQ[0], restmpP[0], tmpModDown)
	ringQ.AddLvl(level, cipherParts[0], tmpModDown, res[0])

	for i := uint64(1); i <= k; i++ {

		ringQ.AddLvl(level, cipherParts[i], cipherParts[(k+1)*i], res[i])

		baseconverter.ModDownSplitNTTPQ(level, restmpQ[i], restmpP[i], tmpModDown)
		ringQ.AddLvl(level, res[i], tmpModDown, res[i])

	}

	ciphertexts.ciphertexts.SetValue(res)
}

// prepare evaluation key for operations in split crt basis
func prepareEvalKey(i, level, modulus, beta uint64, evaluationKeys []*MKEvaluationKey) (d0Q, d1Q, d2Q, d0P, d1P, d2P *MKDecomposedPoly) {

	di0 := evaluationKeys[i-1].key[0]
	di1 := evaluationKeys[i-1].key[1]
	di2 := evaluationKeys[i-1].key[2]

	d0Q = new(MKDecomposedPoly)
	d0Q.poly = make([]*ring.Poly, beta)
	d1Q = new(MKDecomposedPoly)
	d1Q.poly = make([]*ring.Poly, beta)
	d2Q = new(MKDecomposedPoly)
	d2Q.poly = make([]*ring.Poly, beta)
	d0P = new(MKDecomposedPoly)
	d0P.poly = make([]*ring.Poly, beta)
	d1P = new(MKDecomposedPoly)
	d1P.poly = make([]*ring.Poly, beta)
	d2P = new(MKDecomposedPoly)
	d2P.poly = make([]*ring.Poly, beta)

	for u := uint64(0); u < beta; u++ {
		d0Q.poly[u] = di0.poly[u].CopyNew()
		d0Q.poly[u].Coeffs = d0Q.poly[u].Coeffs[:level+1]
		d1Q.poly[u] = di1.poly[u].CopyNew()
		d1Q.poly[u].Coeffs = d1Q.poly[u].Coeffs[:level+1]
		d2Q.poly[u] = di2.poly[u].CopyNew()
		d2Q.poly[u].Coeffs = d2Q.poly[u].Coeffs[:level+1]

		d0P.poly[u] = di0.poly[u].CopyNew()
		d0P.poly[u].Coeffs = d0P.poly[u].Coeffs[modulus:]
		d1P.poly[u] = di1.poly[u].CopyNew()
		d1P.poly[u].Coeffs = d1P.poly[u].Coeffs[modulus:]
		d2P.poly[u] = di2.poly[u].CopyNew()
		d2P.poly[u].Coeffs = d2P.poly[u].Coeffs[modulus:]
	}

	return
}

// prepare public key for operations in split crt basis
func preparePublicKey(j, level, modulus, beta uint64, publicKeys []*MKPublicKey) (pkQ, pkP *MKDecomposedPoly) {

	pkQ = new(MKDecomposedPoly)
	pkQ.poly = make([]*ring.Poly, beta)
	pkP = new(MKDecomposedPoly)
	pkP.poly = make([]*ring.Poly, beta)

	for u := uint64(0); u < beta; u++ {
		pkQ.poly[u] = publicKeys[j-1].key[0].poly[u].CopyNew()
		pkQ.poly[u].Coeffs = pkQ.poly[u].Coeffs[:level+1]

		pkP.poly[u] = publicKeys[j-1].key[0].poly[u].CopyNew()
		pkP.poly[u].Coeffs = pkP.poly[u].Coeffs[modulus:]

	}

	return
}
