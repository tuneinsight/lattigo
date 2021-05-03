package mkbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
)

// Relinearization implements the algorithm 3 in the appendix of the Chen paper
// It does relin directly by linearizing each entry of the extended ciphertext and stores it in cPrime (of size k+1)
// There are (k+1)**2 ciphertexts, and k pairs of (evaluation keys Di,bi)
func Relinearization(evaluationKeys []*MKEvaluationKey, publicKeys []*MKPublicKey, ct *MKCiphertext, params *bfv.Parameters) {

	ringQ := GetRingQ(params)
	ringP := GetRingP(params)

	baseconverter := ring.NewFastBasisExtender(ringQ, ringP)
	level := uint64(len(ringQ.Modulus)) - 1

	k := uint64(len(evaluationKeys))
	restmpQ := make([]*ring.Poly, k+1)
	res := make([]*ring.Poly, k+1)
	restmpP := make([]*ring.Poly, k+1)

	for i := uint64(0); i < k+1; i++ {
		restmpQ[i] = ringQ.NewPoly()
		restmpP[i] = ringP.NewPoly()

		res[i] = ringQ.NewPoly()
	}

	cipherParts := ct.ciphertexts.Value()

	for _, v := range cipherParts {
		ringQ.NTT(v, v)
	}

	for i := uint64(1); i <= k; i++ {

		d0Q, d1Q, d2Q, d0P, d1P, d2P := prepareEvalKey(i, level, params.Beta(), evaluationKeys)

		for j := uint64(1); j <= k; j++ {

			pkQ, pkP := preparePublicKey(j, level, params.Beta(), publicKeys)

			decomposedIJQ, decomposedIJP := GInverse(cipherParts[i*(k+1)+j], params) // line 3

			cIJtmpQ := DotLvl(level, decomposedIJQ, pkQ, ringQ)
			cIJtmpP := Dot(decomposedIJP, pkP, ringP)

			cIJPrime := ringQ.NewPoly()

			baseconverter.ModDownSplitNTTPQ(level, cIJtmpQ, cIJtmpP, cIJPrime) // line 4

			decomposedTmpQ, decomposedTmpP := GInverse(cIJPrime, params) // inverse and matrix mult (line 5)

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

	for _, v := range res {
		ringQ.InvNTT(v, v)
	}

	ct.ciphertexts.SetValue(res)
}

// prepare evaluation key for operations in split crt basis
func prepareEvalKey(i, level, beta uint64, evaluationKeys []*MKEvaluationKey) (d0Q, d1Q, d2Q, d0P, d1P, d2P *MKDecomposedPoly) {

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
		d0P.poly[u].Coeffs = d0P.poly[u].Coeffs[level+1:]
		d1P.poly[u] = di1.poly[u].CopyNew()
		d1P.poly[u].Coeffs = d1P.poly[u].Coeffs[level+1:]
		d2P.poly[u] = di2.poly[u].CopyNew()
		d2P.poly[u].Coeffs = d2P.poly[u].Coeffs[level+1:]
	}

	return
}

// prepare public key for operations in split crt basis
func preparePublicKey(j, level, beta uint64, publicKeys []*MKPublicKey) (pkQ, pkP *MKDecomposedPoly) {

	pkQ = new(MKDecomposedPoly)
	pkQ.poly = make([]*ring.Poly, beta)
	pkP = new(MKDecomposedPoly)
	pkP.poly = make([]*ring.Poly, beta)

	for u := uint64(0); u < beta; u++ {
		pkQ.poly[u] = publicKeys[j-1].key[0].poly[u].CopyNew()
		pkQ.poly[u].Coeffs = pkQ.poly[u].Coeffs[:level+1]

		pkP.poly[u] = publicKeys[j-1].key[0].poly[u].CopyNew()
		pkP.poly[u].Coeffs = pkP.poly[u].Coeffs[level+1:]

	}

	return
}

// GInverse is a method that returns the decomposition of a polynomial from R_qp to R_qp^beta
func GInverse(p *ring.Poly, params *bfv.Parameters) (*MKDecomposedPoly, *MKDecomposedPoly) {

	beta := params.Beta()
	ringQ := GetRingQ(params)
	ringP := GetRingP(params)

	resQ := new(MKDecomposedPoly)
	resP := new(MKDecomposedPoly)

	polynomialsQ := make([]*ring.Poly, beta)
	polynomialsP := make([]*ring.Poly, beta)

	invPoly := ringQ.NewPoly()
	ringQ.InvNTT(p, invPoly)

	level := uint64(len(ringQ.Modulus)) - 1

	// generate each poly decomposed in the base
	for i := uint64(0); i < beta; i++ {

		polynomialsQ[i] = ringQ.NewPoly()
		polynomialsP[i] = ringP.NewPoly()

		decomposeAndSplitNTT(level, i, p, invPoly, polynomialsQ[i], polynomialsP[i], params, ringQ, ringP)

	}

	resQ.poly = polynomialsQ
	resP.poly = polynomialsP

	return resQ, resP
}

// decomposeAndSplitNTT decomposes the input polynomial into the target CRT basis.
// this function was copied from ckks evaluator.go in order not to break the encapsulation
func decomposeAndSplitNTT(level, beta uint64, c2NTT, c2InvNTT, c2QiQ, c2QiP *ring.Poly, params *bfv.Parameters, ringQ, ringP *ring.Ring) {

	decomposer := ring.NewDecomposer(ringQ.Modulus, ringP.Modulus)

	decomposer.DecomposeAndSplit(level, beta, c2InvNTT, c2QiQ, c2QiP)

	p0idxst := beta * params.Alpha()
	p0idxed := p0idxst + decomposer.Xalpha()[beta]

	// c2_qi = cx mod qi mod qi
	for x := uint64(0); x < level+1; x++ {

		qi := ringQ.Modulus[x]
		nttPsi := ringQ.GetNttPsi()[x]
		bredParams := ringQ.GetBredParams()[x]
		mredParams := ringQ.GetMredParams()[x]

		if p0idxst <= x && x < p0idxed {
			p0tmp := c2NTT.Coeffs[x]
			p1tmp := c2QiQ.Coeffs[x]
			for j := uint64(0); j < ringQ.N; j++ {
				p1tmp[j] = p0tmp[j]
			}
		} else {
			ring.NTTLazy(c2QiQ.Coeffs[x], c2QiQ.Coeffs[x], ringQ.N, nttPsi, qi, mredParams, bredParams)
		}
	}
	// c2QiP = c2 mod qi mod pj
	ringP.NTTLazy(c2QiP, c2QiP)
}
