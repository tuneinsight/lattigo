package mkrlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// Relinearization implements the algorithm 3 in the appendix of the Chen paper
// It relinearize each entry of the extended ciphertext and stores it in cPrime (of size k+1)
// There are (k+1)**2 ciphertexts, and k pairs of evaluation and public keys as input
func Relinearization(relinKeys []*MKRelinearizationKey, publicKeys []*MKPublicKey, ct *[]*ring.Poly, params *rlwe.Parameters, level uint64) {

	ringQ := GetRingQ(params)
	ringP := GetRingP(params)

	baseconverter := ring.NewFastBasisExtender(ringQ, ringP)

	k := uint64(len(relinKeys))
	restmpQ := make([]*ring.Poly, k+1)
	res := make([]*ring.Poly, k+1)
	restmpP := make([]*ring.Poly, k+1)

	for i := uint64(0); i < k+1; i++ {
		restmpQ[i] = ringQ.NewPoly()
		restmpP[i] = ringP.NewPoly()

		res[i] = ringQ.NewPoly()
	}
	d0QP := rlwe.NewSwitchingKey(*params)
	d1QP := rlwe.NewSwitchingKey(*params)
	d2QP := rlwe.NewSwitchingKey(*params)

	pkQ := NewDecomposedPoly(ringQ, params.Beta())
	pkP := NewDecomposedPoly(ringP, params.Beta())

	for i := uint64(1); i <= k; i++ {

		prepareRelinKeyKeySwitch(i, level, uint64(len(ringQ.Modulus)), params.Beta(), relinKeys, ringQ, ringP, d0QP, d1QP, d2QP)

		for j := uint64(1); j <= k; j++ {

			if (*ct)[i*(k+1)+j] != nil {

				preparePublicKey(j, level, uint64(len(ringQ.Modulus)), params.Beta(), publicKeys, ringQ, ringP, pkQ, pkP)
				decomposedIJQP := GInverseKeySwitch((*ct)[i*(k+1)+j], params, level) // line 3

				cIJtmpQ, cIJtmpP := DotQP(level, pkQ, pkP, decomposedIJQP, ringQ, ringP)

				cIJPrime := ringQ.NewPoly()

				baseconverter.ModDownSplitNTTPQ(level, cIJtmpQ, cIJtmpP, cIJPrime)  // line 4
				decomposedTmpQ, decomposedTmpP := GInverse(cIJPrime, params, level) // inverse and matrix mult (line 5)

				tmpC0Q, tmpC0P := DotQP(level, decomposedTmpQ, decomposedTmpP, d0QP, ringQ, ringP)
				tmpCiQ, tmpCiP := DotQP(level, decomposedTmpQ, decomposedTmpP, d1QP, ringQ, ringP)

				ringQ.AddLvl(level, restmpQ[0], tmpC0Q, restmpQ[0])
				ringQ.AddLvl(level, restmpQ[i], tmpCiQ, restmpQ[i])

				ringP.Add(restmpP[0], tmpC0P, restmpP[0])
				ringP.Add(restmpP[i], tmpCiP, restmpP[i])

				tmpIJQ, tmpIJP := DotSwk(level, d2QP, decomposedIJQP, ringQ, ringP, params.Beta()) // line 6 of algorithm
				ringQ.AddLvl(level, restmpQ[j], tmpIJQ, restmpQ[j])
				ringP.Add(restmpP[j], tmpIJP, restmpP[j])
			}

		}
	}
	tmpModDown := ringQ.NewPoly()

	baseconverter.ModDownSplitNTTPQ(level, restmpQ[0], restmpP[0], tmpModDown)
	ringQ.AddLvl(level, (*ct)[0], tmpModDown, res[0])
	for i := uint64(1); i <= k; i++ {

		if (*ct)[i] != nil && (*ct)[(k+1)*i] != nil {
			ringQ.AddLvl(level, (*ct)[i], (*ct)[(k+1)*i], res[i])
		} else if (*ct)[i] != nil && (*ct)[(k+1)*i] == nil {
			res[i] = (*ct)[i]
		} else {
			res[i] = (*ct)[(k+1)*i]
		}

		baseconverter.ModDownSplitNTTPQ(level, restmpQ[i], restmpP[i], tmpModDown)
		ringQ.AddLvl(level, res[i], tmpModDown, res[i])

	}

	*ct = res
}

// prepare relinearization key for operations in split crt basis
func prepareRelinKey(i, level, modulus, beta uint64, evaluationKeys []*MKRelinearizationKey, ringQ, ringP *ring.Ring, d0Q, d1Q, d2Q, d0P, d1P, d2P *MKDecomposedPoly) {

	di01 := evaluationKeys[i-1].Key01
	di2 := evaluationKeys[i-1].Key2

	for u := uint64(0); u < beta; u++ {

		d0Q.Poly[u].Coeffs = di01.Value[u][0].Coeffs[:level+1]
		d1Q.Poly[u].Coeffs = di01.Value[u][1].Coeffs[:level+1]
		d2Q.Poly[u].Coeffs = di2.Poly[u].Coeffs[:level+1]

		d0P.Poly[u].Coeffs = di01.Value[u][0].Coeffs[modulus:]
		d1P.Poly[u].Coeffs = di01.Value[u][1].Coeffs[modulus:]
		d2P.Poly[u].Coeffs = di2.Poly[u].Coeffs[modulus:]

	}

}

// prepare relinearization key for operations in split crt basis as three switching keys
func prepareRelinKeyKeySwitch(i, level, modulus, beta uint64, evaluationKeys []*MKRelinearizationKey, ringQ, ringP *ring.Ring, d0QP, d1QP, d2QP *rlwe.SwitchingKey) {

	di01 := evaluationKeys[i-1].Key01
	di2 := evaluationKeys[i-1].Key2

	for u := uint64(0); u < beta; u++ {

		d0QP.Value[u][0].Coeffs = di01.Value[u][0].Coeffs[:level+1]
		d0QP.Value[u][1].Coeffs = di01.Value[u][0].Coeffs[modulus:]

		d1QP.Value[u][0].Coeffs = di01.Value[u][1].Coeffs[:level+1]
		d1QP.Value[u][1].Coeffs = di01.Value[u][1].Coeffs[modulus:]

		d2QP.Value[u][0].Coeffs = di2.Poly[u].Coeffs[:level+1]
		d2QP.Value[u][1].Coeffs = di2.Poly[u].Coeffs[modulus:]

	}

}

// prepare public key for operations in split crt basis
func preparePublicKey(j, level, modulus, beta uint64, publicKeys []*MKPublicKey, ringQ, ringP *ring.Ring, pkQ, pkP *MKDecomposedPoly) {

	for u := uint64(0); u < beta; u++ {

		pkQ.Poly[u].Coeffs = publicKeys[j-1].Key[0].Poly[u].Coeffs[:level+1]
		pkP.Poly[u].Coeffs = publicKeys[j-1].Key[0].Poly[u].Coeffs[modulus:]
	}

}

// GInverse is a method that returns the decomposition of a polynomial from R_qp to R_qp^beta
// polynomials are returned in MForm
func GInverse(p *ring.Poly, params *rlwe.Parameters, level uint64) (*MKDecomposedPoly, *MKDecomposedPoly) {

	beta := params.Beta()
	ringQ := GetRingQ(params)
	ringP := GetRingP(params)

	resQ := new(MKDecomposedPoly)
	resP := new(MKDecomposedPoly)

	polynomialsQ := make([]*ring.Poly, beta)
	polynomialsP := make([]*ring.Poly, beta)
	invPoly := ringQ.NewPoly()

	ringQ.InvNTTLvl(level, p, invPoly)

	// generate each poly decomposed in the base
	for i := uint64(0); i < beta; i++ {

		polynomialsQ[i] = ringQ.NewPoly()
		polynomialsP[i] = ringP.NewPoly()
		decomposeAndSplitNTT(level, i, p, invPoly, polynomialsQ[i], polynomialsP[i], params, ringQ, ringP)
		//pass polynomials in MForm
		ringQ.MFormLvl(level, polynomialsQ[i], polynomialsQ[i])
		ringP.MForm(polynomialsP[i], polynomialsP[i])

	}
	resQ.Poly = polynomialsQ
	resP.Poly = polynomialsP

	return resQ, resP
}

// GInverseKeySwitch is a method that returns the decomposition of a polynomial from R_qp to R_qp^beta as a switching key
// polynomials are returned in MForm
func GInverseKeySwitch(p *ring.Poly, params *rlwe.Parameters, level uint64) *rlwe.SwitchingKey {

	beta := params.Beta()
	ringQ := GetRingQ(params)
	ringP := GetRingP(params)

	swk := new(rlwe.SwitchingKey)
	swk.Value = make([][2]*ring.Poly, int(beta))
	res := new(rlwe.SwitchingKey)

	polynomialsQ := make([]*ring.Poly, beta)
	polynomialsP := make([]*ring.Poly, beta)
	invPoly := ringQ.NewPoly()

	ringQ.InvNTTLvl(level, p, invPoly)

	// generate each poly decomposed in the base
	for i := uint64(0); i < beta; i++ {

		swk.Value[i][0] = ringQ.NewPoly()
		swk.Value[i][1] = ringP.NewPoly()

		polynomialsQ[i] = ringQ.NewPoly()
		polynomialsP[i] = ringP.NewPoly()
		decomposeAndSplitNTT(level, i, p, invPoly, swk.Value[i][0], swk.Value[i][1], params, ringQ, ringP)
		//pass polynomials in MForm
		ringQ.MFormLvl(level, swk.Value[i][0], swk.Value[i][0])
		ringP.MForm(swk.Value[i][1], swk.Value[i][1])

	}
	res.Value = swk.Value

	return res
}

// decomposeAndSplitNTT decomposes the input polynomial into the target CRT basis.
func decomposeAndSplitNTT(level, beta uint64, c2NTT, c2InvNTT, c2QiQ, c2QiP *ring.Poly, params *rlwe.Parameters, ringQ, ringP *ring.Ring) {

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
