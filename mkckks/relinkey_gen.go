package mkckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
)

/*
// Convert creates a switching key K_ij from the evaluation key of the i-th peer and the public key of the j-th peer.
func Convert(D *MKEvaluationKey, publicKey *MKPublicKey, params *ckks.Parameters) *MKSwitchingKey {

	ringQP := GetRingQP(params)
	res := NewMKSwitchingKey(ringQP, params)

	d0 := D.key[0]
	d1 := D.key[1]
	beta := params.Beta() // size of the decomposition

	k0 := NewDecomposedPoly(ringQP, beta)
	k1 := NewDecomposedPoly(ringQP, beta)
	k2 := NewDecomposedPoly(ringQP, beta)

	for l := uint64(0); l < beta; l++ {

		gInv := GInverse(publicKey.key[0].poly[l], params)

		Dot(gInv, d0, k0.poly[l], ringQP) // g^-1 ( b_j[l]) dot d_i[0]
		Dot(gInv, d1, k1.poly[l], ringQP) // g^-1 ( b_j[l]) dot d_i[1]

	}

	for i := uint64(0); i < uint64(len(k2.poly)); i++ {
		ringQP.Copy(D.key[2].poly[i], k2.poly[i])
	}

	res.key[0] = k0
	res.key[1] = k1
	res.key[2] = k2

	return res
}

// GenSharedRelinearizationKey generates a shared relinearization key containing the switching key for all pair of participants.
func GenSharedRelinearizationKey(params *ckks.Parameters, pubKeys []*MKPublicKey, evalKeys []*MKEvaluationKey) *MKRelinearizationKey {

	if len(evalKeys) != len(pubKeys) {
		panic("Number of evaluation keys should be the same as the number of public keys")
	}

	res := new(MKRelinearizationKey)
	nbrParticipants := uint64(len(pubKeys))

	res.key = make([][]*MKSwitchingKey, nbrParticipants)

	for i := uint64(0); i < nbrParticipants; i++ {

		res.key[i] = make([]*MKSwitchingKey, nbrParticipants)

		for j := uint64(0); j < nbrParticipants; j++ {

			res.key[i][j] = Convert(evalKeys[i], pubKeys[j], params)
		}
	}

	return res
}*/

// GInverse is a method that returns the decomposition of a polynomial from R_qp to R_qp^beta
func GInverse(p *ring.Poly, params *ckks.Parameters) (*MKDecomposedPoly, *MKDecomposedPoly) {

	beta := params.Beta()
	ringQ := GetRingQ(params)
	ringP := GetRingP(params)

	level := uint64(len(ringQ.Modulus)) - 1
	resQ := new(MKDecomposedPoly)
	resP := new(MKDecomposedPoly)

	polynomialsQ := make([]*ring.Poly, beta)
	polynomialsP := make([]*ring.Poly, beta)

	invPoly := params.NewPolyQ()
	ringQ.InvNTT(p, invPoly)

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
func decomposeAndSplitNTT(level, beta uint64, c2NTT, c2InvNTT, c2QiQ, c2QiP *ring.Poly, params *ckks.Parameters, ringQ, ringP *ring.Ring) {

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
