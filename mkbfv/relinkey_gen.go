package mkbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
)

// Convert creates a switching key K_ij from the evaluation key of the i-th peer and the public key of the j-th peer.
func Convert(D *MKEvaluationKey, publicKey *MKPublicKey, params *bfv.Parameters) *MKSwitchingKey {
	res := new(MKSwitchingKey)
	d0 := D.key[0]
	d1 := D.key[1]
	beta := params.Beta() // size of the decomposition
	k0 := new(MKDecomposedPoly)
	k1 := new(MKDecomposedPoly)
	k2 := new(MKDecomposedPoly)
	ring := GetRingQP(params)

	for l := uint64(0); l < beta; l++ {
		gInv := GInverse(publicKey.key[0].poly[l], params)
		Dot(gInv, d0, k0.poly[l], ring, beta) // g^-1 ( b_j[l]) dot d_i[0]
		Dot(gInv, d1, k1.poly[l], ring, beta) // g^-1 ( b_j[l]) dot d_i[1]

	}
	copy(k2.poly, D.key[2].poly) // TODO : We could also use the ring_poly.Copy method, if this does not do a deep copy
	res.key = [3]*MKDecomposedPoly{k0, k1, k2}

	return res
}

// CreateSharedRelinearizationKey generates a shared relinearization key containing the switching key for all pair of participants.
func CreateSharedRelinearizationKey(params *bfv.Parameters, pubKeys []*MKPublicKey, evalKeys []*MKEvaluationKey) *MKRelinearizationKey {

	res := new(MKRelinearizationKey)

	nbrParticipants := uint64(len(pubKeys))

	tmpArray := make([][]*MKSwitchingKey, nbrParticipants)

	for i := uint64(0); i < nbrParticipants; i++ {

		for j := uint64(0); j < nbrParticipants; j++ {

			tmpArray[i][j] = Convert(evalKeys[i], pubKeys[j], params)
		}
	}

	res.key = tmpArray

	return res
}

// GInverse is a method that returns the decomposition of a polynomial from R_qp to R_qp^beta
func GInverse(p *ring.Poly, params *bfv.Parameters) *MKDecomposedPoly {

	beta := params.Beta()
	ringQ := GetRingQ(params)
	ringP := GetRingP(params)
	ringQP := GetRingQP(params)

	level := uint64(len(ringQ.Modulus)) - 1
	res := new(MKDecomposedPoly)

	polynomials := make([]*ring.Poly, beta)

	convertor := ring.NewFastBasisExtender(ringQ, ringP)

	c2QiQ := params.NewPolyQ()
	c2QiP := params.NewPolyP()
	invPoly := params.NewPolyQ()
	ringQP.InvNTT(p, invPoly)

	// generate each poly decomposed in the base
	for i := uint64(0); i < beta; i++ {

		decomposeAndSplitNTT(level, i, p, invPoly, c2QiQ, c2QiP, params, ringQ, ringP) // TODO: ask if this indeed decompose from R_q to R_q^beta

		currPoly := params.NewPolyQP()

		convertor.ModUpSplitQP(level, c2QiQ, currPoly)

		polynomials[i] = currPoly
	}

	res.poly = polynomials

	return res
}

// decomposeAndSplitNTT decomposes the input polynomial into the target CRT basis.
// this function was copied from bfv evaluator.go in order not to break the encypsulation
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
