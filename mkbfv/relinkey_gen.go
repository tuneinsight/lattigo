package mkbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
)

// GInverse is a method that returns the decomposition of a polynomial from R_qp to R_qp^beta
func GInverse(p *ring.Poly, params *bfv.Parameters) *MKDecomposedPoly {

	beta := params.Beta()
	ringQ := GetRingQ(params)
	ringP := GetRingP(params)
	ringQP := GetRingQP(params)

	level := uint64(len(ringQ.Modulus)) - 1
	res := new(MKDecomposedPoly)

	polynomials := make([]*ring.Poly, beta)

	c2QiQ := params.NewPolyQ()
	c2QiP := params.NewPolyP()
	invPoly := params.NewPolyQ()
	ringQ.InvNTT(p, invPoly)

	// generate each poly decomposed in the base
	for i := uint64(0); i < beta; i++ {

		decomposeAndSplitNTT(level, i, p, invPoly, c2QiQ, c2QiP, params, ringQ, ringP)

		polynomials[i] = toRingQP(c2QiQ, c2QiP, ringQP)
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

// reassemble two elements of Rq and Rp in Rqp
func toRingQP(p1, p2 *ring.Poly, ringQp *ring.Ring) *ring.Poly {

	res := ringQp.NewPoly()

	// copy coefficients
	res.SetCoefficients(append(p1.Coeffs, p2.Coeffs...))

	return res
}
