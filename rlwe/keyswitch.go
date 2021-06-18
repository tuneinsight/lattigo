package rlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
)

// KeySwitcher is a struct for RLWE key-switching
type KeySwitcher struct {
	decomposer *ring.Decomposer
}

// DecomposeAndSplitNTT decomposes the input polynomial into the target CRT basis.
func DecomposeAndSplitNTT(level, beta int, ringQ, ringP *ring.Ring, decomposer *ring.Decomposer, c2NTT, c2InvNTT, c2QiQ, c2QiP *ring.Poly) {

	decomposer.DecomposeAndSplit(level, beta, c2InvNTT, c2QiQ, c2QiP)

	p0idxst := beta * len(ringP.Modulus)
	p0idxed := p0idxst + decomposer.Xalpha()[beta]

	// c2_qi = cx mod qi mod qi
	for x := 0; x < level+1; x++ {

		qi := ringQ.Modulus[x]
		nttPsi := ringQ.NttPsi[x]
		bredParams := ringQ.BredParams[x]
		mredParams := ringQ.MredParams[x]

		if p0idxst <= x && x < p0idxed {
			p0tmp := c2NTT.Coeffs[x]
			p1tmp := c2QiQ.Coeffs[x]
			for j := 0; j < ringQ.N; j++ {
				p1tmp[j] = p0tmp[j]
			}
		} else {
			ring.NTTLazy(c2QiQ.Coeffs[x], c2QiQ.Coeffs[x], ringQ.N, nttPsi, qi, mredParams, bredParams)
		}
	}
	// c2QiP = c2 mod qi mod pj
	ringP.NTTLazy(c2QiP, c2QiP)
}
