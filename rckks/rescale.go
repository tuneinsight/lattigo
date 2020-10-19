package rckks

import (
	"github.com/ldsec/lattigo/v2/ring"
	"unsafe"
)

// DivRoundByLastModulusNTTRCKKS divides (rounded) the polynomial by its last modulus. The input must be in the NTTRCKKS domain.
func DivRoundByLastModulusNTTRCKKS(r *ring.Ring, p0 *ring.Poly) {

	var pHalf, pHalfNegQi uint64

	level := len(p0.Coeffs) - 1

	pTmp := make([]uint64, r.N)

	invnttrckks(p0.Coeffs[level], p0.Coeffs[level], r.N, r.NttPsiInv[level], r.NttNInv[level], r.Modulus[level], r.MredParams[level])

	// Center by (p-1)/2
	pHalf = (r.Modulus[level] - 1) >> 1
	p0tmp := p0.Coeffs[level]
	pj := r.Modulus[level]

	for i := uint64(0); i < r.N; i = i + 8 {

		z := (*[8]uint64)(unsafe.Pointer(&p0tmp[i]))

		z[0] = ring.CRed(z[0]+pHalf, pj)
		z[1] = ring.CRed(z[1]+pHalf, pj)
		z[2] = ring.CRed(z[2]+pHalf, pj)
		z[3] = ring.CRed(z[3]+pHalf, pj)
		z[4] = ring.CRed(z[4]+pHalf, pj)
		z[5] = ring.CRed(z[5]+pHalf, pj)
		z[6] = ring.CRed(z[6]+pHalf, pj)
		z[7] = ring.CRed(z[7]+pHalf, pj)
	}

	for i := 0; i < level; i++ {

		p1tmp := p0.Coeffs[i]

		qi := r.Modulus[i]
		bredParams := r.BredParams[i]
		mredParams := r.MredParams[i]
		rescaleParams := r.RescaleParams[level-1][i]

		pHalfNegQi = r.Modulus[i] - ring.BRedAdd(pHalf, qi, bredParams)

		for j := uint64(0); j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p0tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&pTmp[j]))

			z[0] = x[0] + pHalfNegQi
			z[1] = x[1] + pHalfNegQi
			z[2] = x[2] + pHalfNegQi
			z[3] = x[3] + pHalfNegQi
			z[4] = x[4] + pHalfNegQi
			z[5] = x[5] + pHalfNegQi
			z[6] = x[6] + pHalfNegQi
			z[7] = x[7] + pHalfNegQi
		}

		nttrckks(pTmp, pTmp, r.N, r.NttPsi[i], qi, mredParams, bredParams)

		// (x[i] - x[-1]) * InvQ
		for j := uint64(0); j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&pTmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))

			z[0] = ring.MRed(z[0]+(qi-x[0]), rescaleParams, qi, mredParams)
			z[1] = ring.MRed(z[1]+(qi-x[1]), rescaleParams, qi, mredParams)
			z[2] = ring.MRed(z[2]+(qi-x[2]), rescaleParams, qi, mredParams)
			z[3] = ring.MRed(z[3]+(qi-x[3]), rescaleParams, qi, mredParams)
			z[4] = ring.MRed(z[4]+(qi-x[4]), rescaleParams, qi, mredParams)
			z[5] = ring.MRed(z[5]+(qi-x[5]), rescaleParams, qi, mredParams)
			z[6] = ring.MRed(z[6]+(qi-x[6]), rescaleParams, qi, mredParams)
			z[7] = ring.MRed(z[7]+(qi-x[7]), rescaleParams, qi, mredParams)
		}
	}

	p0.Coeffs = p0.Coeffs[:level]
}

// DivRoundByLastModulusManyNTTRCKKS applies DivRoundByLastModulusNTTRCKKS "nbRescales" times.
func DivRoundByLastModulusManyNTTRCKKS(r *ring.Ring, p0 *ring.Poly, nbRescales uint64) {
	level := uint64(len(p0.Coeffs) - 1)
	InvNTTRCKKSLvl(r, level, p0, p0)
	r.DivRoundByLastModulusMany(p0, nbRescales)
	NTTRCKKSLvl(r, level-nbRescales, p0, p0)
}

// ModDownSplitNTTPQRCKKS reduces the basis of a polynomial.
// Given a polynomial with coefficients in basis {Q0,Q1....Qi} and {P0,P1...Pj},
// it reduces its basis from {Q0,Q1....Qi} and {P0,P1...Pj} to {Q0,Q1....Qi}
// and does a rounded integer division of the result by P.
// Inputs must be in the NTTRCKKS domain.
func ModDownSplitNTTPQRCKKS(basisextender *ring.FastBasisExtender, level uint64, p1Q, p1P, p2 *ring.Poly) {

	ringQ := basisextender.RingQ
	ringP := basisextender.RingP
	modDownParams := basisextender.ModDownParamsPQ
	polypool := basisextender.PolyPoolQ

	// First we get the P basis part of p1 out of the NTT domain
	InvNTTRCKKS(ringP, p1P, p1P)

	// Then we target this P basis of p1 and convert it to a Q basis (at the "level" of p1) and copy it on polypool
	// polypool is now the representation of the P basis of p1 but in basis Q (at the "level" of p1)
	ring.ModUpExact(p1P.Coeffs, polypool.Coeffs[:level+1], basisextender.ParamsPQ)

	// Finally, for each level of p1 (and polypool since they now share the same basis) we compute p2 = (P^-1) * (p1 - polypool) mod Q
	for i := uint64(0); i < level+1; i++ {

		qi := ringQ.Modulus[i]
		p1tmp := p1Q.Coeffs[i]
		p2tmp := p2.Coeffs[i]
		p3tmp := polypool.Coeffs[i]
		params := modDownParams[i]
		mredParams := ringQ.MredParams[i]
		bredParams := ringQ.BredParams[i]

		// First we switch back the relevant polypool CRT array back to the NTT domain
		nttrckks(p3tmp, p3tmp, ringQ.N, ringQ.GetNttPsi()[i], ringQ.Modulus[i], mredParams, bredParams)

		// Then for each coefficient we compute (P^-1) * (p1[i][j] - polypool[i][j]) mod qi
		for j := uint64(0); j < ringQ.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = ring.MRed(x[0]+(qi-y[0]), params, qi, mredParams)
			z[1] = ring.MRed(x[1]+(qi-y[1]), params, qi, mredParams)
			z[2] = ring.MRed(x[2]+(qi-y[2]), params, qi, mredParams)
			z[3] = ring.MRed(x[3]+(qi-y[3]), params, qi, mredParams)
			z[4] = ring.MRed(x[4]+(qi-y[4]), params, qi, mredParams)
			z[5] = ring.MRed(x[5]+(qi-y[5]), params, qi, mredParams)
			z[6] = ring.MRed(x[6]+(qi-y[6]), params, qi, mredParams)
			z[7] = ring.MRed(x[7]+(qi-y[7]), params, qi, mredParams)
		}
	}

	// In total we do len(P) + len(Q) NTT, which is optimal (linear in the number of moduli of P and Q)
}
