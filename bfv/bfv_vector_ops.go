package bfv

import (
	"math/big"
	"unsafe"

	"github.com/tuneinsight/lattigo/v3/ring"
)

// ScaleUpVec takes a Poly pIn in ringT, scales its coefficients up by (Q/T) mod Q, and writes the result in a
// Poly pOut in ringQ.
func ScaleUpVec(ringQ, ringT *ring.Ring, rescaleParams, tmp []uint64, pIn, pOut *ring.Poly) {

	qModTmontgomery := ring.MForm(new(big.Int).Mod(ringQ.ModulusBigint, ringT.ModulusBigint).Uint64(), ringT.Modulus[0], ringT.BredParams[0])

	t := ringT.Modulus[0]
	tHalf := t >> 1
	tInv := ringT.MredParams[0]

	// (x * Q + T/2) mod T
	for i := 0; i < ringQ.N; i = i + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&pIn.Coeffs[0][i]))
		z := (*[8]uint64)(unsafe.Pointer(&tmp[i]))

		z[0] = ring.CRed(ring.MRed(x[0], qModTmontgomery, t, tInv)+tHalf, t)
		z[1] = ring.CRed(ring.MRed(x[1], qModTmontgomery, t, tInv)+tHalf, t)
		z[2] = ring.CRed(ring.MRed(x[2], qModTmontgomery, t, tInv)+tHalf, t)
		z[3] = ring.CRed(ring.MRed(x[3], qModTmontgomery, t, tInv)+tHalf, t)
		z[4] = ring.CRed(ring.MRed(x[4], qModTmontgomery, t, tInv)+tHalf, t)
		z[5] = ring.CRed(ring.MRed(x[5], qModTmontgomery, t, tInv)+tHalf, t)
		z[6] = ring.CRed(ring.MRed(x[6], qModTmontgomery, t, tInv)+tHalf, t)
		z[7] = ring.CRed(ring.MRed(x[7], qModTmontgomery, t, tInv)+tHalf, t)
	}

	// (x * T^-1 - T/2) mod Qi
	for i := 0; i < len(pOut.Coeffs); i++ {
		p0tmp := tmp
		p1tmp := pOut.Coeffs[i]
		qi := ringQ.Modulus[i]
		bredParams := ringQ.BredParams[i]
		mredParams := ringQ.MredParams[i]
		rescaleParams := qi - rescaleParams[i]

		tHalfNegQi := qi - ring.BRedAdd(tHalf, qi, bredParams)

		for j := 0; j < ringQ.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p0tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))

			z[0] = ring.MRed(x[0]+tHalfNegQi, rescaleParams, qi, mredParams)
			z[1] = ring.MRed(x[1]+tHalfNegQi, rescaleParams, qi, mredParams)
			z[2] = ring.MRed(x[2]+tHalfNegQi, rescaleParams, qi, mredParams)
			z[3] = ring.MRed(x[3]+tHalfNegQi, rescaleParams, qi, mredParams)
			z[4] = ring.MRed(x[4]+tHalfNegQi, rescaleParams, qi, mredParams)
			z[5] = ring.MRed(x[5]+tHalfNegQi, rescaleParams, qi, mredParams)
			z[6] = ring.MRed(x[6]+tHalfNegQi, rescaleParams, qi, mredParams)
			z[7] = ring.MRed(x[7]+tHalfNegQi, rescaleParams, qi, mredParams)
		}
	}
}
