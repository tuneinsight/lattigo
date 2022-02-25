package rlwe

import (
	"github.com/tuneinsight/lattigo/v3/ring"

	"math"
)

// MulRGSWSingleModulus multiplies ct0 (RLWE) with rgsw (RGSW) and writes the result on ctOut (RLWE).
// Assumes that the level of the modulus Q of all the inputs and outputs is zero.
func (ks *KeySwitcher) MulRGSWSingleModulus(ct0 *Ciphertext, rgsw *RGSWCiphertext, ctOut *Ciphertext) {

	// rgsw = [(-as + P*w*m1 + e, a), (-bs + e, b + P*w*m1)]
	// ct = [-cs + m0 + e, c]
	// ctOut = [<ct, rgsw[0]>, <ct, rgsw[1]>] = [ct[0] * rgsw[0][0] + ct[1] * rgsw[0][1], ct[0] * rgsw[1][0] + ct[1] * rgsw[1][1]]
	ringQ := ks.RingQ()
	ringP := ks.RingP()
	ringQP := ks.RingQP()

	c0QP, c1QP := ks.Pool[1], ks.Pool[2]

	c2QP := PolyQP{ct0.Value[0], ks.Pool[0].P}
	ringQ.InvNTTLazyLvl(0, ct0.Value[0], c2QP.P)
	ringP.NTTLazyLvl(0, c2QP.P, c2QP.P)

	// (a, b) + (c0 * rgsw[0][0], c0 * rgsw[0][1])
	ringQP.MulCoeffsMontgomeryConstantLvl(0, 0, rgsw.Value[0][0][0], c2QP, c0QP)
	ringQP.MulCoeffsMontgomeryConstantLvl(0, 0, rgsw.Value[0][0][1], c2QP, c1QP)

	c2QP.Q = ct0.Value[1]
	ringQ.InvNTTLazyLvl(0, ct0.Value[1], c2QP.P)
	ringP.NTTLazyLvl(0, c2QP.P, c2QP.P)

	// (a, b) + (c1 * rgsw[1][0], c1 * rgsw[1][1])
	ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(0, 0, rgsw.Value[0][1][0], c2QP, c0QP)
	ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(0, 0, rgsw.Value[0][1][1], c2QP, c1QP)

	ks.BasisExtender.ModDownQPtoQNTT(0, 0, c0QP.Q, c0QP.P, ctOut.Value[0])
	ks.BasisExtender.ModDownQPtoQNTT(0, 0, c1QP.Q, c1QP.P, ctOut.Value[1])
}

// MulRGSW multiplies ct0 (RLWE) with rgsw (RGSW) and writes the result on ctOut (RLWE).
func (ks *KeySwitcher) MulRGSW(ct0 *Ciphertext, rgsw *RGSWCiphertext, ctOut *Ciphertext) {

	levelQ, levelP := ct0.Level(), len(rgsw.Value[0][0][0].P.Coeffs)-1

	ks.MulRGSWNoModDown(levelQ, levelP, ct0, rgsw, ks.Pool[1].Q, ks.Pool[1].P, ks.Pool[2].Q, ks.Pool[2].P)

	ks.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, ks.Pool[1].Q, ks.Pool[1].P, ctOut.Value[0])
	ks.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, ks.Pool[2].Q, ks.Pool[2].P, ctOut.Value[1])
}

// MulRGSWNoModDown multiplies ct0 (RLWE) with rgsw (RGSW) without the division by P and returns the result in modulus Q and modulus P
// in separate poly.
func (ks *KeySwitcher) MulRGSWNoModDown(levelQ, levelP int, ct0 *Ciphertext, rgsw *RGSWCiphertext, c0OutQ, c0OutP, c1OutQ, c1OutP *ring.Poly) {
	// rgsw = [(-as + P*w*m1 + e, a), (-bs + e, b + P*w*m1)]
	// ct = [-cs + m0 + e, c]
	// ctOut = [<ct, rgsw[0]>, <ct, rgsw[1]>] = [ct[0] * rgsw[0][0] + ct[1] * rgsw[0][1], ct[0] * rgsw[1][0] + ct[1] * rgsw[1][1]]

	var reduce int

	ringQ := ks.RingQ()
	ringP := ks.RingP()
	ringQP := ks.RingQP()

	c2QP := ks.Pool[0]

	c0QP := PolyQP{c0OutQ, c0OutP}
	c1QP := PolyQP{c1OutQ, c1OutP}

	alpha := levelP + 1
	beta := int(math.Ceil(float64(levelQ+1) / float64(levelP+1)))

	QiOverF := ks.Parameters.QiOverflowMargin(levelQ) >> 1
	PiOverF := ks.Parameters.PiOverflowMargin(levelP) >> 1

	var c2NTT, c2InvNTT *ring.Poly
	if ct0.Value[0].IsNTT {
		c2NTT = ct0.Value[0]
		c2InvNTT = ks.PoolInvNTT
		ringQ.InvNTTLvl(levelQ, c2NTT, c2InvNTT)
	} else {
		c2NTT = ks.PoolInvNTT
		c2InvNTT = ct0.Value[0]
		ringQ.NTTLvl(levelQ, c2InvNTT, c2NTT)
	}

	// (a, b) + (c0 * rgsw[0][0], c0 * rgsw[0][1])
	for i := 0; i < beta; i++ {

		ks.DecomposeSingleNTT(levelQ, levelP, alpha, i, c2NTT, c2InvNTT, c2QP.Q, c2QP.P)

		if i == 0 {
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, rgsw.Value[i][0][0], c2QP, c0QP)
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, rgsw.Value[i][0][1], c2QP, c1QP)
		} else {
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, rgsw.Value[i][0][0], c2QP, c0QP)
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, rgsw.Value[i][0][1], c2QP, c1QP)
		}

		if reduce%QiOverF == QiOverF-1 {
			ringQ.ReduceLvl(levelQ, c0QP.Q, c0QP.Q)
			ringQ.ReduceLvl(levelQ, c1QP.Q, c1QP.Q)
		}

		if reduce%PiOverF == PiOverF-1 {
			ringP.ReduceLvl(levelP, c0QP.P, c0QP.P)
			ringP.ReduceLvl(levelP, c1QP.P, c1QP.P)
		}

		reduce++
	}

	if ct0.Value[0].IsNTT {
		c2NTT = ct0.Value[1]
		c2InvNTT = ks.PoolInvNTT
		ringQ.InvNTTLvl(levelQ, c2NTT, c2InvNTT)
	} else {
		c2NTT = ks.PoolInvNTT
		c2InvNTT = ct0.Value[1]
		ringQ.NTTLvl(levelQ, c2InvNTT, c2NTT)
	}

	// (a, b) + (c1 * rgsw[1][0], c1 * rgsw[1][1])
	for i := 0; i < beta; i++ {

		ks.DecomposeSingleNTT(levelQ, levelP, alpha, i, c2NTT, c2InvNTT, c2QP.Q, c2QP.P)

		if i == 0 {
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, rgsw.Value[i][1][0], c2QP, c0QP)
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, rgsw.Value[i][1][1], c2QP, c1QP)
		} else {
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, rgsw.Value[i][1][0], c2QP, c0QP)
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, rgsw.Value[i][1][1], c2QP, c1QP)
		}

		if reduce%QiOverF == QiOverF-1 {
			ringQ.ReduceLvl(levelQ, c0QP.Q, c0QP.Q)
			ringQ.ReduceLvl(levelQ, c1QP.Q, c1QP.Q)
		}

		if reduce%PiOverF == PiOverF-1 {
			ringP.ReduceLvl(levelP, c0QP.P, c0QP.P)
			ringP.ReduceLvl(levelP, c1QP.P, c1QP.P)
		}

		reduce++
	}

	if reduce%QiOverF != 0 {
		ringQ.ReduceLvl(levelQ, c0QP.Q, c0QP.Q)
		ringQ.ReduceLvl(levelQ, c1QP.Q, c1QP.Q)
	}

	if reduce%PiOverF != 0 {
		ringP.ReduceLvl(levelP, c0QP.P, c0QP.P)
		ringP.ReduceLvl(levelP, c1QP.P, c1QP.P)
	}
}
