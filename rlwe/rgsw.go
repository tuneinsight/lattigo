package rlwe

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"math"
)

// ExternalProduct computes RLWE x RGSW -> RLWE
// RLWE : (-as + m + e, a)
//  x
// RGSW : [(-as + P*w*m1 + e, a), (-bs + e, b + P*w*m1)]
//  =
// RLWE : (<RLWE, RGSW[0]>, <RLWE, RGSW[1]>)
func (ks *KeySwitcher) ExternalProduct(op0 *Ciphertext, op1 *RGSWCiphertext, op2 *Ciphertext) {

	levelQ, levelP := op1.LevelQ(), op1.LevelP()

	var c0QP, c1QP PolyQP
	if op0 == op2 {
		c0QP, c1QP = ks.Pool[1], ks.Pool[2]
	} else {
		c0QP, c1QP = PolyQP{op2.Value[0], ks.Pool[1].P}, PolyQP{op2.Value[1], ks.Pool[2].P}
	}

	if levelP < 1 {
		ks.externalProductInPlaceSinglePAndBitDecomp(op0, op1, ks.Pool[1], ks.Pool[2])
	} else {
		ks.externalProductInPlaceMultipleP(levelQ, levelP, op0, op1, ks.Pool[1].Q, ks.Pool[1].P, ks.Pool[2].Q, ks.Pool[2].P)
	}

	if levelP != -1 {
		ks.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c0QP.Q, c0QP.P, op2.Value[0])
		ks.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c1QP.Q, c1QP.P, op2.Value[1])
	}
}

func (ks *KeySwitcher) externalProductInPlaceSinglePAndBitDecomp(ct0 *Ciphertext, rgsw *RGSWCiphertext, c0QP, c1QP PolyQP) {

	// rgsw = [(-as + P*w*m1 + e, a), (-bs + e, b + P*w*m1)]
	// ct = [-cs + m0 + e, c]
	// ctOut = [<ct, rgsw[0]>, <ct, rgsw[1]>] = [ct[0] * rgsw[0][0] + ct[1] * rgsw[0][1], ct[0] * rgsw[1][0] + ct[1] * rgsw[1][1]]
	ringQ := ks.RingQ()
	ringP := ks.RingP()

	levelQ := rgsw.LevelQ()
	levelP := rgsw.LevelP()

	lb2 := ks.logbase2
	mask := uint64(((1 << lb2) - 1))
	if mask == 0 {
		mask = 0xFFFFFFFFFFFFFFFF
	}

	decompRNS := ks.DecompRNS(levelQ, levelP)
	decompBIT := ks.DecompBIT(levelQ, levelP)

	QiOverF := ks.QiOverflowMargin(levelQ) >> 1
	PiOverF := ks.PiOverflowMargin(levelP) >> 1

	ringQ.InvNTTLvl(levelQ, ct0.Value[0], ks.PoolInvNTT)
	cw := ks.Pool[0].Q.Coeffs[0]
	cwNTT := ks.PoolBitDecomp
	var reduce int
	for i := 0; i < decompRNS; i++ {
		for j := 0; j < decompBIT; j++ {

			ring.MaskVec(ks.PoolInvNTT.Coeffs[i], cw, j*lb2, mask)

			// (a, b) + (c0 * rgsw[0][0], c0 * rgsw[0][1])
			if i == 0 && j == 0 {
				for u := 0; u < levelQ+1; u++ {
					ringQ.NTTSingleLazy(u, cw, cwNTT)
					ring.MulCoeffsMontgomeryVec(rgsw.Value[i][j][0][0].Q.Coeffs[u], cwNTT, c0QP.Q.Coeffs[u], ringQ.Modulus[u], ringQ.MredParams[u])
					ring.MulCoeffsMontgomeryVec(rgsw.Value[i][j][0][1].Q.Coeffs[u], cwNTT, c1QP.Q.Coeffs[u], ringQ.Modulus[u], ringQ.MredParams[u])
				}

				for u := 0; u < levelP+1; u++ {
					ringP.NTTSingleLazy(u, cw, cwNTT)
					ring.MulCoeffsMontgomeryVec(rgsw.Value[i][j][0][0].P.Coeffs[u], cwNTT, c0QP.P.Coeffs[u], ringP.Modulus[u], ringP.MredParams[u])
					ring.MulCoeffsMontgomeryVec(rgsw.Value[i][j][0][1].P.Coeffs[u], cwNTT, c1QP.P.Coeffs[u], ringP.Modulus[u], ringP.MredParams[u])
				}
			} else {
				for u := 0; u < levelQ+1; u++ {
					ringQ.NTTSingleLazy(u, cw, cwNTT)
					ring.MulCoeffsMontgomeryAndAddVec(rgsw.Value[i][j][0][0].Q.Coeffs[u], cwNTT, c0QP.Q.Coeffs[u], ringQ.Modulus[u], ringQ.MredParams[u])
					ring.MulCoeffsMontgomeryAndAddVec(rgsw.Value[i][j][0][1].Q.Coeffs[u], cwNTT, c1QP.Q.Coeffs[u], ringQ.Modulus[u], ringQ.MredParams[u])
				}

				for u := 0; u < levelP+1; u++ {
					ringP.NTTSingleLazy(u, cw, cwNTT)
					ring.MulCoeffsMontgomeryAndAddVec(rgsw.Value[i][j][0][0].P.Coeffs[u], cwNTT, c0QP.P.Coeffs[u], ringP.Modulus[u], ringP.MredParams[u])
					ring.MulCoeffsMontgomeryAndAddVec(rgsw.Value[i][j][0][1].P.Coeffs[u], cwNTT, c1QP.P.Coeffs[u], ringP.Modulus[u], ringP.MredParams[u])
				}
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
	}

	ringQ.InvNTTLvl(levelQ, ct0.Value[1], ks.PoolInvNTT)
	for i := 0; i < decompRNS; i++ {
		for j := 0; j < decompBIT; j++ {
			ring.MaskVec(ks.PoolInvNTT.Coeffs[i], cw, j*lb2, mask)

			// (a, b) + (c1 * rgsw[1][0], c1 * rgsw[1][1])

			for u := 0; u < levelQ+1; u++ {
				ringQ.NTTSingleLazy(u, cw, cwNTT)
				ring.MulCoeffsMontgomeryAndAddVec(rgsw.Value[i][j][1][0].Q.Coeffs[u], cwNTT, c0QP.Q.Coeffs[u], ringQ.Modulus[u], ringQ.MredParams[u])
				ring.MulCoeffsMontgomeryAndAddVec(rgsw.Value[i][j][1][1].Q.Coeffs[u], cwNTT, c1QP.Q.Coeffs[u], ringQ.Modulus[u], ringQ.MredParams[u])
			}

			for u := 0; u < levelP+1; u++ {
				ringP.NTTSingleLazy(u, cw, cwNTT)
				ring.MulCoeffsMontgomeryAndAddVec(rgsw.Value[i][j][1][0].P.Coeffs[u], cwNTT, c0QP.P.Coeffs[u], ringP.Modulus[u], ringP.MredParams[u])
				ring.MulCoeffsMontgomeryAndAddVec(rgsw.Value[i][j][1][1].P.Coeffs[u], cwNTT, c1QP.P.Coeffs[u], ringP.Modulus[u], ringP.MredParams[u])
			}
		}
	}
}

func (ks *KeySwitcher) externalProductInPlaceMultipleP(levelQ, levelP int, ct0 *Ciphertext, rgsw *RGSWCiphertext, c0OutQ, c0OutP, c1OutQ, c1OutP *ring.Poly) {
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
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, rgsw.Value[i][0][0][0], c2QP, c0QP)
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, rgsw.Value[i][0][0][1], c2QP, c1QP)
		} else {
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, rgsw.Value[i][0][0][0], c2QP, c0QP)
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, rgsw.Value[i][0][0][1], c2QP, c1QP)
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
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, rgsw.Value[i][0][1][0], c2QP, c0QP)
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, rgsw.Value[i][0][1][1], c2QP, c1QP)
		} else {
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, rgsw.Value[i][0][1][0], c2QP, c0QP)
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, rgsw.Value[i][0][1][1], c2QP, c1QP)
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
