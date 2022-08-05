package rlwe

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
)

// GadgetProduct evaluates poly x Gadget -> RLWE where
//
// ctOut = [dot(decomp(cx) * gadget[0]), dot(decomp(cx) * gadget[1])] mod Q
//
// Expects the flag IsNTT of cx to correctly reflect the domain of cx.
func (eval *Evaluator) GadgetProduct(levelQ int, cx *ring.Poly, gadgetCt GadgetCiphertext, ctOut Ciphertext) {

	levelP := gadgetCt.LevelP()

	ctOutQP := CiphertextQP{Value: []ringqp.Poly{{Q: ctOut.Value[0], P: eval.BuffQP[0].P}, {Q: ctOut.Value[1], P: eval.BuffQP[1].P}}}

	if levelP > 0 {
		eval.GadgetProductNoModDown(levelQ, cx, gadgetCt, ctOutQP)
	} else {
		eval.GadgetProductSinglePAndBitDecompNoModDown(levelQ, cx, gadgetCt, ctOutQP)
	}

	if cx.IsNTT && levelP != -1 {
		eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, ctOutQP.Value[0].Q, ctOutQP.Value[0].P, ctOutQP.Value[0].Q)
		eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, ctOutQP.Value[1].Q, ctOutQP.Value[1].P, ctOutQP.Value[1].Q)
	} else if !cx.IsNTT {
		eval.params.RingQ().InvNTTLazyLvl(levelQ, ctOutQP.Value[0].Q, ctOutQP.Value[0].Q)
		eval.params.RingQ().InvNTTLazyLvl(levelQ, ctOutQP.Value[1].Q, ctOutQP.Value[1].Q)

		if levelP != -1 {
			eval.params.RingP().InvNTTLazyLvl(levelP, ctOutQP.Value[0].P, ctOutQP.Value[0].P)
			eval.params.RingP().InvNTTLazyLvl(levelP, ctOutQP.Value[1].P, ctOutQP.Value[1].P)
			eval.BasisExtender.ModDownQPtoQ(levelQ, levelP, ctOutQP.Value[0].Q, ctOutQP.Value[0].P, ctOutQP.Value[0].Q)
			eval.BasisExtender.ModDownQPtoQ(levelQ, levelP, ctOutQP.Value[1].Q, ctOutQP.Value[1].P, ctOutQP.Value[1].Q)
		}
	}
}

// GadgetProductNoModDown applies the gadget prodcut to the polynomial cx :
//
// p0QP = dot(decomp(cx) * gadget[0]) mod QP (encrypted input is multiplied by P factor)
// p1QP = dot(decomp(cx) * gadget[1]) mod QP (encrypted input is multiplied by P factor)
//
// Expects the flag IsNTT of cx to correctly reflect the domain of cx.
func (eval *Evaluator) GadgetProductNoModDown(levelQ int, cx *ring.Poly, gadgetCt GadgetCiphertext, ctOut CiphertextQP) {

	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()
	ringQP := eval.params.RingQP()

	c2QP := eval.BuffDecompQP[0]

	var cxNTT, cxInvNTT *ring.Poly
	if cx.IsNTT {
		cxNTT = cx
		cxInvNTT = eval.BuffInvNTT
		ringQ.InvNTTLvl(levelQ, cxNTT, cxInvNTT)
	} else {
		cxNTT = eval.BuffInvNTT
		cxInvNTT = cx
		ringQ.NTTLvl(levelQ, cxInvNTT, cxNTT)
	}

	levelP := gadgetCt.LevelP()

	decompRNS := eval.params.DecompRNS(levelQ, levelP)

	QiOverF := eval.params.QiOverflowMargin(levelQ) >> 1
	PiOverF := eval.params.PiOverflowMargin(levelP) >> 1

	el := gadgetCt.Value

	// Key switching with CRT decomposition for the Qi
	var reduce int
	for i := 0; i < decompRNS; i++ {

		eval.DecomposeSingleNTT(levelQ, levelP, levelP+1, i, cxNTT, cxInvNTT, c2QP.Q, c2QP.P)

		if i == 0 {
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, el[i][0].Value[0], c2QP, ctOut.Value[0])
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, el[i][0].Value[1], c2QP, ctOut.Value[1])
		} else {
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, el[i][0].Value[0], c2QP, ctOut.Value[0])
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, el[i][0].Value[1], c2QP, ctOut.Value[1])
		}

		if reduce%QiOverF == QiOverF-1 {
			ringQ.ReduceLvl(levelQ, ctOut.Value[0].Q, ctOut.Value[0].Q)
			ringQ.ReduceLvl(levelQ, ctOut.Value[1].Q, ctOut.Value[1].Q)
		}

		if reduce%PiOverF == PiOverF-1 {
			ringP.ReduceLvl(levelP, ctOut.Value[0].P, ctOut.Value[0].P)
			ringP.ReduceLvl(levelP, ctOut.Value[1].P, ctOut.Value[1].P)
		}

		reduce++
	}

	if reduce%QiOverF != 0 {
		ringQ.ReduceLvl(levelQ, ctOut.Value[0].Q, ctOut.Value[0].Q)
		ringQ.ReduceLvl(levelQ, ctOut.Value[1].Q, ctOut.Value[1].Q)
	}

	if reduce%PiOverF != 0 {
		ringP.ReduceLvl(levelP, ctOut.Value[0].P, ctOut.Value[0].P)
		ringP.ReduceLvl(levelP, ctOut.Value[1].P, ctOut.Value[1].P)
	}
}

// GadgetProductHoisted applies the gadget-product to the decomposed polynomial c2 mod QP (BuffQPDecompQP)
// and divides the result by P, reducing the basis from QP to Q.
//
// ctOut = [P^-1 * <BuffQPDecompQP, evakey[0]>, P^-1 * <BuffQPDecompQP, evakey[1]>] mod Q
func (eval *Evaluator) GadgetProductHoisted(levelQ int, BuffQPDecompQP []ringqp.Poly, gadgetCt GadgetCiphertext, ctOut Ciphertext) {

	ctOutQP := CiphertextQP{Value: []ringqp.Poly{{Q: ctOut.Value[0], P: eval.BuffQP[0].P}, {Q: ctOut.Value[1], P: eval.BuffQP[1].P}}}

	eval.GadgetProductHoistedNoModDown(levelQ, BuffQPDecompQP, gadgetCt, ctOutQP)

	levelP := gadgetCt.LevelP()

	// Computes c0Q = c0Q/c0P and c1Q = c1Q/c1P
	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, ctOutQP.Value[0].Q, ctOutQP.Value[0].P, ctOutQP.Value[0].Q)
	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, ctOutQP.Value[1].Q, ctOutQP.Value[1].P, ctOutQP.Value[1].Q)
}

// GadgetProductHoistedNoModDown applies the gadget-product to the decomposed polynomial c2 mod QP (BuffQPDecompQP)
//
// ctOut = [<BuffQPDecompQP, evakey[0]>, <BuffQPDecompQP, evakey[1]>] mod QP
func (eval *Evaluator) GadgetProductHoistedNoModDown(levelQ int, BuffQPDecompQP []ringqp.Poly, gadgetCt GadgetCiphertext, ctOut CiphertextQP) {

	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()
	ringQP := eval.params.RingQP()

	levelP := gadgetCt.LevelP()

	decompRNS := eval.params.DecompRNS(levelQ, levelP)

	QiOverF := eval.params.QiOverflowMargin(levelQ) >> 1
	PiOverF := eval.params.PiOverflowMargin(levelP) >> 1

	// Key switching with CRT decomposition for the Qi
	var reduce int
	for i := 0; i < decompRNS; i++ {

		if i == 0 {
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, gadgetCt.Value[i][0].Value[0], BuffQPDecompQP[i], ctOut.Value[0])
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, gadgetCt.Value[i][0].Value[1], BuffQPDecompQP[i], ctOut.Value[1])
		} else {
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, gadgetCt.Value[i][0].Value[0], BuffQPDecompQP[i], ctOut.Value[0])
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, gadgetCt.Value[i][0].Value[1], BuffQPDecompQP[i], ctOut.Value[1])
		}

		if reduce%QiOverF == QiOverF-1 {
			ringQ.ReduceLvl(levelQ, ctOut.Value[0].Q, ctOut.Value[0].Q)
			ringQ.ReduceLvl(levelQ, ctOut.Value[1].Q, ctOut.Value[1].Q)
		}

		if reduce%PiOverF == PiOverF-1 {
			ringP.ReduceLvl(levelP, ctOut.Value[0].P, ctOut.Value[0].P)
			ringP.ReduceLvl(levelP, ctOut.Value[1].P, ctOut.Value[1].P)
		}

		reduce++
	}

	if reduce%QiOverF != 0 {
		ringQ.ReduceLvl(levelQ, ctOut.Value[0].Q, ctOut.Value[0].Q)
		ringQ.ReduceLvl(levelQ, ctOut.Value[1].Q, ctOut.Value[1].Q)
	}

	if reduce%PiOverF != 0 {
		ringP.ReduceLvl(levelP, ctOut.Value[0].P, ctOut.Value[0].P)
		ringP.ReduceLvl(levelP, ctOut.Value[1].P, ctOut.Value[1].P)
	}
}

// GadgetProductSinglePAndBitDecompNoModDown applies the key-switch to the polynomial cx :
//
// p0QP = dot(decomp(cx) * evakey[0]) mod QP (encrypted input is multiplied by P factor)
// p1QP = dot(decomp(cx) * evakey[1]) mod QP (encrypted input is multiplied by P factor)
//
// Expects the flag IsNTT of cx to correctly reflect the domain of cx.
func (eval *Evaluator) GadgetProductSinglePAndBitDecompNoModDown(levelQ int, cx *ring.Poly, gadgetCt GadgetCiphertext, ctOut CiphertextQP) {

	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()

	var cxInvNTT *ring.Poly
	if cx.IsNTT {
		cxInvNTT = eval.BuffInvNTT
		ringQ.InvNTTLvl(levelQ, cx, cxInvNTT)
	} else {
		cxInvNTT = cx
	}

	levelP := gadgetCt.LevelP()

	decompRNS := eval.params.DecompRNS(levelQ, levelP)
	decompPw2 := eval.params.DecompPw2(levelQ, levelP)

	pw2 := eval.params.pow2Base

	mask := uint64(((1 << pw2) - 1))

	if mask == 0 {
		mask = 0xFFFFFFFFFFFFFFFF
	}

	cw := eval.BuffDecompQP[0].Q.Coeffs[0]
	cwNTT := eval.BuffBitDecomp

	QiOverF := eval.params.QiOverflowMargin(levelQ) >> 1
	PiOverF := eval.params.PiOverflowMargin(levelP) >> 1

	el := gadgetCt.Value

	// Key switching with CRT decomposition for the Qi
	var reduce int
	for i := 0; i < decompRNS; i++ {
		for j := 0; j < decompPw2; j++ {

			ring.MaskVec(cxInvNTT.Coeffs[i], cw, j*pw2, mask)

			if i == 0 && j == 0 {
				for u := 0; u < levelQ+1; u++ {
					ringQ.NTTSingleLazy(u, cw, cwNTT)
					ring.MulCoeffsMontgomeryConstantVec(el[i][j].Value[0].Q.Coeffs[u], cwNTT, ctOut.Value[0].Q.Coeffs[u], ringQ.Modulus[u], ringQ.MredParams[u])
					ring.MulCoeffsMontgomeryConstantVec(el[i][j].Value[1].Q.Coeffs[u], cwNTT, ctOut.Value[1].Q.Coeffs[u], ringQ.Modulus[u], ringQ.MredParams[u])
				}

				for u := 0; u < levelP+1; u++ {
					ringP.NTTSingleLazy(u, cw, cwNTT)
					ring.MulCoeffsMontgomeryConstantVec(el[i][j].Value[0].P.Coeffs[u], cwNTT, ctOut.Value[0].P.Coeffs[u], ringP.Modulus[u], ringP.MredParams[u])
					ring.MulCoeffsMontgomeryConstantVec(el[i][j].Value[1].P.Coeffs[u], cwNTT, ctOut.Value[1].P.Coeffs[u], ringP.Modulus[u], ringP.MredParams[u])
				}
			} else {
				for u := 0; u < levelQ+1; u++ {
					ringQ.NTTSingleLazy(u, cw, cwNTT)
					ring.MulCoeffsMontgomeryConstantAndAddNoModVec(el[i][j].Value[0].Q.Coeffs[u], cwNTT, ctOut.Value[0].Q.Coeffs[u], ringQ.Modulus[u], ringQ.MredParams[u])
					ring.MulCoeffsMontgomeryConstantAndAddNoModVec(el[i][j].Value[1].Q.Coeffs[u], cwNTT, ctOut.Value[1].Q.Coeffs[u], ringQ.Modulus[u], ringQ.MredParams[u])
				}

				for u := 0; u < levelP+1; u++ {
					ringP.NTTSingleLazy(u, cw, cwNTT)
					ring.MulCoeffsMontgomeryConstantAndAddNoModVec(el[i][j].Value[0].P.Coeffs[u], cwNTT, ctOut.Value[0].P.Coeffs[u], ringP.Modulus[u], ringP.MredParams[u])
					ring.MulCoeffsMontgomeryConstantAndAddNoModVec(el[i][j].Value[1].P.Coeffs[u], cwNTT, ctOut.Value[1].P.Coeffs[u], ringP.Modulus[u], ringP.MredParams[u])
				}
			}

			if reduce%QiOverF == QiOverF-1 {
				ringQ.ReduceLvl(levelQ, ctOut.Value[0].Q, ctOut.Value[0].Q)
				ringQ.ReduceLvl(levelQ, ctOut.Value[1].Q, ctOut.Value[1].Q)
			}

			if reduce%PiOverF == PiOverF-1 {
				ringP.ReduceLvl(levelP, ctOut.Value[0].P, ctOut.Value[0].P)
				ringP.ReduceLvl(levelP, ctOut.Value[1].P, ctOut.Value[1].P)
			}

			reduce++
		}
	}

	if reduce%QiOverF != 0 {
		ringQ.ReduceLvl(levelQ, ctOut.Value[0].Q, ctOut.Value[0].Q)
		ringQ.ReduceLvl(levelQ, ctOut.Value[1].Q, ctOut.Value[1].Q)
	}

	if reduce%PiOverF != 0 {
		ringP.ReduceLvl(levelP, ctOut.Value[0].P, ctOut.Value[0].P)
		ringP.ReduceLvl(levelP, ctOut.Value[1].P, ctOut.Value[1].P)
	}
}
