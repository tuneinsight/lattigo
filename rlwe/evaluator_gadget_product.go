package rlwe

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
)

// GadgetProduct evaluates poly x Gadget -> RLWE where
//
// p0 = dot(decomp(cx) * gadget[0]) mod Q
// p1 = dot(decomp(cx) * gadget[1]) mod Q
//
// Expects the flag IsNTT of cx to correctly reflect the domain of cx.
func (eval *Evaluator) GadgetProduct(levelQ int, cx *ring.Poly, gadgetCt GadgetCiphertext, p0, p1 *ring.Poly) {

	levelP := gadgetCt.LevelP()

	p0QP := ringqp.Poly{Q: p0, P: eval.BuffQP[1].P}
	p1QP := ringqp.Poly{Q: p1, P: eval.BuffQP[2].P}

	if levelP > 0 {
		eval.GadgetProductNoModDown(levelQ, cx, gadgetCt, p0QP, p1QP)
	} else {
		eval.GadgetProductSinglePAndBitDecompNoModDown(levelQ, cx, gadgetCt, p0QP, p1QP)
	}

	if cx.IsNTT && levelP != -1 {
		eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, p0QP.Q, p0QP.P, p0QP.Q)
		eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, p1QP.Q, p1QP.P, p1QP.Q)
	} else if !cx.IsNTT {
		eval.params.RingQ().InvNTTLazyLvl(levelQ, p0QP.Q, p0QP.Q)
		eval.params.RingQ().InvNTTLazyLvl(levelQ, p1QP.Q, p1QP.Q)

		if levelP != -1 {
			eval.params.RingP().InvNTTLazyLvl(levelP, p0QP.P, p0QP.P)
			eval.params.RingP().InvNTTLazyLvl(levelP, p1QP.P, p1QP.P)
			eval.BasisExtender.ModDownQPtoQ(levelQ, levelP, p0QP.Q, p0QP.P, p0QP.Q)
			eval.BasisExtender.ModDownQPtoQ(levelQ, levelP, p1QP.Q, p1QP.P, p1QP.Q)
		}
	}
}

// GadgetProductNoModDown applies the gadget prodcut to the polynomial cx :
//
// p0QP = dot(decomp(cx) * gadget[0]) mod QP (encrypted input is multiplied by P factor)
// p1QP = dot(decomp(cx) * gadget[1]) mod QP (encrypted input is multiplied by P factor)
//
// Expects the flag IsNTT of cx to correctly reflect the domain of cx.
func (eval *Evaluator) GadgetProductNoModDown(levelQ int, cx *ring.Poly, gadgetCt GadgetCiphertext, p0QP, p1QP ringqp.Poly) {

	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()
	ringQP := eval.params.RingQP()

	c2QP := eval.BuffQP[0]

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
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, el[i][0].Value[0], c2QP, p0QP)
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, el[i][0].Value[1], c2QP, p1QP)
		} else {
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, el[i][0].Value[0], c2QP, p0QP)
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, el[i][0].Value[1], c2QP, p1QP)
		}

		if reduce%QiOverF == QiOverF-1 {
			ringQ.ReduceLvl(levelQ, p0QP.Q, p0QP.Q)
			ringQ.ReduceLvl(levelQ, p1QP.Q, p1QP.Q)
		}

		if reduce%PiOverF == PiOverF-1 {
			ringP.ReduceLvl(levelP, p0QP.P, p0QP.P)
			ringP.ReduceLvl(levelP, p1QP.P, p1QP.P)
		}

		reduce++
	}

	if reduce%QiOverF != 0 {
		ringQ.ReduceLvl(levelQ, p0QP.Q, p0QP.Q)
		ringQ.ReduceLvl(levelQ, p1QP.Q, p1QP.Q)
	}

	if reduce%PiOverF != 0 {
		ringP.ReduceLvl(levelP, p0QP.P, p0QP.P)
		ringP.ReduceLvl(levelP, p1QP.P, p1QP.P)
	}
}

// GadgetProductSinglePAndBitDecompNoModDown applies the key-switch to the polynomial cx :
//
// p0QP = dot(decomp(cx) * evakey[0]) mod QP (encrypted input is multiplied by P factor)
// p1QP = dot(decomp(cx) * evakey[1]) mod QP (encrypted input is multiplied by P factor)
//
// Expects the flag IsNTT of cx to correctly reflect the domain of cx.
func (eval *Evaluator) GadgetProductSinglePAndBitDecompNoModDown(levelQ int, cx *ring.Poly, gadgetCt GadgetCiphertext, p0QP, p1QP ringqp.Poly) {

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
	decompBIT := eval.params.DecompBIT(levelQ, levelP)

	lb2 := eval.params.logbase2

	mask := uint64(((1 << lb2) - 1))

	if mask == 0 {
		mask = 0xFFFFFFFFFFFFFFFF
	}

	cw := eval.BuffQP[0].Q.Coeffs[0]
	cwNTT := eval.BuffBitDecomp

	QiOverF := eval.params.QiOverflowMargin(levelQ) >> 1
	PiOverF := eval.params.PiOverflowMargin(levelP) >> 1

	el := gadgetCt.Value

	// Key switching with CRT decomposition for the Qi
	var reduce int
	for i := 0; i < decompRNS; i++ {
		for j := 0; j < decompBIT; j++ {

			ring.MaskVec(cxInvNTT.Coeffs[i], cw, j*lb2, mask)

			if i == 0 && j == 0 {
				for u := 0; u < levelQ+1; u++ {
					ringQ.NTTSingleLazy(u, cw, cwNTT)
					ring.MulCoeffsMontgomeryConstantVec(el[i][j].Value[0].Q.Coeffs[u], cwNTT, p0QP.Q.Coeffs[u], ringQ.Modulus[u], ringQ.MredParams[u])
					ring.MulCoeffsMontgomeryConstantVec(el[i][j].Value[1].Q.Coeffs[u], cwNTT, p1QP.Q.Coeffs[u], ringQ.Modulus[u], ringQ.MredParams[u])
				}

				for u := 0; u < levelP+1; u++ {
					ringP.NTTSingleLazy(u, cw, cwNTT)
					ring.MulCoeffsMontgomeryConstantVec(el[i][j].Value[0].P.Coeffs[u], cwNTT, p0QP.P.Coeffs[u], ringP.Modulus[u], ringP.MredParams[u])
					ring.MulCoeffsMontgomeryConstantVec(el[i][j].Value[1].P.Coeffs[u], cwNTT, p1QP.P.Coeffs[u], ringP.Modulus[u], ringP.MredParams[u])
				}
			} else {
				for u := 0; u < levelQ+1; u++ {
					ringQ.NTTSingleLazy(u, cw, cwNTT)
					ring.MulCoeffsMontgomeryConstantAndAddNoModVec(el[i][j].Value[0].Q.Coeffs[u], cwNTT, p0QP.Q.Coeffs[u], ringQ.Modulus[u], ringQ.MredParams[u])
					ring.MulCoeffsMontgomeryConstantAndAddNoModVec(el[i][j].Value[1].Q.Coeffs[u], cwNTT, p1QP.Q.Coeffs[u], ringQ.Modulus[u], ringQ.MredParams[u])
				}

				for u := 0; u < levelP+1; u++ {
					ringP.NTTSingleLazy(u, cw, cwNTT)
					ring.MulCoeffsMontgomeryConstantAndAddNoModVec(el[i][j].Value[0].P.Coeffs[u], cwNTT, p0QP.P.Coeffs[u], ringP.Modulus[u], ringP.MredParams[u])
					ring.MulCoeffsMontgomeryConstantAndAddNoModVec(el[i][j].Value[1].P.Coeffs[u], cwNTT, p1QP.P.Coeffs[u], ringP.Modulus[u], ringP.MredParams[u])
				}
			}

			if reduce%QiOverF == QiOverF-1 {
				ringQ.ReduceLvl(levelQ, p0QP.Q, p0QP.Q)
				ringQ.ReduceLvl(levelQ, p1QP.Q, p1QP.Q)
			}

			if reduce%PiOverF == PiOverF-1 {
				ringP.ReduceLvl(levelP, p0QP.P, p0QP.P)
				ringP.ReduceLvl(levelP, p1QP.P, p1QP.P)
			}

			reduce++
		}
	}

	if reduce%QiOverF != 0 {
		ringQ.ReduceLvl(levelQ, p0QP.Q, p0QP.Q)
		ringQ.ReduceLvl(levelQ, p1QP.Q, p1QP.Q)
	}

	if reduce%PiOverF != 0 {
		ringP.ReduceLvl(levelP, p0QP.P, p0QP.P)
		ringP.ReduceLvl(levelP, p1QP.P, p1QP.P)
	}
}
