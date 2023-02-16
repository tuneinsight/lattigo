package rlwe

import (
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// GadgetProduct evaluates poly x Gadget -> RLWE where
//
// p0 = dot(decomp(cx) * gadget[0]) mod Q
// p1 = dot(decomp(cx) * gadget[1]) mod Q
//
// Expects the flag IsNTT of ct to correctly reflect the domain of cx.
func (eval *Evaluator) GadgetProduct(levelQ int, cx *ring.Poly, gadgetCt GadgetCiphertext, ct *Ciphertext) {

	levelQ = utils.MinInt(levelQ, gadgetCt.LevelQ())
	levelP := gadgetCt.LevelP()

	ctTmp := CiphertextQP{}
	ctTmp.Value = [2]ringqp.Poly{{Q: ct.Value[0], P: eval.BuffQP[1].P}, {Q: ct.Value[1], P: eval.BuffQP[2].P}}
	ctTmp.IsNTT = ct.IsNTT

	if levelP > 0 {
		eval.GadgetProductLazy(levelQ, cx, gadgetCt, ctTmp)
	} else {
		eval.GadgetProductSinglePAndBitDecompLazy(levelQ, cx, gadgetCt, ctTmp)
	}

	if ct.IsNTT && levelP != -1 {
		eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, ct.Value[0], ctTmp.Value[0].P, ct.Value[0])
		eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, ct.Value[1], ctTmp.Value[1].P, ct.Value[1])
	} else if !ct.IsNTT {

		ringQ := eval.params.RingQ().AtLevel(levelQ)

		if levelP != -1 {

			ringQ.INTTLazy(ct.Value[0], ct.Value[0])
			ringQ.INTTLazy(ct.Value[1], ct.Value[1])

			ringP := eval.params.RingP().AtLevel(levelP)

			ringP.INTTLazy(ctTmp.Value[0].P, ctTmp.Value[0].P)
			ringP.INTTLazy(ctTmp.Value[1].P, ctTmp.Value[1].P)

			eval.BasisExtender.ModDownQPtoQ(levelQ, levelP, ct.Value[0], ctTmp.Value[0].P, ct.Value[0])
			eval.BasisExtender.ModDownQPtoQ(levelQ, levelP, ct.Value[1], ctTmp.Value[1].P, ct.Value[1])
		} else {
			ringQ.INTT(ct.Value[0], ct.Value[0])
			ringQ.INTT(ct.Value[1], ct.Value[1])
		}
	}
}

// GadgetProductLazy applies the gadget prodcut to the polynomial cx:
//
// ct.Value[0] = dot(decomp(cx) * gadget[0]) mod QP (encrypted input is multiplied by P factor)
// ct.Value[1] = dot(decomp(cx) * gadget[1]) mod QP (encrypted input is multiplied by P factor)
//
// Expects the flag IsNTT of ct to correctly reflect the domain of cx.
func (eval *Evaluator) GadgetProductLazy(levelQ int, cx *ring.Poly, gadgetCt GadgetCiphertext, ct CiphertextQP) {

	levelP := gadgetCt.LevelP()

	ringQP := eval.params.RingQP().AtLevel(levelQ, levelP)

	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	c2QP := eval.BuffQP[0]

	var cxNTT, cxInvNTT *ring.Poly
	if ct.IsNTT {
		cxNTT = cx
		cxInvNTT = eval.BuffInvNTT
		ringQ.INTT(cxNTT, cxInvNTT)
	} else {
		cxNTT = eval.BuffInvNTT
		cxInvNTT = cx
		ringQ.NTT(cxInvNTT, cxNTT)
	}

	decompRNS := eval.params.DecompRNS(levelQ, levelP)

	QiOverF := eval.params.QiOverflowMargin(levelQ) >> 1
	PiOverF := eval.params.PiOverflowMargin(levelP) >> 1

	el := gadgetCt.Value

	// Key switching with CRT decomposition for the Qi
	var reduce int
	for i := 0; i < decompRNS; i++ {

		eval.DecomposeSingleNTT(levelQ, levelP, levelP+1, i, cxNTT, cxInvNTT, c2QP.Q, c2QP.P)

		if i == 0 {
			ringQP.MulCoeffsMontgomeryLazy(el[i][0].Value[0], c2QP, ct.Value[0])
			ringQP.MulCoeffsMontgomeryLazy(el[i][0].Value[1], c2QP, ct.Value[1])
		} else {
			ringQP.MulCoeffsMontgomeryLazyThenAddLazy(el[i][0].Value[0], c2QP, ct.Value[0])
			ringQP.MulCoeffsMontgomeryLazyThenAddLazy(el[i][0].Value[1], c2QP, ct.Value[1])
		}

		if reduce%QiOverF == QiOverF-1 {
			ringQ.Reduce(ct.Value[0].Q, ct.Value[0].Q)
			ringQ.Reduce(ct.Value[1].Q, ct.Value[1].Q)
		}

		if reduce%PiOverF == PiOverF-1 {
			ringP.Reduce(ct.Value[0].P, ct.Value[0].P)
			ringP.Reduce(ct.Value[1].P, ct.Value[1].P)
		}

		reduce++
	}

	if reduce%QiOverF != 0 {
		ringQ.Reduce(ct.Value[0].Q, ct.Value[0].Q)
		ringQ.Reduce(ct.Value[1].Q, ct.Value[1].Q)
	}

	if reduce%PiOverF != 0 {
		ringP.Reduce(ct.Value[0].P, ct.Value[0].P)
		ringP.Reduce(ct.Value[1].P, ct.Value[1].P)
	}
}

// GadgetProductSinglePAndBitDecompLazy applies the key-switch to the polynomial cx:
//
// ct.Value[0] = dot(decomp(cx) * evakey[0]) mod QP (encrypted input is multiplied by P factor)
// ct.Value[1] = dot(decomp(cx) * evakey[1]) mod QP (encrypted input is multiplied by P factor)
//
// Expects the flag IsNTT of ct to correctly reflect the domain of cx.
func (eval *Evaluator) GadgetProductSinglePAndBitDecompLazy(levelQ int, cx *ring.Poly, gadgetCt GadgetCiphertext, ct CiphertextQP) {

	levelP := gadgetCt.LevelP()

	ringQP := eval.params.RingQP().AtLevel(levelQ, levelP)

	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	var cxInvNTT *ring.Poly
	if ct.IsNTT {
		cxInvNTT = eval.BuffInvNTT
		ringQ.INTT(cx, cxInvNTT)
	} else {
		cxInvNTT = cx
	}

	decompRNS := eval.params.DecompRNS(levelQ, levelP)
	decompPw2 := eval.params.DecompPw2(levelQ, levelP)

	pw2 := eval.params.pow2Base

	mask := uint64(((1 << pw2) - 1))

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
		for j := 0; j < decompPw2; j++ {

			ring.MaskVec(cxInvNTT.Coeffs[i], j*pw2, mask, cw)

			if i == 0 && j == 0 {
				for u, s := range ringQ.SubRings[:levelQ+1] {
					s.NTTLazy(cw, cwNTT)
					s.MulCoeffsMontgomeryLazy(el[i][j].Value[0].Q.Coeffs[u], cwNTT, ct.Value[0].Q.Coeffs[u])
					s.MulCoeffsMontgomeryLazy(el[i][j].Value[1].Q.Coeffs[u], cwNTT, ct.Value[1].Q.Coeffs[u])
				}

				if ringP != nil {
					for u, s := range ringP.SubRings[:levelP+1] {
						s.NTTLazy(cw, cwNTT)
						s.MulCoeffsMontgomeryLazy(el[i][j].Value[0].P.Coeffs[u], cwNTT, ct.Value[0].P.Coeffs[u])
						s.MulCoeffsMontgomeryLazy(el[i][j].Value[1].P.Coeffs[u], cwNTT, ct.Value[1].P.Coeffs[u])
					}
				}

			} else {
				for u, s := range ringQ.SubRings[:levelQ+1] {
					s.NTTLazy(cw, cwNTT)
					s.MulCoeffsMontgomeryLazyThenAddLazy(el[i][j].Value[0].Q.Coeffs[u], cwNTT, ct.Value[0].Q.Coeffs[u])
					s.MulCoeffsMontgomeryLazyThenAddLazy(el[i][j].Value[1].Q.Coeffs[u], cwNTT, ct.Value[1].Q.Coeffs[u])
				}

				if ringP != nil {
					for u, s := range ringP.SubRings[:levelP+1] {
						s.NTTLazy(cw, cwNTT)
						s.MulCoeffsMontgomeryLazyThenAddLazy(el[i][j].Value[0].P.Coeffs[u], cwNTT, ct.Value[0].P.Coeffs[u])
						s.MulCoeffsMontgomeryLazyThenAddLazy(el[i][j].Value[1].P.Coeffs[u], cwNTT, ct.Value[1].P.Coeffs[u])
					}
				}
			}

			if reduce%QiOverF == QiOverF-1 {
				ringQ.Reduce(ct.Value[0].Q, ct.Value[0].Q)
				ringQ.Reduce(ct.Value[1].Q, ct.Value[1].Q)
			}

			if reduce%PiOverF == PiOverF-1 {
				ringP.Reduce(ct.Value[0].P, ct.Value[0].P)
				ringP.Reduce(ct.Value[1].P, ct.Value[1].P)
			}

			reduce++
		}
	}

	if reduce%QiOverF != 0 {
		ringQ.Reduce(ct.Value[0].Q, ct.Value[0].Q)
		ringQ.Reduce(ct.Value[1].Q, ct.Value[1].Q)
	}

	if ringP != nil {
		if reduce%PiOverF != 0 {
			ringP.Reduce(ct.Value[0].P, ct.Value[0].P)
			ringP.Reduce(ct.Value[1].P, ct.Value[1].P)
		}
	}
}
