package rlwe

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// SwitchKeys re-encrypts ctIn under a different key and returns the result in ctOut.
// It requires a SwitchingKey, which is computed from the key under which the Ciphertext is currently encrypted,
// and the key under which the Ciphertext will be re-encrypted.
// The method will panic if either ctIn or ctOut degree isn't 1.
func (eval *Evaluator) SwitchKeys(ctIn *Ciphertext, switchingKey *SwitchingKey, ctOut *Ciphertext) {

	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		panic("cannot SwitchKeys: input and output Ciphertext must be of degree 1")
	}

	level := utils.MinInt(ctIn.Level(), ctOut.Level())
	ringQ := eval.params.RingQ()

	eval.SwitchKeysInPlace(level, ctIn.Value[1], switchingKey, eval.BuffQP[1].Q, eval.BuffQP[2].Q)

	ringQ.AddLvl(level, ctIn.Value[0], eval.BuffQP[1].Q, ctOut.Value[0])
	ring.CopyValuesLvl(level, eval.BuffQP[2].Q, ctOut.Value[1])
}

// SwitchKeysInPlace applies the general key-switching procedure of the form [c0 + cx*evakey[0], c1 + cx*evakey[1]]
// Will return the result in the same NTT domain as the input cx.
func (eval *Evaluator) SwitchKeysInPlace(levelQ int, cx *ring.Poly, evakey *SwitchingKey, p0, p1 *ring.Poly) {

	levelP := evakey.LevelP()

	if levelP > 0 {
		eval.SwitchKeysInPlaceNoModDown(levelQ, cx, evakey, p0, eval.BuffQP[1].P, p1, eval.BuffQP[2].P)
	} else {
		eval.SwitchKeyInPlaceSinglePAndBitDecomp(levelQ, cx, evakey, p0, eval.BuffQP[1].P, p1, eval.BuffQP[2].P)
	}

	if cx.IsNTT && levelP != -1 {
		eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, p0, eval.BuffQP[1].P, p0)
		eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, p1, eval.BuffQP[2].P, p1)
	} else if !cx.IsNTT {
		eval.params.RingQ().InvNTTLazyLvl(levelQ, p0, p0)
		eval.params.RingQ().InvNTTLazyLvl(levelQ, p1, p1)

		if levelP != -1 {
			eval.params.RingP().InvNTTLazyLvl(levelP, eval.BuffQP[1].P, eval.BuffQP[1].P)
			eval.params.RingP().InvNTTLazyLvl(levelP, eval.BuffQP[2].P, eval.BuffQP[2].P)
			eval.BasisExtender.ModDownQPtoQ(levelQ, levelP, p0, eval.BuffQP[1].P, p0)
			eval.BasisExtender.ModDownQPtoQ(levelQ, levelP, p1, eval.BuffQP[2].P, p1)
		}
	}
}

// Relinearize applies the relinearization procedure on ct0 and returns the result in ctOut.
// The method will panic if the corresponding relinearization key to the ciphertext degree
// is missing.
func (eval *Evaluator) Relinearize(ctIn *Ciphertext, ctOut *Ciphertext) {
	if eval.Rlk == nil || ctIn.Degree()-1 > len(eval.Rlk.Keys) {
		panic("cannot Relinearize: relinearization key missing (or ciphertext degree is too large)")
	}

	level := utils.MinInt(ctIn.Level(), ctOut.Level())

	ringQ := eval.params.RingQ()

	eval.SwitchKeysInPlace(level, ctIn.Value[2], eval.Rlk.Keys[0], eval.BuffQP[1].Q, eval.BuffQP[2].Q)
	ringQ.AddLvl(level, ctIn.Value[0], eval.BuffQP[1].Q, ctOut.Value[0])
	ringQ.AddLvl(level, ctIn.Value[1], eval.BuffQP[2].Q, ctOut.Value[1])

	for deg := ctIn.Degree() - 1; deg > 1; deg-- {
		eval.SwitchKeysInPlace(level, ctIn.Value[deg], eval.Rlk.Keys[deg-2], eval.BuffQP[1].Q, eval.BuffQP[2].Q)
		ringQ.AddLvl(level, ctOut.Value[0], eval.BuffQP[1].Q, ctOut.Value[0])
		ringQ.AddLvl(level, ctOut.Value[1], eval.BuffQP[2].Q, ctOut.Value[1])
	}

	ctOut.Value = ctOut.Value[:2]

	ctOut.Value[0].Coeffs = ctOut.Value[0].Coeffs[:level+1]
	ctOut.Value[1].Coeffs = ctOut.Value[1].Coeffs[:level+1]
}

// DecomposeNTT applies the full RNS basis decomposition for all q_alpha_i on c2.
// Expects the IsNTT flag of c2 to correctly reflect the domain of c2.
// BuffQPDecompQ and BuffQPDecompQ are vectors of polynomials (mod Q and mod P) that store the
// special RNS decomposition of c2 (in the NTT domain)
func (eval *Evaluator) DecomposeNTT(levelQ, levelP, alpha int, c2 *ring.Poly, BuffDecompQP []ringqp.Poly) {

	ringQ := eval.params.RingQ()

	var polyNTT, polyInvNTT *ring.Poly

	if c2.IsNTT {
		polyNTT = c2
		polyInvNTT = eval.BuffInvNTT
		ringQ.InvNTTLvl(levelQ, polyNTT, polyInvNTT)
	} else {
		polyNTT = eval.BuffInvNTT
		polyInvNTT = c2
		ringQ.NTTLvl(levelQ, polyInvNTT, polyNTT)
	}

	beta := (levelQ + 1 + levelP) / (levelP + 1)

	for i := 0; i < beta; i++ {
		eval.DecomposeSingleNTT(levelQ, levelP, alpha, i, polyNTT, polyInvNTT, BuffDecompQP[i].Q, BuffDecompQP[i].P)
	}
}

// DecomposeSingleNTT takes the input polynomial c2 (c2NTT and c2InvNTT, respectively in the NTT and out of the NTT domain)
// modulo q_alpha_beta, and returns the result on c2QiQ are c2QiP the receiver polynomials
// respectively mod Q and mod P (in the NTT domain)
func (eval *Evaluator) DecomposeSingleNTT(levelQ, levelP, alpha, beta int, c2NTT, c2InvNTT, c2QiQ, c2QiP *ring.Poly) {

	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()

	eval.Decomposer.DecomposeAndSplit(levelQ, levelP, alpha, beta, c2InvNTT, c2QiQ, c2QiP)

	p0idxst := beta * alpha
	p0idxed := p0idxst + 1

	// c2_qi = cx mod qi mod qi
	for x := 0; x < levelQ+1; x++ {
		if p0idxst <= x && x < p0idxed {
			copy(c2QiQ.Coeffs[x], c2NTT.Coeffs[x])
		} else {
			ringQ.NTTSingle(x, c2QiQ.Coeffs[x], c2QiQ.Coeffs[x])
		}
	}

	if ringP != nil {
		// c2QiP = c2 mod qi mod pj
		ringP.NTTLvl(levelP, c2QiP, c2QiP)
	}
}

// SwitchKeysInPlaceNoModDown applies the key-switch to the polynomial cx :
//
// BuffQP2 = dot(decomp(cx) * evakey[0]) mod QP (encrypted input is multiplied by P factor)
// BuffQP3 = dot(decomp(cx) * evakey[1]) mod QP (encrypted input is multiplied by P factor)
//
// Expects the flag IsNTT of cx to correctly reflect the domain of cx.
func (eval *Evaluator) SwitchKeysInPlaceNoModDown(levelQ int, cx *ring.Poly, evakey *SwitchingKey, c0Q, c0P, c1Q, c1P *ring.Poly) {

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

	c0QP := ringqp.Poly{Q: c0Q, P: c0P}
	c1QP := ringqp.Poly{Q: c1Q, P: c1P}

	levelP := evakey.Value[0][0][0].P.Level()
	decompRNS := (levelQ + 1 + levelP) / (levelP + 1)

	QiOverF := eval.params.QiOverflowMargin(levelQ) >> 1
	PiOverF := eval.params.PiOverflowMargin(levelP) >> 1

	// Key switching with CRT decomposition for the Qi
	var reduce int
	for i := 0; i < decompRNS; i++ {

		eval.DecomposeSingleNTT(levelQ, levelP, levelP+1, i, cxNTT, cxInvNTT, c2QP.Q, c2QP.P)

		if i == 0 {
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, evakey.Value[i][0][0], c2QP, c0QP)
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, evakey.Value[i][0][1], c2QP, c1QP)
		} else {
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, evakey.Value[i][0][0], c2QP, c0QP)
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, evakey.Value[i][0][1], c2QP, c1QP)
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

// SwitchKeyInPlaceSinglePAndBitDecomp applies the key-switch to the polynomial cx :
//
// BuffQP2 = dot(decomp(cx) * evakey[0]) mod QP (encrypted input is multiplied by P factor)
// BuffQP3 = dot(decomp(cx) * evakey[1]) mod QP (encrypted input is multiplied by P factor)
//
// Expects the flag IsNTT of cx to correctly reflect the domain of cx.
func (eval *Evaluator) SwitchKeyInPlaceSinglePAndBitDecomp(levelQ int, cx *ring.Poly, evakey *SwitchingKey, c0Q, c0P, c1Q, c1P *ring.Poly) {

	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()

	var cxInvNTT *ring.Poly
	if cx.IsNTT {
		cxInvNTT = eval.BuffInvNTT
		ringQ.InvNTTLvl(levelQ, cx, cxInvNTT)
	} else {
		cxInvNTT = cx
	}

	c0QP := ringqp.Poly{Q: c0Q, P: c0P}
	c1QP := ringqp.Poly{Q: c1Q, P: c1P}

	var levelP int
	if evakey.Value[0][0][0].P != nil {
		levelP = evakey.Value[0][0][0].P.Level()
	} else {
		levelP = -1
	}

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

	// Key switching with CRT decomposition for the Qi
	var reduce int
	for i := 0; i < decompRNS; i++ {
		for j := 0; j < decompBIT; j++ {

			ring.MaskVec(cxInvNTT.Coeffs[i], cw, j*lb2, mask)

			if i == 0 && j == 0 {
				for u := 0; u < levelQ+1; u++ {
					ringQ.NTTSingleLazy(u, cw, cwNTT)
					ring.MulCoeffsMontgomeryConstantVec(evakey.Value[i][j][0].Q.Coeffs[u], cwNTT, c0QP.Q.Coeffs[u], ringQ.Modulus[u], ringQ.MredParams[u])
					ring.MulCoeffsMontgomeryConstantVec(evakey.Value[i][j][1].Q.Coeffs[u], cwNTT, c1QP.Q.Coeffs[u], ringQ.Modulus[u], ringQ.MredParams[u])
				}

				for u := 0; u < levelP+1; u++ {
					ringP.NTTSingleLazy(u, cw, cwNTT)
					ring.MulCoeffsMontgomeryConstantVec(evakey.Value[i][j][0].P.Coeffs[u], cwNTT, c0QP.P.Coeffs[u], ringP.Modulus[u], ringP.MredParams[u])
					ring.MulCoeffsMontgomeryConstantVec(evakey.Value[i][j][1].P.Coeffs[u], cwNTT, c1QP.P.Coeffs[u], ringP.Modulus[u], ringP.MredParams[u])
				}
			} else {
				for u := 0; u < levelQ+1; u++ {
					ringQ.NTTSingleLazy(u, cw, cwNTT)
					ring.MulCoeffsMontgomeryConstantAndAddNoModVec(evakey.Value[i][j][0].Q.Coeffs[u], cwNTT, c0QP.Q.Coeffs[u], ringQ.Modulus[u], ringQ.MredParams[u])
					ring.MulCoeffsMontgomeryConstantAndAddNoModVec(evakey.Value[i][j][1].Q.Coeffs[u], cwNTT, c1QP.Q.Coeffs[u], ringQ.Modulus[u], ringQ.MredParams[u])
				}

				for u := 0; u < levelP+1; u++ {
					ringP.NTTSingleLazy(u, cw, cwNTT)
					ring.MulCoeffsMontgomeryConstantAndAddNoModVec(evakey.Value[i][j][0].P.Coeffs[u], cwNTT, c0QP.P.Coeffs[u], ringP.Modulus[u], ringP.MredParams[u])
					ring.MulCoeffsMontgomeryConstantAndAddNoModVec(evakey.Value[i][j][1].P.Coeffs[u], cwNTT, c1QP.P.Coeffs[u], ringP.Modulus[u], ringP.MredParams[u])
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

	if reduce%QiOverF != 0 {
		ringQ.ReduceLvl(levelQ, c0QP.Q, c0QP.Q)
		ringQ.ReduceLvl(levelQ, c1QP.Q, c1QP.Q)
	}

	if reduce%PiOverF != 0 {
		ringP.ReduceLvl(levelP, c0QP.P, c0QP.P)
		ringP.ReduceLvl(levelP, c1QP.P, c1QP.P)
	}
}

// KeyswitchHoisted applies the key-switch to the decomposed polynomial c2 mod QP (BuffQPDecompQ and BuffQPDecompP)
// and divides the result by P, reducing the basis from QP to Q.
//
// BuffQP2 = dot(BuffQPDecompQ||BuffQPDecompP * evakey[0]) mod Q
// BuffQP3 = dot(BuffQPDecompQ||BuffQPDecompP * evakey[1]) mod Q
func (eval *Evaluator) KeyswitchHoisted(levelQ int, BuffQPDecompQP []ringqp.Poly, evakey *SwitchingKey, c0Q, c1Q, c0P, c1P *ring.Poly) {

	eval.KeyswitchHoistedNoModDown(levelQ, BuffQPDecompQP, evakey, c0Q, c1Q, c0P, c1P)

	levelP := evakey.Value[0][0][0].P.Level()

	// Computes c0Q = c0Q/c0P and c1Q = c1Q/c1P
	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c0Q, c0P, c0Q)
	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c1Q, c1P, c1Q)
}

// KeyswitchHoistedNoModDown applies the key-switch to the decomposed polynomial c2 mod QP (BuffQPDecompQ and BuffQPDecompP)
//
// BuffQP2 = dot(BuffQPDecompQ||BuffQPDecompP * evakey[0]) mod QP
// BuffQP3 = dot(BuffQPDecompQ||BuffQPDecompP * evakey[1]) mod QP
func (eval *Evaluator) KeyswitchHoistedNoModDown(levelQ int, BuffQPDecompQP []ringqp.Poly, evakey *SwitchingKey, c0Q, c1Q, c0P, c1P *ring.Poly) {

	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()
	ringQP := eval.params.RingQP()

	c0QP := ringqp.Poly{Q: c0Q, P: c0P}
	c1QP := ringqp.Poly{Q: c1Q, P: c1P}

	levelP := evakey.Value[0][0][0].P.Level()
	decompRNS := (levelQ + 1 + levelP) / (levelP + 1)

	QiOverF := eval.params.QiOverflowMargin(levelQ) >> 1
	PiOverF := eval.params.PiOverflowMargin(levelP) >> 1

	// Key switching with CRT decomposition for the Qi
	var reduce int
	for i := 0; i < decompRNS; i++ {

		if i == 0 {
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, evakey.Value[i][0][0], BuffQPDecompQP[i], c0QP)
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, evakey.Value[i][0][1], BuffQPDecompQP[i], c1QP)
		} else {
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, evakey.Value[i][0][0], BuffQPDecompQP[i], c0QP)
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, evakey.Value[i][0][1], BuffQPDecompQP[i], c1QP)
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
