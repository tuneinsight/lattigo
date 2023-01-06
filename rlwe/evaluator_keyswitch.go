package rlwe

import (
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// SwitchKeys re-encrypts ctIn under a different key and returns the result in ctOut.
// It requires a SwitchingKey, which is computed from the key under which the Ciphertext is currently encrypted
// and the key under which the Ciphertext will be re-encrypted.
// The method will panic if either ctIn or ctOut degree isn't 1.
func (eval *Evaluator) SwitchKeys(ctIn *Ciphertext, switchingKey *SwitchingKey, ctOut *Ciphertext) {

	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		panic("cannot SwitchKeys: input and output Ciphertext must be of degree 1")
	}

	level := utils.MinInt(ctIn.Level(), ctOut.Level())
	ringQ := eval.params.RingQ().AtLevel(level)

	ctTmp := &Ciphertext{}
	ctTmp.Value = []*ring.Poly{eval.BuffQP[1].Q, eval.BuffQP[2].Q}
	ctTmp.IsNTT = ctIn.IsNTT

	eval.GadgetProduct(level, ctIn.Value[1], switchingKey.GadgetCiphertext, ctTmp)

	ringQ.Add(ctIn.Value[0], ctTmp.Value[0], ctOut.Value[0])
	ring.CopyLvl(level, ctTmp.Value[1], ctOut.Value[1])

	ctOut.MetaData = ctIn.MetaData
}

// Relinearize applies the relinearization procedure on ct0 and returns the result in ctOut.
// The method will panic if the corresponding relinearization key to the ciphertext degree
// is missing.
func (eval *Evaluator) Relinearize(ctIn *Ciphertext, ctOut *Ciphertext) {
	if eval.Rlk == nil || ctIn.Degree()-1 > len(eval.Rlk.Keys) {
		panic("cannot Relinearize: relinearization key missing (or ciphertext degree is too large)")
	}

	level := utils.MinInt(ctIn.Level(), ctOut.Level())

	ringQ := eval.params.RingQ().AtLevel(level)

	ctTmp := &Ciphertext{}
	ctTmp.Value = []*ring.Poly{eval.BuffQP[1].Q, eval.BuffQP[2].Q}
	ctTmp.IsNTT = ctIn.IsNTT

	eval.GadgetProduct(level, ctIn.Value[2], eval.Rlk.Keys[0].GadgetCiphertext, ctTmp)
	ringQ.Add(ctIn.Value[0], ctTmp.Value[0], ctOut.Value[0])
	ringQ.Add(ctIn.Value[1], ctTmp.Value[1], ctOut.Value[1])

	for deg := ctIn.Degree() - 1; deg > 1; deg-- {
		eval.GadgetProduct(level, ctIn.Value[deg], eval.Rlk.Keys[deg-2].GadgetCiphertext, ctTmp)
		ringQ.Add(ctOut.Value[0], ctTmp.Value[0], ctOut.Value[0])
		ringQ.Add(ctOut.Value[1], ctTmp.Value[1], ctOut.Value[1])
	}

	ctOut.Resize(1, level)

	ctOut.MetaData = ctIn.MetaData
}

// DecomposeNTT applies the full RNS basis decomposition on c2.
// Expects the IsNTT flag of c2 to correctly reflect the domain of c2.
// BuffQPDecompQ and BuffQPDecompQ are vectors of polynomials (mod Q and mod P) that store the
// special RNS decomposition of c2 (in the NTT domain)
func (eval *Evaluator) DecomposeNTT(levelQ, levelP, nbPi int, c2 *ring.Poly, c2IsNTT bool, BuffDecompQP []ringqp.Poly) {

	ringQ := eval.params.RingQ().AtLevel(levelQ)

	var polyNTT, polyInvNTT *ring.Poly

	if c2IsNTT {
		polyNTT = c2
		polyInvNTT = eval.BuffInvNTT
		ringQ.INTT(polyNTT, polyInvNTT)
	} else {
		polyNTT = eval.BuffInvNTT
		polyInvNTT = c2
		ringQ.NTT(polyInvNTT, polyNTT)
	}

	decompRNS := eval.params.DecompRNS(levelQ, levelP)
	for i := 0; i < decompRNS; i++ {
		eval.DecomposeSingleNTT(levelQ, levelP, nbPi, i, polyNTT, polyInvNTT, BuffDecompQP[i].Q, BuffDecompQP[i].P)
	}
}

// DecomposeSingleNTT takes the input polynomial c2 (c2NTT and c2InvNTT, respectively in the NTT and out of the NTT domain)
// modulo the RNS basis, and returns the result on c2QiQ and c2QiP, the receiver polynomials respectively mod Q and mod P (in the NTT domain)
func (eval *Evaluator) DecomposeSingleNTT(levelQ, levelP, nbPi, decompRNS int, c2NTT, c2InvNTT, c2QiQ, c2QiP *ring.Poly) {

	ringQ := eval.params.RingQ().AtLevel(levelQ)
	ringP := eval.params.RingP().AtLevel(levelP)

	eval.Decomposer.DecomposeAndSplit(levelQ, levelP, nbPi, decompRNS, c2InvNTT, c2QiQ, c2QiP)

	p0idxst := decompRNS * nbPi
	p0idxed := p0idxst + nbPi

	// c2_qi = cx mod qi mod qi
	for x := 0; x < levelQ+1; x++ {
		if p0idxst <= x && x < p0idxed {
			copy(c2QiQ.Coeffs[x], c2NTT.Coeffs[x])
		} else {
			ringQ.SubRings[x].NTT(c2QiQ.Coeffs[x], c2QiQ.Coeffs[x])
		}
	}

	if ringP != nil {
		// c2QiP = c2 mod qi mod pj
		ringP.NTT(c2QiP, c2QiP)
	}
}

// KeyswitchHoisted applies the key-switch to the decomposed polynomial c2 mod QP (BuffQPDecompQ and BuffQPDecompP)
// and divides the result by P, reducing the basis from QP to Q.
//
// BuffQP2 = dot(BuffQPDecompQ||BuffQPDecompP * evakey[0]) mod Q
// BuffQP3 = dot(BuffQPDecompQ||BuffQPDecompP * evakey[1]) mod Q
func (eval *Evaluator) KeyswitchHoisted(levelQ int, BuffQPDecompQP []ringqp.Poly, evakey *SwitchingKey, c0Q, c1Q, c0P, c1P *ring.Poly) {

	eval.KeyswitchHoistedLazy(levelQ, BuffQPDecompQP, evakey, c0Q, c1Q, c0P, c1P)

	levelP := evakey.Value[0][0].Value[0].P.Level()

	// Computes c0Q = c0Q/c0P and c1Q = c1Q/c1P
	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c0Q, c0P, c0Q)
	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c1Q, c1P, c1Q)
}

// KeyswitchHoistedLazy applies the key-switch to the decomposed polynomial c2 mod QP (BuffQPDecompQ and BuffQPDecompP)
//
// BuffQP2 = dot(BuffQPDecompQ||BuffQPDecompP * evakey[0]) mod QP
// BuffQP3 = dot(BuffQPDecompQ||BuffQPDecompP * evakey[1]) mod QP
func (eval *Evaluator) KeyswitchHoistedLazy(levelQ int, BuffQPDecompQP []ringqp.Poly, evakey *SwitchingKey, c0Q, c1Q, c0P, c1P *ring.Poly) {

	levelP := evakey.LevelP()

	ringQP := eval.params.RingQP().AtLevel(levelQ, levelP)

	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	c0QP := ringqp.Poly{Q: c0Q, P: c0P}
	c1QP := ringqp.Poly{Q: c1Q, P: c1P}

	decompRNS := (levelQ + 1 + levelP) / (levelP + 1)

	QiOverF := eval.params.QiOverflowMargin(levelQ) >> 1
	PiOverF := eval.params.PiOverflowMargin(levelP) >> 1

	// Key switching with CRT decomposition for the Qi
	var reduce int
	for i := 0; i < decompRNS; i++ {

		if i == 0 {
			ringQP.MulCoeffsMontgomeryLazy(evakey.Value[i][0].Value[0], BuffQPDecompQP[i], c0QP)
			ringQP.MulCoeffsMontgomeryLazy(evakey.Value[i][0].Value[1], BuffQPDecompQP[i], c1QP)
		} else {
			ringQP.MulCoeffsMontgomeryLazyThenAddLazy(evakey.Value[i][0].Value[0], BuffQPDecompQP[i], c0QP)
			ringQP.MulCoeffsMontgomeryLazyThenAddLazy(evakey.Value[i][0].Value[1], BuffQPDecompQP[i], c1QP)
		}

		if reduce%QiOverF == QiOverF-1 {
			ringQ.Reduce(c0QP.Q, c0QP.Q)
			ringQ.Reduce(c1QP.Q, c1QP.Q)
		}

		if reduce%PiOverF == PiOverF-1 {
			ringP.Reduce(c0QP.P, c0QP.P)
			ringP.Reduce(c1QP.P, c1QP.P)
		}

		reduce++
	}

	if reduce%QiOverF != 0 {
		ringQ.Reduce(c0QP.Q, c0QP.Q)
		ringQ.Reduce(c1QP.Q, c1QP.Q)
	}

	if reduce%PiOverF != 0 {
		ringP.Reduce(c0QP.P, c0QP.P)
		ringP.Reduce(c1QP.P, c1QP.P)
	}
}
