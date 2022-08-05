package rlwe

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// Automorphism computes phi(ct), where phi is the map X -> X^galEl. The method requires
// that the corresponding RotationKey has been added to the Evaluator. The method will
// panic if either ctIn or ctOut degree is not equal to 1.
func (eval *Evaluator) Automorphism(ctIn *Ciphertext, galEl uint64, ctOut *Ciphertext) {

	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		panic("cannot apply Automorphism: input and output Ciphertext must be of degree 1")
	}

	if galEl == 1 {
		if ctOut != ctIn {
			ctOut.Copy(ctIn)
		}
		return
	}

	rtk, generated := eval.Rtks.GetRotationKey(galEl)
	if !generated {
		panic(fmt.Sprintf("galEl key 5^%d missing", eval.params.InverseGaloisElement(galEl)))
	}

	level := utils.MinInt(ctIn.Level(), ctOut.Level())

	ringQ := eval.params.RingQ()

	eval.GadgetProduct(level, ctIn.Value[1], rtk.GadgetCiphertext, Ciphertext{Value: []*ring.Poly{eval.BuffQP[0].Q, eval.BuffQP[1].Q}})
	ringQ.AddLvl(level, eval.BuffQP[0].Q, ctIn.Value[0], eval.BuffQP[0].Q)

	if ctIn.Value[0].IsNTT {
		ringQ.PermuteNTTWithIndexLvl(level, eval.BuffQP[0].Q, eval.PermuteNTTIndex[galEl], ctOut.Value[0])
		ringQ.PermuteNTTWithIndexLvl(level, eval.BuffQP[1].Q, eval.PermuteNTTIndex[galEl], ctOut.Value[1])
	} else {
		ringQ.PermuteLvl(level, eval.BuffQP[0].Q, galEl, ctOut.Value[0])
		ringQ.PermuteLvl(level, eval.BuffQP[1].Q, galEl, ctOut.Value[1])
	}

	ctOut.Resize(ctOut.Degree(), level)
}

// AutomorphismHoisted is similar to Automorphism, except that it takes as input ctIn and c1DecompQP, where c1DecompQP is the RNS
// decomposition of its element of degree 1. This decomposition can be obtained with DecomposeNTT.
// The method requires that the corresponding RotationKey  has been added to the Evaluator.
// The method will panic if either ctIn or ctOut degree is not equal to 1.
func (eval *Evaluator) AutomorphismHoisted(level int, ctIn *Ciphertext, c1DecompQP []ringqp.Poly, galEl uint64, ctOut *Ciphertext) {

	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		panic("cannot apply AutomorphismHoisted: input and output Ciphertext must be of degree 1")
	}

	if galEl == 1 {
		if ctIn != ctOut {
			ctOut.Copy(ctIn)
		}
		return
	}

	rtk, generated := eval.Rtks.GetRotationKey(galEl)
	if !generated {
		panic(fmt.Sprintf("galEl key 5^%d missing", eval.params.InverseGaloisElement(galEl)))
	}

	ringQ := eval.params.RingQ()

	eval.GadgetProductHoisted(level, c1DecompQP, rtk.GadgetCiphertext, Ciphertext{Value: []*ring.Poly{eval.BuffQP[0].Q, eval.BuffQP[1].Q}})
	ringQ.AddLvl(level, eval.BuffQP[0].Q, ctIn.Value[0], eval.BuffQP[0].Q)

	if ctIn.Value[0].IsNTT {
		ringQ.PermuteNTTWithIndexLvl(level, eval.BuffQP[0].Q, eval.PermuteNTTIndex[galEl], ctOut.Value[0])
		ringQ.PermuteNTTWithIndexLvl(level, eval.BuffQP[1].Q, eval.PermuteNTTIndex[galEl], ctOut.Value[1])
	} else {
		ringQ.PermuteLvl(level, eval.BuffQP[0].Q, galEl, ctOut.Value[0])
		ringQ.PermuteLvl(level, eval.BuffQP[1].Q, galEl, ctOut.Value[1])
	}

	ctOut.Resize(ctOut.Degree(), level)
}

// AutomorphismHoistedNoModDown is similar to AutomorphismHoisted, except that it returns a ciphertext modulo QP and scaled by P.
// The method requires that the corresponding RotationKey has been added to the Evaluator.The method will panic if either ctIn or ctOut degree is not equal to 1.
func (eval *Evaluator) AutomorphismHoistedNoModDown(levelQ int, c0 *ring.Poly, c1DecompQP []ringqp.Poly, galEl uint64, ctOut CiphertextQP) {

	rtk, generated := eval.Rtks.GetRotationKey(galEl)
	if !generated {
		panic(fmt.Sprintf("galEl key 5^%d missing", eval.params.InverseGaloisElement(galEl)))
	}

	levelP := rtk.LevelP()

	eval.GadgetProductHoistedNoModDown(levelQ, c1DecompQP, rtk.GadgetCiphertext, CiphertextQP{Value: []ringqp.Poly{eval.BuffQP[0], eval.BuffQP[1]}})

	ringQ := eval.params.RingQ()

	if c0.IsNTT {

		index := eval.PermuteNTTIndex[galEl]

		ringQ.PermuteNTTWithIndexLvl(levelQ, eval.BuffQP[1].Q, index, ctOut.Value[1].Q)
		ringQ.PermuteNTTWithIndexLvl(levelP, eval.BuffQP[1].P, index, ctOut.Value[1].P)

		if levelP > -1 {
			ringQ.MulScalarBigintLvl(levelQ, c0, eval.params.RingP().ModulusAtLevel[levelP], eval.BuffQP[1].Q)
		}

		ringQ.AddLvl(levelQ, eval.BuffQP[0].Q, eval.BuffQP[1].Q, eval.BuffQP[0].Q)

		ringQ.PermuteNTTWithIndexLvl(levelQ, eval.BuffQP[0].Q, index, ctOut.Value[0].Q)
		ringQ.PermuteNTTWithIndexLvl(levelP, eval.BuffQP[0].P, index, ctOut.Value[0].P)
	} else {
		ringQ.PermuteLvl(levelQ, eval.BuffQP[1].Q, galEl, ctOut.Value[1].Q)
		ringQ.PermuteLvl(levelP, eval.BuffQP[1].P, galEl, ctOut.Value[1].P)

		if levelP > -1 {
			ringQ.MulScalarBigintLvl(levelQ, c0, eval.params.RingP().ModulusAtLevel[levelP], eval.BuffQP[1].Q)
		}

		ringQ.AddLvl(levelQ, eval.BuffQP[0].Q, eval.BuffQP[1].Q, eval.BuffQP[0].Q)

		ringQ.PermuteLvl(levelQ, eval.BuffQP[0].Q, galEl, ctOut.Value[1].Q)
		ringQ.PermuteLvl(levelP, eval.BuffQP[0].P, galEl, ctOut.Value[1].P)
	}
}

// Trace maps X -> sum((-1)^i * X^{i*n+1}) for n <= i < N
// Monomial X^k vanishes if k is not divisible by (N/n), else it is multiplied by (N/n).
// Ciphertext is pre-multiplied by (N/n)^-1 to remove the (N/n) factor.
// Examples of full Trace for [0 + 1X + 2X^2 + 3X^3 + 4X^4 + 5X^5 + 6X^6 + 7X^7]
//
// 1)	   [1 + 2X + 3X^2 + 4X^3 + 5X^4 + 6X^5 + 7X^6 + 8X^7]
//       + [1 - 6X - 3X^2 + 8X^3 + 5X^4 + 2X^5 - 7X^6 - 4X^7]  {X-> X^(i * 5^1)}
//   	 = [2 - 4X + 0X^2 +12X^3 +10X^4 + 8X^5 - 0X^6 + 4X^7]
//
// 2)      [2 - 4X + 0X^2 +12X^3 +10X^4 + 8X^5 - 0X^6 + 4X^7]
//       + [2 + 4X + 0X^2 -12X^3 +10X^4 - 8X^5 + 0X^6 - 4X^7]  {X-> X^(i * 5^2)}
//       = [4 + 0X + 0X^2 - 0X^3 +20X^4 + 0X^5 + 0X^6 - 0X^7]
//
// 3)      [4 + 0X + 0X^2 - 0X^3 +20X^4 + 0X^5 + 0X^6 - 0X^7]
//       + [4 + 0X + 0X^2 - 0X^3 -20X^4 + 0X^5 + 0X^6 - 0X^7]  {X-> X^(i * -1)}
//       = [8 + 0X + 0X^2 - 0X^3 + 0X^4 + 0X^5 + 0X^6 - 0X^7]
func (eval *Evaluator) Trace(ctIn *Ciphertext, logN int, ctOut *Ciphertext) {

	levelQ := utils.MinInt(ctIn.Level(), ctOut.Level())

	ctOut.Resize(ctOut.Degree(), levelQ)

	gap := 1 << (eval.params.LogN() - logN - 1)

	if logN == 0 {
		gap <<= 1
	}

	if gap > 1 {

		ringQ := eval.params.RingQ()

		// pre-multiplication by (N/n)^-1
		for i := 0; i < levelQ+1; i++ {
			Q := ringQ.Modulus[i]
			bredParams := ringQ.BredParams[i]
			mredparams := ringQ.MredParams[i]
			invN := ring.ModExp(uint64(gap), Q-2, Q)
			invN = ring.MForm(invN, Q, bredParams)

			ring.MulScalarMontgomeryVec(ctIn.Value[0].Coeffs[i], ctOut.Value[0].Coeffs[i], invN, Q, mredparams)
			ring.MulScalarMontgomeryVec(ctIn.Value[1].Coeffs[i], ctOut.Value[1].Coeffs[i], invN, Q, mredparams)
		}

		buff := NewCiphertextAtLevelFromPoly(levelQ, [2]*ring.Poly{eval.BuffQP[3].Q, eval.BuffQP[4].Q})
		buff.Value[0].IsNTT = true
		buff.Value[1].IsNTT = true

		for i := logN; i < eval.params.LogN()-1; i++ {
			eval.Automorphism(ctOut, eval.params.GaloisElementForColumnRotationBy(1<<i), buff)
			ringQ.AddLvl(levelQ, ctOut.Value[0], buff.Value[0], ctOut.Value[0])
			ringQ.AddLvl(levelQ, ctOut.Value[1], buff.Value[1], ctOut.Value[1])
		}

		if logN == 0 {
			eval.Automorphism(ctOut, ringQ.NthRoot-1, buff)
			ringQ.AddLvl(levelQ, ctOut.Value[0], buff.Value[0], ctOut.Value[0])
			ringQ.AddLvl(levelQ, ctOut.Value[1], buff.Value[1], ctOut.Value[1])
		}

	} else {
		if ctIn != ctOut {
			ctOut.Copy(ctIn)
		}
	}
}
