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

	eval.SwitchKeysInPlace(level, ctIn.Value[1], rtk, eval.Pool[1].Q, eval.Pool[2].Q)
	ringQ.AddLvl(level, eval.Pool[1].Q, ctIn.Value[0], eval.Pool[1].Q)

	if ctIn.Value[0].IsNTT {
		ringQ.PermuteNTTWithIndexLvl(level, eval.Pool[1].Q, eval.PermuteNTTIndex[galEl], ctOut.Value[0])
		ringQ.PermuteNTTWithIndexLvl(level, eval.Pool[2].Q, eval.PermuteNTTIndex[galEl], ctOut.Value[1])
	} else {
		ringQ.PermuteLvl(level, eval.Pool[1].Q, galEl, ctOut.Value[0])
		ringQ.PermuteLvl(level, eval.Pool[2].Q, galEl, ctOut.Value[1])
	}

	ctOut.Value[0].Coeffs = ctOut.Value[0].Coeffs[:level+1]
	ctOut.Value[1].Coeffs = ctOut.Value[1].Coeffs[:level+1]
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

	eval.KeyswitchHoisted(level, c1DecompQP, rtk, eval.Pool[0].Q, eval.Pool[1].Q, eval.Pool[0].P, eval.Pool[1].P)
	ringQ.AddLvl(level, eval.Pool[0].Q, ctIn.Value[0], eval.Pool[0].Q)

	if ctIn.Value[0].IsNTT {
		ringQ.PermuteNTTWithIndexLvl(level, eval.Pool[0].Q, eval.PermuteNTTIndex[galEl], ctOut.Value[0])
		ringQ.PermuteNTTWithIndexLvl(level, eval.Pool[1].Q, eval.PermuteNTTIndex[galEl], ctOut.Value[1])
	} else {
		ringQ.PermuteLvl(level, eval.Pool[0].Q, galEl, ctOut.Value[0])
		ringQ.PermuteLvl(level, eval.Pool[1].Q, galEl, ctOut.Value[1])
	}

	ctOut.Value[0].Coeffs = ctOut.Value[0].Coeffs[:level+1]
	ctOut.Value[1].Coeffs = ctOut.Value[1].Coeffs[:level+1]
}

// AutomorphismHoistedNoModDown is similar to AutomorphismHoisted, except that it returns a ciphertext modulo QP and scaled by P.
// The method requires that the corresponding RotationKey has been added to the Evaluator.The method will panic if either ctIn or ctOut degree is not equal to 1.
func (eval *Evaluator) AutomorphismHoistedNoModDown(levelQ int, c0 *ring.Poly, c1DecompQP []ringqp.Poly, galEl uint64, ct0OutQ, ct1OutQ, ct0OutP, ct1OutP *ring.Poly) {

	levelP := eval.params.PCount() - 1

	rtk, generated := eval.Rtks.GetRotationKey(galEl)
	if !generated {
		panic(fmt.Sprintf("galEl key 5^%d missing", eval.params.InverseGaloisElement(galEl)))
	}

	eval.KeyswitchHoistedNoModDown(levelQ, c1DecompQP, rtk, eval.Pool[0].Q, eval.Pool[1].Q, eval.Pool[0].P, eval.Pool[1].P)

	ringQ := eval.params.RingQ()

	if c0.IsNTT {

		index := eval.PermuteNTTIndex[galEl]

		ringQ.PermuteNTTWithIndexLvl(levelQ, eval.Pool[1].Q, index, ct1OutQ)
		ringQ.PermuteNTTWithIndexLvl(levelP, eval.Pool[1].P, index, ct1OutP)

		ringQ.MulScalarBigintLvl(levelQ, c0, eval.params.RingP().ModulusBigint, eval.Pool[1].Q)
		ringQ.AddLvl(levelQ, eval.Pool[0].Q, eval.Pool[1].Q, eval.Pool[0].Q)

		ringQ.PermuteNTTWithIndexLvl(levelQ, eval.Pool[0].Q, index, ct0OutQ)
		ringQ.PermuteNTTWithIndexLvl(levelP, eval.Pool[0].P, index, ct0OutP)
	} else {
		ringQ.PermuteLvl(levelQ, eval.Pool[1].Q, galEl, ct1OutQ)
		ringQ.PermuteLvl(levelP, eval.Pool[1].P, galEl, ct1OutP)

		ringQ.MulScalarBigintLvl(levelQ, c0, eval.params.RingP().ModulusBigint, eval.Pool[1].Q)
		ringQ.AddLvl(levelQ, eval.Pool[0].Q, eval.Pool[1].Q, eval.Pool[0].Q)

		ringQ.PermuteLvl(levelQ, eval.Pool[0].Q, galEl, ct0OutQ)
		ringQ.PermuteLvl(levelP, eval.Pool[0].P, galEl, ct0OutP)
	}
}
