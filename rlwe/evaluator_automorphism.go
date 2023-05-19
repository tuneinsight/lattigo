package rlwe

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
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

	var evk *GaloisKey
	var err error
	if evk, err = eval.CheckAndGetGaloisKey(galEl); err != nil {
		panic(fmt.Errorf("cannot apply Automorphism: %w", err))
	}

	level := utils.Min(ctIn.Level(), ctOut.Level())

	ctOut.Resize(ctOut.Degree(), level)

	ringQ := eval.params.RingQ().AtLevel(level)

	ctTmp := &Ciphertext{OperandQ{Value: []*ring.Poly{eval.BuffQP[0].Q, eval.BuffQP[1].Q}}}
	ctTmp.IsNTT = ctIn.IsNTT

	eval.GadgetProduct(level, ctIn.Value[1], &evk.GadgetCiphertext, ctTmp)

	ringQ.Add(ctTmp.Value[0], ctIn.Value[0], ctTmp.Value[0])

	if ctIn.IsNTT {
		ringQ.AutomorphismNTTWithIndex(ctTmp.Value[0], eval.AutomorphismIndex[galEl], ctOut.Value[0])
		ringQ.AutomorphismNTTWithIndex(ctTmp.Value[1], eval.AutomorphismIndex[galEl], ctOut.Value[1])
	} else {
		ringQ.Automorphism(ctTmp.Value[0], galEl, ctOut.Value[0])
		ringQ.Automorphism(ctTmp.Value[1], galEl, ctOut.Value[1])
	}

	ctOut.MetaData = ctIn.MetaData
}

// AutomorphismHoisted is similar to Automorphism, except that it takes as input ctIn and c1DecompQP, where c1DecompQP is the RNS
// decomposition of its element of degree 1. This decomposition can be obtained with DecomposeNTT.
// The method requires that the corresponding RotationKey has been added to the Evaluator.
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

	var evk *GaloisKey
	var err error
	if evk, err = eval.CheckAndGetGaloisKey(galEl); err != nil {
		panic(fmt.Errorf("cannot apply AutomorphismHoisted: %w", err))
	}

	ctOut.Resize(ctOut.Degree(), level)

	ringQ := eval.params.RingQ().AtLevel(level)

	ctTmp := &Ciphertext{}
	ctTmp.Value = []*ring.Poly{eval.BuffQP[0].Q, eval.BuffQP[1].Q} // GadgetProductHoisted uses the same buffers for its ciphertext QP
	ctTmp.IsNTT = ctIn.IsNTT

	eval.GadgetProductHoisted(level, c1DecompQP, &evk.EvaluationKey.GadgetCiphertext, ctTmp)
	ringQ.Add(ctTmp.Value[0], ctIn.Value[0], ctTmp.Value[0])

	if ctIn.IsNTT {
		ringQ.AutomorphismNTTWithIndex(ctTmp.Value[0], eval.AutomorphismIndex[galEl], ctOut.Value[0])
		ringQ.AutomorphismNTTWithIndex(ctTmp.Value[1], eval.AutomorphismIndex[galEl], ctOut.Value[1])
	} else {
		ringQ.Automorphism(ctTmp.Value[0], galEl, ctOut.Value[0])
		ringQ.Automorphism(ctTmp.Value[1], galEl, ctOut.Value[1])
	}

	ctOut.MetaData = ctIn.MetaData
}

// AutomorphismHoistedLazy is similar to AutomorphismHoisted, except that it returns a ciphertext modulo QP and scaled by P.
// The method requires that the corresponding RotationKey has been added to the Evaluator.
// Result NTT domain is returned according to the NTT flag of ctQP.
func (eval *Evaluator) AutomorphismHoistedLazy(levelQ int, ctIn *Ciphertext, c1DecompQP []ringqp.Poly, galEl uint64, ctQP *OperandQP) {

	var evk *GaloisKey
	var err error
	if evk, err = eval.CheckAndGetGaloisKey(galEl); err != nil {
		panic(fmt.Errorf("cannot apply AutomorphismHoistedLazy: %w", err))
	}

	levelP := evk.LevelP()

	ctTmp := &OperandQP{}
	ctTmp.Value = []*ringqp.Poly{&eval.BuffQP[0], &eval.BuffQP[1]}
	ctTmp.IsNTT = ctQP.IsNTT

	eval.GadgetProductHoistedLazy(levelQ, c1DecompQP, &evk.GadgetCiphertext, ctTmp)

	ringQP := eval.params.RingQP().AtLevel(levelQ, levelP)

	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	index := eval.AutomorphismIndex[galEl]

	if ctQP.IsNTT {

		ringQP.AutomorphismNTTWithIndex(ctTmp.Value[1], index, ctQP.Value[1])

		if levelP > -1 {
			ringQ.MulScalarBigint(ctIn.Value[0], ringP.ModulusAtLevel[levelP], ctTmp.Value[1].Q)
		}

		ringQ.Add(ctTmp.Value[0].Q, ctTmp.Value[1].Q, ctTmp.Value[0].Q)

		ringQP.AutomorphismNTTWithIndex(ctTmp.Value[0], index, ctQP.Value[0])
	} else {

		ringQP.Automorphism(ctTmp.Value[1], galEl, ctQP.Value[1])

		if levelP > -1 {
			ringQ.MulScalarBigint(ctIn.Value[0], ringP.ModulusAtLevel[levelP], ctTmp.Value[1].Q)
		}

		ringQ.Add(ctTmp.Value[0].Q, ctTmp.Value[1].Q, ctTmp.Value[0].Q)

		ringQP.Automorphism(ctTmp.Value[0], galEl, ctQP.Value[0])
	}
}
