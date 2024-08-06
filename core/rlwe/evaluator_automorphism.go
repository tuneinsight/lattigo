package rlwe

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
	"github.com/tuneinsight/lattigo/v6/utils"
)

// Automorphism computes phi(ct), where phi is the map X -> X^galEl. The method requires
// that the corresponding RotationKey has been added to the [Evaluator]. The method will
// return an error if either ctIn or opOut degree is not equal to 1.
func (eval Evaluator) Automorphism(ctIn *Ciphertext, galEl uint64, opOut *Ciphertext) (err error) {

	if ctIn.Degree() != 1 || opOut.Degree() != 1 {
		return fmt.Errorf("cannot apply Automorphism: input and output Ciphertext must be of degree 1")
	}

	if galEl == 1 {
		if opOut != ctIn {
			opOut.Copy(ctIn)
		}
		return
	}

	var evk *GaloisKey
	if evk, err = eval.CheckAndGetGaloisKey(galEl); err != nil {
		return fmt.Errorf("cannot apply Automorphism: %w", err)
	}

	level := utils.Min(ctIn.Level(), opOut.Level())

	opOut.Resize(opOut.Degree(), level)

	ringQ := eval.params.RingQ().AtLevel(level)

	ctTmp := &Ciphertext{Element: Element[ring.Poly]{Value: []ring.Poly{eval.BuffQP[0].Q, eval.BuffQP[1].Q}}}
	ctTmp.MetaData = ctIn.MetaData

	eval.GadgetProduct(level, ctIn.Value[1], &evk.GadgetCiphertext, ctTmp)

	ringQ.Add(ctTmp.Value[0], ctIn.Value[0], ctTmp.Value[0])

	if ctIn.IsNTT {
		ringQ.AutomorphismNTTWithIndex(ctTmp.Value[0], eval.automorphismIndex[galEl], opOut.Value[0])
		ringQ.AutomorphismNTTWithIndex(ctTmp.Value[1], eval.automorphismIndex[galEl], opOut.Value[1])
	} else {
		ringQ.Automorphism(ctTmp.Value[0], galEl, opOut.Value[0])
		ringQ.Automorphism(ctTmp.Value[1], galEl, opOut.Value[1])
	}

	*opOut.MetaData = *ctIn.MetaData

	return
}

// AutomorphismHoisted is similar to [Evaluator.Automorphism], except that it takes as input ctIn and c1DecompQP, where c1DecompQP is the RNS
// decomposition of its element of degree 1. This decomposition can be obtained with [Evaluator.DecomposeNTT].
// The method requires that the corresponding RotationKey has been added to the [Evaluator].
// The method will return an error if either ctIn or opOut degree is not equal to 1.
func (eval Evaluator) AutomorphismHoisted(level int, ctIn *Ciphertext, c1DecompQP []ringqp.Poly, galEl uint64, opOut *Ciphertext) (err error) {

	if ctIn.Degree() != 1 || opOut.Degree() != 1 {
		return fmt.Errorf("cannot apply AutomorphismHoisted: input and output Ciphertext must be of degree 1")
	}

	if galEl == 1 {
		if ctIn != opOut {
			opOut.Copy(ctIn)
		}
		return
	}

	var evk *GaloisKey
	if evk, err = eval.CheckAndGetGaloisKey(galEl); err != nil {
		return fmt.Errorf("cannot apply AutomorphismHoisted: %w", err)
	}

	opOut.Resize(opOut.Degree(), level)

	ringQ := eval.params.RingQ().AtLevel(level)

	ctTmp := &Ciphertext{}
	ctTmp.Value = []ring.Poly{eval.BuffQP[0].Q, eval.BuffQP[1].Q} // GadgetProductHoisted uses the same buffers for its ciphertext QP
	ctTmp.MetaData = ctIn.MetaData

	eval.GadgetProductHoisted(level, c1DecompQP, &evk.EvaluationKey.GadgetCiphertext, ctTmp)
	ringQ.Add(ctTmp.Value[0], ctIn.Value[0], ctTmp.Value[0])

	if ctIn.IsNTT {
		ringQ.AutomorphismNTTWithIndex(ctTmp.Value[0], eval.automorphismIndex[galEl], opOut.Value[0])
		ringQ.AutomorphismNTTWithIndex(ctTmp.Value[1], eval.automorphismIndex[galEl], opOut.Value[1])
	} else {
		ringQ.Automorphism(ctTmp.Value[0], galEl, opOut.Value[0])
		ringQ.Automorphism(ctTmp.Value[1], galEl, opOut.Value[1])
	}

	*opOut.MetaData = *ctIn.MetaData

	return
}

// AutomorphismHoistedLazy is similar to [Evaluator.AutomorphismHoisted], except that it returns a ciphertext modulo QP and scaled by P.
// The method requires that the corresponding RotationKey has been added to the [Evaluator].
// Result NTT domain is returned according to the NTT flag of ctQP.
func (eval Evaluator) AutomorphismHoistedLazy(levelQ int, ctIn *Ciphertext, c1DecompQP []ringqp.Poly, galEl uint64, ctQP *Element[ringqp.Poly]) (err error) {

	var evk *GaloisKey
	if evk, err = eval.CheckAndGetGaloisKey(galEl); err != nil {
		return fmt.Errorf("cannot apply AutomorphismHoistedLazy: %w", err)
	}

	levelP := evk.LevelP()

	if ctQP.LevelP() < levelP {
		return fmt.Errorf("ctQP.LevelP()=%d < GaloisKey[%d].LevelP()=%d", ctQP.LevelP(), galEl, levelP)
	}

	ctTmp := &Element[ringqp.Poly]{}
	ctTmp.Value = []ringqp.Poly{eval.BuffQP[0], eval.BuffQP[1]}
	ctTmp.MetaData = ctIn.MetaData

	if err = eval.GadgetProductHoistedLazy(levelQ, c1DecompQP, &evk.GadgetCiphertext, ctTmp); err != nil {
		panic(fmt.Errorf("eval.GadgetProductHoistedLazy: %w", err))
	}

	ringQP := eval.params.RingQP().AtLevel(levelQ, levelP)

	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	index := eval.automorphismIndex[galEl]

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

	return
}
