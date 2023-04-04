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

	ctTmp := &Ciphertext{Value: []*ring.Poly{eval.BuffQP[0].Q, eval.BuffQP[1].Q}}
	ctTmp.IsNTT = ctIn.IsNTT

	eval.GadgetProduct(level, ctIn.Value[1], evk.GadgetCiphertext, ctTmp)

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

	eval.GadgetProductHoisted(level, c1DecompQP, evk.EvaluationKey.GadgetCiphertext, ctTmp)
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
func (eval *Evaluator) AutomorphismHoistedLazy(levelQ int, ctIn *Ciphertext, c1DecompQP []ringqp.Poly, galEl uint64, ctQP CiphertextQP) {

	var evk *GaloisKey
	var err error
	if evk, err = eval.CheckAndGetGaloisKey(galEl); err != nil {
		panic(fmt.Errorf("cannot apply AutomorphismHoistedLazy: %w", err))
	}

	levelP := evk.LevelP()

	ctTmp := CiphertextQP{}
	ctTmp.Value = [2]ringqp.Poly{eval.BuffQP[0], eval.BuffQP[1]}
	ctTmp.IsNTT = ctQP.IsNTT

	eval.GadgetProductHoistedLazy(levelQ, c1DecompQP, evk.GadgetCiphertext, ctTmp)

	ringQP := eval.params.RingQP().AtLevel(levelQ, levelP)

	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	index := eval.AutomorphismIndex[galEl]

	if ctQP.IsNTT {

		ringQP.AutomorphismNTTWithIndex(&ctTmp.Value[1], index, &ctQP.Value[1])

		if levelP > -1 {
			ringQ.MulScalarBigint(ctIn.Value[0], ringP.ModulusAtLevel[levelP], ctTmp.Value[1].Q)
		}

		ringQ.Add(ctTmp.Value[0].Q, ctTmp.Value[1].Q, ctTmp.Value[0].Q)

		ringQP.AutomorphismNTTWithIndex(&ctTmp.Value[0], index, &ctQP.Value[0])
	} else {

		ringQP.Automorphism(&ctTmp.Value[1], galEl, &ctQP.Value[1])

		if levelP > -1 {
			ringQ.MulScalarBigint(ctIn.Value[0], ringP.ModulusAtLevel[levelP], ctTmp.Value[1].Q)
		}

		ringQ.Add(ctTmp.Value[0].Q, ctTmp.Value[1].Q, ctTmp.Value[0].Q)

		ringQP.Automorphism(&ctTmp.Value[0], galEl, &ctQP.Value[0])
	}
}

// Trace maps X -> sum((-1)^i * X^{i*n+1}) for n <= i < N
// Monomial X^k vanishes if k is not divisible by (N/n), otherwise it is multiplied by (N/n).
// Ciphertext is pre-multiplied by (N/n)^-1 to remove the (N/n) factor.
// Examples of full Trace for [0 + 1X + 2X^2 + 3X^3 + 4X^4 + 5X^5 + 6X^6 + 7X^7]
//
// 1.
//
//	  [1 + 2X + 3X^2 + 4X^3 + 5X^4 + 6X^5 + 7X^6 + 8X^7]
//	+ [1 - 6X - 3X^2 + 8X^3 + 5X^4 + 2X^5 - 7X^6 - 4X^7]  {X-> X^(i * 5^1)}
//	= [2 - 4X + 0X^2 +12X^3 +10X^4 + 8X^5 - 0X^6 + 4X^7]
//
// 2.
//
//	  [2 - 4X + 0X^2 +12X^3 +10X^4 + 8X^5 - 0X^6 + 4X^7]
//	+ [2 + 4X + 0X^2 -12X^3 +10X^4 - 8X^5 + 0X^6 - 4X^7]  {X-> X^(i * 5^2)}
//	= [4 + 0X + 0X^2 - 0X^3 +20X^4 + 0X^5 + 0X^6 - 0X^7]
//
// 3.
//
//	  [4 + 0X + 0X^2 - 0X^3 +20X^4 + 0X^5 + 0X^6 - 0X^7]
//	+ [4 + 0X + 0X^2 - 0X^3 -20X^4 + 0X^5 + 0X^6 - 0X^7]  {X-> X^(i * -1)}
//	= [8 + 0X + 0X^2 - 0X^3 + 0X^4 + 0X^5 + 0X^6 - 0X^7]
func (eval *Evaluator) Trace(ctIn *Ciphertext, logN int, ctOut *Ciphertext) {

	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		panic("ctIn.Degree() != 1 or ctOut.Degree() != 1")
	}

	levelQ := utils.Min(ctIn.Level(), ctOut.Level())

	ctOut.Resize(ctOut.Degree(), levelQ)

	ctOut.MetaData = ctIn.MetaData

	gap := 1 << (eval.params.LogN() - logN - 1)

	if logN == 0 {
		gap <<= 1
	}

	if gap > 1 {

		ringQ := eval.params.RingQ().AtLevel(levelQ)

		// pre-multiplication by (N/n)^-1
		for i, s := range ringQ.SubRings[:levelQ+1] {

			NInv := ring.MForm(ring.ModExp(uint64(gap), s.Modulus-2, s.Modulus), s.Modulus, s.BRedConstant)

			s.MulScalarMontgomery(ctIn.Value[0].Coeffs[i], NInv, ctOut.Value[0].Coeffs[i])
			s.MulScalarMontgomery(ctIn.Value[1].Coeffs[i], NInv, ctOut.Value[1].Coeffs[i])
		}

		buff := NewCiphertextAtLevelFromPoly(levelQ, []*ring.Poly{eval.BuffQP[3].Q, eval.BuffQP[4].Q})
		buff.IsNTT = ctIn.IsNTT

		for i := logN; i < eval.params.LogN()-1; i++ {
			eval.Automorphism(ctOut, eval.params.GaloisElementForColumnRotationBy(1<<i), buff)
			ringQ.Add(ctOut.Value[0], buff.Value[0], ctOut.Value[0])
			ringQ.Add(ctOut.Value[1], buff.Value[1], ctOut.Value[1])
		}

		if logN == 0 {
			eval.Automorphism(ctOut, ringQ.NthRoot()-1, buff)
			ringQ.Add(ctOut.Value[0], buff.Value[0], ctOut.Value[0])
			ringQ.Add(ctOut.Value[1], buff.Value[1], ctOut.Value[1])
		}

	} else {
		if ctIn != ctOut {
			ctOut.Copy(ctIn)
		}
	}
}
