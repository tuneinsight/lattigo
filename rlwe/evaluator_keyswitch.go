package rlwe

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/utils"
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
	ringQ := eval.params.RingQ()

	eval.GadgetProduct(level, ctIn.Value[1], switchingKey.GadgetCiphertext, Ciphertext{Value: []*ring.Poly{eval.BuffQP[1].Q, eval.BuffQP[2].Q}})

	ringQ.AddLvl(level, ctIn.Value[0], eval.BuffQP[1].Q, ctOut.Value[0])
	ring.CopyValuesLvl(level, eval.BuffQP[2].Q, ctOut.Value[1])

	if ctIn.Scale != nil{
		ctOut.Scale = ctIn.Scale.CopyNew()
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

	eval.GadgetProduct(level, ctIn.Value[2], eval.Rlk.Keys[0].GadgetCiphertext, Ciphertext{Value: []*ring.Poly{eval.BuffQP[1].Q, eval.BuffQP[2].Q}})
	ringQ.AddLvl(level, ctIn.Value[0], eval.BuffQP[1].Q, ctOut.Value[0])
	ringQ.AddLvl(level, ctIn.Value[1], eval.BuffQP[2].Q, ctOut.Value[1])

	for deg := ctIn.Degree() - 1; deg > 1; deg-- {
		eval.GadgetProduct(level, ctIn.Value[deg], eval.Rlk.Keys[deg-2].GadgetCiphertext, Ciphertext{Value: []*ring.Poly{eval.BuffQP[1].Q, eval.BuffQP[2].Q}})
		ringQ.AddLvl(level, ctOut.Value[0], eval.BuffQP[1].Q, ctOut.Value[0])
		ringQ.AddLvl(level, ctOut.Value[1], eval.BuffQP[2].Q, ctOut.Value[1])
	}

	ctOut.Value = ctOut.Value[:2]

	ctOut.Resize(ctOut.Degree(), level)

	if ctIn.Scale != nil{
		ctOut.Scale = ctIn.Scale.CopyNew()
	}
}
