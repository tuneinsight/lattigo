package rlwe

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils"
)

// ApplyEvaluationKey is a generic method to apply an [EvaluationKey] on a ciphertext.
// An EvaluationKey is a type of public key that is be used during the evaluation of
// a homomorphic circuit to provide additional functionalities, like relinearization
// or rotations.
//
// In a nutshell, an [EvaluationKey] encrypts a secret skIn under a secret skOut and
// enables the public and non interactive re-encryption of any ciphertext encrypted
// under skIn to a new ciphertext encrypted under skOut.
//
// The method will return an error if either ctIn or opOut degree isn't 1.
//
// This method can also be used to switch a ciphertext to one with a different ring degree.
// Note that the parameters of the smaller ring degree must be the same or a subset of the
// moduli Q and P of the one for the larger ring degree.
//
// To do so, it must be provided with the appropriate [EvaluationKey], and have the operands
// matching the target ring degrees.
//
// To switch a ciphertext to a smaller ring degree:
//   - ctIn ring degree must match the evaluator's ring degree.
//   - opOut ring degree must match the smaller ring degree.
//   - evk must have been generated using the key-generator of the large ring degree with as input large-key -> small-key.
//
// To switch a ciphertext to a smaller ring degree:
//   - ctIn ring degree must match the smaller ring degree.
//   - opOut ring degree must match the evaluator's ring degree.
//   - evk must have been generated using the key-generator of the large ring degree with as input small-key -> large-key.
func (eval Evaluator) ApplyEvaluationKey(ctIn *Ciphertext, evk *EvaluationKey, opOut *Ciphertext) (err error) {

	if ctIn.Degree() != 1 || opOut.Degree() != 1 {
		return fmt.Errorf("cannot ApplyEvaluationKey: input and output Ciphertext must be of degree 1")
	}

	level := utils.Min(ctIn.Level(), opOut.Level())
	ringQ := eval.params.RingQ().AtLevel(level)

	NIn := ctIn.Value[0].N()
	NOut := opOut.Value[0].N()

	// Re-encryption to a larger ring degree.
	if NIn < NOut {

		if NOut != ringQ.N() {
			return fmt.Errorf("cannot ApplyEvaluationKey: opOut ring degree does not match evaluator params ring degree")
		}

		// Maps to larger ring degree Y = X^{N/n} -> X
		if ctIn.IsNTT {
			SwitchCiphertextRingDegreeNTT(ctIn.El(), nil, opOut.El())
		} else {
			SwitchCiphertextRingDegree(ctIn.El(), opOut.El())
		}

		// Re-encrypt opOut from the key from small to larger ring degree
		eval.applyEvaluationKey(level, opOut, evk, opOut)

		// Re-encryption to a smaller ring degree.
	} else if NIn > NOut {

		if NIn != ringQ.N() {
			return fmt.Errorf("cannot ApplyEvaluationKey: ctIn ring degree does not match evaluator params ring degree")
		}

		level := utils.Min(ctIn.Level(), opOut.Level())

		ctTmp, err := NewCiphertextAtLevelFromPoly(level, eval.BuffCt.Value)

		// Sanity check, this error should not happen unless the
		// evaluator's buffer have been improperly tempered with.
		if err != nil {
			panic(err)
		}

		ctTmp.MetaData = ctIn.MetaData

		// Switches key from large to small degree
		eval.applyEvaluationKey(level, ctIn, evk, ctTmp)

		// Maps to smaller ring degree X -> Y = X^{N/n}
		if ctIn.IsNTT {
			SwitchCiphertextRingDegreeNTT(ctTmp.El(), ringQ, opOut.El())
		} else {
			SwitchCiphertextRingDegree(ctTmp.El(), opOut.El())
		}

		// Re-encryption to the same ring degree.
	} else {
		eval.applyEvaluationKey(level, ctIn, evk, opOut)
	}

	*opOut.MetaData = *ctIn.MetaData

	return
}

func (eval Evaluator) applyEvaluationKey(level int, ctIn *Ciphertext, evk *EvaluationKey, opOut *Ciphertext) {
	ctTmp := &Ciphertext{}
	ctTmp.Value = []ring.Poly{eval.BuffQP[0].Q, eval.BuffQP[1].Q}
	ctTmp.MetaData = ctIn.MetaData
	eval.GadgetProduct(level, ctIn.Value[1], &evk.GadgetCiphertext, ctTmp)
	eval.params.RingQ().AtLevel(level).Add(ctIn.Value[0], ctTmp.Value[0], opOut.Value[0])
	opOut.Value[1].CopyLvl(level, ctTmp.Value[1])
}

// Relinearize applies the relinearization procedure on ct0 and returns the result in opOut.
// Relinearization is a special procedure required to ensure ciphertext compactness.
// It takes as input a quadratic ciphertext, that decrypts with the key (1, sk, sk^2) and
// outputs a linear ciphertext that decrypts with the key (1, sk).
// In a nutshell, the relinearization re-encrypt the term that decrypts using sk^2 to one
// that decrypts using sk.
// The method will return an error if:
//   - The input ciphertext degree isn't 2.
//   - The corresponding relinearization key to the ciphertext degree
//
// is missing.
func (eval Evaluator) Relinearize(ctIn *Ciphertext, opOut *Ciphertext) (err error) {

	if ctIn.Degree() != 2 {
		return fmt.Errorf("cannot relinearize: ctIn.Degree() should be 2 but is %d", ctIn.Degree())
	}

	var rlk *RelinearizationKey
	if rlk, err = eval.CheckAndGetRelinearizationKey(); err != nil {
		return fmt.Errorf("cannot relinearize: %w", err)
	}

	level := utils.Min(ctIn.Level(), opOut.Level())

	ringQ := eval.params.RingQ().AtLevel(level)

	ctTmp := &Ciphertext{}
	ctTmp.Value = []ring.Poly{eval.BuffQP[0].Q, eval.BuffQP[1].Q}
	ctTmp.MetaData = ctIn.MetaData

	eval.GadgetProduct(level, ctIn.Value[2], &rlk.GadgetCiphertext, ctTmp)
	ringQ.Add(ctIn.Value[0], ctTmp.Value[0], opOut.Value[0])
	ringQ.Add(ctIn.Value[1], ctTmp.Value[1], opOut.Value[1])

	opOut.Resize(1, level)

	*opOut.MetaData = *ctIn.MetaData

	return
}
