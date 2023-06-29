package lut

import (
	"github.com/tuneinsight/lattigo/v4/rgsw"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// EvaluationKey is a struct storing the encryption
// of the bits of the LWE key.
type EvaluationKey struct {
	SkPos []*rgsw.Ciphertext
	SkNeg []*rgsw.Ciphertext
}

func (evk EvaluationKey) Base2Decomposition() int {
	return evk.SkPos[0].Value[0].BaseTwoDecomposition
}

// GenEvaluationKeyNew generates a new LUT evaluation key
func GenEvaluationKeyNew(paramsRLWE rlwe.Parameters, skRLWE *rlwe.SecretKey, paramsLWE rlwe.Parameters, skLWE *rlwe.SecretKey, Base2Decomposition int) (key EvaluationKey) {

	skLWEInvNTT := paramsLWE.RingQ().NewPoly()

	paramsLWE.RingQ().INTT(skLWE.Value.Q, skLWEInvNTT)

	plaintextRGSWOne := rlwe.NewPlaintext(paramsRLWE, paramsRLWE.MaxLevel())
	plaintextRGSWOne.IsNTT = true
	NRLWE := paramsRLWE.N()
	for j := 0; j < paramsRLWE.QCount(); j++ {
		for i := 0; i < NRLWE; i++ {
			plaintextRGSWOne.Value.Coeffs[j][i] = 1
		}
	}

	encryptor := rgsw.NewEncryptor(paramsRLWE, skRLWE)

	levelQ := paramsRLWE.QCount() - 1
	levelP := paramsRLWE.PCount() - 1

	skRGSWPos := make([]*rgsw.Ciphertext, paramsLWE.N())
	skRGSWNeg := make([]*rgsw.Ciphertext, paramsLWE.N())

	ringQ := paramsLWE.RingQ()
	Q := ringQ.SubRings[0].Modulus
	OneMForm := ring.MForm(1, Q, ringQ.SubRings[0].BRedConstant)
	MinusOneMform := ring.MForm(Q-1, Q, ringQ.SubRings[0].BRedConstant)

	for i, si := range skLWEInvNTT.Coeffs[0] {

		skRGSWPos[i] = rgsw.NewCiphertext(paramsRLWE, levelQ, levelP, Base2Decomposition)
		skRGSWNeg[i] = rgsw.NewCiphertext(paramsRLWE, levelQ, levelP, Base2Decomposition)

		// sk_i =  1 -> [RGSW(1), RGSW(0)]
		if si == OneMForm {
			encryptor.Encrypt(plaintextRGSWOne, skRGSWPos[i])
			encryptor.EncryptZero(skRGSWNeg[i])
			// sk_i = -1 -> [RGSW(0), RGSW(1)]
		} else if si == MinusOneMform {
			encryptor.EncryptZero(skRGSWPos[i])
			encryptor.Encrypt(plaintextRGSWOne, skRGSWNeg[i])
			// sk_i =  0 -> [RGSW(0), RGSW(0)]
		} else {
			encryptor.EncryptZero(skRGSWPos[i])
			encryptor.EncryptZero(skRGSWNeg[i])
		}
	}

	return EvaluationKey{SkPos: skRGSWPos, SkNeg: skRGSWNeg}
}
