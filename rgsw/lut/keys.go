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

// GenEvaluationKey generates the LUT evaluation key
func GenEvaluationKey(paramsRLWE rlwe.Parameters, skRLWE *rlwe.SecretKey, paramsLWE rlwe.Parameters, skLWE *rlwe.SecretKey) (key EvaluationKey) {

	skLWEInvNTT := paramsLWE.RingQ().NewPoly()

	paramsLWE.RingQ().InvNTT(skLWE.Value.Q, skLWEInvNTT)

	plaintextRGSWOne := rlwe.NewPlaintext(paramsRLWE, paramsRLWE.MaxLevel())
	plaintextRGSWOne.Value.IsNTT = true
	for j := 0; j < paramsRLWE.QCount(); j++ {
		for i := 0; i < paramsRLWE.N(); i++ {
			plaintextRGSWOne.Value.Coeffs[j][i] = 1
		}
	}

	encryptor := rgsw.NewEncryptor(paramsRLWE, skRLWE)

	levelQ := paramsRLWE.QCount() - 1
	levelP := paramsRLWE.PCount() - 1

	skRGSWPos := make([]*rgsw.Ciphertext, paramsLWE.N())
	skRGSWNeg := make([]*rgsw.Ciphertext, paramsLWE.N())

	ringQ := paramsLWE.RingQ()
	Q := ringQ.Modulus[0]
	OneMForm := ring.MForm(1, Q, ringQ.BredParams[0])
	MinusOneMform := ring.MForm(Q-1, Q, ringQ.BredParams[0])

	decompRNS := paramsRLWE.DecompRNS(levelQ, levelP)
	decompPw2 := paramsRLWE.DecompPw2(levelQ, levelP)
	ringQP := *paramsRLWE.RingQP()

	for i, si := range skLWEInvNTT.Coeffs[0] {

		skRGSWPos[i] = rgsw.NewCiphertext(levelQ, levelP, decompRNS, decompPw2, ringQP)
		skRGSWNeg[i] = rgsw.NewCiphertext(levelQ, levelP, decompRNS, decompPw2, ringQP)

		// sk_i =  1 -> [RGSW(1), RGSW(0)]
		if si == OneMForm {
			encryptor.Encrypt(plaintextRGSWOne, skRGSWPos[i])
			encryptor.Encrypt(nil, skRGSWNeg[i])
			// sk_i = -1 -> [RGSW(0), RGSW(1)]
		} else if si == MinusOneMform {
			encryptor.Encrypt(nil, skRGSWPos[i])
			encryptor.Encrypt(plaintextRGSWOne, skRGSWNeg[i])
			// sk_i =  0 -> [RGSW(0), RGSW(0)]
		} else {
			encryptor.Encrypt(nil, skRGSWPos[i])
			encryptor.Encrypt(nil, skRGSWNeg[i])
		}
	}

	return EvaluationKey{SkPos: skRGSWPos, SkNeg: skRGSWNeg}
}
