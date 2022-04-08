package lut

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/rlwe/rgsw"
)

// Key is a struct storing the encryption
// of the bits of the LWE key.
type Key struct {
	SkPos []*rgsw.Ciphertext
	SkNeg []*rgsw.Ciphertext
	One   *rgsw.Plaintext
}

// GenKey generates the LUT evaluation key
func GenKey(paramsRLWE rlwe.Parameters, skRLWE *rlwe.SecretKey, paramsLWE rlwe.Parameters, skLWE *rlwe.SecretKey) (key Key) {

	skLWEInvNTT := paramsLWE.RingQ().NewPoly()

	paramsLWE.RingQ().InvNTT(skLWE.Value.Q, skLWEInvNTT)

	plaintextRGSWOne := rlwe.NewPlaintext(paramsRLWE, paramsRLWE.MaxLevel())
	plaintextRGSWOne.Value.IsNTT = true
	for j := 0; j < paramsRLWE.QCount(); j++ {
		for i := 0; i < paramsRLWE.N(); i++ {
			plaintextRGSWOne.Value.Coeffs[j][i] = 1
		}
	}

	encryptor := rlwe.NewEncryptor(paramsRLWE, skRLWE)

	levelQ := paramsRLWE.QCount() - 1
	levelP := paramsRLWE.PCount() - 1

	skRGSWPos := make([]*rgsw.Ciphertext, paramsLWE.N())
	skRGSWNeg := make([]*rgsw.Ciphertext, paramsLWE.N())

	ringQ := paramsLWE.RingQ()
	Q := ringQ.Modulus[0]
	OneMForm := ring.MForm(1, Q, ringQ.BredParams[0])
	MinusOneMform := ring.MForm(Q-1, Q, ringQ.BredParams[0])

	decompRNS := paramsRLWE.DecompRNS(levelQ, levelP)
	decompBIT := paramsRLWE.DecompBIT(levelQ, levelP)
	ringQP := *paramsRLWE.RingQP()

	for i, si := range skLWEInvNTT.Coeffs[0] {

		skRGSWPos[i] = rgsw.NewCiphertextNTT(levelQ, levelP, decompRNS, decompBIT, ringQP)
		skRGSWNeg[i] = rgsw.NewCiphertextNTT(levelQ, levelP, decompRNS, decompBIT, ringQP)

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

	return Key{SkPos: skRGSWPos, SkNeg: skRGSWNeg, One: rgsw.NewPlaintext(uint64(1), levelQ, levelP, paramsRLWE.LogBase2(), decompBIT, ringQP)}
}
