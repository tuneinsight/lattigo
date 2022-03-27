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

// GenLUTKey generates the LUT evaluation key
func (eval *Evaluator) GenLUTKey(skRLWE, skLWE *rlwe.SecretKey) (key Key) {

	paramsLUT := eval.paramsLUT
	paramsLWE := eval.paramsLWE

	skLWEInvNTT := eval.paramsLWE.RingQ().NewPoly()

	paramsLWE.RingQ().InvNTT(skLWE.Value.Q, skLWEInvNTT)

	plaintextRGSWOne := rlwe.NewPlaintext(paramsLUT, paramsLUT.MaxLevel())
	plaintextRGSWOne.Value.IsNTT = true
	for j := 0; j < paramsLUT.QCount(); j++ {
		for i := 0; i < paramsLUT.N(); i++ {
			plaintextRGSWOne.Value.Coeffs[j][i] = 1
		}
	}

	encryptor := rlwe.NewEncryptor(paramsLUT, skRLWE)

	levelQ := paramsLUT.QCount() - 1
	levelP := paramsLUT.PCount() - 1

	skRGSWPos := make([]*rgsw.Ciphertext, paramsLWE.N())
	skRGSWNeg := make([]*rgsw.Ciphertext, paramsLWE.N())

	ringQ := paramsLWE.RingQ()
	Q := ringQ.Modulus[0]
	OneMForm := ring.MForm(1, Q, ringQ.BredParams[0])
	MinusOneMform := ring.MForm(Q-1, Q, ringQ.BredParams[0])

	decompRNS := paramsLUT.DecompRNS(levelQ, levelP)
	decompBIT := paramsLUT.DecompBIT(levelQ, levelP)
	ringQP := paramsLUT.RingQP()

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

	return Key{SkPos: skRGSWPos, SkNeg: skRGSWNeg, One: rgsw.NewPlaintext(uint64(1), levelQ, levelP, paramsLUT.LogBase2(), decompBIT, *ringQP)}
}
