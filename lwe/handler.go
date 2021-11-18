package lwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

type Handler struct {
	*rlwe.KeySwitcher
	paramsLUT rlwe.Parameters
	paramsLWE rlwe.Parameters
	rtks      *rlwe.RotationKeySet

	nPowInv      [][]uint64
	xPow         []*ring.Poly
	xPowMinusOne []rlwe.PolyQP

	permuteNTTIndex map[uint64][]uint64

	poolMod2N [2]*ring.Poly

	accumulator *rlwe.Ciphertext
}

func NewHandler(paramsLUT, paramsLWE rlwe.Parameters, rtks *rlwe.RotationKeySet) (h *Handler) {
	h = new(Handler)
	h.KeySwitcher = rlwe.NewKeySwitcher(paramsLUT)
	h.paramsLUT = paramsLUT
	h.paramsLWE = paramsLWE

	ringQ := paramsLUT.RingQ()
	ringP := paramsLUT.RingP()

	h.poolMod2N = [2]*ring.Poly{paramsLWE.RingQ().NewPolyLvl(0), paramsLWE.RingQ().NewPolyLvl(0)}
	h.accumulator = rlwe.NewCiphertextNTT(paramsLUT, 1, paramsLUT.MaxLevel())

	h.nPowInv = make([][]uint64, paramsLUT.LogN())
	h.xPow = make([]*ring.Poly, paramsLUT.LogN())

	for i := 0; i < paramsLUT.LogN(); i++ {
		h.nPowInv[i] = make([]uint64, paramsLUT.MaxLevel()+1)
		var nTimesN uint64 = 1 << (paramsLUT.LogN() + i)
		for j := 0; j < paramsLUT.MaxLevel()+1; j++ {
			h.nPowInv[i][j] = ring.MForm(ring.ModExp(nTimesN, ringQ.Modulus[j]-2, ringQ.Modulus[j]), ringQ.Modulus[j], ringQ.BredParams[j])
		}

		h.xPow[i] = ringQ.NewPoly()
		if i == 0 {
			for j := 0; j < paramsLUT.MaxLevel()+1; j++ {
				h.xPow[i].Coeffs[j][1<<i] = ring.MForm(1, ringQ.Modulus[j], ringQ.BredParams[j])
			}
			ringQ.NTT(h.xPow[i], h.xPow[i])
		} else {
			ringQ.MulCoeffsMontgomery(h.xPow[i-1], h.xPow[i-1], h.xPow[i]) // X^{n} = X^{1} * X^{n-1}
		}
	}

	oneNTTMFormQ := ringQ.NewPoly()

	for i := range ringQ.Modulus {
		for j := 0; j < ringQ.N; j++ {
			oneNTTMFormQ.Coeffs[i][j] = ring.MForm(1, ringQ.Modulus[i], ringQ.BredParams[i])
		}
	}

	oneNTTMFormP := ringP.NewPoly()
	for i := range ringP.Modulus {
		for j := 0; j < ringP.N; j++ {
			oneNTTMFormP.Coeffs[i][j] = ring.MForm(1, ringP.Modulus[i], ringP.BredParams[i])
		}
	}

	N := ringQ.N

	h.xPowMinusOne = make([]rlwe.PolyQP, 2*N)
	for i := 0; i < N; i++ {
		h.xPowMinusOne[i].Q = ringQ.NewPoly()
		h.xPowMinusOne[i].P = ringP.NewPoly()
		h.xPowMinusOne[i+N].Q = ringQ.NewPoly()
		h.xPowMinusOne[i+N].P = ringP.NewPoly()
		if i == 0 || i == 1 {
			for j := range ringQ.Modulus {
				h.xPowMinusOne[i].Q.Coeffs[j][i] = ring.MForm(1, ringQ.Modulus[j], ringQ.BredParams[j])
			}

			for j := range ringP.Modulus {
				h.xPowMinusOne[i].P.Coeffs[j][i] = ring.MForm(1, ringP.Modulus[j], ringP.BredParams[j])
			}

			ringQ.NTT(h.xPowMinusOne[i].Q, h.xPowMinusOne[i].Q)
			ringP.NTT(h.xPowMinusOne[i].P, h.xPowMinusOne[i].P)

			ringQ.Neg(h.xPowMinusOne[i].Q, h.xPowMinusOne[i+N].Q)
			ringP.Neg(h.xPowMinusOne[i].P, h.xPowMinusOne[i+N].P)

		} else {
			ringQ.MulCoeffsMontgomery(h.xPowMinusOne[1].Q, h.xPowMinusOne[i-1].Q, h.xPowMinusOne[i].Q) // X^{n} = X^{1} * X^{n-1}
			ringP.MulCoeffsMontgomery(h.xPowMinusOne[1].P, h.xPowMinusOne[i-1].P, h.xPowMinusOne[i].P)
			ringQ.Neg(h.xPowMinusOne[i].Q, h.xPowMinusOne[i+N].Q) // X^{2n} = -X^{1} * X^{n-1}
			ringP.Neg(h.xPowMinusOne[i].P, h.xPowMinusOne[i+N].P)
		}
	}

	for i := 0; i < 2*N; i++ {
		ringQ.Sub(h.xPowMinusOne[i].Q, oneNTTMFormQ, h.xPowMinusOne[i].Q) // X^{n} - 1
		ringP.Sub(h.xPowMinusOne[i].P, oneNTTMFormP, h.xPowMinusOne[i].P)
	}

	h.rtks = rtks
	h.permuteNTTIndex = *h.permuteNTTIndexesForKey(rtks)

	return
}

func (h *Handler) permuteNTTIndexesForKey(rtks *rlwe.RotationKeySet) *map[uint64][]uint64 {
	if rtks == nil {
		return &map[uint64][]uint64{}
	}
	permuteNTTIndex := make(map[uint64][]uint64, len(rtks.Keys))
	for galEl := range rtks.Keys {
		permuteNTTIndex[galEl] = h.paramsLUT.RingQ().PermuteNTTIndex(galEl)
	}
	return &permuteNTTIndex
}

type LUTKey struct {
	SkPos  []*rlwe.RGSWCiphertext
	SkNeg  []*rlwe.RGSWCiphertext
	EncOne *rlwe.RGSWCiphertext
}

func (h *Handler) GenLUTKey(skRLWE, skLWE *rlwe.SecretKey) (lutkey *LUTKey) {

	paramsLUT := h.paramsLUT
	paramsLWE := h.paramsLWE

	skLWEInvNTT := h.paramsLWE.RingQ().NewPoly()

	paramsLWE.RingQ().InvNTT(skLWE.Value.Q, skLWEInvNTT)

	plaintextRGSWOne := rlwe.NewPlaintext(paramsLUT, paramsLUT.MaxLevel())
	plaintextRGSWOne.Value.IsNTT = true
	for j := 0; j < paramsLUT.QCount(); j++{
		for i := 0; i < paramsLUT.N(); i++ {
			plaintextRGSWOne.Value.Coeffs[j][i] = 1
		}
	}

	encryptor := rlwe.NewEncryptor(paramsLUT, skRLWE)

	EncOneRGSW := rlwe.NewCiphertextRGSWNTT(paramsLUT, paramsLUT.MaxLevel())
	encryptor.EncryptRGSW(plaintextRGSWOne, EncOneRGSW)

	skRGSWPos := make([]*rlwe.RGSWCiphertext, paramsLWE.N())
	skRGSWNeg := make([]*rlwe.RGSWCiphertext, paramsLWE.N())

	ringQ := paramsLWE.RingQ()
	Q := ringQ.Modulus[0]
	OneMForm := ring.MForm(1, Q, ringQ.BredParams[0])
	MinusOneMform := ring.MForm(Q-1, Q, ringQ.BredParams[0])

	for i, si := range skLWEInvNTT.Coeffs[0] {

		skRGSWPos[i] = rlwe.NewCiphertextRGSWNTT(paramsLUT, paramsLUT.MaxLevel())
		skRGSWNeg[i] = rlwe.NewCiphertextRGSWNTT(paramsLUT, paramsLUT.MaxLevel())

		if si == OneMForm {
			encryptor.EncryptRGSW(plaintextRGSWOne, skRGSWPos[i])
			encryptor.EncryptRGSW(nil, skRGSWNeg[i])
		} else if si == MinusOneMform {
			encryptor.EncryptRGSW(nil, skRGSWPos[i])
			encryptor.EncryptRGSW(plaintextRGSWOne, skRGSWNeg[i])
		} else {
			encryptor.EncryptRGSW(nil, skRGSWPos[i])
			encryptor.EncryptRGSW(nil, skRGSWNeg[i])
		}
	}

	return &LUTKey{SkPos: skRGSWPos, SkNeg: skRGSWNeg, EncOne: EncOneRGSW}
}

func ReduceRGSW(rgsw *rlwe.RGSWCiphertext, ringQP *rlwe.RingQP, res *rlwe.RGSWCiphertext) {

	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	for i := range rgsw.Value {
		ringQ.Reduce(rgsw.Value[i][0][0].Q, res.Value[i][0][0].Q)
		ringP.Reduce(rgsw.Value[i][0][0].P, res.Value[i][0][0].P)

		ringQ.Reduce(rgsw.Value[i][0][1].Q, res.Value[i][0][1].Q)
		ringP.Reduce(rgsw.Value[i][0][1].P, res.Value[i][0][1].P)

		ringQ.Reduce(rgsw.Value[i][1][0].Q, res.Value[i][1][0].Q)
		ringP.Reduce(rgsw.Value[i][1][0].P, res.Value[i][1][0].P)

		ringQ.Reduce(rgsw.Value[i][1][1].Q, res.Value[i][1][1].Q)
		ringP.Reduce(rgsw.Value[i][1][1].P, res.Value[i][1][1].P)
	}
}

func AddRGSW(rgsw *rlwe.RGSWCiphertext, ringQP *rlwe.RingQP, res *rlwe.RGSWCiphertext) {

	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	for i := range rgsw.Value {
		ringQ.AddNoMod(res.Value[i][0][0].Q, rgsw.Value[i][0][0].Q, res.Value[i][0][0].Q)
		ringP.AddNoMod(res.Value[i][0][0].P, rgsw.Value[i][0][0].P, res.Value[i][0][0].P)

		ringQ.AddNoMod(res.Value[i][0][1].Q, rgsw.Value[i][0][1].Q, res.Value[i][0][1].Q)
		ringP.AddNoMod(res.Value[i][0][1].P, rgsw.Value[i][0][1].P, res.Value[i][0][1].P)

		ringQ.AddNoMod(res.Value[i][1][0].Q, rgsw.Value[i][1][0].Q, res.Value[i][1][0].Q)
		ringP.AddNoMod(res.Value[i][1][0].P, rgsw.Value[i][1][0].P, res.Value[i][1][0].P)

		ringQ.AddNoMod(res.Value[i][1][1].Q, rgsw.Value[i][1][1].Q, res.Value[i][1][1].Q)
		ringP.AddNoMod(res.Value[i][1][1].P, rgsw.Value[i][1][1].P, res.Value[i][1][1].P)
	}
}

func MulRGSWByXPowAlphaMinusOne(rgsw *rlwe.RGSWCiphertext, powXMinusOne rlwe.PolyQP, ringQP *rlwe.RingQP, res *rlwe.RGSWCiphertext) {

	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	for i := range rgsw.Value {
		ringQ.MulCoeffsMontgomeryConstant(rgsw.Value[i][0][0].Q, powXMinusOne.Q, res.Value[i][0][0].Q)
		ringP.MulCoeffsMontgomeryConstant(rgsw.Value[i][0][0].P, powXMinusOne.P, res.Value[i][0][0].P)

		ringQ.MulCoeffsMontgomeryConstant(rgsw.Value[i][0][1].Q, powXMinusOne.Q, res.Value[i][0][1].Q)
		ringP.MulCoeffsMontgomeryConstant(rgsw.Value[i][0][1].P, powXMinusOne.P, res.Value[i][0][1].P)

		ringQ.MulCoeffsMontgomeryConstant(rgsw.Value[i][1][0].Q, powXMinusOne.Q, res.Value[i][1][0].Q)
		ringP.MulCoeffsMontgomeryConstant(rgsw.Value[i][1][0].P, powXMinusOne.P, res.Value[i][1][0].P)

		ringQ.MulCoeffsMontgomeryConstant(rgsw.Value[i][1][1].Q, powXMinusOne.Q, res.Value[i][1][1].Q)
		ringP.MulCoeffsMontgomeryConstant(rgsw.Value[i][1][1].P, powXMinusOne.P, res.Value[i][1][1].P)
	}
}

func MulRGSWByXPowAlphaMinusOneAndAdd(rgsw *rlwe.RGSWCiphertext, powXMinusOne rlwe.PolyQP, ringQP *rlwe.RingQP, res *rlwe.RGSWCiphertext) {

	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	for i := range rgsw.Value {
		ringQ.MulCoeffsMontgomeryConstantAndAddNoMod(rgsw.Value[i][0][0].Q, powXMinusOne.Q, res.Value[i][0][0].Q)
		ringP.MulCoeffsMontgomeryConstantAndAddNoMod(rgsw.Value[i][0][0].P, powXMinusOne.P, res.Value[i][0][0].P)

		ringQ.MulCoeffsMontgomeryConstantAndAddNoMod(rgsw.Value[i][0][1].Q, powXMinusOne.Q, res.Value[i][0][1].Q)
		ringP.MulCoeffsMontgomeryConstantAndAddNoMod(rgsw.Value[i][0][1].P, powXMinusOne.P, res.Value[i][0][1].P)

		ringQ.MulCoeffsMontgomeryConstantAndAddNoMod(rgsw.Value[i][1][0].Q, powXMinusOne.Q, res.Value[i][1][0].Q)
		ringP.MulCoeffsMontgomeryConstantAndAddNoMod(rgsw.Value[i][1][0].P, powXMinusOne.P, res.Value[i][1][0].P)

		ringQ.MulCoeffsMontgomeryConstantAndAddNoMod(rgsw.Value[i][1][1].Q, powXMinusOne.Q, res.Value[i][1][1].Q)
		ringP.MulCoeffsMontgomeryConstantAndAddNoMod(rgsw.Value[i][1][1].P, powXMinusOne.P, res.Value[i][1][1].P)
	}
}