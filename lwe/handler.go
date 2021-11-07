package lwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
	"math/big"
	"math/bits"
)

type Handler struct {
	*rlwe.KeySwitcher
	paramsLUT rlwe.Parameters
	paramsLWE  rlwe.Parameters
	rtks       *rlwe.RotationKeySet

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
	for i := 0; i < paramsLUT.N(); i++ {
		plaintextRGSWOne.Value.Coeffs[0][i] = 1
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

// ModSwitchRLWETo2N applys round(x * 2N / Q) to the coefficients of polQ and returns the
// result on pol2N.
func (h *Handler) ModSwitchRLWETo2N(polQ *ring.Poly, pol2N *ring.Poly) {
	coeffsBigint := make([]*big.Int, len(polQ.Coeffs[0]))

	level := polQ.Level()

	ringQ := h.paramsLWE.RingQ()

	ringQ.PolyToBigintLvl(polQ.Level(), polQ, coeffsBigint)

	QBig := ring.NewUint(1)
	for i := 0; i < level+1; i++ {
		QBig.Mul(QBig, ring.NewUint(ringQ.Modulus[i]))
	}

	twoN := uint64(h.paramsLUT.N() << 1)
	twoNBig := ring.NewUint(twoN)
	tmp := pol2N.Coeffs[0]
	for i := 0; i < ringQ.N; i++ {
		coeffsBigint[i].Mul(coeffsBigint[i], twoNBig)
		ring.DivRound(coeffsBigint[i], QBig, coeffsBigint[i])
		tmp[i] = coeffsBigint[i].Uint64() & (twoN - 1)
	}
}

func (h *Handler) ExtractAndEvaluateLUT(ct *rlwe.Ciphertext, lut *ring.Poly, lutKey *LUTKey) (lwe []*rlwe.Ciphertext) {

	ks := h.KeySwitcher

	bRLWEMod2N := h.poolMod2N[0]
	aRLWEMod2N := h.poolMod2N[1]

	acc := h.accumulator

	ringQLUT := h.paramsLUT.RingQ()
	ringQLWE := h.paramsLWE.RingQ()
	ringQPLUT := h.paramsLUT.RingQP()

	mask := uint64(ringQLUT.N<<1) - 1

	ringQLWE.InvNTT(ct.Value[0], acc.Value[0])
	ringQLWE.InvNTT(ct.Value[1], acc.Value[1])

	h.ModSwitchRLWETo2N(acc.Value[1], acc.Value[1])

	// Copy coefficients multiplied by X^{N-1} in reverse order:
	// a_{0} -a_{N-1} -a2_{N-2} ... -a_{1}
	tmp0 := aRLWEMod2N.Coeffs[0]
	tmp1 := acc.Value[1].Coeffs[0]
	tmp0[0] = tmp1[0]
	for j := 1; j < ringQLWE.N; j++ {
		tmp0[j] = -tmp1[ringQLWE.N-j] & mask
	}

	h.ModSwitchRLWETo2N(acc.Value[0], bRLWEMod2N)

	lwe = make([]*rlwe.Ciphertext, ringQLWE.N)

	for i := 0; i < 32; i++ {

		a := aRLWEMod2N.Coeffs[0]
		b := bRLWEMod2N.Coeffs[0][i]

		ringQLUT.MulCoeffsMontgomery(lut, h.xPowMinusOne[b].Q, acc.Value[0])
		ringQLUT.Add(acc.Value[0], lut, acc.Value[0])
		acc.Value[1].Zero() // TODO remove

		tmpRGSW := rlwe.NewCiphertextRGSWNTT(h.paramsLUT, h.paramsLUT.MaxLevel())

		levelQ, levelP := acc.Level(), len(tmpRGSW.Value[0][0][0].P.Coeffs)-1

		for j := 0; j < ringQLWE.N; j++ {
			MulRGSWByXPowAlphaMinusOne(lutKey.SkPos[j], h.xPowMinusOne[a[j]], ringQPLUT, tmpRGSW)
			MulRGSWByXPowAlphaMinusOneAndAdd(lutKey.SkNeg[j], h.xPowMinusOne[-a[j]&mask], ringQPLUT, tmpRGSW)
			AddRGSW(lutKey.EncOne, ringQPLUT, tmpRGSW) // TODO : add 1 in plaintext

			
			ks.MulRGSWNoModDown(levelQ, levelP, acc, tmpRGSW, ks.Pool[1].Q, ks.Pool[1].P, ks.Pool[2].Q, ks.Pool[2].P)

			ks.Baseconverter.ModDownQPtoQNTT(levelQ, levelP, ks.Pool[1].Q, ks.Pool[1].P, acc.Value[0])
			ks.Baseconverter.ModDownQPtoQNTT(levelQ, levelP, ks.Pool[2].Q, ks.Pool[2].P, acc.Value[1])

		}

		// Extracts the first coefficient
		lwe[i] = acc.CopyNew() //ExtractLWEFromRLWESingle(acc, ringQRLWE)

		MulBySmallMonomialMod2N(mask, aRLWEMod2N, 1)
	}

	return lwe
}

func AddRGSW(rgsw *rlwe.RGSWCiphertext, ringQP *rlwe.RingQP, res *rlwe.RGSWCiphertext) {

	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	for i := range rgsw.Value {
		ringQ.Add(res.Value[i][0][0].Q, rgsw.Value[i][0][0].Q, res.Value[i][0][0].Q)
		ringP.Add(res.Value[i][0][0].P, rgsw.Value[i][0][0].P, res.Value[i][0][0].P)

		ringQ.Add(res.Value[i][0][1].Q, rgsw.Value[i][0][1].Q, res.Value[i][0][1].Q)
		ringP.Add(res.Value[i][0][1].P, rgsw.Value[i][0][1].P, res.Value[i][0][1].P)

		ringQ.Add(res.Value[i][1][0].Q, rgsw.Value[i][1][0].Q, res.Value[i][1][0].Q)
		ringP.Add(res.Value[i][1][0].P, rgsw.Value[i][1][0].P, res.Value[i][1][0].P)

		ringQ.Add(res.Value[i][1][1].Q, rgsw.Value[i][1][1].Q, res.Value[i][1][1].Q)
		ringP.Add(res.Value[i][1][1].P, rgsw.Value[i][1][1].P, res.Value[i][1][1].P)
	}
}

func MulRGSWByXPowAlphaMinusOne(rgsw *rlwe.RGSWCiphertext, powXMinusOne rlwe.PolyQP, ringQP *rlwe.RingQP, res *rlwe.RGSWCiphertext) {

	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	for i := range rgsw.Value {
		ringQ.MulCoeffsMontgomery(rgsw.Value[i][0][0].Q, powXMinusOne.Q, res.Value[i][0][0].Q)
		ringP.MulCoeffsMontgomery(rgsw.Value[i][0][0].P, powXMinusOne.P, res.Value[i][0][0].P)

		ringQ.MulCoeffsMontgomery(rgsw.Value[i][0][1].Q, powXMinusOne.Q, res.Value[i][0][1].Q)
		ringP.MulCoeffsMontgomery(rgsw.Value[i][0][1].P, powXMinusOne.P, res.Value[i][0][1].P)

		ringQ.MulCoeffsMontgomery(rgsw.Value[i][1][0].Q, powXMinusOne.Q, res.Value[i][1][0].Q)
		ringP.MulCoeffsMontgomery(rgsw.Value[i][1][0].P, powXMinusOne.P, res.Value[i][1][0].P)

		ringQ.MulCoeffsMontgomery(rgsw.Value[i][1][1].Q, powXMinusOne.Q, res.Value[i][1][1].Q)
		ringP.MulCoeffsMontgomery(rgsw.Value[i][1][1].P, powXMinusOne.P, res.Value[i][1][1].P)
	}
}

func MulRGSWByXPowAlphaMinusOneAndAdd(rgsw *rlwe.RGSWCiphertext, powXMinusOne rlwe.PolyQP, ringQP *rlwe.RingQP, res *rlwe.RGSWCiphertext) {

	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	for i := range rgsw.Value {
		ringQ.MulCoeffsMontgomeryAndAdd(rgsw.Value[i][0][0].Q, powXMinusOne.Q, res.Value[i][0][0].Q)
		ringP.MulCoeffsMontgomeryAndAdd(rgsw.Value[i][0][0].P, powXMinusOne.P, res.Value[i][0][0].P)

		ringQ.MulCoeffsMontgomeryAndAdd(rgsw.Value[i][0][1].Q, powXMinusOne.Q, res.Value[i][0][1].Q)
		ringP.MulCoeffsMontgomeryAndAdd(rgsw.Value[i][0][1].P, powXMinusOne.P, res.Value[i][0][1].P)

		ringQ.MulCoeffsMontgomeryAndAdd(rgsw.Value[i][1][0].Q, powXMinusOne.Q, res.Value[i][1][0].Q)
		ringP.MulCoeffsMontgomeryAndAdd(rgsw.Value[i][1][0].P, powXMinusOne.P, res.Value[i][1][0].P)

		ringQ.MulCoeffsMontgomeryAndAdd(rgsw.Value[i][1][1].Q, powXMinusOne.Q, res.Value[i][1][1].Q)
		ringP.MulCoeffsMontgomeryAndAdd(rgsw.Value[i][1][1].P, powXMinusOne.P, res.Value[i][1][1].P)
	}
}

//MulBySmallMonomial multiplies pol by x^n
func MulBySmallMonomialMod2N(mask uint64, pol *ring.Poly, n int) {
	N := len(pol.Coeffs[0])
	pol.Coeffs[0] = append(pol.Coeffs[0][N-n:], pol.Coeffs[0][:N-n]...)
	tmp := pol.Coeffs[0]
	for j := 0; j < n; j++ {
		tmp[j] = -tmp[j] & mask
	}

}

func (h *Handler) MergeRLWE(ciphertexts []*rlwe.Ciphertext) (ciphertext *rlwe.Ciphertext) {

	slots := len(ciphertexts)

	if slots&(slots-1) != 0 {
		panic("len(ciphertext) must be a power of two smaller or equal to the ring degree")
	}

	logSlots := bits.Len64(uint64(len(ciphertexts))) - 1

	level := ciphertexts[0].Level()

	ringQ := h.paramsLUT.RingQ()

	nPowInv := h.nPowInv[h.paramsLUT.LogN()-logSlots]
	Q := ringQ.Modulus
	mredParams := ringQ.MredParams

	// Multiplies by (Slots * N) ^-1 mod Q
	for i := range ciphertexts {
		v0, v1 := ciphertexts[i].Value[0], ciphertexts[i].Value[1]
		for j := 0; j < ciphertexts[0].Level()+1; j++ {
			ring.MulScalarMontgomeryVec(v0.Coeffs[j], v0.Coeffs[j], nPowInv[j], Q[j], mredParams[j])
			ring.MulScalarMontgomeryVec(v1.Coeffs[j], v1.Coeffs[j], nPowInv[j], Q[j], mredParams[j])
		}
	}

	// Padds for the repacking algorithm
	if slots != h.paramsLUT.N() {
		ciphertexts = append(ciphertexts, make([]*rlwe.Ciphertext, h.paramsLUT.N()-len(ciphertexts))...)
		N := ringQ.N
		gap := N / slots
		for i := 0; i < slots; i++ {
			ciphertexts[N-(i+1)*gap], ciphertexts[slots-i-1] = ciphertexts[slots-i-1], ciphertexts[N-(i+1)*gap]
		}
	}

	ciphertext = h.mergeRLWERecurse(ciphertexts)

	tmp := rlwe.NewCiphertextNTT(h.paramsLUT, 1, ciphertext.Level())
	for i := logSlots - 1; i < h.paramsLUT.LogN()-1; i++ {
		Rotate(ciphertext, h.paramsLUT.GaloisElementForColumnRotationBy(1<<i), h.permuteNTTIndex, h.paramsLUT, h.KeySwitcher, h.rtks, tmp)
		ringQ.AddLvl(level, ciphertext.Value[0], tmp.Value[0], ciphertext.Value[0])
		ringQ.AddLvl(level, ciphertext.Value[1], tmp.Value[1], ciphertext.Value[1])
	}

	return
}

// PackLWEs repacks LWE ciphertexts into a RLWE ciphertext
func (h *Handler) mergeRLWERecurse(ciphertexts []*rlwe.Ciphertext) *rlwe.Ciphertext {

	ringQ := h.paramsLUT.RingQ()

	L := bits.Len64(uint64(len(ciphertexts))) - 1

	if L == 0 {
		return ciphertexts[0]
	}

	odd := make([]*rlwe.Ciphertext, len(ciphertexts)>>1)
	even := make([]*rlwe.Ciphertext, len(ciphertexts)>>1)

	for i := 0; i < len(ciphertexts)>>1; i++ {
		odd[i] = ciphertexts[2*i]
		even[i] = ciphertexts[2*i+1]
	}

	ctEven := h.mergeRLWERecurse(odd)
	ctOdd := h.mergeRLWERecurse(even)

	if ctEven == nil && ctOdd == nil {
		return nil
	}

	var tmpEven *rlwe.Ciphertext
	if ctEven != nil {
		tmpEven = ctEven.CopyNew()
	}

	// ctOdd * X^(N/2^L)
	if ctOdd != nil {

		level := ctOdd.Level()

		//X^(N/2^L)
		ringQ.MulCoeffsMontgomeryLvl(level, ctOdd.Value[0], h.xPow[len(h.xPow)-L], ctOdd.Value[0])
		ringQ.MulCoeffsMontgomeryLvl(level, ctOdd.Value[1], h.xPow[len(h.xPow)-L], ctOdd.Value[1])

		// ctEven + ctOdd * X^(N/2^L)
		ringQ.AddLvl(level, ctEven.Value[0], ctOdd.Value[0], ctEven.Value[0])
		ringQ.AddLvl(level, ctEven.Value[1], ctOdd.Value[1], ctEven.Value[1])

		// phi(ctEven - ctOdd * X^(N/2^L), 2^(L-2))
		ringQ.SubLvl(level, tmpEven.Value[0], ctOdd.Value[0], tmpEven.Value[0])
		ringQ.SubLvl(level, tmpEven.Value[1], ctOdd.Value[1], tmpEven.Value[1])
	}

	if ctEven != nil {

		level := ctEven.Level()

		// if L-2 == -1, then gal = 2N-1
		if L == 1 {
			Rotate(tmpEven, uint64(2*ringQ.N-1), h.permuteNTTIndex, h.paramsLUT, h.KeySwitcher, h.rtks, tmpEven)
		} else {
			Rotate(tmpEven, h.paramsLUT.GaloisElementForColumnRotationBy(1<<(L-2)), h.permuteNTTIndex, h.paramsLUT, h.KeySwitcher, h.rtks, tmpEven)
		}

		// ctEven + ctOdd * X^(N/2^L) + phi(ctEven - ctOdd * X^(N/2^L), 2^(L-2))
		ringQ.AddLvl(level, ctEven.Value[0], tmpEven.Value[0], ctEven.Value[0])
		ringQ.AddLvl(level, ctEven.Value[1], tmpEven.Value[1], ctEven.Value[1])
	}

	return ctEven
}

// Rotate rotates a ciphertext
func Rotate(ctIn *rlwe.Ciphertext, galEl uint64, permuteNTTindex map[uint64][]uint64, paramsLUT rlwe.Parameters, ks *rlwe.KeySwitcher, rtks *rlwe.RotationKeySet, ctOut *rlwe.Ciphertext) {

	ringQ := paramsLUT.RingQ()
	rtk, _ := rtks.GetRotationKey(galEl)
	level := utils.MinInt(ctIn.Level(), ctOut.Level())
	index := permuteNTTindex[galEl]
	ks.SwitchKeysInPlace(level, ctIn.Value[1], rtk, ks.Pool[1].Q, ks.Pool[2].Q)
	ringQ.AddLvl(level, ks.Pool[1].Q, ctIn.Value[0], ks.Pool[1].Q)
	ringQ.PermuteNTTWithIndexLvl(level, ks.Pool[1].Q, index, ctOut.Value[0])
	ringQ.PermuteNTTWithIndexLvl(level, ks.Pool[2].Q, index, ctOut.Value[1])
}

func (h *Handler) Add(ct0, ct1, ct2 *Ciphertext) {

	level := utils.MinInt(utils.MinInt(ct0.Level(), ct1.Level()), ct2.Level())

	for i := 0; i < level+1; i++ {
		Q := h.paramsLUT.RingQ().Modulus[i]
		ring.AddVec(ct0.Value[i][1:], ct1.Value[i][1:], ct2.Value[i][1:], Q)
		ct2.Value[i][0] = ring.CRed(ct0.Value[i][0]+ct1.Value[i][0], Q)
	}

	ct2.Value = ct2.Value[:level+1]
}

func (h *Handler) Sub(ct0, ct1, ct2 *Ciphertext) {

	level := utils.MinInt(utils.MinInt(ct0.Level(), ct1.Level()), ct2.Level())

	Q := h.paramsLUT.RingQ().Modulus
	for i := 0; i < level+1; i++ {
		ring.SubVec(ct0.Value[i][1:], ct1.Value[i][1:], ct2.Value[i][1:], Q[i])
		ct2.Value[i][0] = ring.CRed(Q[i]+ct0.Value[i][0]-ct1.Value[i][0], Q[i])
	}

	ct2.Value = ct2.Value[:level+1]
}

func (h *Handler) MulScalar(ct0 *Ciphertext, scalar uint64, ct1 *Ciphertext) {

	level := utils.MinInt(ct0.Level(), ct1.Level())

	ringQ := h.paramsLUT.RingQ()
	for i := 0; i < level+1; i++ {
		Q := ringQ.Modulus[i]
		mredParams := ringQ.MredParams[i]
		scalarMont := ring.MForm(scalar, Q, ringQ.BredParams[i])
		ring.MulScalarMontgomeryVec(ct0.Value[i][1:], ct1.Value[i][1:], scalarMont, Q, mredParams)
		ct1.Value[i][0] = ring.MRed(ct0.Value[i][0], scalarMont, Q, mredParams)
	}

	ct1.Value = ct1.Value[:level+1]
}
