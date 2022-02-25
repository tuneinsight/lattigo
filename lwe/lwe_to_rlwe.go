package lwe

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
	"math/bits"
)

// LWEToRLWE transforms a set of LWE samples into their respective RLWE ciphertext such that decrypt(RLWE)[0] = decrypt(LWE)
func (h *Handler) LWEToRLWE(ctLWE map[int]*Ciphertext) (ctRLWE map[int]*rlwe.Ciphertext) {

	var level int
	for i := range ctLWE {
		level = ctLWE[i].Level()
		break
	}

	ringQ := h.paramsLUT.RingQ()
	acc := ringQ.NewPolyLvl(level)
	ctRLWE = make(map[int]*rlwe.Ciphertext)
	for i := range ctLWE {

		// Alocates ciphertext
		ctRLWE[i] = rlwe.NewCiphertextNTT(h.paramsLUT, 1, level)

		for u := 0; u < level+1; u++ {

			ctRLWE[i].Value[0].Coeffs[u][0] = ctLWE[i].Value[u][0]

			// Copy coefficients multiplied by X^{N-1} in reverse order:
			// a_{0} -a_{N-1} -a2_{N-2} ... -a_{1}
			tmp0, tmp1 := acc.Coeffs[u], ctLWE[i].Value[u][1:]
			tmp0[0] = tmp1[0]
			for k := 1; k < ringQ.N; k++ {
				tmp0[k] = ringQ.Modulus[u] - tmp1[ringQ.N-k]
			}

			copy(ctRLWE[i].Value[1].Coeffs[u], acc.Coeffs[u])
		}

		// Switches to NTT domain
		ringQ.NTTLvl(level, ctRLWE[i].Value[0], ctRLWE[i].Value[0])
		ringQ.NTTLvl(level, ctRLWE[i].Value[1], ctRLWE[i].Value[1])
	}

	return
}

// MergeRLWE merges a series of RLWE where decrypt(RLWE)[0] = decrypt(LWE) into a single RLWE.
// inputs RLWE should have been created using the method LWEToRLWE.
func (h *Handler) MergeRLWE(ciphertexts map[int]*rlwe.Ciphertext) (ciphertext *rlwe.Ciphertext) {

	var level int
	for i := range ciphertexts {
		level = ciphertexts[i].Level()
		break
	}

	ringQ := h.paramsLUT.RingQ()

	NInv := ringQ.NttNInv
	Q := ringQ.Modulus
	mredParams := ringQ.MredParams

	// Multiplies by (Slots * N) ^-1 mod Q
	for i := range ciphertexts {
		if ciphertexts[i] != nil {
			v0, v1 := ciphertexts[i].Value[0], ciphertexts[i].Value[1]
			for j := 0; j < level+1; j++ {
				ring.MulScalarMontgomeryVec(v0.Coeffs[j], v0.Coeffs[j], NInv[j], Q[j], mredParams[j])
				ring.MulScalarMontgomeryVec(v1.Coeffs[j], v1.Coeffs[j], NInv[j], Q[j], mredParams[j])
			}
		}
	}

	ciphertextslist := make([]*rlwe.Ciphertext, ringQ.N)

	for i := range ciphertexts {
		ciphertextslist[i] = ciphertexts[i]
	}

	if ciphertextslist[0] == nil {
		ciphertextslist[0] = rlwe.NewCiphertextNTT(h.paramsLUT, 1, level)
	}

	return h.mergeRLWERecurse(ciphertextslist)
}

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

		if ctEven != nil {
			// ctEven + ctOdd * X^(N/2^L)
			ringQ.AddLvl(level, ctEven.Value[0], ctOdd.Value[0], ctEven.Value[0])
			ringQ.AddLvl(level, ctEven.Value[1], ctOdd.Value[1], ctEven.Value[1])

			// phi(ctEven - ctOdd * X^(N/2^L), 2^(L-2))
			ringQ.SubLvl(level, tmpEven.Value[0], ctOdd.Value[0], tmpEven.Value[0])
			ringQ.SubLvl(level, tmpEven.Value[1], ctOdd.Value[1], tmpEven.Value[1])
		}
	}

	if ctEven != nil {

		level := ctEven.Level()

		// if L-2 == -1, then gal = 2N-1
		if L == 1 {
			rotate(tmpEven, uint64(2*ringQ.N-1), h.permuteNTTIndex, h.paramsLUT, h.KeySwitcher, h.rtks, tmpEven)
		} else {
			rotate(tmpEven, h.paramsLUT.GaloisElementForColumnRotationBy(1<<(L-2)), h.permuteNTTIndex, h.paramsLUT, h.KeySwitcher, h.rtks, tmpEven)
		}

		// ctEven + ctOdd * X^(N/2^L) + phi(ctEven - ctOdd * X^(N/2^L), 2^(L-2))
		ringQ.AddLvl(level, ctEven.Value[0], tmpEven.Value[0], ctEven.Value[0])
		ringQ.AddLvl(level, ctEven.Value[1], tmpEven.Value[1], ctEven.Value[1])
	}

	return ctEven
}

func rotate(ctIn *rlwe.Ciphertext, galEl uint64, permuteNTTindex map[uint64][]uint64, paramsLUT rlwe.Parameters, ks *rlwe.KeySwitcher, rtks *rlwe.RotationKeySet, ctOut *rlwe.Ciphertext) {

	ringQ := paramsLUT.RingQ()
	rtk, _ := rtks.GetRotationKey(galEl)
	level := utils.MinInt(ctIn.Level(), ctOut.Level())
	index := permuteNTTindex[galEl]
	ks.SwitchKeysInPlace(level, ctIn.Value[1], rtk, ks.Pool[1].Q, ks.Pool[2].Q)
	ringQ.AddLvl(level, ks.Pool[1].Q, ctIn.Value[0], ks.Pool[1].Q)
	ringQ.PermuteNTTWithIndexLvl(level, ks.Pool[1].Q, index, ctOut.Value[0])
	ringQ.PermuteNTTWithIndexLvl(level, ks.Pool[2].Q, index, ctOut.Value[1])
}
