package lwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
	"math/bits"
)

type Handler struct {
	*rlwe.KeySwitcher
	params rlwe.Parameters
	rtks   *rlwe.RotationKeySet

	nPowInv      [][]uint64
	xPow         []*ring.Poly
	xPowMinusOne []*ring.Poly

	permuteNTTIndex map[uint64][]uint64
}

func NewHandler(params rlwe.Parameters, rtks *rlwe.RotationKeySet) (h *Handler) {
	h = new(Handler)
	h.KeySwitcher = rlwe.NewKeySwitcher(params)
	h.params = params

	ringQ := params.RingQ()

	h.nPowInv = make([][]uint64, params.LogN())
	h.xPow = make([]*ring.Poly, params.LogN())

	for i := 0; i < params.LogN(); i++ {
		h.nPowInv[i] = make([]uint64, params.MaxLevel()+1)
		var nTimesN uint64 = 1 << (params.LogN() + i)
		for j := 0; j < params.MaxLevel()+1; j++ {
			h.nPowInv[i][j] = ring.MForm(ring.ModExp(nTimesN, ringQ.Modulus[j]-2, ringQ.Modulus[j]), ringQ.Modulus[j], ringQ.BredParams[j])
		}

		h.xPow[i] = ringQ.NewPoly()
		if i == 0 {
			for j := 0; j < params.MaxLevel()+1; j++ {
				h.xPow[i].Coeffs[j][1<<i] = ring.MForm(1, ringQ.Modulus[j], ringQ.BredParams[j])
			}
			ringQ.NTT(h.xPow[i], h.xPow[i])
		} else {
			ringQ.MulCoeffsMontgomery(h.xPow[i-1], h.xPow[i-1], h.xPow[i]) // X^{n} = X^{1} * X^{n-1}
		}
	}

	one := ringQ.NewPoly()
	for i := 0; i < params.MaxLevel()+1; i++ {
		one.Coeffs[i][0] = 1
	}
	ringQ.NTT(one, one)

	N := ringQ.N

	h.xPowMinusOne = make([]*ring.Poly, 2*N)
	for i := 0; i < N; i++ {
		h.xPowMinusOne[i] = ringQ.NewPoly()
		h.xPowMinusOne[i+N] = ringQ.NewPoly()
		if i == 0 {
			for j := 0; j < params.MaxLevel()+1; j++ {
				h.xPowMinusOne[i].Coeffs[j][i] = ring.MForm(1, ringQ.Modulus[j], ringQ.BredParams[j])
			}
			ringQ.NTT(h.xPowMinusOne[i], h.xPowMinusOne[i])
		} else {
			ringQ.MulCoeffsMontgomery(h.xPowMinusOne[0], h.xPowMinusOne[i-1], h.xPowMinusOne[i]) // X^{n} = X^{1} * X^{n-1}
			ringQ.Neg(h.xPowMinusOne[i], h.xPowMinusOne[i+N])                                    // X^{2n} = -X^{1} * X^{n-1}
		}
	}
	for i := 0; i < 2*N; i++ {
		ringQ.Sub(h.xPowMinusOne[i], one, h.xPowMinusOne[i]) // X^{n} - 1
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
		permuteNTTIndex[galEl] = h.params.RingQ().PermuteNTTIndex(galEl)
	}
	return &permuteNTTIndex
}

func (h *Handler) MergeRLWE(ciphertexts []*rlwe.Ciphertext) (ciphertext *rlwe.Ciphertext) {

	slots := len(ciphertexts)

	if slots&(slots-1) != 0 {
		panic("len(ciphertext) must be a power of two smaller or equal to the ring degree")
	}

	logSlots := bits.Len64(uint64(len(ciphertexts))) - 1

	level := ciphertexts[0].Level()

	ringQ := h.params.RingQ()

	nPowInv := h.nPowInv[h.params.LogN()-logSlots]
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
	if slots != h.params.N() {
		ciphertexts = append(ciphertexts, make([]*rlwe.Ciphertext, h.params.N()-len(ciphertexts))...)
		N := ringQ.N
		gap := N / slots
		for i := 0; i < slots; i++ {
			ciphertexts[N-(i+1)*gap], ciphertexts[slots-i-1] = ciphertexts[slots-i-1], ciphertexts[N-(i+1)*gap]
		}
	}

	ciphertext = h.mergeRLWERecurse(ciphertexts)

	tmp := rlwe.NewCiphertextNTT(h.params, 1, ciphertext.Level())
	for i := logSlots - 1; i < h.params.LogN()-1; i++ {
		Rotate(ciphertext, h.params.GaloisElementForColumnRotationBy(1<<i), h.permuteNTTIndex, h.params, h.KeySwitcher, h.rtks, tmp)
		ringQ.AddLvl(level, ciphertext.Value[0], tmp.Value[0], ciphertext.Value[0])
		ringQ.AddLvl(level, ciphertext.Value[1], tmp.Value[1], ciphertext.Value[1])
	}

	return
}

// PackLWEs repacks LWE ciphertexts into a RLWE ciphertext
func (h *Handler) mergeRLWERecurse(ciphertexts []*rlwe.Ciphertext) *rlwe.Ciphertext {

	ringQ := h.params.RingQ()

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
			Rotate(tmpEven, uint64(2*ringQ.N-1), h.permuteNTTIndex, h.params, h.KeySwitcher, h.rtks, tmpEven)
		} else {
			Rotate(tmpEven, h.params.GaloisElementForColumnRotationBy(1<<(L-2)), h.permuteNTTIndex, h.params, h.KeySwitcher, h.rtks, tmpEven)
		}

		// ctEven + ctOdd * X^(N/2^L) + phi(ctEven - ctOdd * X^(N/2^L), 2^(L-2))
		ringQ.AddLvl(level, ctEven.Value[0], tmpEven.Value[0], ctEven.Value[0])
		ringQ.AddLvl(level, ctEven.Value[1], tmpEven.Value[1], ctEven.Value[1])
	}

	return ctEven
}

// Rotate rotates a ciphertext
func Rotate(ctIn *rlwe.Ciphertext, galEl uint64, permuteNTTindex map[uint64][]uint64, params rlwe.Parameters, ks *rlwe.KeySwitcher, rtks *rlwe.RotationKeySet, ctOut *rlwe.Ciphertext) {

	ringQ := params.RingQ()
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
		Q := h.params.RingQ().Modulus[i]
		ring.AddVec(ct0.Value[i][1:], ct1.Value[i][1:], ct2.Value[i][1:], Q)
		ct2.Value[i][0] = ring.CRed(ct0.Value[i][0]+ct1.Value[i][0], Q)
	}

	ct2.Value = ct2.Value[:level+1]
}

func (h *Handler) Sub(ct0, ct1, ct2 *Ciphertext) {

	level := utils.MinInt(utils.MinInt(ct0.Level(), ct1.Level()), ct2.Level())

	Q := h.params.RingQ().Modulus
	for i := 0; i < level+1; i++ {
		ring.SubVec(ct0.Value[i][1:], ct1.Value[i][1:], ct2.Value[i][1:], Q[i])
		ct2.Value[i][0] = ring.CRed(Q[i]+ct0.Value[i][0]-ct1.Value[i][0], Q[i])
	}

	ct2.Value = ct2.Value[:level+1]
}

func (h *Handler) MulScalar(ct0 *Ciphertext, scalar uint64, ct1 *Ciphertext) {

	level := utils.MinInt(ct0.Level(), ct1.Level())

	ringQ := h.params.RingQ()
	for i := 0; i < level+1; i++ {
		Q := ringQ.Modulus[i]
		mredParams := ringQ.MredParams[i]
		scalarMont := ring.MForm(scalar, Q, ringQ.BredParams[i])
		ring.MulScalarMontgomeryVec(ct0.Value[i][1:], ct1.Value[i][1:], scalarMont, Q, mredParams)
		ct1.Value[i][0] = ring.MRed(ct0.Value[i][0], scalarMont, Q, mredParams)
	}

	ct1.Value = ct1.Value[:level+1]
}
