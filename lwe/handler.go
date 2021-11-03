package lwe


import(
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
	"math/bits"
)


type Handler struct{
	*rlwe.KeySwitcher
	params rlwe.Parameters
	rtks *rlwe.RotationKeySet

	nPowInv [][]uint64
	xPow []*ring.Poly

	permuteNTTIndex map[uint64][]uint64
}

func NewHandler(params rlwe.Parameters, rtks *rlwe.RotationKeySet) (h *Handler){
	h = new(Handler)
	h.KeySwitcher = rlwe.NewKeySwitcher(params)
	h.params = params

	h.nPowInv = make([][]uint64, params.LogN())
	h.xPow = make([]*ring.Poly, params.LogN())

	ringQ := params.RingQ()
	for i := 0; i < params.LogN(); i++ {

		h.nPowInv[i] = make([]uint64, params.MaxLevel()+1)
		for j := 0; j < params.MaxLevel()+1; j++{
			h.nPowInv[i][j] = ring.MForm(ring.ModExp(2<<i, ringQ.Modulus[j]-2, ringQ.Modulus[j]), ringQ.Modulus[j], ringQ.BredParams[j])
		}

		h.xPow[i] = ringQ.NewPoly()
		if i == 0{
			for j := 0; j < params.MaxLevel()+1; j++{
				h.xPow[i].Coeffs[j][1<<i] = ring.MForm(1, ringQ.Modulus[j], ringQ.BredParams[j])
			}
			ringQ.NTT(h.xPow[i], h.xPow[i])
		}else{
			ringQ.MulCoeffsMontgomery(h.xPow[i-1], h.xPow[i-1], h.xPow[i]) // X^{2^(n+1)} = X^{2^n} * X^{2^n}
		}		
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


func (h *Handler) MergeRLWE(ciphertexts []*rlwe.Ciphertext) (ciphertext *rlwe.Ciphertext){

	slots := len(ciphertexts)

	if slots & (slots-1) != 0 {
		panic("len(ciphertext) must be a power of two smaller or equal to the ring degree")
	}

	logSlots := bits.Len64(uint64(len(ciphertexts))) - 1

	level := ciphertexts[0].Level()

	ringQ := h.params.RingQ()

	nPowInv := h.nPowInv[logSlots-1]
	Q := ringQ.Modulus
	mredParams := ringQ.MredParams

	// Multiplies by (Slots * N) ^-1 mod Q
	for i := range ciphertexts{
		v0, v1 := ciphertexts[i].Value[0], ciphertexts[i].Value[1]
		for j := 0; j < ciphertexts[0].Level()+1; j++{
			ring.MulScalarMontgomeryVec(v0.Coeffs[j], v0.Coeffs[j], nPowInv[j], Q[j], mredParams[j])
			ring.MulScalarMontgomeryVec(v1.Coeffs[j], v1.Coeffs[j], nPowInv[j], Q[j], mredParams[j])
		}
	}

	// Padds for the repacking algorithm
	if slots != h.params.N(){
		ciphertexts = append(ciphertexts, make([]*rlwe.Ciphertext, h.params.N()-len(ciphertexts))...)
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
func (h *Handler) mergeRLWERecurse(ciphertexts []*rlwe.Ciphertext) (*rlwe.Ciphertext){

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