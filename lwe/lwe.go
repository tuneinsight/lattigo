package lwe

import(
	"github.com/ldsec/lattigo/v2/utils"
	"github.com/ldsec/lattigo/v2/ring"
)

type Plaintext struct {
	Value []uint64
}

type Ciphertext struct {
	Value [][]uint64
}

func (pt *Plaintext) Level() int {
	return len(pt.Value) - 1
}

func (ct *Ciphertext) Level() int {
	return len(ct.Value) - 1
}

func NewCiphertext(N, level int) (ct *Ciphertext) {
	ct = new(Ciphertext)
	ct.Value = make([][]uint64, level+1)
	for i := 0; i < level+1; i++ {
		ct.Value[i] = make([]uint64, N+1)
	}
	return ct
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