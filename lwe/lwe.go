package lwe

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// Plaintext is CRT representation of an
// integer m.
type Plaintext struct {
	Value []uint64
}

// Ciphertext is a CRT representation of
// the LWE sample of size N+1, for N
// the degree of the LWE sample.
// The first element of the slice stores the CRT
// representation of <-a, s> + m + e and the N
// next element the CRT representation of a.
type Ciphertext struct {
	Value [][]uint64
}

// Level returns the CRT level of the target.
func (pt *Plaintext) Level() int {
	return len(pt.Value) - 1
}

// Level returns the CRT level of the target.
func (ct *Ciphertext) Level() int {
	return len(ct.Value) - 1
}

// NewCiphertext allocates a new LWE sample of degree N
// and level level.
func NewCiphertext(N, level int) (ct *Ciphertext) {
	ct = new(Ciphertext)
	ct.Value = make([][]uint64, level+1)
	for i := 0; i < level+1; i++ {
		ct.Value[i] = make([]uint64, N+1)
	}
	return ct
}

// Add adds ct0 to ct1 and returns the result on ct2.
func (h *Handler) Add(ct0, ct1, ct2 *Ciphertext) {

	level := utils.MinInt(utils.MinInt(ct0.Level(), ct1.Level()), ct2.Level())

	for i := 0; i < level+1; i++ {
		Q := h.paramsLUT.RingQ().Modulus[i]
		ring.AddVec(ct0.Value[i][1:], ct1.Value[i][1:], ct2.Value[i][1:], Q)
		ct2.Value[i][0] = ring.CRed(ct0.Value[i][0]+ct1.Value[i][0], Q)
	}

	ct2.Value = ct2.Value[:level+1]
}

// Sub subtracts ct1 to ct0 and returns the result on ct2.
func (h *Handler) Sub(ct0, ct1, ct2 *Ciphertext) {

	level := utils.MinInt(utils.MinInt(ct0.Level(), ct1.Level()), ct2.Level())

	Q := h.paramsLUT.RingQ().Modulus
	for i := 0; i < level+1; i++ {
		ring.SubVec(ct0.Value[i][1:], ct1.Value[i][1:], ct2.Value[i][1:], Q[i])
		ct2.Value[i][0] = ring.CRed(Q[i]+ct0.Value[i][0]-ct1.Value[i][0], Q[i])
	}

	ct2.Value = ct2.Value[:level+1]
}

// MulScalar multiplies ct0 by the provided scalar and returns the result on ct1.
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
