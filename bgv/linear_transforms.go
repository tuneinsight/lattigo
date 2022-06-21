package bgv

import (
	"github.com/tuneinsight/lattigo/v3/ring"
)

// InnerSumLog applies an optimized inner sum on the ciphertext (log2(n) + HW(n) rotations with double hoisting).
// The operation assumes that `ctIn` encrypts SlotCount/`batchSize` sub-vectors of size `batchSize` which it adds together (in parallel) by groups of `n`.
// It outputs in ctOut a ciphertext for which the "leftmost" sub-vector of each group is equal to the sum of the group.
// This method is faster than InnerSum when the number of rotations is large and uses log2(n) + HW(n) instead of 'n' keys.
func (eval *evaluator) InnerSumLog(ctIn *Ciphertext, batchSize, n int, ctOut *Ciphertext) {

	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()
	ringQP := eval.params.RingQP()

	levelQ := ctIn.Level()
	levelP := len(ringP.Modulus) - 1

	ctOut.Resize(ctOut.Degree(), levelQ)
	ctOut.Scale = ctIn.Scale

	if n == 1 {
		if ctIn != ctOut {
			ring.CopyValuesLvl(levelQ, ctIn.Value[0], ctOut.Value[0])
			ring.CopyValuesLvl(levelQ, ctIn.Value[1], ctOut.Value[1])
		}
	} else {

		// Memory buffer for ctIn = ctIn + rot(ctIn, 2^i) in Q
		tmpc0 := eval.buffQ[0] // unused memory buffer from evaluator
		tmpc1 := eval.buffQ[1] // unused memory buffer from evaluator
		tmpc2 := eval.buffQ[2] // unused memory buffer from evaluator

		c0OutQP := eval.BuffQP[2]
		c1OutQP := eval.BuffQP[3]
		c0QP := eval.BuffQP[4]
		c1QP := eval.BuffQP[5]

		tmpc0.IsNTT = true
		tmpc1.IsNTT = true
		c0QP.Q.IsNTT = true
		c1QP.Q.IsNTT = true

		tmpct := NewCiphertextAtLevelFromPoly(levelQ, [2]*ring.Poly{tmpc0, tmpc1})
		ctqp := NewCiphertextAtLevelFromPoly(levelQ, [2]*ring.Poly{c0QP.Q, c1QP.Q})

		state := false
		copy := true
		// Binary reading of the input n
		for i, j := 0, n; j > 0; i, j = i+1, j>>1 {

			// Starts by decomposing the input ciphertext
			if i == 0 {
				// If first iteration, then copies directly from the input ciphertext that hasn't been rotated
				ringQ.MulScalarBigintLvl(levelQ, ctIn.Value[1], eval.tInvModQ[levelQ], tmpc2)
				eval.DecomposeNTT(levelQ, levelP, levelP+1, tmpc2, eval.BuffDecompQP)
			} else {
				// Else copies from the rotated input ciphertext
				ringQ.MulScalarBigintLvl(levelQ, tmpc1, eval.tInvModQ[levelQ], tmpc2)
				eval.DecomposeNTT(levelQ, levelP, levelP+1, tmpc2, eval.BuffDecompQP)
			}

			// If the binary reading scans a 1
			if j&1 == 1 {

				k := n - (n & ((2 << i) - 1))
				k *= batchSize

				// If the rotation is not zero
				if k != 0 {

					// Rotate((tmpc0, tmpc1), k)
					if i == 0 {
						eval.AutomorphismHoistedNoModDown(levelQ, ctIn.Value[0], eval.BuffDecompQP, eval.params.GaloisElementForColumnRotationBy(k), c0QP.Q, c1QP.Q, c0QP.P, c1QP.P)
					} else {
						eval.AutomorphismHoistedNoModDown(levelQ, tmpc0, eval.BuffDecompQP, eval.params.GaloisElementForColumnRotationBy(k), c0QP.Q, c1QP.Q, c0QP.P, c1QP.P)
					}

					// ctOut += Rotate((tmpc0, tmpc1), k)
					if copy {
						ringQP.CopyValuesLvl(levelQ, levelP, c0QP, c0OutQP)
						ringQP.CopyValuesLvl(levelQ, levelP, c1QP, c1OutQP)
						copy = false
					} else {
						ringQP.AddLvl(levelQ, levelP, c0OutQP, c0QP, c0OutQP)
						ringQP.AddLvl(levelQ, levelP, c1OutQP, c1QP, c1OutQP)
					}
				} else {

					state = true

					// if n is not a power of two
					if n&(n-1) != 0 {
						eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c0OutQP.Q, c0OutQP.P, c0OutQP.Q) // Division by P
						eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c1OutQP.Q, c1OutQP.P, c1OutQP.Q) // Division by P

						ringQ.MulScalarLvl(levelQ, c0OutQP.Q, eval.params.T(), c0OutQP.Q)
						ringQ.MulScalarLvl(levelQ, c1OutQP.Q, eval.params.T(), c1OutQP.Q)

						// ctOut += (tmpc0, tmpc1)
						ringQ.AddLvl(levelQ, c0OutQP.Q, tmpc0, ctOut.Value[0])
						ringQ.AddLvl(levelQ, c1OutQP.Q, tmpc1, ctOut.Value[1])

					} else {

						ring.CopyValuesLvl(levelQ, tmpc0, ctOut.Value[0])
						ring.CopyValuesLvl(levelQ, tmpc1, ctOut.Value[1])

						ctOut.Value[0].IsNTT = true
						ctOut.Value[1].IsNTT = true
					}
				}
			}

			if !state {
				rot := eval.params.GaloisElementForColumnRotationBy((1 << i) * batchSize)
				if i == 0 {
					eval.automorphismHoisted(levelQ, ctIn, eval.BuffDecompQP, rot, tmpct)
					ringQ.AddLvl(levelQ, tmpc0, ctIn.Value[0], tmpc0)
					ringQ.AddLvl(levelQ, tmpc1, ctIn.Value[1], tmpc1)
				} else {
					// (tmpc0, tmpc1) = Rotate((tmpc0, tmpc1), 2^i)
					eval.automorphismHoisted(levelQ, tmpct, eval.BuffDecompQP, rot, ctqp)
					ringQ.AddLvl(levelQ, tmpc0, c0QP.Q, tmpc0)
					ringQ.AddLvl(levelQ, tmpc1, c1QP.Q, tmpc1)
				}
			}
		}
	}
}
