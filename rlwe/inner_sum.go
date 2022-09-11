package rlwe

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
)

// InnerSumLog applies an optimized inner sum on the ciphertext (log2(n) + HW(n) rotations with double hoisting).
// The operation assumes that `ctIn` encrypts SlotCount/`batchSize` sub-vectors of size `batchSize` which it adds together (in parallel) by groups of `n`.
// It outputs in ctOut a ciphertext for which the "leftmost" sub-vector of each group is equal to the sum of the group.
// This method is faster than InnerSum when the number of rotations is large and uses log2(n) + HW(n) instead of 'n' keys.
func (eval *Evaluator) InnerSumLog(ctIn *Ciphertext, batchSize, n int, ctOut *Ciphertext) {

	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		panic("ctIn.Degree() != 1 or ctOut.Degree() != 1")
	}

	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()
	ringQP := eval.params.RingQP()

	levelQ := ctIn.Level()
	levelP := len(ringP.Modulus) - 1

	ctOut.Resize(ctOut.Degree(), levelQ)

	if ctIn.Scale != nil {
		ctOut.Scale = ctIn.Scale.CopyNew()
	}

	if n == 1 {
		if ctIn != ctOut {
			ring.CopyValuesLvl(levelQ, ctIn.Value[0], ctOut.Value[0])
			ring.CopyValuesLvl(levelQ, ctIn.Value[1], ctOut.Value[1])
		}
	} else {

		// Memory buffer for ctIn = ctIn + rot(ctIn, 2^i) in Q
		tmpct := &eval.BuffCt
		tmpct.Value[0].IsNTT = true
		tmpct.Value[1].IsNTT = true

		c0OutQP := eval.BuffQP[2]
		c1OutQP := eval.BuffQP[3]
		c0QP := eval.BuffQP[4]
		c1QP := eval.BuffQP[5]

		cQP := CiphertextQP{Value: []ringqp.Poly{c0QP, c1QP}}

		c0QP.Q.IsNTT = true
		c1QP.Q.IsNTT = true

		ctqp := NewCiphertextAtLevelFromPoly(levelQ, [2]*ring.Poly{c0QP.Q, c1QP.Q})

		state := false
		copy := true
		// Binary reading of the input n
		for i, j := 0, n; j > 0; i, j = i+1, j>>1 {

			// Starts by decomposing the input ciphertext
			if i == 0 {
				// If first iteration, then copies directly from the input ciphertext that hasn't been rotated
				eval.DecomposeNTT(levelQ, levelP, levelP+1, ctIn.Value[1], eval.BuffDecompQP)
			} else {
				// Else copies from the rotated input ciphertext
				//ringQ.MulScalarBigintLvl(levelQ, tmpc1, eval.tInvModQ[levelQ], tmpc2)
				eval.DecomposeNTT(levelQ, levelP, levelP+1, tmpct.Value[1], eval.BuffDecompQP)
			}

			// If the binary reading scans a 1
			if j&1 == 1 {

				k := n - (n & ((2 << i) - 1))
				k *= batchSize

				// If the rotation is not zero
				if k != 0 {

					// Rotate(tmpct, k)
					if i == 0 {
						eval.AutomorphismHoistedNoModDown(levelQ, ctIn.Value[0], eval.BuffDecompQP, eval.params.GaloisElementForColumnRotationBy(k), cQP)
					} else {
						eval.AutomorphismHoistedNoModDown(levelQ, tmpct.Value[0], eval.BuffDecompQP, eval.params.GaloisElementForColumnRotationBy(k), cQP)
					}

					// ctOut += Rotate(tmpct, k)
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

						// ctOut += (tmpc0, tmpc1)
						ringQ.AddLvl(levelQ, c0OutQP.Q, tmpct.Value[0], ctOut.Value[0])
						ringQ.AddLvl(levelQ, c1OutQP.Q, tmpct.Value[1], ctOut.Value[1])

					} else {

						ring.CopyValuesLvl(levelQ, tmpct.Value[0], ctOut.Value[0])
						ring.CopyValuesLvl(levelQ, tmpct.Value[1], ctOut.Value[1])

						ctOut.Value[0].IsNTT = true
						ctOut.Value[1].IsNTT = true
					}
				}
			}

			if !state {
				rot := eval.params.GaloisElementForColumnRotationBy((1 << i) * batchSize)
				if i == 0 {

					eval.AutomorphismHoisted(levelQ, ctIn, eval.BuffDecompQP, rot, tmpct)

					ringQ.AddLvl(levelQ, tmpct.Value[0], ctIn.Value[0], tmpct.Value[0])
					ringQ.AddLvl(levelQ, tmpct.Value[1], ctIn.Value[1], tmpct.Value[1])

				} else {
					// (tmpc0, tmpc1) = Rotate((tmpc0, tmpc1), 2^i)

					eval.AutomorphismHoisted(levelQ, tmpct, eval.BuffDecompQP, rot, ctqp)

					ringQ.AddLvl(levelQ, tmpct.Value[0], ctqp.Value[0], tmpct.Value[0])
					ringQ.AddLvl(levelQ, tmpct.Value[1], ctqp.Value[1], tmpct.Value[1])
				}
			}
		}
	}
}
