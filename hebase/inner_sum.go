package hebase

import (
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
)

// InnerSum applies an optimized inner sum on the Ciphertext (log2(n) + HW(n) rotations with double hoisting).
// The operation assumes that `ctIn` encrypts SlotCount/`batchSize` sub-vectors of size `batchSize` which it adds together (in parallel) in groups of `n`.
// It outputs in opOut a Ciphertext for which the "leftmost" sub-vector of each group is equal to the sum of the group.
func (eval Evaluator) InnerSum(ctIn *rlwe.Ciphertext, batchSize, n int, opOut *rlwe.Ciphertext) (err error) {

	params := eval.GetRLWEParameters()

	levelQ := ctIn.Level()
	levelP := params.PCount() - 1

	ringQP := params.RingQP().AtLevel(ctIn.Level(), levelP)

	ringQ := ringQP.RingQ

	opOut.Resize(opOut.Degree(), levelQ)
	*opOut.MetaData = *ctIn.MetaData

	ctInNTT, err := rlwe.NewCiphertextAtLevelFromPoly(levelQ, eval.BuffCt.Value[:2])

	if err != nil {
		panic(err)
	}

	ctInNTT.MetaData = &rlwe.MetaData{}
	ctInNTT.IsNTT = true

	if !ctIn.IsNTT {
		ringQ.NTT(ctIn.Value[0], ctInNTT.Value[0])
		ringQ.NTT(ctIn.Value[1], ctInNTT.Value[1])
	} else {
		ring.CopyLvl(levelQ, ctIn.Value[0], ctInNTT.Value[0])
		ring.CopyLvl(levelQ, ctIn.Value[1], ctInNTT.Value[1])
	}

	if n == 1 {
		if ctIn != opOut {
			ring.CopyLvl(levelQ, ctIn.Value[0], opOut.Value[0])
			ring.CopyLvl(levelQ, ctIn.Value[1], opOut.Value[1])
		}
	} else {

		// BuffQP[0:2] are used by AutomorphismHoistedLazy

		// Accumulator mod QP (i.e. opOut Mod QP)
		accQP := &rlwe.Operand[ringqp.Poly]{Value: []ringqp.Poly{eval.BuffQP[2], eval.BuffQP[3]}}
		accQP.MetaData = ctInNTT.MetaData

		// Buffer mod QP (i.e. to store the result of lazy gadget products)
		cQP := &rlwe.Operand[ringqp.Poly]{Value: []ringqp.Poly{eval.BuffQP[4], eval.BuffQP[5]}}
		cQP.MetaData = ctInNTT.MetaData

		// Buffer mod Q (i.e. to store the result of gadget products)
		cQ, err := rlwe.NewCiphertextAtLevelFromPoly(levelQ, []ring.Poly{cQP.Value[0].Q, cQP.Value[1].Q})

		if err != nil {
			panic(err)
		}

		cQ.MetaData = ctInNTT.MetaData

		state := false
		copy := true
		// Binary reading of the input n
		for i, j := 0, n; j > 0; i, j = i+1, j>>1 {

			// Starts by decomposing the input ciphertext
			eval.DecomposeNTT(levelQ, levelP, levelP+1, ctInNTT.Value[1], true, eval.BuffDecompQP)

			// If the binary reading scans a 1 (j is odd)
			if j&1 == 1 {

				k := n - (n & ((2 << i) - 1))
				k *= batchSize

				// If the rotation is not zero
				if k != 0 {

					rot := params.GaloisElement(k)

					// opOutQP = opOutQP + Rotate(ctInNTT, k)
					if copy {
						if err = eval.AutomorphismHoistedLazy(levelQ, ctInNTT, eval.BuffDecompQP, rot, accQP); err != nil {
							return err
						}
						copy = false
					} else {
						if err = eval.AutomorphismHoistedLazy(levelQ, ctInNTT, eval.BuffDecompQP, rot, cQP); err != nil {
							return err
						}
						ringQP.Add(accQP.Value[0], cQP.Value[0], accQP.Value[0])
						ringQP.Add(accQP.Value[1], cQP.Value[1], accQP.Value[1])
					}

					// j is even
				} else {

					state = true

					// if n is not a power of two, then at least one j was odd, and thus the buffer opOutQP is not empty
					if n&(n-1) != 0 {

						// opOut = opOutQP/P + ctInNTT
						eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, accQP.Value[0].Q, accQP.Value[0].P, opOut.Value[0]) // Division by P
						eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, accQP.Value[1].Q, accQP.Value[1].P, opOut.Value[1]) // Division by P

						ringQ.Add(opOut.Value[0], ctInNTT.Value[0], opOut.Value[0])
						ringQ.Add(opOut.Value[1], ctInNTT.Value[1], opOut.Value[1])

					} else {
						ring.CopyLvl(levelQ, ctInNTT.Value[0], opOut.Value[0])
						ring.CopyLvl(levelQ, ctInNTT.Value[1], opOut.Value[1])
					}
				}
			}

			if !state {

				rot := params.GaloisElement((1 << i) * batchSize)

				// ctInNTT = ctInNTT + Rotate(ctInNTT, 2^i)
				if err = eval.AutomorphismHoisted(levelQ, ctInNTT, eval.BuffDecompQP, rot, cQ); err != nil {
					return err
				}
				ringQ.Add(ctInNTT.Value[0], cQ.Value[0], ctInNTT.Value[0])
				ringQ.Add(ctInNTT.Value[1], cQ.Value[1], ctInNTT.Value[1])
			}
		}
	}

	if !ctIn.IsNTT {
		ringQ.INTT(opOut.Value[0], opOut.Value[0])
		ringQ.INTT(opOut.Value[1], opOut.Value[1])
	}

	return
}

// GaloisElementsForInnerSum returns the list of Galois elements necessary to apply the method
// `InnerSum` operation with parameters `batch` and `n`.
func GaloisElementsForInnerSum(params rlwe.ParameterProvider, batch, n int) (galEls []uint64) {

	rotIndex := make(map[int]bool)

	var k int
	for i := 1; i < n; i <<= 1 {

		k = i
		k *= batch
		rotIndex[k] = true

		k = n - (n & ((i << 1) - 1))
		k *= batch
		rotIndex[k] = true
	}

	rotations := make([]int, len(rotIndex))
	var i int
	for j := range rotIndex {
		rotations[i] = j
		i++
	}

	return params.GetRLWEParameters().GaloisElements(rotations)
}

// Replicate applies an optimized replication on the Ciphertext (log2(n) + HW(n) rotations with double hoisting).
// It acts as the inverse of a inner sum (summing elements from left to right).
// The replication is parameterized by the size of the sub-vectors to replicate "batchSize" and
// the number of times 'n' they need to be replicated.
// To ensure correctness, a gap of zero values of size batchSize * (n-1) must exist between
// two consecutive sub-vectors to replicate.
// This method is faster than Replicate when the number of rotations is large and it uses log2(n) + HW(n) instead of 'n'.
func (eval Evaluator) Replicate(ctIn *rlwe.Ciphertext, batchSize, n int, opOut *rlwe.Ciphertext) (err error) {
	return eval.InnerSum(ctIn, -batchSize, n, opOut)
}

// GaloisElementsForReplicate returns the list of Galois elements necessary to perform the
// `Replicate` operation with parameters `batch` and `n`.
func GaloisElementsForReplicate(params rlwe.ParameterProvider, batch, n int) (galEls []uint64) {
	return GaloisElementsForInnerSum(params, -batch, n)
}
