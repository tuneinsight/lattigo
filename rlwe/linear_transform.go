package rlwe

import (
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// Expand expands a RLWE Ciphertext encrypting sum ai * X^i to 2^logN ciphertexts,
// each encrypting ai * X^0 for 0 <= i < 2^LogN. That is, it extracts the first 2^logN
// coefficients, whose degree is a multiple of 2^logGap, of ctIn and returns an RLWE
// Ciphertext for each coefficient extracted.
func (eval *Evaluator) Expand(ctIn *Ciphertext, logN, logGap int) (ctOut []*Ciphertext) {

	if ctIn.Degree() != 1 {
		panic("ctIn.Degree() != 1")
	}

	if eval.params.RingType() != ring.Standard {
		panic("Expand is only supported for ring.Type = ring.Standard (X^{-2^{i}} does not exist in the sub-ring Z[X + X^{-1}])")
	}

	params := eval.params

	level := ctIn.Level()

	ringQ := params.RingQ().AtLevel(level)

	// Compute X^{-2^{i}} from 1 to LogN
	xPow2 := genXPow2(ringQ, logN, true)

	ctOut = make([]*Ciphertext, 1<<(logN-logGap))
	ctOut[0] = ctIn.CopyNew()

	if ct := ctOut[0]; !ctIn.IsNTT {
		ringQ.NTT(ct.Value[0], ct.Value[0])
		ringQ.NTT(ct.Value[1], ct.Value[1])
		ct.IsNTT = true
	}

	// Multiplies by 2^{-logN} mod Q
	NInv := new(big.Int).SetUint64(1 << logN)
	NInv.ModInverse(NInv, ringQ.ModulusAtLevel[level])

	ringQ.MulScalarBigint(ctOut[0].Value[0], NInv, ctOut[0].Value[0])
	ringQ.MulScalarBigint(ctOut[0].Value[1], NInv, ctOut[0].Value[1])

	gap := 1 << logGap

	tmp := NewCiphertextAtLevelFromPoly(level, []*ring.Poly{eval.BuffCt.Value[0], eval.BuffCt.Value[1]})
	tmp.MetaData = ctIn.MetaData

	for i := 0; i < logN; i++ {

		n := 1 << i

		galEl := uint64(ringQ.N()/n + 1)

		half := n / gap

		for j := 0; j < (n+gap-1)/gap; j++ {

			c0 := ctOut[j]

			// X -> X^{N/n + 1}
			//[a, b, c, d] -> [a, -b, c, -d]
			eval.Automorphism(c0, galEl, tmp)

			if j+half > 0 {

				c1 := ctOut[j].CopyNew()

				// Zeroes odd coeffs: [a, b, c, d] + [a, -b, c, -d] -> [2a, 0, 2b, 0]
				ringQ.Add(c0.Value[0], tmp.Value[0], c0.Value[0])
				ringQ.Add(c0.Value[1], tmp.Value[1], c0.Value[1])

				// Zeroes even coeffs: [a, b, c, d] - [a, -b, c, -d] -> [0, 2b, 0, 2d]
				ringQ.Sub(c1.Value[0], tmp.Value[0], c1.Value[0])
				ringQ.Sub(c1.Value[1], tmp.Value[1], c1.Value[1])

				// c1 * X^{-2^{i}}: [0, 2b, 0, 2d] * X^{-n} -> [2b, 0, 2d, 0]
				ringQ.MulCoeffsMontgomery(c1.Value[0], xPow2[i], c1.Value[0])
				ringQ.MulCoeffsMontgomery(c1.Value[1], xPow2[i], c1.Value[1])

				ctOut[j+half] = c1

			} else {

				// Zeroes odd coeffs: [a, b, c, d] + [a, -b, c, -d] -> [2a, 0, 2b, 0]
				ringQ.Add(c0.Value[0], tmp.Value[0], c0.Value[0])
				ringQ.Add(c0.Value[1], tmp.Value[1], c0.Value[1])
			}
		}
	}

	for _, ct := range ctOut {
		if ct != nil && !ctIn.IsNTT {
			ringQ.INTT(ct.Value[0], ct.Value[0])
			ringQ.INTT(ct.Value[1], ct.Value[1])
			ct.IsNTT = false
		}
	}
	return
}

// Merge merges a batch of RLWE, packing the first coefficient of each RLWE into a single RLWE.
// The operation will require N/gap + log(gap) key-switches, where gap is the minimum gap between
// two non-zero coefficients of the final Ciphertext.
// The method takes as input a map of Ciphertext, indexing in which coefficient of the final
// Ciphertext the first coefficient of each Ciphertext of the map must be packed.
// This method accepts ciphertexts both in and out of the NTT domain, but the result
// is always returned in the NTT domain.
func (eval *Evaluator) Merge(ctIn map[int]*Ciphertext) (ctOut *Ciphertext) {

	if eval.params.RingType() != ring.Standard {
		panic("Merge is only supported for ring.Type = ring.Standard (X^{2^{i}} does not exist in the sub-ring Z[X + X^{-1}])")
	}

	params := eval.params
	ringQ := params.RingQ()

	var levelQ int
	for i := range ctIn {
		levelQ = ctIn[i].Level()
		break
	}

	for i := range ctIn {
		levelQ = utils.Min(levelQ, ctIn[i].Level())
	}

	xPow2 := genXPow2(ringQ.AtLevel(levelQ), params.LogN(), false)

	// Multiplies by (Slots * N) ^-1 mod Q
	for i := range ctIn {
		if ctIn[i] != nil {

			if ctIn[i].Degree() != 1 {
				panic("cannot Merge: ctIn.Degree() != 1")
			}

			v0, v1 := ctIn[i].Value[0], ctIn[i].Value[1]
			for j, s := range ringQ.SubRings[:levelQ+1] {
				s.MulScalarMontgomery(v0.Coeffs[j], s.NInv, v0.Coeffs[j])
				s.MulScalarMontgomery(v1.Coeffs[j], s.NInv, v1.Coeffs[j])
			}
		}
	}

	ciphertextslist := make([]*Ciphertext, ringQ.N())

	for i := range ctIn {
		ciphertextslist[i] = ctIn[i]
	}

	if ciphertextslist[0] == nil {
		ciphertextslist[0] = NewCiphertext(params, 1, levelQ)
		ciphertextslist[0].IsNTT = true
	}

	return eval.mergeRLWERecurse(ciphertextslist, xPow2)
}

func (eval *Evaluator) mergeRLWERecurse(ciphertexts []*Ciphertext, xPow []*ring.Poly) *Ciphertext {

	L := bits.Len64(uint64(len(ciphertexts))) - 1

	if L == 0 {
		return ciphertexts[0]
	}

	odd := make([]*Ciphertext, len(ciphertexts)>>1)
	even := make([]*Ciphertext, len(ciphertexts)>>1)

	for i := 0; i < len(ciphertexts)>>1; i++ {
		odd[i] = ciphertexts[2*i]
		even[i] = ciphertexts[2*i+1]
	}

	ctEven := eval.mergeRLWERecurse(odd, xPow)
	ctOdd := eval.mergeRLWERecurse(even, xPow)

	if ctEven == nil && ctOdd == nil {
		return nil
	}

	var level = 0xFFFF // Case if ctOdd == nil

	if ctOdd != nil {
		level = ctOdd.Level()
	}

	if ctEven != nil {
		level = utils.Min(level, ctEven.Level())
	}

	ringQ := eval.params.RingQ().AtLevel(level)

	if ctOdd != nil {
		if !ctOdd.IsNTT {
			ringQ.NTT(ctOdd.Value[0], ctOdd.Value[0])
			ringQ.NTT(ctOdd.Value[1], ctOdd.Value[1])
			ctOdd.IsNTT = true
		}
	}

	if ctEven != nil {
		if !ctEven.IsNTT {
			ringQ.NTT(ctEven.Value[0], ctEven.Value[0])
			ringQ.NTT(ctEven.Value[1], ctEven.Value[1])
			ctEven.IsNTT = true
		}
	}

	var tmpEven *Ciphertext
	if ctEven != nil {
		tmpEven = ctEven.CopyNew()
	}

	// ctOdd * X^(N/2^L)
	if ctOdd != nil {

		//X^(N/2^L)
		ringQ.MulCoeffsMontgomery(ctOdd.Value[0], xPow[len(xPow)-L], ctOdd.Value[0])
		ringQ.MulCoeffsMontgomery(ctOdd.Value[1], xPow[len(xPow)-L], ctOdd.Value[1])

		if ctEven != nil {
			// ctEven + ctOdd * X^(N/2^L)
			ringQ.Add(ctEven.Value[0], ctOdd.Value[0], ctEven.Value[0])
			ringQ.Add(ctEven.Value[1], ctOdd.Value[1], ctEven.Value[1])

			// phi(ctEven - ctOdd * X^(N/2^L), 2^(L-2))
			ringQ.Sub(tmpEven.Value[0], ctOdd.Value[0], tmpEven.Value[0])
			ringQ.Sub(tmpEven.Value[1], ctOdd.Value[1], tmpEven.Value[1])
		}
	}

	if ctEven != nil {

		// if L-2 == -1, then gal = -1
		if L == 1 {
			eval.Automorphism(tmpEven, ringQ.NthRoot()-1, tmpEven)
		} else {
			eval.Automorphism(tmpEven, eval.params.GaloisElementForColumnRotationBy(1<<(L-2)), tmpEven)
		}

		// ctEven + ctOdd * X^(N/2^L) + phi(ctEven - ctOdd * X^(N/2^L), 2^(L-2))
		ringQ.Add(ctEven.Value[0], tmpEven.Value[0], ctEven.Value[0])
		ringQ.Add(ctEven.Value[1], tmpEven.Value[1], ctEven.Value[1])
	}

	return ctEven
}

func genXPow2(r *ring.Ring, logN int, div bool) (xPow []*ring.Poly) {

	// Compute X^{-n} from 0 to LogN
	xPow = make([]*ring.Poly, logN)

	moduli := r.ModuliChain()[:r.Level()+1]
	BRC := r.BRedConstants()

	var idx int
	for i := 0; i < logN; i++ {

		idx = 1 << i

		if div {
			idx = r.N() - idx
		}

		xPow[i] = r.NewPoly()

		if i == 0 {

			for j := range moduli {
				xPow[i].Coeffs[j][idx] = ring.MForm(1, moduli[j], BRC[j])
			}

			r.NTT(xPow[i], xPow[i])

		} else {
			r.MulCoeffsMontgomery(xPow[i-1], xPow[i-1], xPow[i]) // X^{n} = X^{1} * X^{n-1}
		}
	}

	if div {
		r.Neg(xPow[0], xPow[0])
	}

	return
}

// InnerSum applies an optimized inner sum on the Ciphertext (log2(n) + HW(n) rotations with double hoisting).
// The operation assumes that `ctIn` encrypts SlotCount/`batchSize` sub-vectors of size `batchSize` which it adds together (in parallel) in groups of `n`.
// It outputs in ctOut a Ciphertext for which the "leftmost" sub-vector of each group is equal to the sum of the group.
func (eval *Evaluator) InnerSum(ctIn *Ciphertext, batchSize, n int, ctOut *Ciphertext) {

	levelQ := ctIn.Level()
	levelP := eval.params.PCount() - 1

	ringQP := eval.params.RingQP().AtLevel(ctIn.Level(), levelP)

	ringQ := ringQP.RingQ

	ctOut.Resize(ctOut.Degree(), levelQ)
	ctOut.MetaData = ctIn.MetaData

	ctInNTT := NewCiphertextAtLevelFromPoly(levelQ, eval.BuffCt.Value[:2])
	ctInNTT.IsNTT = true

	if !ctIn.IsNTT {
		ringQ.NTT(ctIn.Value[0], ctInNTT.Value[0])
		ringQ.NTT(ctIn.Value[1], ctInNTT.Value[1])
	} else {
		ring.CopyLvl(levelQ, ctIn.Value[0], ctInNTT.Value[0])
		ring.CopyLvl(levelQ, ctIn.Value[1], ctInNTT.Value[1])
	}

	if n == 1 {
		if ctIn != ctOut {
			ring.CopyLvl(levelQ, ctIn.Value[0], ctOut.Value[0])
			ring.CopyLvl(levelQ, ctIn.Value[1], ctOut.Value[1])
		}
	} else {

		// BuffQP[0:2] are used by AutomorphismHoistedLazy

		// Accumulator mod QP (i.e. ctOut Mod QP)
		accQP := CiphertextQP{Value: [2]ringqp.Poly{eval.BuffQP[2], eval.BuffQP[3]}}
		accQP.IsNTT = true

		// Buffer mod QP (i.e. to store the result of lazy gadget products)
		cQP := CiphertextQP{Value: [2]ringqp.Poly{eval.BuffQP[4], eval.BuffQP[5]}}
		cQP.IsNTT = true

		// Buffer mod Q (i.e. to store the result of gadget products)
		cQ := NewCiphertextAtLevelFromPoly(levelQ, []*ring.Poly{cQP.Value[0].Q, cQP.Value[1].Q})
		cQ.IsNTT = true

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

					rot := eval.params.GaloisElementForColumnRotationBy(k)

					// ctOutQP = ctOutQP + Rotate(ctInNTT, k)
					if copy {
						eval.AutomorphismHoistedLazy(levelQ, ctInNTT, eval.BuffDecompQP, rot, accQP)
						copy = false
					} else {
						eval.AutomorphismHoistedLazy(levelQ, ctInNTT, eval.BuffDecompQP, rot, cQP)
						ringQP.Add(accQP.Value[0], cQP.Value[0], accQP.Value[0])
						ringQP.Add(accQP.Value[1], cQP.Value[1], accQP.Value[1])
					}

					// j is even
				} else {

					state = true

					// if n is not a power of two, then at least one j was odd, and thus the buffer ctOutQP is not empty
					if n&(n-1) != 0 {

						// ctOut = ctOutQP/P + ctInNTT
						eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, accQP.Value[0].Q, accQP.Value[0].P, ctOut.Value[0]) // Division by P
						eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, accQP.Value[1].Q, accQP.Value[1].P, ctOut.Value[1]) // Division by P

						ringQ.Add(ctOut.Value[0], ctInNTT.Value[0], ctOut.Value[0])
						ringQ.Add(ctOut.Value[1], ctInNTT.Value[1], ctOut.Value[1])

					} else {
						ring.CopyLvl(levelQ, ctInNTT.Value[0], ctOut.Value[0])
						ring.CopyLvl(levelQ, ctInNTT.Value[1], ctOut.Value[1])
					}
				}
			}

			if !state {

				rot := eval.params.GaloisElementForColumnRotationBy((1 << i) * batchSize)

				// ctInNTT = ctInNTT + Rotate(ctInNTT, 2^i)
				eval.AutomorphismHoisted(levelQ, ctInNTT, eval.BuffDecompQP, rot, cQ)
				ringQ.Add(ctInNTT.Value[0], cQ.Value[0], ctInNTT.Value[0])
				ringQ.Add(ctInNTT.Value[1], cQ.Value[1], ctInNTT.Value[1])
			}
		}
	}

	if !ctIn.IsNTT {
		ringQ.INTT(ctOut.Value[0], ctOut.Value[0])
		ringQ.INTT(ctOut.Value[1], ctOut.Value[1])
	}
}

// Replicate applies an optimized replication on the Ciphertext (log2(n) + HW(n) rotations with double hoisting).
// It acts as the inverse of a inner sum (summing elements from left to right).
// The replication is parameterized by the size of the sub-vectors to replicate "batchSize" and
// the number of times 'n' they need to be replicated.
// To ensure correctness, a gap of zero values of size batchSize * (n-1) must exist between
// two consecutive sub-vectors to replicate.
// This method is faster than Replicate when the number of rotations is large and it uses log2(n) + HW(n) instead of 'n'.
func (eval *Evaluator) Replicate(ctIn *Ciphertext, batchSize, n int, ctOut *Ciphertext) {
	eval.InnerSum(ctIn, -batchSize, n, ctOut)
}
