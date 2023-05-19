package rlwe

import (
	"fmt"
	"math/big"

	//"math/bits"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// Trace maps X -> sum((-1)^i * X^{i*n+1}) for n <= i < N
// Monomial X^k vanishes if k is not divisible by (N/n), otherwise it is multiplied by (N/n).
// Ciphertext is pre-multiplied by (N/n)^-1 to remove the (N/n) factor.
// Examples of full Trace for [0 + 1X + 2X^2 + 3X^3 + 4X^4 + 5X^5 + 6X^6 + 7X^7]
//
// 1.
//
//	  [1 + 2X + 3X^2 + 4X^3 + 5X^4 + 6X^5 + 7X^6 + 8X^7]
//	+ [1 - 6X - 3X^2 + 8X^3 + 5X^4 + 2X^5 - 7X^6 - 4X^7]  {X-> X^(i * 5^1)}
//	= [2 - 4X + 0X^2 +12X^3 +10X^4 + 8X^5 - 0X^6 + 4X^7]
//
// 2.
//
//	  [2 - 4X + 0X^2 +12X^3 +10X^4 + 8X^5 - 0X^6 + 4X^7]
//	+ [2 + 4X + 0X^2 -12X^3 +10X^4 - 8X^5 + 0X^6 - 4X^7]  {X-> X^(i * 5^2)}
//	= [4 + 0X + 0X^2 - 0X^3 +20X^4 + 0X^5 + 0X^6 - 0X^7]
//
// 3.
//
//	  [4 + 0X + 0X^2 - 0X^3 +20X^4 + 0X^5 + 0X^6 - 0X^7]
//	+ [4 + 0X + 0X^2 - 0X^3 -20X^4 + 0X^5 + 0X^6 - 0X^7]  {X-> X^(i * -1)}
//	= [8 + 0X + 0X^2 - 0X^3 + 0X^4 + 0X^5 + 0X^6 - 0X^7]
func (eval *Evaluator) Trace(ctIn *Ciphertext, logN int, ctOut *Ciphertext) {

	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		panic("ctIn.Degree() != 1 or ctOut.Degree() != 1")
	}

	level := utils.Min(ctIn.Level(), ctOut.Level())

	ctOut.Resize(ctOut.Degree(), level)

	ctOut.MetaData = ctIn.MetaData

	gap := 1 << (eval.params.LogN() - logN - 1)

	if logN == 0 {
		gap <<= 1
	}

	if gap > 1 {

		ringQ := eval.params.RingQ().AtLevel(level)

		if ringQ.Type() == ring.ConjugateInvariant {
			gap >>= 1 // We skip the last step that applies phi(5^{-1})
		}

		NInv := new(big.Int).SetUint64(uint64(gap))
		NInv.ModInverse(NInv, ringQ.ModulusAtLevel[level])

		// pre-multiplication by (N/n)^-1
		ringQ.MulScalarBigint(ctIn.Value[0], NInv, ctOut.Value[0])
		ringQ.MulScalarBigint(ctIn.Value[1], NInv, ctOut.Value[1])

		if !ctIn.IsNTT {
			ringQ.NTT(ctOut.Value[0], ctOut.Value[0])
			ringQ.NTT(ctOut.Value[1], ctOut.Value[1])
			ctOut.IsNTT = true
		}

		buff := NewCiphertextAtLevelFromPoly(level, []*ring.Poly{eval.BuffQP[3].Q, eval.BuffQP[4].Q})
		buff.IsNTT = true

		for i := logN; i < eval.params.LogN()-1; i++ {
			eval.Automorphism(ctOut, eval.params.GaloisElementForColumnRotationBy(1<<i), buff)
			ringQ.Add(ctOut.Value[0], buff.Value[0], ctOut.Value[0])
			ringQ.Add(ctOut.Value[1], buff.Value[1], ctOut.Value[1])
		}

		if logN == 0 && ringQ.Type() == ring.Standard {
			eval.Automorphism(ctOut, ringQ.NthRoot()-1, buff)
			ringQ.Add(ctOut.Value[0], buff.Value[0], ctOut.Value[0])
			ringQ.Add(ctOut.Value[1], buff.Value[1], ctOut.Value[1])
		}

		if !ctIn.IsNTT {
			ringQ.INTT(ctOut.Value[0], ctOut.Value[0])
			ringQ.INTT(ctOut.Value[1], ctOut.Value[1])
			ctOut.IsNTT = false
		}

	} else {
		if ctIn != ctOut {
			ctOut.Copy(ctIn)
		}
	}
}

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

// Pack packs a batch of RLWE ciphertexts, packing the batch of ciphertexts into a single ciphertext.
// The number of key-switching operations is log(N/len(cts)) + len(cts)
//
// Input:
//
//	cts: a map of rlwe.Ciphertext, where the index in the map is the future position of the first coefficient
//	      of the indexed ciphertext in the final ciphertext (see example). Ciphertexts can be in or out of the NTT domain.
//	logGap: all coefficients of the input ciphertexts that are not a multiple of X^{2^{logGap}} will be zeroed
//	        during the merging (see example). This is equivalent to skipping the first 2^{logGap} steps of the
//	        algorithm, i.e. having as input ciphertexts that are already partially packed.
//
// Output: a ciphertext packing all input ciphertexts
//
// Example: we want to pack 4 ciphertexts into one, and keep only coefficients which are a multiple of X^{4}.
//
//	To do so, we must set logGap = 2.
//	Here the `X` slots are treated as garbage slots that we want to discard during the procedure.
//
//	input: map[int]{
//	   0: [x00, X, X, X, x01, X, X, X],   with logGap = 2
//	   1: [x10, X, X, X, x11, X, X, X],
//	   2: [x20, X, X, X, x21, X, X, X],
//	   3: [x30, X, X, X, x31, X, X, X],
//		}
//
//	 Step 1:
//	         map[0]: 2^{-1} * (map[0] + X^2 * map[2] + phi_{5^2}(map[0] - X^2 * map[2]) = [x00, X, x20, X, x01, X, x21, X]
//	         map[1]: 2^{-1} * (map[1] + X^2 * map[3] + phi_{5^2}(map[1] - X^2 * map[3]) = [x10, X, x30, X, x11, X, x31, X]
//	 Step 2:
//	         map[0]: 2^{-1} * (map[0] + X^1 * map[1] + phi_{5^4}(map[0] - X^1 * map[1]) = [x00, x10, x20, x30, x01, x11, x21, x22]
//
// Note: any usued coefficient in the output ciphertext will be zeroed during the procedure.
func (eval *Evaluator) Pack(cts map[int]*Ciphertext, inputLogGap int) (ct *Ciphertext) {

	params := eval.Parameters()

	if params.RingType() != ring.Standard {
		panic(fmt.Errorf("cannot Pack: procedure is only supported for ring.Type = ring.Standard (X^{2^{i}} does not exist in the sub-ring Z[X + X^{-1}])"))
	}

	if len(cts) < 2 {
		panic(fmt.Errorf("cannot Pack: #cts must be at least 2"))
	}

	keys := utils.GetSortedKeys(cts)

	gap := keys[1] - keys[0]
	level := cts[keys[0]].Level()

	for i, key := range keys[1:] {
		level = utils.Min(level, cts[key].Level())

		if i < len(keys)-1 {
			gap = utils.Min(gap, keys[i+1]-keys[i])
		}
	}

	logN := params.LogN()
	ringQ := params.RingQ().AtLevel(level)

	logStart := logN - inputLogGap
	logEnd := logN // Forces the trace to clean unused slots

	if logStart >= logEnd {
		panic(fmt.Errorf("cannot PackRLWE: gaps between ciphertexts is smaller than inputLogGap > N"))
	}

	xPow2 := genXPow2(ringQ.AtLevel(level), params.LogN(), false) // log(N) polynomial to generate, quick

	NInv := new(big.Int).SetUint64(uint64(1 << (logEnd - logStart)))
	NInv.ModInverse(NInv, ringQ.ModulusAtLevel[level])

	for _, key := range keys {

		ct := cts[key]

		if ct.Degree() != 1 {
			panic(fmt.Errorf("cannot PackRLWE: cts[%d].Degree() != 1", key))
		}

		if !ct.IsNTT {
			ringQ.NTT(ct.Value[0], ct.Value[0])
			ringQ.NTT(ct.Value[1], ct.Value[1])
			ct.IsNTT = true
		}

		ringQ.MulScalarBigint(ct.Value[0], NInv, ct.Value[0])
		ringQ.MulScalarBigint(ct.Value[1], NInv, ct.Value[1])
	}

	tmpa := &Ciphertext{}
	tmpa.Value = []*ring.Poly{ringQ.NewPoly(), ringQ.NewPoly()}
	tmpa.IsNTT = true

	for i := logStart; i < logEnd; i++ {

		t := 1 << (logN - 1 - i)

		for jx, jy := 0, t; jx < t; jx, jy = jx+1, jy+1 {

			a := cts[jx]
			b := cts[jy]

			if b != nil {

				//X^(N/2^L)
				ringQ.MulCoeffsMontgomery(b.Value[0], xPow2[len(xPow2)-i-1], b.Value[0])
				ringQ.MulCoeffsMontgomery(b.Value[1], xPow2[len(xPow2)-i-1], b.Value[1])

				if a != nil {

					// tmpa = phi(a - b * X^{N/2^{i}}, 2^{i-1})
					ringQ.Sub(a.Value[0], b.Value[0], tmpa.Value[0])
					ringQ.Sub(a.Value[1], b.Value[1], tmpa.Value[1])

					// a = a + b * X^{N/2^{i}}
					ringQ.Add(a.Value[0], b.Value[0], a.Value[0])
					ringQ.Add(a.Value[1], b.Value[1], a.Value[1])

				} else {
					// if ct[jx] == nil, then simply re-assigns
					cts[jx] = cts[jy]
				}
			}

			if a != nil {

				var galEl uint64

				if i == 0 {
					galEl = ringQ.NthRoot() - 1
				} else {
					galEl = eval.Parameters().GaloisElementForColumnRotationBy(1 << (i - 1))
				}

				if b != nil {
					eval.Automorphism(tmpa, galEl, tmpa)
				} else {
					eval.Automorphism(a, galEl, tmpa)
				}

				// a + b * X^{N/2^{i}} + phi(a - b * X^{N/2^{i}}, 2^{i-1})
				ringQ.Add(a.Value[0], tmpa.Value[0], a.Value[0])
				ringQ.Add(a.Value[1], tmpa.Value[1], a.Value[1])
			}
		}
	}

	return cts[0]
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
		accQP := &OperandQP{Value: []*ringqp.Poly{&eval.BuffQP[2], &eval.BuffQP[3]}}
		accQP.IsNTT = true

		// Buffer mod QP (i.e. to store the result of lazy gadget products)
		cQP := &OperandQP{Value: []*ringqp.Poly{&eval.BuffQP[4], &eval.BuffQP[5]}}
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
