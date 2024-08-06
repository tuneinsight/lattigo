package rlwe

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
	"github.com/tuneinsight/lattigo/v6/utils"
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
//
// The method will return an error if the input and output ciphertexts degree is not one.
func (eval Evaluator) Trace(ctIn *Ciphertext, logN int, opOut *Ciphertext) (err error) {

	if ctIn.Degree() != 1 || opOut.Degree() != 1 {
		return fmt.Errorf("ctIn.Degree() != 1 or opOut.Degree() != 1")
	}

	params := eval.GetRLWEParameters()

	level := utils.Min(ctIn.Level(), opOut.Level())

	opOut.Resize(opOut.Degree(), level)

	*opOut.MetaData = *ctIn.MetaData

	gap := 1 << (params.LogN() - logN - 1)

	if logN == 0 {
		gap <<= 1
	}

	if gap > 1 {

		ringQ := params.RingQ().AtLevel(level)

		if ringQ.Type() == ring.ConjugateInvariant {
			gap >>= 1 // We skip the last step that applies phi(5^{-1})
		}

		NInv := new(big.Int).SetUint64(uint64(gap))
		NInv.ModInverse(NInv, ringQ.ModulusAtLevel[level])

		// pre-multiplication by (N/n)^-1
		ringQ.MulScalarBigint(ctIn.Value[0], NInv, opOut.Value[0])
		ringQ.MulScalarBigint(ctIn.Value[1], NInv, opOut.Value[1])

		if !ctIn.IsNTT {
			ringQ.NTT(opOut.Value[0], opOut.Value[0])
			ringQ.NTT(opOut.Value[1], opOut.Value[1])
			opOut.IsNTT = true
		}

		buff, err := NewCiphertextAtLevelFromPoly(level, []ring.Poly{eval.BuffQP[3].Q, eval.BuffQP[4].Q})

		// Sanity check, this error should not happen unless the
		// evaluator's buffer has been improperly tempered with.
		if err != nil {
			panic(err)
		}

		buff.IsNTT = true

		for i := logN; i < params.LogN()-1; i++ {

			if err = eval.Automorphism(opOut, params.GaloisElement(1<<i), buff); err != nil {
				return err
			}

			ringQ.Add(opOut.Value[0], buff.Value[0], opOut.Value[0])
			ringQ.Add(opOut.Value[1], buff.Value[1], opOut.Value[1])
		}

		if logN == 0 && ringQ.Type() == ring.Standard {

			if err = eval.Automorphism(opOut, ringQ.NthRoot()-1, buff); err != nil {
				return err
			}

			ringQ.Add(opOut.Value[0], buff.Value[0], opOut.Value[0])
			ringQ.Add(opOut.Value[1], buff.Value[1], opOut.Value[1])
		}

		if !ctIn.IsNTT {
			ringQ.INTT(opOut.Value[0], opOut.Value[0])
			ringQ.INTT(opOut.Value[1], opOut.Value[1])
			opOut.IsNTT = false
		}

	} else {
		if ctIn != opOut {
			opOut.Copy(ctIn)
		}
	}

	return
}

// GaloisElementsForTrace returns the list of Galois elements required for the for the `Trace` operation.
// Trace maps X -> sum((-1)^i * X^{i*n+1}) for 2^{LogN} <= i < N.
func GaloisElementsForTrace(params ParameterProvider, logN int) (galEls []uint64) {

	p := params.GetRLWEParameters()

	galEls = []uint64{}
	for i, j := logN, 0; i < p.LogN()-1; i, j = i+1, j+1 {
		galEls = append(galEls, p.GaloisElement(1<<i))
	}

	if logN == 0 {
		switch p.RingType() {
		case ring.Standard:
			galEls = append(galEls, p.GaloisElementOrderTwoOrthogonalSubgroup())
		case ring.ConjugateInvariant:
			panic("cannot GaloisElementsForTrace: Galois element GaloisGen^-1 is undefined in ConjugateInvariant Ring")
		default:
			panic("cannot GaloisElementsForTrace: invalid ring type")
		}
	}

	return
}

// InnerSum applies an optimized inner sum on the Ciphertext (log2(n) + HW(n) rotations with double hoisting).
// The operation assumes that `ctIn` encrypts Slots/`batchSize` sub-vectors of size `batchSize` and will add them together (in parallel) in groups of `n`.
// It outputs in opOut a [Ciphertext] for which the "leftmost" sub-vector of each group is equal to the sum of the group.
//
// The inner sum is computed in a tree fashion. Example for batchSize=2 & n=4 (garbage slots are marked by 'x'):
//
//  1. [{a, b}, {c, d}, {e, f}, {g, h}, {a, b}, {c, d}, {e, f}, {g, h}]
//
//  2. [{a, b}, {c, d}, {e, f}, {g, h}, {a, b}, {c, d}, {e, f}, {g, h}]
//     +
//     [{c, d}, {e, f}, {g, h}, {x, x}, {c, d}, {e, f}, {g, h}, {x, x}] (rotate batchSize * 2^{0})
//     =
//     [{a+c, b+d}, {x, x}, {e+g, f+h}, {x, x}, {a+c, b+d}, {x, x}, {e+g, f+h}, {x, x}]
//
//  3. [{a+c, b+d}, {x, x}, {e+g, f+h}, {x, x}, {a+c, b+d}, {x, x}, {e+g, f+h}, {x, x}] (rotate batchSize * 2^{1})
//     +
//     [{e+g, f+h}, {x, x}, {x, x}, {x, x}, {e+g, f+h}, {x, x}, {x, x}, {x, x}] =
//     =
//     [{a+c+e+g, b+d+f+h}, {x, x}, {x, x}, {x, x}, {a+c+e+g, b+d+f+h}, {x, x}, {x, x}, {x, x}]
func (eval Evaluator) InnerSum(ctIn *Ciphertext, batchSize, n int, opOut *Ciphertext) (err error) {

	params := eval.GetRLWEParameters()

	levelQ := ctIn.Level()
	levelP := params.PCount() - 1

	ringQP := params.RingQP().AtLevel(ctIn.Level(), levelP)

	ringQ := ringQP.RingQ

	opOut.Resize(opOut.Degree(), levelQ)
	*opOut.MetaData = *ctIn.MetaData

	ctInNTT, err := NewCiphertextAtLevelFromPoly(levelQ, eval.BuffCt.Value[:2])

	// Sanity check, this error should not happen unless the
	// evaluator's buffer thave been improperly tempered with.
	if err != nil {
		panic(err)
	}

	ctInNTT.MetaData = &MetaData{}
	ctInNTT.IsNTT = true

	if !ctIn.IsNTT {
		ringQ.NTT(ctIn.Value[0], ctInNTT.Value[0])
		ringQ.NTT(ctIn.Value[1], ctInNTT.Value[1])
	} else {
		ctInNTT.Value[0].CopyLvl(levelQ, ctIn.Value[0])
		ctInNTT.Value[1].CopyLvl(levelQ, ctIn.Value[1])
	}

	if n == 1 {
		if ctIn != opOut {
			opOut.Value[0].CopyLvl(levelQ, ctIn.Value[0])
			opOut.Value[1].CopyLvl(levelQ, ctIn.Value[1])
		}
	} else {

		// BuffQP[0:2] are used by AutomorphismHoistedLazy

		// Accumulator mod QP (i.e. opOut Mod QP)
		accQP := &Element[ringqp.Poly]{Value: []ringqp.Poly{eval.BuffQP[2], eval.BuffQP[3]}}
		accQP.MetaData = ctInNTT.MetaData

		// Buffer mod QP (i.e. to store the result of lazy gadget products)
		cQP := &Element[ringqp.Poly]{Value: []ringqp.Poly{eval.BuffQP[4], eval.BuffQP[5]}}
		cQP.MetaData = ctInNTT.MetaData

		// Buffer mod Q (i.e. to store the result of gadget products)
		cQ, err := NewCiphertextAtLevelFromPoly(levelQ, []ring.Poly{cQP.Value[0].Q, cQP.Value[1].Q})

		// Sanity check, this error should not happen unless the
		// evaluator's buffer has been improperly tempered with.
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
						opOut.Value[0].CopyLvl(levelQ, ctInNTT.Value[0])
						opOut.Value[1].CopyLvl(levelQ, ctInNTT.Value[1])
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

// InnerFunction applies an user defined function on the [Ciphertext] with a tree-like combination requiring log2(n) + HW(n) rotations.
//
// InnerFunction with f = eval.Add(a, b, c) is equivalent to [Evaluator.InnerSum] (although slightly slower).
//
// The operation assumes that `ctIn` encrypts Slots/`batchSize` sub-vectors of size `batchSize` and will add them together (in parallel) in groups of `n`.
// It outputs in opOut a [Ciphertext] for which the "leftmost" sub-vector of each group is equal to the pair-wise recursive evaluation of
// function over the group.
//
// The inner function is computed in a tree fashion. Example for batchSize=2 & n=4 (garbage slots are marked by 'x'):
//
//  1. [{a, b}, {c, d}, {e, f}, {g, h}, {a, b}, {c, d}, {e, f}, {g, h}]
//
//  2. [{a, b}, {c, d}, {e, f}, {g, h}, {a, b}, {c, d}, {e, f}, {g, h}]
//     f
//     [{c, d}, {e, f}, {g, h}, {x, x}, {c, d}, {e, f}, {g, h}, {x, x}] (rotate batchSize * 2^{0})
//     =
//     [{f(a, c), f(b, d)}, {f(c, e), f(d, f)}, {f(e, g), f(f, h)}, {x, x}, {f(a, c), f(b, d)}, {f(c, e), f(d, f)}, {f(e, g), f(f, h)}, {x, x}]
//
//  3. [{f(a, c), f(b, d)}, {x, x}, {f(e, g), f(f, h)}, {x, x}, {f(a, c), f(b, d)}, {x, x}, {f(e, g), f(f, h)}, {x, x}] (rotate batchSize * 2^{1})
//     +
//     [{f(e, g), f(f, h)}, {x, x}, {x, x}, {x, x}, {f(e, g), f(f, h)}, {x, x}, {x, x}, {x, x}] =
//     =
//     [{f(f(a,c),f(e,g)), f(f(b, d), f(f, h))}, {x, x}, {x, x}, {x, x}, {f(f(a,c),f(e,g)), f(f(b, d), f(f, h))}, {x, x}, {x, x}, {x, x}]
func (eval Evaluator) InnerFunction(ctIn *Ciphertext, batchSize, n int, f func(a, b, c *Ciphertext) (err error), opOut *Ciphertext) (err error) {

	params := eval.GetRLWEParameters()

	levelQ := utils.Min(ctIn.Level(), opOut.Level())

	ringQ := params.RingQ().AtLevel(levelQ)

	opOut.Resize(opOut.Degree(), levelQ)
	*opOut.MetaData = *ctIn.MetaData

	P0 := params.RingQ().NewPoly()
	P1 := params.RingQ().NewPoly()
	P2 := params.RingQ().NewPoly()
	P3 := params.RingQ().NewPoly()

	ctInNTT := NewCiphertext(params, 1, levelQ)

	*ctInNTT.MetaData = *ctIn.MetaData
	ctInNTT.IsNTT = true

	if !ctIn.IsNTT {
		ringQ.NTT(ctIn.Value[0], ctInNTT.Value[0])
		ringQ.NTT(ctIn.Value[1], ctInNTT.Value[1])
	} else {
		ctInNTT.Copy(ctIn)
	}

	if n == 1 {
		opOut.Copy(ctIn)
	} else {

		// Accumulator mod Q
		accQ, err := NewCiphertextAtLevelFromPoly(levelQ, []ring.Poly{P0, P1})
		*accQ.MetaData = *ctInNTT.MetaData

		// Sanity check, this error should not happen unless the
		// evaluator's buffer has been improperly tempered with.
		if err != nil {
			panic(err)
		}

		// Buffer mod Q
		cQ, err := NewCiphertextAtLevelFromPoly(levelQ, []ring.Poly{P2, P3})
		*cQ.MetaData = *ctInNTT.MetaData

		// Sanity check, this error should not happen unless the
		// evaluator's buffer has been improperly tempered with.
		if err != nil {
			panic(err)
		}

		state := false
		copy := true
		// Binary reading of the input n
		for i, j := 0, n; j > 0; i, j = i+1, j>>1 {

			// If the binary reading scans a 1 (j is odd)
			if j&1 == 1 {

				k := n - (n & ((2 << i) - 1))
				k *= batchSize

				// If the rotation is not zero
				if k != 0 {

					rot := params.GaloisElement(k)

					// opOutQ = f(opOutQ, Rotate(ctInNTT, k), opOutQ)
					if copy {
						if err = eval.Automorphism(ctInNTT, rot, accQ); err != nil {
							return err
						}
						copy = false
					} else {
						if err = eval.Automorphism(ctInNTT, rot, cQ); err != nil {
							return err
						}

						if err = f(accQ, cQ, accQ); err != nil {
							return err
						}
					}

					// j is even
				} else {

					state = true

					// if n is not a power of two, then at least one j was odd, and thus the buffer opOutQ is not empty
					if n&(n-1) != 0 {

						opOut.Copy(accQ)

						if err = f(opOut, ctInNTT, opOut); err != nil {
							return err
						}

					} else {
						opOut.Copy(ctInNTT)
					}
				}
			}

			if !state {

				// ctInNTT = f(ctInNTT, Rotate(ctInNTT, 2^i), ctInNTT)
				if err = eval.Automorphism(ctInNTT, params.GaloisElement((1<<i)*batchSize), cQ); err != nil {
					return err
				}

				if err = f(ctInNTT, cQ, ctInNTT); err != nil {
					return err
				}
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
// [Evaluator.InnerSum] operation with parameters batch and n.
func GaloisElementsForInnerSum(params ParameterProvider, batch, n int) (galEls []uint64) {

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

// Replicate applies an optimized replication on the [Ciphertext] (log2(n) + HW(n) rotations with double hoisting).
// It acts as the inverse of a inner sum (summing elements from left to right).
// The replication is parameterized by the size of the sub-vectors to replicate batchSize and
// the number of times n they need to be replicated.
// To ensure correctness, a gap of zero values of size batchSize * (n-1) must exist between
// two consecutive sub-vectors to replicate.
// This method is faster than Replicate when the number of rotations is large and it uses log2(n) + HW(n) instead of n.
func (eval Evaluator) Replicate(ctIn *Ciphertext, batchSize, n int, opOut *Ciphertext) (err error) {
	return eval.InnerSum(ctIn, -batchSize, n, opOut)
}

// GaloisElementsForReplicate returns the list of Galois elements necessary to perform the
// [Evaluator.Replicate] operation with parameters batch and n.
func GaloisElementsForReplicate(params ParameterProvider, batch, n int) (galEls []uint64) {
	return GaloisElementsForInnerSum(params, -batch, n)
}
