package rlwe

import (
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// Operand is a common interface for Ciphertext and Plaintext types.
type Operand interface {
	El() *Ciphertext
	Degree() int
	Level() int
	GetScale() Scale
	SetScale(Scale)
}

// Evaluator is a struct that holds the necessary elements to execute general homomorphic
// operation on RLWE ciphertexts, such as automorphisms, key-switching and relinearization.
type Evaluator struct {
	*evaluatorBase
	*evaluatorBuffers

	Rlk             *RelinearizationKey
	Rtks            *RotationKeySet
	PermuteNTTIndex map[uint64][]uint64

	BasisExtender *ring.BasisExtender
	Decomposer    *ring.Decomposer
}

type evaluatorBase struct {
	params Parameters
}

type evaluatorBuffers struct {
	BuffCt Ciphertext
	// BuffQP[0-0]: Key-Switch on the fly decomp(c2)
	// BuffQP[1-2]: Key-Switch output
	// BuffQP[3-5]: Available
	BuffQP        [6]ringqp.Poly
	BuffInvNTT    *ring.Poly
	BuffDecompQP  []ringqp.Poly // Memory Buff for the basis extension in hoisting
	BuffBitDecomp []uint64
}

func newEvaluatorBase(params Parameters) *evaluatorBase {
	return &evaluatorBase{
		params: params,
	}
}

func newEvaluatorBuffers(params Parameters) *evaluatorBuffers {

	buff := new(evaluatorBuffers)
	decompRNS := params.DecompRNS(params.MaxLevelQ(), params.MaxLevelP())
	ringQP := params.RingQP()

	buff.BuffCt = Ciphertext{Value: []*ring.Poly{ringQP.RingQ.NewPoly(), ringQP.RingQ.NewPoly()}}

	buff.BuffQP = [6]ringqp.Poly{ringQP.NewPoly(), ringQP.NewPoly(), ringQP.NewPoly(), ringQP.NewPoly(), ringQP.NewPoly(), ringQP.NewPoly()}

	buff.BuffInvNTT = params.RingQ().NewPoly()

	buff.BuffDecompQP = make([]ringqp.Poly, decompRNS)
	for i := 0; i < decompRNS; i++ {
		buff.BuffDecompQP[i] = ringQP.NewPoly()
	}

	buff.BuffBitDecomp = make([]uint64, params.RingQ().N())

	return buff
}

// NewEvaluator creates a new Evaluator.
func NewEvaluator(params Parameters, evaluationKey *EvaluationKey) (eval *Evaluator) {
	eval = new(Evaluator)
	eval.evaluatorBase = newEvaluatorBase(params)
	eval.evaluatorBuffers = newEvaluatorBuffers(params)

	if params.RingP() != nil {
		eval.BasisExtender = ring.NewBasisExtender(params.RingQ(), params.RingP())
		eval.Decomposer = ring.NewDecomposer(params.RingQ(), params.RingP())
	}

	if evaluationKey != nil {
		if evaluationKey.Rlk != nil {
			eval.Rlk = evaluationKey.Rlk
		}

		if evaluationKey.Rtks != nil {
			eval.Rtks = evaluationKey.Rtks
			eval.PermuteNTTIndex = *eval.permuteNTTIndexesForKey(eval.Rtks)
		}
	}

	return
}

// Parameters returns the parameters used to instantiate the target evaluator.
func (eval *Evaluator) Parameters() Parameters {
	return eval.params
}

// permuteNTTIndexesForKey generates permutation indexes for automorphisms for ciphertexts
// that are given in the NTT domain.
func (eval *Evaluator) permuteNTTIndexesForKey(rtks *RotationKeySet) *map[uint64][]uint64 {
	if rtks == nil {
		return &map[uint64][]uint64{}
	}
	permuteNTTIndex := make(map[uint64][]uint64, len(rtks.Keys))
	for galEl := range rtks.Keys {
		permuteNTTIndex[galEl] = eval.params.RingQ().PermuteNTTIndex(galEl)
	}
	return &permuteNTTIndex
}

// ShallowCopy creates a shallow copy of this Evaluator in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Evaluators can be used concurrently.
func (eval *Evaluator) ShallowCopy() *Evaluator {
	return &Evaluator{
		evaluatorBase:    eval.evaluatorBase,
		Decomposer:       eval.Decomposer,
		BasisExtender:    eval.BasisExtender.ShallowCopy(),
		evaluatorBuffers: newEvaluatorBuffers(eval.params),
		Rlk:              eval.Rlk,
		Rtks:             eval.Rtks,
		PermuteNTTIndex:  eval.PermuteNTTIndex,
	}
}

// WithKey creates a shallow copy of the receiver Evaluator for which the new EvaluationKey is evaluationKey
// and where the temporary buffers are shared. The receiver and the returned Evaluators cannot be used concurrently.
func (eval *Evaluator) WithKey(evaluationKey *EvaluationKey) *Evaluator {
	var indexes map[uint64][]uint64
	if evaluationKey.Rtks == eval.Rtks {
		indexes = eval.PermuteNTTIndex
	} else {
		indexes = *eval.permuteNTTIndexesForKey(evaluationKey.Rtks)
	}
	return &Evaluator{
		evaluatorBase:    eval.evaluatorBase,
		evaluatorBuffers: eval.evaluatorBuffers,
		Decomposer:       eval.Decomposer,
		BasisExtender:    eval.BasisExtender,
		Rlk:              evaluationKey.Rlk,
		Rtks:             evaluationKey.Rtks,
		PermuteNTTIndex:  indexes,
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

	params := eval.params

	levelQ := ctIn.Level()

	ringQ := params.RingQ().AtLevel(levelQ)

	// Compute X^{-2^{i}} from 1 to LogN
	xPow2 := genXPow2(ringQ.AtLevel(levelQ), logN, true)

	ctOut = make([]*Ciphertext, 1<<(logN-logGap))
	ctOut[0] = ctIn.CopyNew()

	// Multiplies by 2^{-logN} mod Q
	v0, v1 := ctOut[0].Value[0], ctOut[0].Value[1]
	for i, s := range ringQ.SubRings[:levelQ+1] {

		NInv := ring.MForm(ring.ModExp(1<<logN, s.Modulus-2, s.Modulus), s.Modulus, s.BRedConstant)

		s.MulScalarMontgomery(v0.Coeffs[i], NInv, v0.Coeffs[i])
		s.MulScalarMontgomery(v1.Coeffs[i], NInv, v1.Coeffs[i])
	}

	gap := 1 << logGap

	tmp := NewCiphertext(params, 1, levelQ)

	for i := 0; i < logN; i++ {

		galEl := uint64(ringQ.N()/(1<<i) + 1)

		for j := 0; j < (1 << i); j += gap {

			c0 := ctOut[j/gap]

			// X -> X^{N/2^{i} + 1}
			eval.Automorphism(c0, galEl, tmp)

			c1 := c0.CopyNew()

			// Zeroes odd coeffs: [a, b, c, d] -> [2a, 0, 2b, 0]
			ringQ.Add(c0.Value[0], tmp.Value[0], c0.Value[0])
			ringQ.Add(c0.Value[1], tmp.Value[1], c0.Value[1])

			if (j+(1<<i))/gap > 0 {

				// Zeroes even coeffs: [a, b, c, d] -> [0, 2b, 0, 2d]
				ringQ.Sub(c1.Value[0], tmp.Value[0], c1.Value[0])
				ringQ.Sub(c1.Value[1], tmp.Value[1], c1.Value[1])

				// c1 * X^{-2^{i}}: [0, 2b, 0, 2d] -> [2b, 0, 2d, 0]
				ringQ.MulCoeffsMontgomery(c1.Value[0], xPow2[i], c1.Value[0])
				ringQ.MulCoeffsMontgomery(c1.Value[1], xPow2[i], c1.Value[1])

				ctOut[(j+(1<<i))/gap] = c1
			}
		}
	}

	return
}

// Merge merges a batch of RLWE, packing the first coefficient of each RLWE into a single RLWE.
// The operation will require N/gap + log(gap) key-switches, where gap is the minimum gap between
// two non-zero coefficients of the final Ciphertext.
// The method takes as input a map of Ciphertext, indexing in which coefficient of the final
// Ciphertext the first coefficient of each Ciphertext of the map must be packed.
// All input ciphertexts must be in the NTT domain; otherwise, the method will panic.
func (eval *Evaluator) Merge(ctIn map[int]*Ciphertext) (ctOut *Ciphertext) {

	params := eval.params
	ringQ := params.RingQ()

	var levelQ int
	for i := range ctIn {
		levelQ = ctIn[i].Level()
		break
	}

	for i := range ctIn {
		levelQ = utils.MinInt(levelQ, ctIn[i].Level())
	}

	xPow2 := genXPow2(ringQ.AtLevel(levelQ), params.LogN(), false)

	// Multiplies by (Slots * N) ^-1 mod Q
	for i := range ctIn {
		if ctIn[i] != nil {

			if !ctIn[i].IsNTT {
				panic("canot Merge: all ctIn must be in the NTT domain")
			}

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

	var tmpEven *Ciphertext
	if ctEven != nil {
		tmpEven = ctEven.CopyNew()
	}

	var level = 0xFFFF // Case if ctOdd == nil

	if ctOdd != nil {
		level = ctOdd.Level()
	}

	if ctEven != nil {
		level = utils.MinInt(level, ctEven.Level())
	}

	ringQ := eval.params.RingQ().AtLevel(level)

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

	if n == 1 {
		if ctIn != ctOut {
			ring.CopyLvl(levelQ, ctIn.Value[0], ctOut.Value[0])
			ring.CopyLvl(levelQ, ctIn.Value[1], ctOut.Value[1])
		}
	} else {

		c0OutQP := eval.BuffQP[2]
		c1OutQP := eval.BuffQP[3]

		cQP := CiphertextQP{Value: [2]ringqp.Poly{eval.BuffQP[4], eval.BuffQP[5]}}
		cQP.IsNTT = true

		// Memory buffer for ctIn = ctIn + rot(ctIn, 2^i) in Q
		tmpct := NewCiphertextAtLevelFromPoly(levelQ, eval.BuffCt.Value[:2])
		tmpct.IsNTT = true

		ctqp := NewCiphertextAtLevelFromPoly(levelQ, []*ring.Poly{cQP.Value[0].Q, cQP.Value[1].Q})
		ctqp.IsNTT = true

		state := false
		copy := true
		// Binary reading of the input n
		for i, j := 0, n; j > 0; i, j = i+1, j>>1 {

			// Starts by decomposing the input ciphertext
			if i == 0 {
				// If first iteration, then copies directly from the input ciphertext that hasn't been rotated
				eval.DecomposeNTT(levelQ, levelP, levelP+1, ctIn.Value[1], true, eval.BuffDecompQP)
			} else {
				// Else copies from the rotated input ciphertext
				eval.DecomposeNTT(levelQ, levelP, levelP+1, tmpct.Value[1], true, eval.BuffDecompQP)
			}

			// If the binary reading scans a 1
			if j&1 == 1 {

				k := n - (n & ((2 << i) - 1))
				k *= batchSize

				// If the rotation is not zero
				if k != 0 {

					// Rotate((tmpc0, tmpc1), k)
					if i == 0 {
						eval.AutomorphismHoistedLazy(levelQ, ctIn.Value[0], eval.BuffDecompQP, eval.params.GaloisElementForColumnRotationBy(k), cQP)
					} else {
						eval.AutomorphismHoistedLazy(levelQ, tmpct.Value[0], eval.BuffDecompQP, eval.params.GaloisElementForColumnRotationBy(k), cQP)
					}

					// ctOut += Rotate((tmpc0, tmpc1), k)
					if copy {
						ringqp.CopyLvl(levelQ, levelP, cQP.Value[0], c0OutQP)
						ringqp.CopyLvl(levelQ, levelP, cQP.Value[1], c1OutQP)
						copy = false
					} else {
						ringQP.Add(c0OutQP, cQP.Value[0], c0OutQP)
						ringQP.Add(c1OutQP, cQP.Value[1], c1OutQP)
					}
				} else {

					state = true

					// if n is not a power of two
					if n&(n-1) != 0 {

						eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c0OutQP.Q, c0OutQP.P, c0OutQP.Q) // Division by P
						eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c1OutQP.Q, c1OutQP.P, c1OutQP.Q) // Division by P

						// ctOut += (tmpc0, tmpc1)
						ringQ.Add(c0OutQP.Q, tmpct.Value[0], ctOut.Value[0])
						ringQ.Add(c1OutQP.Q, tmpct.Value[1], ctOut.Value[1])

					} else {

						ring.CopyLvl(levelQ, tmpct.Value[0], ctOut.Value[0])
						ring.CopyLvl(levelQ, tmpct.Value[1], ctOut.Value[1])
					}
				}
			}

			if !state {

				rot := eval.params.GaloisElementForColumnRotationBy((1 << i) * batchSize)
				if i == 0 {

					eval.AutomorphismHoisted(levelQ, ctIn, eval.BuffDecompQP, rot, tmpct)

					ringQ.Add(tmpct.Value[0], ctIn.Value[0], tmpct.Value[0])
					ringQ.Add(tmpct.Value[1], ctIn.Value[1], tmpct.Value[1])
				} else {
					// (tmpc0, tmpc1) = Rotate((tmpc0, tmpc1), 2^i)
					eval.AutomorphismHoisted(levelQ, tmpct, eval.BuffDecompQP, rot, ctqp)
					ringQ.Add(tmpct.Value[0], cQP.Value[0].Q, tmpct.Value[0])
					ringQ.Add(tmpct.Value[1], cQP.Value[1].Q, tmpct.Value[1])
				}
			}
		}
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
