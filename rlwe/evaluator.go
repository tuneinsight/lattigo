package rlwe

import (
	"math/bits"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
)

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
	ev := new(evaluatorBase)
	ev.params = params
	return ev
}

func newEvaluatorBuffers(params Parameters) *evaluatorBuffers {

	buff := new(evaluatorBuffers)
	decompRNS := params.DecompRNS(params.QCount()-1, params.PCount()-1)
	ringQP := params.RingQP()

	buff.BuffCt = Ciphertext{Value: []*ring.Poly{ringQP.RingQ.NewPoly(), ringQP.RingQ.NewPoly()}}

	buff.BuffQP = [6]ringqp.Poly{ringQP.NewPoly(), ringQP.NewPoly(), ringQP.NewPoly(), ringQP.NewPoly(), ringQP.NewPoly(), ringQP.NewPoly()}

	buff.BuffInvNTT = params.RingQ().NewPoly()

	buff.BuffDecompQP = make([]ringqp.Poly, decompRNS)
	for i := 0; i < decompRNS; i++ {
		buff.BuffDecompQP[i] = ringQP.NewPoly()
	}

	buff.BuffBitDecomp = make([]uint64, params.RingQ().N)

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

// DecomposeNTT applies the full RNS basis decomposition on c2.
// Expects the IsNTT flag of c2 to correctly reflect the domain of c2.
// BuffQPDecompQ and BuffQPDecompQ are vectors of polynomials (mod Q and mod P) that store the
// special RNS decomposition of c2 (in the NTT domain)
func (eval *Evaluator) DecomposeNTT(levelQ, levelP, nbPi int, c2 *ring.Poly, BuffDecompQP []ringqp.Poly) {

	ringQ := eval.params.RingQ()

	var polyNTT, polyInvNTT *ring.Poly

	if c2.IsNTT {
		polyNTT = c2
		polyInvNTT = eval.BuffInvNTT
		ringQ.InvNTTLvl(levelQ, polyNTT, polyInvNTT)
	} else {
		polyNTT = eval.BuffInvNTT
		polyInvNTT = c2
		ringQ.NTTLvl(levelQ, polyInvNTT, polyNTT)
	}

	decompRNS := eval.params.DecompRNS(levelQ, levelP)
	for i := 0; i < decompRNS; i++ {
		eval.DecomposeSingleNTT(levelQ, levelP, nbPi, i, polyNTT, polyInvNTT, BuffDecompQP[i].Q, BuffDecompQP[i].P)
	}
}

// DecomposeSingleNTT takes the input polynomial c2 (c2NTT and c2InvNTT, respectively in the NTT and out of the NTT domain)
// modulo the RNS basis, and returns the result on c2QiQ are c2QiP the receiver polynomials respectively mod Q and mod P (in the NTT domain)
func (eval *Evaluator) DecomposeSingleNTT(levelQ, levelP, nbPi, decompRNS int, c2NTT, c2InvNTT, c2QiQ, c2QiP *ring.Poly) {

	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()

	eval.Decomposer.DecomposeAndSplit(levelQ, levelP, nbPi, decompRNS, c2InvNTT, c2QiQ, c2QiP)

	p0idxst := decompRNS * nbPi
	p0idxed := p0idxst + nbPi

	// c2_qi = cx mod qi mod qi
	for x := 0; x < levelQ+1; x++ {
		if p0idxst <= x && x < p0idxed {
			copy(c2QiQ.Coeffs[x], c2NTT.Coeffs[x])
		} else {
			ringQ.NTTSingle(x, c2QiQ.Coeffs[x], c2QiQ.Coeffs[x])
		}
	}

	if ringP != nil {
		// c2QiP = c2 mod qi mod pj
		ringP.NTTLvl(levelP, c2QiP, c2QiP)
	}
}

// ExpandRLWE expands a RLWE ciphertext encrypting sum ai * X^i to 2^logN ciphertexts,
// each encrypting ai * X^0 for 0 <= i < 2^LogN. That is, it extracts the first 2^logN
// coefficients of ctIn and returns a RLWE ciphetext for each coefficient extracted.
func (eval *Evaluator) ExpandRLWE(ctIn *Ciphertext, logN int) (ctOut []*Ciphertext) {

	params := eval.params
	ringQ := params.RingQ()

	levelQ := ctIn.Level()

	// Compute X^{-2^{i}} from 1 to LogN
	xPow2 := genXPow2(ringQ, levelQ, logN, true)

	ctOut = make([]*Ciphertext, 1<<logN)
	ctOut[0] = ctIn.CopyNew()

	Q := ringQ.Modulus
	mredParams := ringQ.MredParams
	bredParams := ringQ.BredParams

	// Multiplies by 2^{-logN} mod Q
	v0, v1 := ctOut[0].Value[0], ctOut[0].Value[1]
	for i := 0; i < levelQ+1; i++ {
		NInv := ring.MForm(ring.ModExp(1<<logN, Q[i]-2, Q[i]), Q[i], bredParams[i])
		ring.MulScalarMontgomeryVec(v0.Coeffs[i], v0.Coeffs[i], NInv, Q[i], mredParams[i])
		ring.MulScalarMontgomeryVec(v1.Coeffs[i], v1.Coeffs[i], NInv, Q[i], mredParams[i])
	}

	tmp := NewCiphertextNTT(params, 1, levelQ)

	for i := 0; i < logN; i++ {

		galEl := uint64(ringQ.N/(1<<i) + 1)

		for j := 0; j < (1 << i); j++ {

			c0 := ctOut[j]

			// X -> X^{N/2^{i} + 1}
			eval.Automorphism(c0, galEl, tmp)

			c1 := c0.CopyNew()

			// Zeroes odd coeffs: [a, b, c, d] -> [2a, 0, 2b, 0]
			ringQ.AddLvl(levelQ, c0.Value[0], tmp.Value[0], c0.Value[0])
			ringQ.AddLvl(levelQ, c0.Value[1], tmp.Value[1], c0.Value[1])

			// Zeroes even coeffs: [a, b, c, d] -> [0, 2b, 0, 2d]
			ringQ.SubLvl(levelQ, c1.Value[0], tmp.Value[0], c1.Value[0])
			ringQ.SubLvl(levelQ, c1.Value[1], tmp.Value[1], c1.Value[1])

			// c1 * X^{-2^{i}}: [0, 2b, 0, 2d] -> [2b, 0, 2d, 0]
			ringQ.MulCoeffsMontgomeryLvl(levelQ, c1.Value[0], xPow2[i], c1.Value[0])
			ringQ.MulCoeffsMontgomeryLvl(levelQ, c1.Value[1], xPow2[i], c1.Value[1])

			ctOut[j+(1<<i)] = c1
		}
	}

	return
}

// MergeRLWE merges a batch of RLWE, packing the first coefficient of each RLWE into a single RLWE.
// The operation will require N/gap + log(gap) key-switches, where gap is the minimum gap between
// two non-zero coefficients of the final Ciphertext.
// The method takes as input a map of Ciphertext, indexing in which coefficient of the final
// Ciphertext the first coefficient of each Ciphertext of the map must be packed.
func (eval *Evaluator) MergeRLWE(ctIn map[int]*Ciphertext) (ctOut *Ciphertext) {

	params := eval.params
	ringQ := params.RingQ()

	var levelQ int
	for i := range ctIn {
		levelQ = ctIn[i].Level()
		break
	}

	xPow2 := genXPow2(ringQ, levelQ, params.LogN(), false)

	NInv := ringQ.NttNInv
	Q := ringQ.Modulus
	mredParams := ringQ.MredParams

	// Multiplies by (Slots * N) ^-1 mod Q
	for i := range ctIn {
		if ctIn[i] != nil {
			v0, v1 := ctIn[i].Value[0], ctIn[i].Value[1]
			for j := 0; j < levelQ+1; j++ {
				ring.MulScalarMontgomeryVec(v0.Coeffs[j], v0.Coeffs[j], NInv[j], Q[j], mredParams[j])
				ring.MulScalarMontgomeryVec(v1.Coeffs[j], v1.Coeffs[j], NInv[j], Q[j], mredParams[j])
			}
		}
	}

	ciphertextslist := make([]*Ciphertext, ringQ.N)

	for i := range ctIn {
		ciphertextslist[i] = ctIn[i]
	}

	if ciphertextslist[0] == nil {
		ciphertextslist[0] = NewCiphertextNTT(params, 1, levelQ)
	}

	return eval.mergeRLWERecurse(ciphertextslist, xPow2)
}

func (eval *Evaluator) mergeRLWERecurse(ciphertexts []*Ciphertext, xPow []*ring.Poly) *Ciphertext {

	ringQ := eval.params.RingQ()

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

	// ctOdd * X^(N/2^L)
	if ctOdd != nil {

		level := ctOdd.Level()

		//X^(N/2^L)
		ringQ.MulCoeffsMontgomeryLvl(level, ctOdd.Value[0], xPow[len(xPow)-L], ctOdd.Value[0])
		ringQ.MulCoeffsMontgomeryLvl(level, ctOdd.Value[1], xPow[len(xPow)-L], ctOdd.Value[1])

		if ctEven != nil {
			// ctEven + ctOdd * X^(N/2^L)
			ringQ.AddLvl(level, ctEven.Value[0], ctOdd.Value[0], ctEven.Value[0])
			ringQ.AddLvl(level, ctEven.Value[1], ctOdd.Value[1], ctEven.Value[1])

			// phi(ctEven - ctOdd * X^(N/2^L), 2^(L-2))
			ringQ.SubLvl(level, tmpEven.Value[0], ctOdd.Value[0], tmpEven.Value[0])
			ringQ.SubLvl(level, tmpEven.Value[1], ctOdd.Value[1], tmpEven.Value[1])
		}
	}

	if ctEven != nil {

		level := ctEven.Level()

		// if L-2 == -1, then gal = -1
		if L == 1 {
			eval.Automorphism(tmpEven, uint64(2*ringQ.N-1), tmpEven)
		} else {
			eval.Automorphism(tmpEven, eval.params.GaloisElementForColumnRotationBy(1<<(L-2)), tmpEven)
		}

		// ctEven + ctOdd * X^(N/2^L) + phi(ctEven - ctOdd * X^(N/2^L), 2^(L-2))
		ringQ.AddLvl(level, ctEven.Value[0], tmpEven.Value[0], ctEven.Value[0])
		ringQ.AddLvl(level, ctEven.Value[1], tmpEven.Value[1], ctEven.Value[1])
	}

	return ctEven
}

func genXPow2(r *ring.Ring, levelQ, logN int, div bool) (xPow []*ring.Poly) {

	// Compute X^{-n} from 0 to LogN
	xPow = make([]*ring.Poly, logN)

	var idx int
	for i := 0; i < logN; i++ {

		idx = 1 << i

		if div {
			idx = r.N - idx
		}

		xPow[i] = r.NewPoly()

		if i == 0 {

			for j := 0; j < levelQ+1; j++ {
				xPow[i].Coeffs[j][idx] = ring.MForm(1, r.Modulus[j], r.BredParams[j])
			}

			r.NTTLvl(levelQ, xPow[i], xPow[i])

		} else {
			r.MulCoeffsMontgomeryLvl(levelQ, xPow[i-1], xPow[i-1], xPow[i]) // X^{n} = X^{1} * X^{n-1}
		}
	}

	if div {
		r.NegLvl(levelQ, xPow[0], xPow[0])
	}

	return
}
