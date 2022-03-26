package rlwe

import (
	"fmt"
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v3/utils"
	"math/bits"
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
	// PoolQ[0]/PoolP[0] : on the fly decomp(c2)
	// PoolQ[1-5]/PoolP[1-5] : available
	Pool          [6]ringqp.Poly
	PoolInvNTT    *ring.Poly
	PoolDecompQP  []ringqp.Poly // Memory pool for the basis extension in hoisting
	PoolBitDecomp []uint64
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

	buff.Pool = [6]ringqp.Poly{ringQP.NewPoly(), ringQP.NewPoly(), ringQP.NewPoly(), ringQP.NewPoly(), ringQP.NewPoly(), ringQP.NewPoly()}

	buff.PoolInvNTT = params.RingQ().NewPoly()

	buff.PoolDecompQP = make([]ringqp.Poly, decompRNS)
	for i := 0; i < decompRNS; i++ {
		buff.PoolDecompQP[i] = ringQP.NewPoly()
	}

	buff.PoolBitDecomp = make([]uint64, params.RingQ().N)

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

// permuteNTTIndexesForKey generates pemutation indexes for automorphisms for ciphertext
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

// Automorphism computes phi(ct), where phi is the map X -> X^galEl. The method requires
// that the corresponding RotationKey has been added to the Evaluator. The method will
// panic if either ctIn or ctOut degree is not equal to 1.
func (eval *Evaluator) Automorphism(ctIn *Ciphertext, galEl uint64, ctOut *Ciphertext) {

	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		panic("cannot apply Automorphism: input and output Ciphertext must be of degree 1")
	}

	if galEl == 1 {
		if ctOut != ctIn {
			ctOut.Copy(ctIn)
		}
		return
	}

	rtk, generated := eval.Rtks.GetRotationKey(galEl)
	if !generated {
		panic(fmt.Sprintf("galEl key 5^%d missing", eval.params.InverseGaloisElement(galEl)))
	}

	level := utils.MinInt(ctIn.Level(), ctOut.Level())

	ringQ := eval.params.RingQ()

	eval.SwitchKeysInPlace(level, ctIn.Value[1], rtk, eval.Pool[1].Q, eval.Pool[2].Q)
	ringQ.AddLvl(level, eval.Pool[1].Q, ctIn.Value[0], eval.Pool[1].Q)

	if ctIn.Value[0].IsNTT {
		ringQ.PermuteNTTWithIndexLvl(level, eval.Pool[1].Q, eval.PermuteNTTIndex[galEl], ctOut.Value[0])
		ringQ.PermuteNTTWithIndexLvl(level, eval.Pool[2].Q, eval.PermuteNTTIndex[galEl], ctOut.Value[1])
	} else {
		ringQ.PermuteLvl(level, eval.Pool[1].Q, galEl, ctOut.Value[0])
		ringQ.PermuteLvl(level, eval.Pool[2].Q, galEl, ctOut.Value[1])
	}

	ctOut.Value[0].Coeffs = ctOut.Value[0].Coeffs[:level+1]
	ctOut.Value[1].Coeffs = ctOut.Value[1].Coeffs[:level+1]
}

// AutomorphismHoisted is similar to Automorphism, except that it takes as input ctIn and c1DecompQP, where c1DecompQP is the RNS
// decomposition of its element of degree 1. This decomposition can be obtained with DecomposeNTT.
// The method requires that the corresponding RotationKey  has been added to the Evaluator.
// The method will panic if either ctIn or ctOut degree is not equal to 1.
func (eval *Evaluator) AutomorphismHoisted(level int, ctIn *Ciphertext, c1DecompQP []ringqp.Poly, galEl uint64, ctOut *Ciphertext) {

	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		panic("cannot apply AutomorphismHoisted: input and output Ciphertext must be of degree 1")
	}

	if galEl == 1 {
		if ctIn != ctOut {
			ctOut.Copy(ctIn)
		}
		return
	}

	rtk, generated := eval.Rtks.GetRotationKey(galEl)
	if !generated {
		panic(fmt.Sprintf("galEl key 5^%d missing", eval.params.InverseGaloisElement(galEl)))
	}

	ringQ := eval.params.RingQ()

	eval.KeyswitchHoisted(level, c1DecompQP, rtk, eval.Pool[0].Q, eval.Pool[1].Q, eval.Pool[0].P, eval.Pool[1].P)
	ringQ.AddLvl(level, eval.Pool[0].Q, ctIn.Value[0], eval.Pool[0].Q)

	if ctIn.Value[0].IsNTT {
		ringQ.PermuteNTTWithIndexLvl(level, eval.Pool[0].Q, eval.PermuteNTTIndex[galEl], ctOut.Value[0])
		ringQ.PermuteNTTWithIndexLvl(level, eval.Pool[1].Q, eval.PermuteNTTIndex[galEl], ctOut.Value[1])
	} else {
		ringQ.PermuteLvl(level, eval.Pool[0].Q, galEl, ctOut.Value[0])
		ringQ.PermuteLvl(level, eval.Pool[1].Q, galEl, ctOut.Value[1])
	}

	ctOut.Value[0].Coeffs = ctOut.Value[0].Coeffs[:level+1]
	ctOut.Value[1].Coeffs = ctOut.Value[1].Coeffs[:level+1]
}

// AutomorphismHoistedNoModDown is similar to AutomorphismHoisted, except that it returns a ciphertext modulo QP and scaled by P.
// The method requires that the corresponding RotationKey has been added to the Evaluator.The method will panic if either ctIn or ctOut degree is not equal to 1.
func (eval *Evaluator) AutomorphismHoistedNoModDown(levelQ int, c0 *ring.Poly, c1DecompQP []ringqp.Poly, galEl uint64, ct0OutQ, ct1OutQ, ct0OutP, ct1OutP *ring.Poly) {

	levelP := eval.params.PCount() - 1

	rtk, generated := eval.Rtks.GetRotationKey(galEl)
	if !generated {
		panic(fmt.Sprintf("galEl key 5^%d missing", eval.params.InverseGaloisElement(galEl)))
	}

	eval.KeyswitchHoistedNoModDown(levelQ, c1DecompQP, rtk, eval.Pool[0].Q, eval.Pool[1].Q, eval.Pool[0].P, eval.Pool[1].P)

	ringQ := eval.params.RingQ()

	if c0.IsNTT {

		index := eval.PermuteNTTIndex[galEl]

		ringQ.PermuteNTTWithIndexLvl(levelQ, eval.Pool[1].Q, index, ct1OutQ)
		ringQ.PermuteNTTWithIndexLvl(levelP, eval.Pool[1].P, index, ct1OutP)

		ringQ.MulScalarBigintLvl(levelQ, c0, eval.params.RingP().ModulusBigint, eval.Pool[1].Q)
		ringQ.AddLvl(levelQ, eval.Pool[0].Q, eval.Pool[1].Q, eval.Pool[0].Q)

		ringQ.PermuteNTTWithIndexLvl(levelQ, eval.Pool[0].Q, index, ct0OutQ)
		ringQ.PermuteNTTWithIndexLvl(levelP, eval.Pool[0].P, index, ct0OutP)
	} else {
		ringQ.PermuteLvl(levelQ, eval.Pool[1].Q, galEl, ct1OutQ)
		ringQ.PermuteLvl(levelP, eval.Pool[1].P, galEl, ct1OutP)

		ringQ.MulScalarBigintLvl(levelQ, c0, eval.params.RingP().ModulusBigint, eval.Pool[1].Q)
		ringQ.AddLvl(levelQ, eval.Pool[0].Q, eval.Pool[1].Q, eval.Pool[0].Q)

		ringQ.PermuteLvl(levelQ, eval.Pool[0].Q, galEl, ct0OutQ)
		ringQ.PermuteLvl(levelP, eval.Pool[0].P, galEl, ct0OutP)
	}
}

// SwitchKeys re-encrypts ctIn under a different key and returns the result in ctOut.
// It requires a SwitchingKey, which is computed from the key under which the Ciphertext is currently encrypted,
// and the key under which the Ciphertext will be re-encrypted.
// The method will panic if either ctIn or ctOut degree isn't 1.
func (eval *Evaluator) SwitchKeys(ctIn *Ciphertext, switchingKey *SwitchingKey, ctOut *Ciphertext) {

	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		panic("cannot SwitchKeys: input and output Ciphertext must be of degree 1")
	}

	level := utils.MinInt(ctIn.Level(), ctOut.Level())
	ringQ := eval.params.RingQ()

	eval.SwitchKeysInPlace(level, ctIn.Value[1], switchingKey, eval.Pool[1].Q, eval.Pool[2].Q)

	ringQ.AddLvl(level, ctIn.Value[0], eval.Pool[1].Q, ctOut.Value[0])
	ring.CopyValuesLvl(level, eval.Pool[2].Q, ctOut.Value[1])
}

// Relinearize applies the relinearization procedure on ct0 and returns the result in ctOut.
// The method will panic if the corresponding relinearization key to the ciphertext degree
// is missing.
func (eval *Evaluator) Relinearize(ctIn *Ciphertext, ctOut *Ciphertext) {
	if eval.Rlk == nil || ctIn.Degree()-1 > len(eval.Rlk.Keys) {
		panic("cannot Relinearize: relinearization key missing (or ciphertext degree is too large)")
	}

	level := utils.MinInt(ctIn.Level(), ctOut.Level())

	ringQ := eval.params.RingQ()

	eval.SwitchKeysInPlace(level, ctIn.Value[2], eval.Rlk.Keys[0], eval.Pool[1].Q, eval.Pool[2].Q)
	ringQ.AddLvl(level, ctIn.Value[0], eval.Pool[1].Q, ctOut.Value[0])
	ringQ.AddLvl(level, ctIn.Value[1], eval.Pool[2].Q, ctOut.Value[1])

	for deg := ctIn.Degree() - 1; deg > 1; deg-- {
		eval.SwitchKeysInPlace(level, ctIn.Value[deg], eval.Rlk.Keys[deg-2], eval.Pool[1].Q, eval.Pool[2].Q)
		ringQ.AddLvl(level, ctOut.Value[0], eval.Pool[1].Q, ctOut.Value[0])
		ringQ.AddLvl(level, ctOut.Value[1], eval.Pool[2].Q, ctOut.Value[1])
	}

	ctOut.Value = ctOut.Value[:2]

	ctOut.Value[0].Coeffs = ctOut.Value[0].Coeffs[:level+1]
	ctOut.Value[1].Coeffs = ctOut.Value[1].Coeffs[:level+1]
}

// MergeRLWE merges a batch of RLWE, packing the first coefficient of each RLWE into a single RLWE.
// The operation will require N/gap + log(gap) key-switches, where gap is the minimum gap between
// two non-zero coefficients of the final ciphertext.
// The method takes as input a map of Ciphertexts, indexing in which coefficient, of the final
// ciphertext, the first coefficient of each ciphertext of the map must be packed.
func (eval *Evaluator) MergeRLWE(ciphertexts map[int]*Ciphertext) (ciphertext *Ciphertext) {

	params := eval.params
	ringQ := params.RingQ()

	// Compute X^{n} from 0 to LogN LUT
	xPow := make([]*ring.Poly, params.LogN())
	for i := 0; i < params.LogN(); i++ {
		xPow[i] = ringQ.NewPoly()
		if i == 0 {
			for j := 0; j < params.MaxLevel()+1; j++ {
				xPow[i].Coeffs[j][1<<i] = ring.MForm(1, ringQ.Modulus[j], ringQ.BredParams[j])
			}
			ringQ.NTT(xPow[i], xPow[i])
		} else {
			ringQ.MulCoeffsMontgomery(xPow[i-1], xPow[i-1], xPow[i]) // X^{n} = X^{1} * X^{n-1}
		}
	}

	var level int
	for i := range ciphertexts {
		level = ciphertexts[i].Level()
		break
	}

	NInv := ringQ.NttNInv
	Q := ringQ.Modulus
	mredParams := ringQ.MredParams

	// Multiplies by (Slots * N) ^-1 mod Q
	for i := range ciphertexts {
		if ciphertexts[i] != nil {
			v0, v1 := ciphertexts[i].Value[0], ciphertexts[i].Value[1]
			for j := 0; j < level+1; j++ {
				ring.MulScalarMontgomeryVec(v0.Coeffs[j], v0.Coeffs[j], NInv[j], Q[j], mredParams[j])
				ring.MulScalarMontgomeryVec(v1.Coeffs[j], v1.Coeffs[j], NInv[j], Q[j], mredParams[j])
			}
		}
	}

	ciphertextslist := make([]*Ciphertext, ringQ.N)

	for i := range ciphertexts {
		ciphertextslist[i] = ciphertexts[i]
	}

	if ciphertextslist[0] == nil {
		ciphertextslist[0] = NewCiphertextNTT(params, 1, level)
	}

	return eval.mergeRLWERecurse(ciphertextslist, xPow)
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
			eval.rotate(tmpEven, uint64(2*ringQ.N-1), tmpEven)
		} else {
			eval.rotate(tmpEven, eval.params.GaloisElementForColumnRotationBy(1<<(L-2)), tmpEven)
		}

		// ctEven + ctOdd * X^(N/2^L) + phi(ctEven - ctOdd * X^(N/2^L), 2^(L-2))
		ringQ.AddLvl(level, ctEven.Value[0], tmpEven.Value[0], ctEven.Value[0])
		ringQ.AddLvl(level, ctEven.Value[1], tmpEven.Value[1], ctEven.Value[1])
	}

	return ctEven
}

func (eval *Evaluator) rotate(ctIn *Ciphertext, galEl uint64, ctOut *Ciphertext) {
	ringQ := eval.params.RingQ()
	rtk, _ := eval.Rtks.GetRotationKey(galEl)
	level := utils.MinInt(ctIn.Level(), ctOut.Level())
	index := eval.PermuteNTTIndex[galEl]
	eval.SwitchKeysInPlace(level, ctIn.Value[1], rtk, eval.Pool[1].Q, eval.Pool[2].Q)
	ringQ.AddLvl(level, eval.Pool[1].Q, ctIn.Value[0], eval.Pool[1].Q)
	ringQ.PermuteNTTWithIndexLvl(level, eval.Pool[1].Q, index, ctOut.Value[0])
	ringQ.PermuteNTTWithIndexLvl(level, eval.Pool[2].Q, index, ctOut.Value[1])
}

// SwitchKeysInPlace applies the general key-switching procedure of the form [c0 + cx*evakey[0], c1 + cx*evakey[1]]
// Will return the result in the same NTT domain as the input cx.
func (eval *Evaluator) SwitchKeysInPlace(levelQ int, cx *ring.Poly, evakey *SwitchingKey, p0, p1 *ring.Poly) {

	levelP := evakey.LevelP()

	if levelP > 0 {
		eval.SwitchKeysInPlaceNoModDown(levelQ, cx, evakey, p0, eval.Pool[1].P, p1, eval.Pool[2].P)
	} else {
		eval.SwitchKeyInPlaceSinglePAndBitDecomp(levelQ, cx, evakey, p0, eval.Pool[1].P, p1, eval.Pool[2].P)
	}

	if cx.IsNTT && levelP != -1 {
		eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, p0, eval.Pool[1].P, p0)
		eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, p1, eval.Pool[2].P, p1)
	} else if !cx.IsNTT {
		eval.params.RingQ().InvNTTLazyLvl(levelQ, p0, p0)
		eval.params.RingQ().InvNTTLazyLvl(levelQ, p1, p1)

		if levelP != -1 {
			eval.params.RingP().InvNTTLazyLvl(levelP, eval.Pool[1].P, eval.Pool[1].P)
			eval.params.RingP().InvNTTLazyLvl(levelP, eval.Pool[2].P, eval.Pool[2].P)
			eval.BasisExtender.ModDownQPtoQ(levelQ, levelP, p0, eval.Pool[1].P, p0)
			eval.BasisExtender.ModDownQPtoQ(levelQ, levelP, p1, eval.Pool[2].P, p1)
		}
	}
}

// DecomposeNTT applies the full RNS basis decomposition for all q_alpha_i on c2.
// Expects the IsNTT flag of c2 to correctly reflect the domain of c2.
// PoolDecompQ and PoolDecompQ are vectors of polynomials (mod Q and mod P) that store the
// special RNS decomposition of c2 (in the NTT domain)
func (eval *Evaluator) DecomposeNTT(levelQ, levelP, alpha int, c2 *ring.Poly, PoolDecomp []ringqp.Poly) {

	ringQ := eval.params.RingQ()

	var polyNTT, polyInvNTT *ring.Poly

	if c2.IsNTT {
		polyNTT = c2
		polyInvNTT = eval.PoolInvNTT
		ringQ.InvNTTLvl(levelQ, polyNTT, polyInvNTT)
	} else {
		polyNTT = eval.PoolInvNTT
		polyInvNTT = c2
		ringQ.NTTLvl(levelQ, polyInvNTT, polyNTT)
	}

	beta := (levelQ + 1 + levelP) / (levelP + 1)

	for i := 0; i < beta; i++ {
		eval.DecomposeSingleNTT(levelQ, levelP, alpha, i, polyNTT, polyInvNTT, PoolDecomp[i].Q, PoolDecomp[i].P)
	}
}

// DecomposeSingleNTT takes the input polynomial c2 (c2NTT and c2InvNTT, respectively in the NTT and out of the NTT domain)
// modulo q_alpha_beta, and returns the result on c2QiQ are c2QiP the receiver polynomials
// respectively mod Q and mod P (in the NTT domain)
func (eval *Evaluator) DecomposeSingleNTT(levelQ, levelP, alpha, beta int, c2NTT, c2InvNTT, c2QiQ, c2QiP *ring.Poly) {

	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()

	eval.Decomposer.DecomposeAndSplit(levelQ, levelP, alpha, beta, c2InvNTT, c2QiQ, c2QiP)

	p0idxst := beta * alpha
	p0idxed := p0idxst + 1

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

// SwitchKeysInPlaceNoModDown applies the key-switch to the polynomial cx :
//
// pool2 = dot(decomp(cx) * evakey[0]) mod QP (encrypted input is multiplied by P factor)
// pool3 = dot(decomp(cx) * evakey[1]) mod QP (encrypted input is multiplied by P factor)
//
// Expects the flag IsNTT of cx to correctly reflect the domain of cx.
func (eval *Evaluator) SwitchKeysInPlaceNoModDown(levelQ int, cx *ring.Poly, evakey *SwitchingKey, c0Q, c0P, c1Q, c1P *ring.Poly) {

	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()
	ringQP := eval.params.RingQP()

	c2QP := eval.Pool[0]

	var cxNTT, cxInvNTT *ring.Poly
	if cx.IsNTT {
		cxNTT = cx
		cxInvNTT = eval.PoolInvNTT
		ringQ.InvNTTLvl(levelQ, cxNTT, cxInvNTT)
	} else {
		cxNTT = eval.PoolInvNTT
		cxInvNTT = cx
		ringQ.NTTLvl(levelQ, cxInvNTT, cxNTT)
	}

	c0QP := ringqp.Poly{Q: c0Q, P: c0P}
	c1QP := ringqp.Poly{Q: c1Q, P: c1P}

	levelP := evakey.Value[0][0][0].P.Level()
	decompRNS := (levelQ + 1 + levelP) / (levelP + 1)

	QiOverF := eval.params.QiOverflowMargin(levelQ) >> 1
	PiOverF := eval.params.PiOverflowMargin(levelP) >> 1

	// Key switching with CRT decomposition for the Qi
	var reduce int
	for i := 0; i < decompRNS; i++ {

		eval.DecomposeSingleNTT(levelQ, levelP, levelP+1, i, cxNTT, cxInvNTT, c2QP.Q, c2QP.P)

		if i == 0 {
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, evakey.Value[i][0][0], c2QP, c0QP)
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, evakey.Value[i][0][1], c2QP, c1QP)
		} else {
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, evakey.Value[i][0][0], c2QP, c0QP)
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, evakey.Value[i][0][1], c2QP, c1QP)
		}

		if reduce%QiOverF == QiOverF-1 {
			ringQ.ReduceLvl(levelQ, c0QP.Q, c0QP.Q)
			ringQ.ReduceLvl(levelQ, c1QP.Q, c1QP.Q)
		}

		if reduce%PiOverF == PiOverF-1 {
			ringP.ReduceLvl(levelP, c0QP.P, c0QP.P)
			ringP.ReduceLvl(levelP, c1QP.P, c1QP.P)
		}

		reduce++
	}

	if reduce%QiOverF != 0 {
		ringQ.ReduceLvl(levelQ, c0QP.Q, c0QP.Q)
		ringQ.ReduceLvl(levelQ, c1QP.Q, c1QP.Q)
	}

	if reduce%PiOverF != 0 {
		ringP.ReduceLvl(levelP, c0QP.P, c0QP.P)
		ringP.ReduceLvl(levelP, c1QP.P, c1QP.P)
	}
}

// SwitchKeyInPlaceSinglePAndBitDecomp applies the key-switch to the polynomial cx :
//
// pool2 = dot(decomp(cx) * evakey[0]) mod QP (encrypted input is multiplied by P factor)
// pool3 = dot(decomp(cx) * evakey[1]) mod QP (encrypted input is multiplied by P factor)
//
// Expects the flag IsNTT of cx to correctly reflect the domain of cx.
func (eval *Evaluator) SwitchKeyInPlaceSinglePAndBitDecomp(levelQ int, cx *ring.Poly, evakey *SwitchingKey, c0Q, c0P, c1Q, c1P *ring.Poly) {

	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()

	var cxInvNTT *ring.Poly
	if cx.IsNTT {
		cxInvNTT = eval.PoolInvNTT
		ringQ.InvNTTLvl(levelQ, cx, cxInvNTT)
	} else {
		cxInvNTT = cx
	}

	c0QP := ringqp.Poly{Q: c0Q, P: c0P}
	c1QP := ringqp.Poly{Q: c1Q, P: c1P}

	var levelP int
	if evakey.Value[0][0][0].P != nil {
		levelP = evakey.Value[0][0][0].P.Level()
	} else {
		levelP = -1
	}

	decompRNS := eval.params.DecompRNS(levelQ, levelP)
	decompBIT := eval.params.DecompBIT(levelQ, levelP)

	lb2 := eval.params.logbase2

	mask := uint64(((1 << lb2) - 1))

	if mask == 0 {
		mask = 0xFFFFFFFFFFFFFFFF
	}

	cw := eval.Pool[0].Q.Coeffs[0]
	cwNTT := eval.PoolBitDecomp

	QiOverF := eval.params.QiOverflowMargin(levelQ) >> 1
	PiOverF := eval.params.PiOverflowMargin(levelP) >> 1

	// Key switching with CRT decomposition for the Qi
	var reduce int
	for i := 0; i < decompRNS; i++ {
		for j := 0; j < decompBIT; j++ {

			ring.MaskVec(cxInvNTT.Coeffs[i], cw, j*lb2, mask)

			if i == 0 && j == 0 {
				for u := 0; u < levelQ+1; u++ {
					ringQ.NTTSingleLazy(u, cw, cwNTT)
					ring.MulCoeffsMontgomeryConstantVec(evakey.Value[i][j][0].Q.Coeffs[u], cwNTT, c0QP.Q.Coeffs[u], ringQ.Modulus[u], ringQ.MredParams[u])
					ring.MulCoeffsMontgomeryConstantVec(evakey.Value[i][j][1].Q.Coeffs[u], cwNTT, c1QP.Q.Coeffs[u], ringQ.Modulus[u], ringQ.MredParams[u])
				}

				for u := 0; u < levelP+1; u++ {
					ringP.NTTSingleLazy(u, cw, cwNTT)
					ring.MulCoeffsMontgomeryConstantVec(evakey.Value[i][j][0].P.Coeffs[u], cwNTT, c0QP.P.Coeffs[u], ringP.Modulus[u], ringP.MredParams[u])
					ring.MulCoeffsMontgomeryConstantVec(evakey.Value[i][j][1].P.Coeffs[u], cwNTT, c1QP.P.Coeffs[u], ringP.Modulus[u], ringP.MredParams[u])
				}
			} else {
				for u := 0; u < levelQ+1; u++ {
					ringQ.NTTSingleLazy(u, cw, cwNTT)
					ring.MulCoeffsMontgomeryConstantAndAddNoModVec(evakey.Value[i][j][0].Q.Coeffs[u], cwNTT, c0QP.Q.Coeffs[u], ringQ.Modulus[u], ringQ.MredParams[u])
					ring.MulCoeffsMontgomeryConstantAndAddNoModVec(evakey.Value[i][j][1].Q.Coeffs[u], cwNTT, c1QP.Q.Coeffs[u], ringQ.Modulus[u], ringQ.MredParams[u])
				}

				for u := 0; u < levelP+1; u++ {
					ringP.NTTSingleLazy(u, cw, cwNTT)
					ring.MulCoeffsMontgomeryConstantAndAddNoModVec(evakey.Value[i][j][0].P.Coeffs[u], cwNTT, c0QP.P.Coeffs[u], ringP.Modulus[u], ringP.MredParams[u])
					ring.MulCoeffsMontgomeryConstantAndAddNoModVec(evakey.Value[i][j][1].P.Coeffs[u], cwNTT, c1QP.P.Coeffs[u], ringP.Modulus[u], ringP.MredParams[u])
				}
			}

			if reduce%QiOverF == QiOverF-1 {
				ringQ.ReduceLvl(levelQ, c0QP.Q, c0QP.Q)
				ringQ.ReduceLvl(levelQ, c1QP.Q, c1QP.Q)
			}

			if reduce%PiOverF == PiOverF-1 {
				ringP.ReduceLvl(levelP, c0QP.P, c0QP.P)
				ringP.ReduceLvl(levelP, c1QP.P, c1QP.P)
			}

			reduce++
		}
	}

	if reduce%QiOverF != 0 {
		ringQ.ReduceLvl(levelQ, c0QP.Q, c0QP.Q)
		ringQ.ReduceLvl(levelQ, c1QP.Q, c1QP.Q)
	}

	if reduce%PiOverF != 0 {
		ringP.ReduceLvl(levelP, c0QP.P, c0QP.P)
		ringP.ReduceLvl(levelP, c1QP.P, c1QP.P)
	}
}

// KeyswitchHoisted applies the key-switch to the decomposed polynomial c2 mod QP (PoolDecompQ and PoolDecompP)
// and divides the result by P, reducing the basis from QP to Q.
//
// pool2 = dot(PoolDecompQ||PoolDecompP * evakey[0]) mod Q
// pool3 = dot(PoolDecompQ||PoolDecompP * evakey[1]) mod Q
func (eval *Evaluator) KeyswitchHoisted(levelQ int, PoolDecompQP []ringqp.Poly, evakey *SwitchingKey, c0Q, c1Q, c0P, c1P *ring.Poly) {

	eval.KeyswitchHoistedNoModDown(levelQ, PoolDecompQP, evakey, c0Q, c1Q, c0P, c1P)

	levelP := evakey.Value[0][0][0].P.Level()

	// Computes c0Q = c0Q/c0P and c1Q = c1Q/c1P
	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c0Q, c0P, c0Q)
	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c1Q, c1P, c1Q)
}

// KeyswitchHoistedNoModDown applies the key-switch to the decomposed polynomial c2 mod QP (PoolDecompQ and PoolDecompP)
//
// pool2 = dot(PoolDecompQ||PoolDecompP * evakey[0]) mod QP
// pool3 = dot(PoolDecompQ||PoolDecompP * evakey[1]) mod QP
func (eval *Evaluator) KeyswitchHoistedNoModDown(levelQ int, PoolDecompQP []ringqp.Poly, evakey *SwitchingKey, c0Q, c1Q, c0P, c1P *ring.Poly) {

	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()
	ringQP := eval.params.RingQP()

	c0QP := ringqp.Poly{Q: c0Q, P: c0P}
	c1QP := ringqp.Poly{Q: c1Q, P: c1P}

	levelP := evakey.Value[0][0][0].P.Level()
	decompRNS := (levelQ + 1 + levelP) / (levelP + 1)

	QiOverF := eval.params.QiOverflowMargin(levelQ) >> 1
	PiOverF := eval.params.PiOverflowMargin(levelP) >> 1

	// Key switching with CRT decomposition for the Qi
	var reduce int
	for i := 0; i < decompRNS; i++ {

		if i == 0 {
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, evakey.Value[i][0][0], PoolDecompQP[i], c0QP)
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, evakey.Value[i][0][1], PoolDecompQP[i], c1QP)
		} else {
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, evakey.Value[i][0][0], PoolDecompQP[i], c0QP)
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, evakey.Value[i][0][1], PoolDecompQP[i], c1QP)
		}

		if reduce%QiOverF == QiOverF-1 {
			ringQ.ReduceLvl(levelQ, c0QP.Q, c0QP.Q)
			ringQ.ReduceLvl(levelQ, c1QP.Q, c1QP.Q)
		}

		if reduce%PiOverF == PiOverF-1 {
			ringP.ReduceLvl(levelP, c0QP.P, c0QP.P)
			ringP.ReduceLvl(levelP, c1QP.P, c1QP.P)
		}

		reduce++
	}

	if reduce%QiOverF != 0 {
		ringQ.ReduceLvl(levelQ, c0QP.Q, c0QP.Q)
		ringQ.ReduceLvl(levelQ, c1QP.Q, c1QP.Q)
	}

	if reduce%PiOverF != 0 {
		ringP.ReduceLvl(levelP, c0QP.P, c0QP.P)
		ringP.ReduceLvl(levelP, c1QP.P, c1QP.P)
	}
}
