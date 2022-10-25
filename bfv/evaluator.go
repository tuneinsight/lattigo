package bfv

import (
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// Evaluator is an interface implementing the public methodes of the eval.
type Evaluator interface {
	Add(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)
	AddNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext)
	AddNoMod(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)
	AddNoModNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext)
	Sub(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)
	SubNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext)
	SubNoMod(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)
	SubNoModNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext)
	Neg(ctIn *rlwe.Ciphertext, ctOut *rlwe.Ciphertext)
	NegNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext)
	Reduce(ctIn *rlwe.Ciphertext, ctOut *rlwe.Ciphertext)
	ReduceNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext)
	AddScalar(ctIn *rlwe.Ciphertext, scalar uint64, ctOut *rlwe.Ciphertext)
	MulScalar(ctIn *rlwe.Ciphertext, scalar uint64, ctOut *rlwe.Ciphertext)
	MulScalarAndAdd(ctIn *rlwe.Ciphertext, scalar uint64, ctOut *rlwe.Ciphertext)
	MulScalarNew(ctIn *rlwe.Ciphertext, scalar uint64) (ctOut *rlwe.Ciphertext)
	Rescale(ctIn, ctOut *rlwe.Ciphertext)
	RescaleTo(level int, ctIn, ctOut *rlwe.Ciphertext)
	Mul(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)
	MulNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext)
	MulAndAdd(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)
	Relinearize(ctIn *rlwe.Ciphertext, ctOut *rlwe.Ciphertext)
	RelinearizeNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext)
	SwitchKeys(ctIn *rlwe.Ciphertext, switchKey *rlwe.SwitchingKey, ctOut *rlwe.Ciphertext)
	EvaluatePoly(input interface{}, pol *Polynomial) (opOut *rlwe.Ciphertext, err error)
	EvaluatePolyVector(input interface{}, pols []*Polynomial, encoder Encoder, slotsIndex map[int][]int) (opOut *rlwe.Ciphertext, err error)
	SwitchKeysNew(ctIn *rlwe.Ciphertext, switchkey *rlwe.SwitchingKey) (ctOut *rlwe.Ciphertext)
	RotateColumnsNew(ctIn *rlwe.Ciphertext, k int) (ctOut *rlwe.Ciphertext)
	RotateColumns(ctIn *rlwe.Ciphertext, k int, ctOut *rlwe.Ciphertext)
	RotateRows(ctIn *rlwe.Ciphertext, ctOut *rlwe.Ciphertext)
	RotateRowsNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext)
	InnerSum(ctIn *rlwe.Ciphertext, ctOut *rlwe.Ciphertext)
	ShallowCopy() Evaluator
	WithKey(rlwe.EvaluationKey) Evaluator

	BuffQ() [][]*ring.Poly
	BuffQMul() [][]*ring.Poly
	BuffPt() *rlwe.Plaintext
}

// evaluator is a struct that holds the necessary elements to perform the homomorphic operations between ciphertexts and/or plaintexts.
// It also holds a memory buffer used to store intermediate computations.
type evaluator struct {
	*evaluatorBase
	*evaluatorBuffers
	*rlwe.Evaluator

	basisExtenderQ1toQ2 *ring.BasisExtender
}

type evaluatorBase struct {
	params   Parameters
	ringQ    *ring.Ring
	ringP    *ring.Ring
	ringQMul *ring.Ring

	t         uint64
	tInvModQi []uint64
	levelQMul []int      // optimal #QiMul depending on #Qi (variable level)
	pHalf     []*big.Int // all prod(QiMul) / 2 depending on #Qi

	tDividesQ bool
}

func newEvaluatorPrecomp(params Parameters) *evaluatorBase {
	ev := new(evaluatorBase)

	ev.params = params

	ev.t = params.T()

	ev.ringQ = params.RingQ()
	ev.ringP = params.RingP()
	ev.ringQMul = params.RingQMul()

	ev.levelQMul = make([]int, len(ev.ringQ.Modulus))
	Q := new(big.Int).SetUint64(1)
	for i := range ev.levelQMul {
		Q.Mul(Q, new(big.Int).SetUint64(ev.ringQ.Modulus[i]))
		ev.levelQMul[i] = int(math.Ceil(float64(Q.BitLen()+params.LogN())/61.0)) - 1
	}

	ev.pHalf = make([]*big.Int, len(ev.ringQMul.Modulus))

	QMul := new(big.Int).SetUint64(1)
	for i := range ev.pHalf {
		QMul.Mul(QMul, new(big.Int).SetUint64(ev.ringQMul.Modulus[i]))
		ev.pHalf[i] = new(big.Int).Rsh(QMul, 1)
	}

	return ev
}

type evaluatorBuffers struct {
	buffQ    [][]*ring.Poly
	buffQMul [][]*ring.Poly
	buffPt   *rlwe.Plaintext
}

func newEvaluatorBuffer(eval *evaluatorBase) *evaluatorBuffers {
	evb := new(evaluatorBuffers)
	evb.buffQ = make([][]*ring.Poly, 4)
	evb.buffQMul = make([][]*ring.Poly, 4)
	for i := 0; i < 4; i++ {
		evb.buffQ[i] = make([]*ring.Poly, 6)
		evb.buffQMul[i] = make([]*ring.Poly, 6)
		for j := 0; j < 6; j++ {
			evb.buffQ[i][j] = eval.ringQ.NewPoly()
			evb.buffQMul[i][j] = eval.ringQMul.NewPoly()
		}
	}

	evb.buffPt = NewPlaintext(eval.params, eval.params.MaxLevel())

	return evb
}

// NewEvaluator creates a new Evaluator, that can be used to do homomorphic
// operations on ciphertexts and/or plaintexts. It stores a memory buffer
// and ciphertexts that will be used for intermediate values.
func NewEvaluator(params Parameters, evaluationKey rlwe.EvaluationKey) Evaluator {
	ev := new(evaluator)
	ev.evaluatorBase = newEvaluatorPrecomp(params)
	ev.evaluatorBuffers = newEvaluatorBuffer(ev.evaluatorBase)

	if params.T() != params.Q()[0] {
		ev.tInvModQi = make([]uint64, len(params.RingQ().Modulus))
		for i, qi := range params.RingQ().Modulus {
			ev.tInvModQi[i] = ring.MForm(ring.ModExp(params.T(), qi-2, qi), qi, params.RingQ().BredParams[i])
		}
	} else {
		ev.tDividesQ = true
	}

	ev.basisExtenderQ1toQ2 = ring.NewBasisExtender(ev.ringQ, ev.ringQMul)
	ev.Evaluator = rlwe.NewEvaluator(params.Parameters, &evaluationKey)

	return ev
}

// NewEvaluators creates n evaluators sharing the same read-only data-structures.
func NewEvaluators(params Parameters, evaluationKey rlwe.EvaluationKey, n int) []Evaluator {
	if n <= 0 {
		return []Evaluator{}
	}
	evas := make([]Evaluator, n)
	for i := range evas {
		if i == 0 {
			evas[0] = NewEvaluator(params, evaluationKey)
		} else {
			evas[i] = evas[i-1].ShallowCopy()
		}
	}
	return evas
}

// Add adds ctIn to op1 and returns the result in ctOut.
func (eval *evaluator) Add(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext) {
	el0, el1, elOut := eval.getElemAndCheckBinary(ctIn, op1, ctOut, utils.MaxInt(ctIn.Degree(), op1.Degree()), true)
	eval.evaluateInPlaceBinary(el0, el1, elOut, eval.ringQ.AddLvl)
}

// AddNew adds ctIn to op1 and creates a new element ctOut to store the result.
func (eval *evaluator) AddNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext) {
	ctOut = rlwe.NewCiphertext(eval.params.Parameters, utils.MaxInt(ctIn.Degree(), op1.Degree()), ctIn.Level())
	eval.Add(ctIn, op1, ctOut)
	return
}

// AddNoMod adds ctIn to op1 without modular reduction, and returns the result in cOut.
func (eval *evaluator) AddNoMod(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext) {
	el0, el1, elOut := eval.getElemAndCheckBinary(ctIn, op1, ctOut, utils.MaxInt(ctIn.Degree(), op1.Degree()), true)
	eval.evaluateInPlaceBinary(el0, el1, elOut, eval.ringQ.AddNoModLvl)
}

// AddNoModNew adds ctIn to op1 without modular reduction and creates a new element ctOut to store the result.
func (eval *evaluator) AddNoModNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext) {
	ctOut = rlwe.NewCiphertext(eval.params.Parameters, utils.MaxInt(ctIn.Degree(), op1.Degree()), ctIn.Level())
	eval.AddNoMod(ctIn, op1, ctOut)
	return
}

// Sub subtracts op1 from ctIn and returns the result in cOut.
func (eval *evaluator) Sub(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext) {
	el0, el1, elOut := eval.getElemAndCheckBinary(ctIn, op1, ctOut, utils.MaxInt(ctIn.Degree(), op1.Degree()), true)
	eval.evaluateInPlaceBinary(el0, el1, elOut, eval.ringQ.SubLvl)

	if el0.Degree() < el1.Degree() {
		for i := el0.Degree() + 1; i < el1.Degree()+1; i++ {
			eval.ringQ.Neg(ctOut.Value[i], ctOut.Value[i])
		}
	}
}

// SubNew subtracts op1 from ctIn and creates a new element ctOut to store the result.
func (eval *evaluator) SubNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext) {
	ctOut = rlwe.NewCiphertext(eval.params.Parameters, utils.MaxInt(ctIn.Degree(), op1.Degree()), ctIn.Level())
	eval.Sub(ctIn, op1, ctOut)
	return
}

// SubNoMod subtracts op1 from ctIn without modular reduction and returns the result on ctOut.
func (eval *evaluator) SubNoMod(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext) {
	el0, el1, elOut := eval.getElemAndCheckBinary(ctIn, op1, ctOut, utils.MaxInt(ctIn.Degree(), op1.Degree()), true)

	eval.evaluateInPlaceBinary(el0, el1, elOut, eval.ringQ.SubNoModLvl)

	if el0.Degree() < el1.Degree() {
		for i := el0.Degree() + 1; i < el1.Degree()+1; i++ {
			eval.ringQ.Neg(ctOut.Value[i], ctOut.Value[i])
		}
	}
}

// SubNoModNew subtracts op1 from ctIn without modular reduction and creates a new element ctOut to store the result.
func (eval *evaluator) SubNoModNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext) {
	ctOut = rlwe.NewCiphertext(eval.params.Parameters, utils.MaxInt(ctIn.Degree(), op1.Degree()), ctIn.Level())
	eval.SubNoMod(ctIn, op1, ctOut)
	return
}

// Neg negates ctIn and returns the result in ctOut.
func (eval *evaluator) Neg(ctIn, ctOut *rlwe.Ciphertext) {
	el0, elOut := eval.getElemAndCheckUnary(ctIn, ctOut, ctIn.Degree())
	evaluateInPlaceUnary(el0, elOut, eval.ringQ.NegLvl)
}

// NegNew negates ctIn and creates a new element to store the result.
func (eval *evaluator) NegNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {
	ctOut = rlwe.NewCiphertext(eval.params.Parameters, ctIn.Degree(), ctIn.Level())
	eval.Neg(ctIn, ctOut)
	return ctOut
}

// Reduce applies a modular reduction to op and returns the result in ctOut.
func (eval *evaluator) Reduce(ctIn *rlwe.Ciphertext, ctOut *rlwe.Ciphertext) {
	el0, elOut := eval.getElemAndCheckUnary(ctIn, ctOut, ctIn.Degree())
	evaluateInPlaceUnary(el0, elOut, eval.ringQ.ReduceLvl)
}

// ReduceNew applies a modular reduction to ctIn and creates a new element ctOut to store the result.
func (eval *evaluator) ReduceNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {
	ctOut = rlwe.NewCiphertext(eval.params.Parameters, ctIn.Degree(), ctIn.Level())
	eval.Reduce(ctIn, ctOut)
	return ctOut
}

// MulScalar multiplies ctIn by a uint64 scalar and returns the result in ctOut.
func (eval *evaluator) MulScalar(ctIn *rlwe.Ciphertext, scalar uint64, ctOut *rlwe.Ciphertext) {
	el0, elOut := eval.getElemAndCheckUnary(ctIn, ctOut, ctIn.Degree())
	fun := func(level int, el, elOut *ring.Poly) { eval.ringQ.MulScalarLvl(level, el, scalar, elOut) }
	evaluateInPlaceUnary(el0, elOut, fun)
}

// MulScalarAndAdd multiplies ctIn by a uint64 scalar and adds the result on ctOut.
func (eval *evaluator) MulScalarAndAdd(ctIn *rlwe.Ciphertext, scalar uint64, ctOut *rlwe.Ciphertext) {
	el0, elOut := eval.getElemAndCheckUnary(ctIn, ctOut, ctIn.Degree())
	fun := func(level int, el, elOut *ring.Poly) { eval.ringQ.MulScalarAndAddLvl(level, el, scalar, elOut) }
	evaluateInPlaceUnary(el0, elOut, fun)
}

// AddScalar adds the scalar on ctIn and returns the result on ctOut.
func (eval *evaluator) AddScalar(ctIn *rlwe.Ciphertext, scalar uint64, ctOut *rlwe.Ciphertext) {
	scalarBigint := new(big.Int).SetUint64(scalar)
	scalarBigint.Mul(scalarBigint, eval.ringQ.ModulusAtLevel[ctIn.Level()])
	ring.DivRound(scalarBigint, eval.params.RingT().ModulusAtLevel[0], scalarBigint)
	tmp := new(big.Int)
	for i, qi := range eval.ringQ.Modulus[:ctIn.Level()+1] {
		ctOut.Value[0].Coeffs[i][0] = ring.CRed(ctIn.Value[0].Coeffs[i][0]+tmp.Mod(scalarBigint, new(big.Int).SetUint64(qi)).Uint64(), qi)
	}
}

// MulScalarNew multiplies ctIn by a uint64 scalar and creates a new element ctOut to store the result.
func (eval *evaluator) MulScalarNew(ctIn *rlwe.Ciphertext, scalar uint64) (ctOut *rlwe.Ciphertext) {
	ctOut = rlwe.NewCiphertext(eval.params.Parameters, ctIn.Degree(), ctIn.Level())
	eval.MulScalar(ctIn, scalar, ctOut)
	return
}

// Rescale divides the ciphertext by the last modulus.
func (eval *evaluator) Rescale(ctIn, ctOut *rlwe.Ciphertext) {
	eval.RescaleTo(ctIn.Level()-1, ctIn, ctOut)
}

// RescaleTo divides the ciphertext by the last moduli until it has `level+1` moduli left.
func (eval *evaluator) RescaleTo(level int, ctIn, ctOut *rlwe.Ciphertext) {

	if ctIn.Level() < level || ctOut.Level() < ctIn.Level()-level {
		panic("cannot RescaleTo: (ctIn.Level() || ctOut.Level()) < level")
	}

	eval.ringQ.DivRoundByLastModulusManyLvl(ctIn.Level(), ctIn.Level()-level, ctIn.Value[0], eval.buffQ[0][0], ctOut.Value[0])
	eval.ringQ.DivRoundByLastModulusManyLvl(ctIn.Level(), ctIn.Level()-level, ctIn.Value[1], eval.buffQ[0][0], ctOut.Value[1])

	ctOut.Resize(ctOut.Degree(), level)
}

// tensorAndRescale computes (ct0 x ct1) * (t/Q) and stores the result in ctOut.
func (eval *evaluator) tensorAndRescale(ct0, ct1, ctOut *rlwe.Ciphertext) {

	level := utils.MinInt(utils.MinInt(ct0.Level(), ct1.Level()), ctOut.Level())

	levelQMul := eval.levelQMul[level]

	ctOut.Resize(ctOut.Degree(), level)

	c0Q1 := eval.buffQ[0]
	c0Q2 := eval.buffQMul[0]

	c1Q1 := eval.buffQ[1]
	c1Q2 := eval.buffQMul[1]

	// Prepares the ciphertexts for the Tensoring by extending their
	// basis from Q to QP and transforming them to NTT form
	eval.modUpAndNTTLvl(level, levelQMul, ct0, c0Q1, c0Q2)

	if ct0 != ct1 {
		eval.modUpAndNTTLvl(level, levelQMul, ct1, c1Q1, c1Q2)
	}

	// Tensoring: multiplies each elements of the ciphertexts together
	// and adds them to their corresponding position in the new ciphertext
	// based on their respective degree

	// Case where both Elements are of degree 1
	if ct0.Degree() == 1 && ct1.Degree() == 1 {
		eval.tensoreLowDegLvl(level, levelQMul, ct0, ct1)
		// Case where at least one element is not of degree 1
	} else {
		eval.tensortLargeDegLvl(level, levelQMul, ct0, ct1)
	}

	eval.quantizeLvl(level, levelQMul, ctOut)
}

func (eval *evaluator) modUpAndNTTLvl(level, levelQMul int, ct *rlwe.Ciphertext, cQ, cQMul []*ring.Poly) {
	for i := range ct.Value {
		eval.basisExtenderQ1toQ2.ModUpQtoP(level, levelQMul, ct.Value[i], cQMul[i])
		eval.ringQ.NTTLazyLvl(level, ct.Value[i], cQ[i])
		eval.ringQMul.NTTLazyLvl(levelQMul, cQMul[i], cQMul[i])
	}
}

func (eval *evaluator) tensoreLowDegLvl(level, levelQMul int, ct0, ct1 *rlwe.Ciphertext) {

	c0Q1 := eval.buffQ[0]
	c0Q2 := eval.buffQMul[0]

	c1Q1 := eval.buffQ[1]
	c1Q2 := eval.buffQMul[1]

	c2Q1 := eval.buffQ[2]
	c2Q2 := eval.buffQMul[2]

	c00Q := eval.buffQ[3][0]
	c00Q2 := eval.buffQMul[3][0]
	c01Q := eval.buffQ[3][1]
	c01P := eval.buffQMul[3][1]

	eval.ringQ.MFormLvl(level, c0Q1[0], c00Q)
	eval.ringQMul.MForm(c0Q2[0], c00Q2)

	eval.ringQ.MFormLvl(level, c0Q1[1], c01Q)
	eval.ringQMul.MForm(c0Q2[1], c01P)

	// Squaring case
	if ct0 == ct1 {

		// c0 = c0[0]*c0[0]
		eval.ringQ.MulCoeffsMontgomeryLvl(level, c00Q, c0Q1[0], c2Q1[0])
		eval.ringQMul.MulCoeffsMontgomeryLvl(levelQMul, c00Q2, c0Q2[0], c2Q2[0])

		// c1 = 2*c0[0]*c0[1]
		eval.ringQ.MulCoeffsMontgomeryLvl(level, c00Q, c0Q1[1], c2Q1[1])
		eval.ringQMul.MulCoeffsMontgomeryLvl(levelQMul, c00Q2, c0Q2[1], c2Q2[1])

		eval.ringQ.AddNoModLvl(level, c2Q1[1], c2Q1[1], c2Q1[1])
		eval.ringQMul.AddNoModLvl(levelQMul, c2Q2[1], c2Q2[1], c2Q2[1])

		// c2 = c0[1]*c0[1]
		eval.ringQ.MulCoeffsMontgomeryLvl(level, c01Q, c0Q1[1], c2Q1[2])
		eval.ringQMul.MulCoeffsMontgomeryLvl(levelQMul, c01P, c0Q2[1], c2Q2[2])

		// Normal case
	} else {

		// c0 = c0[0]*c1[0]
		eval.ringQ.MulCoeffsMontgomeryLvl(level, c00Q, c1Q1[0], c2Q1[0])
		eval.ringQMul.MulCoeffsMontgomeryLvl(levelQMul, c00Q2, c1Q2[0], c2Q2[0])

		// c1 = c0[0]*c1[1] + c0[1]*c1[0]
		eval.ringQ.MulCoeffsMontgomeryLvl(level, c00Q, c1Q1[1], c2Q1[1])
		eval.ringQMul.MulCoeffsMontgomeryLvl(levelQMul, c00Q2, c1Q2[1], c2Q2[1])

		eval.ringQ.MulCoeffsMontgomeryAndAddNoModLvl(level, c01Q, c1Q1[0], c2Q1[1])
		eval.ringQMul.MulCoeffsMontgomeryAndAddNoModLvl(levelQMul, c01P, c1Q2[0], c2Q2[1])

		// c2 = c0[1]*c1[1]
		eval.ringQ.MulCoeffsMontgomeryLvl(level, c01Q, c1Q1[1], c2Q1[2])
		eval.ringQMul.MulCoeffsMontgomeryLvl(levelQMul, c01P, c1Q2[1], c2Q2[2])
	}
}

func (eval *evaluator) tensortLargeDegLvl(level, levelQMul int, ct0, ct1 *rlwe.Ciphertext) {

	c0Q1 := eval.buffQ[0]
	c0Q2 := eval.buffQMul[0]

	c1Q1 := eval.buffQ[1]
	c1Q2 := eval.buffQMul[1]

	c2Q1 := eval.buffQ[2]
	c2Q2 := eval.buffQMul[2]

	for i := 0; i < ct0.Degree()+ct1.Degree()+1; i++ {
		c2Q1[i].Zero()
		c2Q2[i].Zero()
	}

	// Squaring case
	if ct0 == ct1 {

		c00Q1 := eval.buffQ[3]
		c00Q2 := eval.buffQMul[3]

		for i := range ct0.Value {
			eval.ringQ.MFormLvl(level, c0Q1[i], c00Q1[i])
			eval.ringQMul.MFormLvl(levelQMul, c0Q2[i], c00Q2[i])
		}

		for i := 0; i < ct0.Degree()+1; i++ {
			for j := i + 1; j < ct0.Degree()+1; j++ {
				eval.ringQ.MulCoeffsMontgomeryLvl(level, c00Q1[i], c0Q1[j], c2Q1[i+j])
				eval.ringQMul.MulCoeffsMontgomeryLvl(levelQMul, c00Q2[i], c0Q2[j], c2Q2[i+j])

				eval.ringQ.AddLvl(level, c2Q1[i+j], c2Q1[i+j], c2Q1[i+j])
				eval.ringQMul.AddLvl(levelQMul, c2Q2[i+j], c2Q2[i+j], c2Q2[i+j])
			}
		}

		for i := 0; i < ct0.Degree()+1; i++ {
			eval.ringQ.MulCoeffsMontgomeryAndAddLvl(level, c00Q1[i], c0Q1[i], c2Q1[i<<1])
			eval.ringQMul.MulCoeffsMontgomeryAndAddLvl(levelQMul, c00Q2[i], c0Q2[i], c2Q2[i<<1])
		}

		// Normal case
	} else {
		for i := range ct0.Value {
			eval.ringQ.MFormLvl(level, c0Q1[i], c0Q1[i])
			eval.ringQMul.MFormLvl(levelQMul, c0Q2[i], c0Q2[i])
			for j := range ct1.Value {
				eval.ringQ.MulCoeffsMontgomeryAndAddLvl(level, c0Q1[i], c1Q1[j], c2Q1[i+j])
				eval.ringQMul.MulCoeffsMontgomeryAndAddLvl(levelQMul, c0Q2[i], c1Q2[j], c2Q2[i+j])
			}
		}
	}
}

func (eval *evaluator) quantizeLvl(level, levelQMul int, ctOut *rlwe.Ciphertext) {

	c2Q1 := eval.buffQ[2]
	c2Q2 := eval.buffQMul[2]

	// Applies the inverse NTT to the ciphertext, scales down the ciphertext
	// by t/q and reduces its basis from QP to Q
	for i := range ctOut.Value {
		eval.ringQ.InvNTTLazyLvl(level, c2Q1[i], c2Q1[i])
		eval.ringQMul.InvNTTLazyLvl(levelQMul, c2Q2[i], c2Q2[i])

		// Extends the basis Q of ct(x) to the basis P and Divides (ct(x)Q -> P) by Q
		eval.basisExtenderQ1toQ2.ModDownQPtoP(level, levelQMul, c2Q1[i], c2Q2[i], c2Q2[i]) // QP / Q -> P

		// Centers ct(x)P by (P-1)/2 and extends ct(x)P to the basis Q
		eval.ringQMul.AddScalarBigintLvl(levelQMul, c2Q2[i], eval.pHalf[levelQMul], c2Q2[i])
		eval.basisExtenderQ1toQ2.ModUpPtoQ(levelQMul, level, c2Q2[i], ctOut.Value[i])
		eval.ringQ.SubScalarBigintLvl(level, ctOut.Value[i], eval.pHalf[levelQMul], ctOut.Value[i])

		// (ct(x)/Q)*T, doing so only requires that Q*P > Q*Q, faster but adds error ~|T|
		eval.ringQ.MulScalarLvl(level, ctOut.Value[i], eval.t, ctOut.Value[i])
	}
}

// Mul multiplies ctIn by op1 and returns the result in ctOut.
func (eval *evaluator) Mul(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext) {
	el0, el1, elOut := eval.getElemAndCheckBinary(ctIn, op1, ctOut, ctIn.Degree()+op1.Degree(), false)
	switch op1 := op1.(type) {
	case *PlaintextMul:
		eval.mulPlaintextMul(ctIn, op1, ctOut)
	case *PlaintextRingT:
		eval.mulPlaintextRingT(ctIn, op1, ctOut)
	case *rlwe.Plaintext, *rlwe.Ciphertext:
		eval.tensorAndRescale(el0, el1, elOut)
	default:
		panic(fmt.Errorf("invalid rlwe.Operand type for Mul: %T", op1))
	}
}

// MulAndAdd multiplies ctIn with op1 and adds the result on ctOut.
func (eval *evaluator) MulAndAdd(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext) {

	level := utils.MinInt(ctIn.Level(), ctOut.Level())

	ct2 := &rlwe.Ciphertext{Value: make([]*ring.Poly, ctIn.Degree()+op1.Degree()+1)}
	for i := range ct2.Value {
		ct2.Value[i] = new(ring.Poly)
		ct2.Value[i].Coeffs = eval.buffQ[2][i].Coeffs[:level+1]
	}

	eval.Mul(ctIn, op1, ct2)
	eval.Add(ctOut, ct2, ctOut)
}

func (eval *evaluator) mulPlaintextMul(ctIn *rlwe.Ciphertext, ptRt *PlaintextMul, ctOut *rlwe.Ciphertext) {

	level := utils.MinInt(ctIn.Level(), ctOut.Level())

	for i := range ctIn.Value {
		eval.ringQ.NTTLazyLvl(level, ctIn.Value[i], ctOut.Value[i])
		eval.ringQ.MulCoeffsMontgomeryConstantLvl(level, ctOut.Value[i], ptRt.Value, ctOut.Value[i])
		eval.ringQ.InvNTTLvl(level, ctOut.Value[i], ctOut.Value[i])
	}
}

func (eval *evaluator) mulPlaintextRingT(ctIn *rlwe.Ciphertext, ptRt *PlaintextRingT, ctOut *rlwe.Ciphertext) {

	level := utils.MinInt(ctIn.Level(), ctOut.Level())

	ctOut.Resize(ctOut.Degree(), level)

	ringQ := eval.ringQ

	coeffs := ptRt.Value.Coeffs[0]
	coeffsNTT := eval.buffQ[0][0].Coeffs[0]

	for i := range ctIn.Value {

		// Copies the inputCT on the outputCT and switches to the NTT domain
		eval.ringQ.NTTLazyLvl(level, ctIn.Value[i], ctOut.Value[i])

		// Switches the outputCT in the Montgomery domain
		eval.ringQ.MFormLvl(level, ctOut.Value[i], ctOut.Value[i])

		// For each qi in Q
		for j := range ringQ.Modulus[:level+1] {

			tmp := ctOut.Value[i].Coeffs[j]
			qi := ringQ.Modulus[j]
			nttPsi := ringQ.NttPsi[j]
			bredParams := ringQ.BredParams[j]
			mredParams := ringQ.MredParams[j]

			// Transforms the plaintext in the NTT domain of that qi
			ring.NTTLazy(coeffs, coeffsNTT, ringQ.N, nttPsi, qi, mredParams, bredParams)

			// Multiplies NTT_qi(pt) * NTT_qi(ct)
			ring.MulCoeffsMontgomeryVec(tmp, coeffsNTT, tmp, qi, mredParams)

		}

		// Switches the ciphertext out of the NTT domain
		eval.ringQ.InvNTTLvl(level, ctOut.Value[i], ctOut.Value[i])
	}
}

// MulNew multiplies ctIn by op1 and creates a new element ctOut to store the result.
func (eval *evaluator) MulNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext) {
	ctOut = rlwe.NewCiphertext(eval.params.Parameters, ctIn.Degree()+op1.Degree(), ctIn.Level())
	eval.Mul(ctIn, op1, ctOut)
	return
}

// RelinearizeNew relinearizes the ciphertext ctIn of degree > 1 until it is of degree 1, and creates a new ciphertext to store the result.
func (eval *evaluator) RelinearizeNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {
	ctOut = rlwe.NewCiphertext(eval.params.Parameters, 1, ctIn.Level())
	eval.Relinearize(ctIn, ctOut)
	return
}

// SwitchKeysNew applies the key-switching procedure to the ciphertext ct0 and creates a new ciphertext to store the result. It requires as an additional input a valid switching-key:
// it must encrypt the target key under the public key under which ct0 is currently encrypted.
func (eval *evaluator) SwitchKeysNew(ctIn *rlwe.Ciphertext, switchkey *rlwe.SwitchingKey) (ctOut *rlwe.Ciphertext) {
	ctOut = rlwe.NewCiphertext(eval.params.Parameters, 1, ctIn.Level())
	eval.SwitchKeys(ctIn, switchkey, ctOut)
	return
}

// RotateColumns rotates the columns of ct0 by k positions to the left and returns the result in ctOut. As an additional input it requires a RotationKeys struct.
func (eval *evaluator) RotateColumns(ct0 *rlwe.Ciphertext, k int, ctOut *rlwe.Ciphertext) {
	eval.Automorphism(ct0, eval.params.GaloisElementForColumnRotationBy(k), ctOut)
}

// RotateColumnsNew applies RotateColumns and returns the result in a new Ciphertext.
func (eval *evaluator) RotateColumnsNew(ctIn *rlwe.Ciphertext, k int) (ctOut *rlwe.Ciphertext) {
	ctOut = rlwe.NewCiphertext(eval.params.Parameters, 1, ctIn.Level())
	eval.RotateColumns(ctIn, k, ctOut)
	return
}

// RotateRows rotates the rows of ct0 and returns the result in ctOut.
func (eval *evaluator) RotateRows(ct0 *rlwe.Ciphertext, ctOut *rlwe.Ciphertext) {
	eval.Automorphism(ct0, eval.params.GaloisElementForRowRotation(), ctOut)
}

// RotateRowsNew rotates the rows of ctIn and returns the result a new Ciphertext.
func (eval *evaluator) RotateRowsNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {
	ctOut = rlwe.NewCiphertext(eval.params.Parameters, 1, ctIn.Level())
	eval.RotateRows(ctIn, ctOut)
	return
}

// InnerSum computes the inner sum of ctIn and returns the result in ctOut. It requires a rotation key that stores all the left powers of two rotations.
// The resulting vector will be of the form [sum, sum, .., sum, sum].
func (eval *evaluator) InnerSum(ctIn *rlwe.Ciphertext, ctOut *rlwe.Ciphertext) {
	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		panic("cannot InnerSum: input and output must be of degree 1")
	}
	cTmp := rlwe.NewCiphertext(eval.params.Parameters, 1, ctIn.Level())
	ctOut.Copy(ctIn.El())

	for i := 1; i < int(eval.ringQ.N>>1); i <<= 1 {
		eval.RotateColumns(ctOut, i, cTmp)
		eval.Add(cTmp, ctOut, ctOut)
	}

	eval.RotateRows(ctOut, cTmp)
	eval.Add(ctOut, cTmp, ctOut)
}

// ShallowCopy creates a shallow copy of this evaluator in which the read-only data-structures are
// shared with the receiver.
func (eval *evaluator) ShallowCopy() Evaluator {
	return &evaluator{
		evaluatorBase:       eval.evaluatorBase,
		Evaluator:           eval.Evaluator.ShallowCopy(),
		evaluatorBuffers:    newEvaluatorBuffer(eval.evaluatorBase),
		basisExtenderQ1toQ2: eval.basisExtenderQ1toQ2.ShallowCopy(),
	}
}

// WithKey creates a shallow copy of this evaluator in which the read-only data-structures are
// shared with the receiver but the EvaluationKey is evaluationKey.
func (eval *evaluator) WithKey(evaluationKey rlwe.EvaluationKey) Evaluator {
	return &evaluator{
		evaluatorBase:       eval.evaluatorBase,
		Evaluator:           eval.Evaluator.WithKey(&evaluationKey),
		evaluatorBuffers:    eval.evaluatorBuffers,
		basisExtenderQ1toQ2: eval.basisExtenderQ1toQ2,
	}
}

// BuffQ returns the internal evaluator buffQ buffer.
func (eval *evaluator) BuffQ() [][]*ring.Poly {
	return eval.buffQ
}

// BuffQMul returns the internal evaluator buffQMul buffer.
func (eval *evaluator) BuffQMul() [][]*ring.Poly {
	return eval.buffQMul
}

// BuffPt returns the internal evaluator plaintext buffer.
func (eval *evaluator) BuffPt() *rlwe.Plaintext {
	return eval.buffPt
}

func (eval *evaluator) getRingQElem(level int, op rlwe.Operand) *rlwe.Ciphertext {
	switch o := op.(type) {
	case *rlwe.Ciphertext, *rlwe.Plaintext:
		return o.El()
	case *PlaintextRingT:
		if eval.tDividesQ {
			ScaleUpTIsQ0VecLvl(level, eval.params.RingQ(), o.Value, eval.buffPt.Value)
		} else {
			ScaleUpTCoprimeWithQVecLvl(level, eval.params.RingQ(), eval.params.RingT(), eval.tInvModQi, eval.BuffQP[0].Q.Coeffs[0], o.Value, eval.buffPt.Value)
		}
		return eval.buffPt.El()
	default:
		panic(fmt.Errorf("cannot getRingQElem: invalid rlwe.Operand type for operation: %T", o))
	}
}

// getElemAndCheckBinary unwraps the elements from the rlwe.Operands and checks that the receiver has sufficiently large degree.
func (eval *evaluator) getElemAndCheckBinary(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext, opOutMinDegree int, ensureRingQ bool) (el0, el1, elOut *rlwe.Ciphertext) {
	if ctIn == nil || op1 == nil || ctOut == nil {
		panic("cannot getElemAndCheckBinary: ctIn, op1 or ctOut cannot be nil")
	}

	if ctOut.Degree() < opOutMinDegree {
		panic("cannot getElemAndCheckBinary: ctOut.Degree() degree is too small")
	}

	if ctIn.Level() != ctOut.Level() {
		panic("cannot getElemAndCheckBinary: ctIn and ctOut must be at the same level")
	}

	if _, isPtRingT := op1.(*PlaintextRingT); !isPtRingT && ctIn.Level() != op1.Level() {
		panic("cannot getElemAndCheckBinary: ctIn & op1 of type *bfv.Plaintext or *bfv.PlaintextMul must be at the same level")
	}

	level := ctIn.Level()

	if ensureRingQ {
		return ctIn.El(), eval.getRingQElem(level, op1), ctOut.El()
	}

	return ctIn.El(), op1.El(), ctOut.El()
}

func (eval *evaluator) getElemAndCheckUnary(ctIn, ctOut *rlwe.Ciphertext, opOutMinDegree int) (el0, elOut *rlwe.Ciphertext) {

	if ctIn == nil || ctOut == nil {
		panic("cannot getElemAndCheckUnary: ctIn or ctOut cannot be nil")
	}

	if ctIn.Level() != ctOut.Level() {
		panic("cannot getElemAndCheckUnary: ctIn.Level() is not equal to ctOut.Level()")
	}

	if ctOut.Degree() < opOutMinDegree {
		panic("cannot getElemAndCheckUnary: ctOut.Degree() is too small")
	}

	return ctIn.El(), ctOut.El()
}

// evaluateInPlaceBinary applies the provided function in place on el0 and el1 and returns the result in elOut.
func (eval *evaluator) evaluateInPlaceBinary(el0, el1, elOut *rlwe.Ciphertext, evaluate func(int, *ring.Poly, *ring.Poly, *ring.Poly)) {

	smallest, largest, _ := rlwe.GetSmallestLargest(el0, el1)

	level := utils.MinInt(utils.MinInt(el0.Level(), el1.Level()), elOut.Level())

	elOut.Resize(elOut.Degree(), level)

	for i := 0; i < smallest.Degree()+1; i++ {
		evaluate(level, el0.Value[i], el1.Value[i], elOut.Value[i])
	}

	// If the inputs degrees differ, it copies the remaining degree on the receiver.
	if largest != nil && largest != elOut { // checks to avoid unnecessary work.
		for i := smallest.Degree() + 1; i < largest.Degree()+1; i++ {
			elOut.Value[i].Copy(largest.Value[i])
		}
	}
}

// evaluateInPlaceUnary applies the provided function in place on el0 and returns the result in elOut.
func evaluateInPlaceUnary(el0, elOut *rlwe.Ciphertext, evaluate func(int, *ring.Poly, *ring.Poly)) {

	level := utils.MinInt(el0.Level(), elOut.Level())

	elOut.Resize(elOut.Degree(), level)

	for i := range el0.Value {
		evaluate(level, el0.Value[i], elOut.Value[i])
	}
}
