package bfv

import (
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// Operand is a common interface for Ciphertext and Plaintext.
type Operand interface {
	El() *rlwe.Ciphertext
	Level() int
	Degree() int
}

// Evaluator is an interface implementing the public methodes of the eval.
type Evaluator interface {
	Add(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext)
	AddNew(ctIn *Ciphertext, op1 Operand) (ctOut *Ciphertext)
	AddNoMod(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext)
	AddNoModNew(ctIn *Ciphertext, op1 Operand) (ctOut *Ciphertext)
	Sub(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext)
	SubNew(ctIn *Ciphertext, op1 Operand) (ctOut *Ciphertext)
	SubNoMod(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext)
	SubNoModNew(ctIn *Ciphertext, op1 Operand) (ctOut *Ciphertext)
	Neg(ctIn *Ciphertext, ctOut *Ciphertext)
	NegNew(ctIn *Ciphertext) (ctOut *Ciphertext)
	Reduce(ctIn *Ciphertext, ctOut *Ciphertext)
	ReduceNew(ctIn *Ciphertext) (ctOut *Ciphertext)
	AddScalar(ctIn *Ciphertext, scalar uint64, ctOut *Ciphertext)
	MulScalar(ctIn *Ciphertext, scalar uint64, ctOut *Ciphertext)
	MulScalarAndAdd(ctIn *Ciphertext, scalar uint64, ctOut *Ciphertext)
	MulScalarNew(ctIn *Ciphertext, scalar uint64) (ctOut *Ciphertext)
	Rescale(ctIn, ctOut *Ciphertext)
	RescaleTo(level int, ctIn, ctOut *Ciphertext)
	Mul(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext)
	MulNew(ctIn *Ciphertext, op1 Operand) (ctOut *Ciphertext)
	MulAndAdd(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext)
	Relinearize(ctIn *Ciphertext, ctOut *Ciphertext)
	RelinearizeNew(ctIn *Ciphertext) (ctOut *Ciphertext)
	SwitchKeys(ctIn *Ciphertext, switchKey *rlwe.SwitchingKey, ctOut *Ciphertext)
	EvaluatePoly(input interface{}, pol *Polynomial) (opOut *Ciphertext, err error)
	EvaluatePolyVector(input interface{}, pols []*Polynomial, encoder Encoder, slotsIndex map[int][]int) (opOut *Ciphertext, err error)
	SwitchKeysNew(ctIn *Ciphertext, switchkey *rlwe.SwitchingKey) (ctOut *Ciphertext)
	RotateColumnsNew(ctIn *Ciphertext, k int) (ctOut *Ciphertext)
	RotateColumns(ctIn *Ciphertext, k int, ctOut *Ciphertext)
	RotateRows(ctIn *Ciphertext, ctOut *Ciphertext)
	RotateRowsNew(ctIn *Ciphertext) (ctOut *Ciphertext)
	InnerSum(ctIn *Ciphertext, ctOut *Ciphertext)
	ShallowCopy() Evaluator
	WithKey(rlwe.EvaluationKey) Evaluator

	BuffQ() [][]*ring.Poly
	BuffQMul() [][]*ring.Poly
	BuffPt() *Plaintext
}

// evaluator is a struct that holds the necessary elements to perform the homomorphic operations between ciphertexts and/or plaintexts.
// It also holds a memory buffer used to store intermediate computations.
type evaluator struct {
	*evaluatorBase
	*evaluatorBuffers
	*rlwe.KeySwitcher

	rlk  *rlwe.RelinearizationKey
	rtks *rlwe.RotationKeySet

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
	buffPt   *Plaintext
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

	evb.buffPt = NewPlaintext(eval.params)

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
	if params.PCount() != 0 {
		ev.KeySwitcher = rlwe.NewKeySwitcher(params.Parameters)
	}
	ev.rlk = evaluationKey.Rlk
	ev.rtks = evaluationKey.Rtks
	return ev
}

// NewEvaluators creates n evaluators sharing the same read-only data-structures.
func NewEvaluators(params Parameters, evaluationKey rlwe.EvaluationKey, n int) []Evaluator {
	if n <= 0 {
		return []Evaluator{}
	}
	evas := make([]Evaluator, n, n)
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
func (eval *evaluator) Add(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext) {
	el0, el1, elOut := eval.getElemAndCheckBinary(ctIn, op1, ctOut, utils.MaxInt(ctIn.Degree(), op1.Degree()), true)
	eval.evaluateInPlaceBinary(el0, el1, elOut, eval.ringQ.AddLvl)
}

// AddNew adds ctIn to op1 and creates a new element ctOut to store the result.
func (eval *evaluator) AddNew(ctIn *Ciphertext, op1 Operand) (ctOut *Ciphertext) {
	ctOut = NewCiphertextLvl(eval.params, utils.MaxInt(ctIn.Degree(), op1.Degree()), ctIn.Level())
	eval.Add(ctIn, op1, ctOut)
	return
}

// AddNoMod adds ctIn to op1 without modular reduction, and returns the result in cOut.
func (eval *evaluator) AddNoMod(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext) {
	el0, el1, elOut := eval.getElemAndCheckBinary(ctIn, op1, ctOut, utils.MaxInt(ctIn.Degree(), op1.Degree()), true)
	eval.evaluateInPlaceBinary(el0, el1, elOut, eval.ringQ.AddNoModLvl)
}

// AddNoModNew adds ctIn to op1 without modular reduction and creates a new element ctOut to store the result.
func (eval *evaluator) AddNoModNew(ctIn *Ciphertext, op1 Operand) (ctOut *Ciphertext) {
	ctOut = NewCiphertextLvl(eval.params, utils.MaxInt(ctIn.Degree(), op1.Degree()), ctIn.Level())
	eval.AddNoMod(ctIn, op1, ctOut)
	return
}

// Sub subtracts op1 from ctIn and returns the result in cOut.
func (eval *evaluator) Sub(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext) {
	el0, el1, elOut := eval.getElemAndCheckBinary(ctIn, op1, ctOut, utils.MaxInt(ctIn.Degree(), op1.Degree()), true)
	eval.evaluateInPlaceBinary(el0, el1, elOut, eval.ringQ.SubLvl)

	if el0.Degree() < el1.Degree() {
		for i := el0.Degree() + 1; i < el1.Degree()+1; i++ {
			eval.ringQ.Neg(ctOut.Value[i], ctOut.Value[i])
		}
	}
}

// SubNew subtracts op1 from ctIn and creates a new element ctOut to store the result.
func (eval *evaluator) SubNew(ctIn *Ciphertext, op1 Operand) (ctOut *Ciphertext) {
	ctOut = NewCiphertextLvl(eval.params, utils.MaxInt(ctIn.Degree(), op1.Degree()), ctIn.Level())
	eval.Sub(ctIn, op1, ctOut)
	return
}

// SubNoMod subtracts op1 from ctIn without modular reduction and returns the result on ctOut.
func (eval *evaluator) SubNoMod(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext) {
	el0, el1, elOut := eval.getElemAndCheckBinary(ctIn, op1, ctOut, utils.MaxInt(ctIn.Degree(), op1.Degree()), true)

	eval.evaluateInPlaceBinary(el0, el1, elOut, eval.ringQ.SubNoModLvl)

	if el0.Degree() < el1.Degree() {
		for i := el0.Degree() + 1; i < el1.Degree()+1; i++ {
			eval.ringQ.Neg(ctOut.Value[i], ctOut.Value[i])
		}
	}
}

// SubNoModNew subtracts op1 from ctIn without modular reduction and creates a new element ctOut to store the result.
func (eval *evaluator) SubNoModNew(ctIn *Ciphertext, op1 Operand) (ctOut *Ciphertext) {
	ctOut = NewCiphertextLvl(eval.params, utils.MaxInt(ctIn.Degree(), op1.Degree()), ctIn.Level())
	eval.SubNoMod(ctIn, op1, ctOut)
	return
}

// Neg negates ctIn and returns the result in ctOut.
func (eval *evaluator) Neg(ctIn, ctOut *Ciphertext) {
	el0, elOut := eval.getElemAndCheckUnary(ctIn, ctOut, ctIn.Degree())
	evaluateInPlaceUnary(el0, elOut, eval.ringQ.NegLvl)
}

// NegNew negates ctIn and creates a new element to store the result.
func (eval *evaluator) NegNew(ctIn *Ciphertext) (ctOut *Ciphertext) {
	ctOut = NewCiphertextLvl(eval.params, ctIn.Degree(), ctIn.Level())
	eval.Neg(ctIn, ctOut)
	return ctOut
}

// Reduce applies a modular reduction to op and returns the result in ctOut.
func (eval *evaluator) Reduce(ctIn *Ciphertext, ctOut *Ciphertext) {
	el0, elOut := eval.getElemAndCheckUnary(ctIn, ctOut, ctIn.Degree())
	evaluateInPlaceUnary(el0, elOut, eval.ringQ.ReduceLvl)
}

// ReduceNew applies a modular reduction to ctIn and creates a new element ctOut to store the result.
func (eval *evaluator) ReduceNew(ctIn *Ciphertext) (ctOut *Ciphertext) {
	ctOut = NewCiphertextLvl(eval.params, ctIn.Degree(), ctIn.Level())
	eval.Reduce(ctIn, ctOut)
	return ctOut
}

// MulScalar multiplies ctIn by a uint64 scalar and returns the result in ctOut.
func (eval *evaluator) MulScalar(ctIn *Ciphertext, scalar uint64, ctOut *Ciphertext) {
	el0, elOut := eval.getElemAndCheckUnary(ctIn, ctOut, ctIn.Degree())
	fun := func(level int, el, elOut *ring.Poly) { eval.ringQ.MulScalarLvl(level, el, scalar, elOut) }
	evaluateInPlaceUnary(el0, elOut, fun)
}

// MulScalarAndAdd multiplies ctIn by a uint64 scalar and adds the result on ctOut.
func (eval *evaluator) MulScalarAndAdd(ctIn *Ciphertext, scalar uint64, ctOut *Ciphertext) {
	el0, elOut := eval.getElemAndCheckUnary(ctIn, ctOut, ctIn.Degree())
	fun := func(level int, el, elOut *ring.Poly) { eval.ringQ.MulScalarAndAddLvl(level, el, scalar, elOut) }
	evaluateInPlaceUnary(el0, elOut, fun)
}

// AddScalar adds the scalar on ctIn and returns the result on ctOut.
func (eval *evaluator) AddScalar(ctIn *Ciphertext, scalar uint64, ctOut *Ciphertext) {
	scalarBigint := new(big.Int).SetUint64(scalar)
	scalarBigint.Mul(scalarBigint, eval.ringQ.ModulusAtLevel[ctIn.Level()])
	ring.DivRound(scalarBigint, eval.params.RingT().ModulusAtLevel[0], scalarBigint)
	tmp := new(big.Int)
	for i, qi := range eval.ringQ.Modulus[:ctIn.Level()+1] {
		ctOut.Value[0].Coeffs[i][0] = ring.CRed(ctIn.Value[0].Coeffs[i][0]+tmp.Mod(scalarBigint, new(big.Int).SetUint64(qi)).Uint64(), qi)
	}
}

// MulScalarNew multiplies ctIn by a uint64 scalar and creates a new element ctOut to store the result.
func (eval *evaluator) MulScalarNew(ctIn *Ciphertext, scalar uint64) (ctOut *Ciphertext) {
	ctOut = NewCiphertextLvl(eval.params, ctIn.Degree(), ctIn.Level())
	eval.MulScalar(ctIn, scalar, ctOut)
	return
}

// Rescale divides the ciphertext by the last modulus.
func (eval *evaluator) Rescale(ctIn, ctOut *Ciphertext) {
	eval.RescaleTo(ctIn.Level()-1, ctIn, ctOut)
}

// RescaleTo divides the ciphertext by the last moduli until it has `level+1` moduli left.
func (eval *evaluator) RescaleTo(level int, ctIn, ctOut *Ciphertext) {

	if ctIn.Level() < level || ctOut.Level() < ctIn.Level()-level {
		panic("cannot RescaleTo: (ctIn.Level() || ctOut.Level()) < level")
	}

	eval.ringQ.DivRoundByLastModulusManyLvl(ctIn.Level(), ctIn.Level()-level, ctIn.Value[0], eval.buffQ[0][0], ctOut.Value[0])
	eval.ringQ.DivRoundByLastModulusManyLvl(ctIn.Level(), ctIn.Level()-level, ctIn.Value[1], eval.buffQ[0][0], ctOut.Value[1])

	ctOut.Value[0].Coeffs = ctOut.Value[0].Coeffs[:level+1]
	ctOut.Value[1].Coeffs = ctOut.Value[1].Coeffs[:level+1]
}

// tensorAndRescale computes (ct0 x ct1) * (t/Q) and stores the result in ctOut.
func (eval *evaluator) tensorAndRescale(ct0, ct1, ctOut *rlwe.Ciphertext) {

	level := utils.MinInt(utils.MinInt(ct0.Level(), ct1.Level()), ctOut.Level())

	levelQMul := eval.levelQMul[level]

	ctOut.Value[0].Coeffs = ctOut.Value[0].Coeffs[:level+1]
	ctOut.Value[1].Coeffs = ctOut.Value[1].Coeffs[:level+1]

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
func (eval *evaluator) Mul(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext) {
	el0, el1, elOut := eval.getElemAndCheckBinary(ctIn, op1, ctOut, ctIn.Degree()+op1.Degree(), false)
	switch op1 := op1.(type) {
	case *PlaintextMul:
		eval.mulPlaintextMul(ctIn, op1, ctOut)
	case *PlaintextRingT:
		eval.mulPlaintextRingT(ctIn, op1, ctOut)
	case *Plaintext, *Ciphertext:
		eval.tensorAndRescale(el0, el1, elOut)
	default:
		panic(fmt.Errorf("invalid operand type for Mul: %T", op1))
	}
}

// MulAndAdd multiplies ctIn with op1 and adds the result on ctOut.
func (eval *evaluator) MulAndAdd(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext) {

	level := utils.MinInt(ctIn.Level(), ctOut.Level())

	ct2 := &Ciphertext{&rlwe.Ciphertext{Value: make([]*ring.Poly, ctIn.Degree()+op1.Degree()+1)}}
	for i := range ct2.Value {
		ct2.Value[i] = new(ring.Poly)
		ct2.Value[i].Coeffs = eval.buffQ[2][i].Coeffs[:level+1]
	}

	eval.Mul(ctIn, op1, ct2)
	eval.Add(ctOut, ct2, ctOut)
}

func (eval *evaluator) mulPlaintextMul(ctIn *Ciphertext, ptRt *PlaintextMul, ctOut *Ciphertext) {

	level := utils.MinInt(ctIn.Level(), ctOut.Level())

	for i := range ctIn.Value {
		eval.ringQ.NTTLazyLvl(level, ctIn.Value[i], ctOut.Value[i])
		eval.ringQ.MulCoeffsMontgomeryConstantLvl(level, ctOut.Value[i], ptRt.Value, ctOut.Value[i])
		eval.ringQ.InvNTTLvl(level, ctOut.Value[i], ctOut.Value[i])
	}
}

func (eval *evaluator) mulPlaintextRingT(ctIn *Ciphertext, ptRt *PlaintextRingT, ctOut *Ciphertext) {

	level := utils.MinInt(ctIn.Level(), ctOut.Level())

	ctOut.Value[0].Coeffs = ctOut.Value[0].Coeffs[:level+1]
	ctOut.Value[1].Coeffs = ctOut.Value[1].Coeffs[:level+1]

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
func (eval *evaluator) MulNew(ctIn *Ciphertext, op1 Operand) (ctOut *Ciphertext) {
	ctOut = NewCiphertextLvl(eval.params, ctIn.Degree()+op1.Degree(), ctIn.Level())
	eval.Mul(ctIn, op1, ctOut)
	return
}

// relinearize is a method common to Relinearize and RelinearizeNew. It switches ctIn to the NTT domain, applies the keyswitch, and returns the result out of the NTT domain.
func (eval *evaluator) relinearize(ctIn *Ciphertext, ctOut *Ciphertext) {

	level := utils.MinInt(ctIn.Level(), ctOut.Level())

	ctOut.Value[0].Coeffs = ctOut.Value[0].Coeffs[:level+1]
	ctOut.Value[1].Coeffs = ctOut.Value[1].Coeffs[:level+1]

	if ctOut != ctIn {
		ring.CopyValues(ctIn.Value[0], ctOut.Value[0])
		ring.CopyValues(ctIn.Value[1], ctOut.Value[1])
	}

	for deg := uint64(ctIn.Degree()); deg > 1; deg-- {
		eval.SwitchKeysInPlace(level, ctIn.Value[deg], eval.rlk.Keys[deg-2], eval.BuffQP[1].Q, eval.BuffQP[2].Q)
		eval.ringQ.AddLvl(level, ctOut.Value[0], eval.BuffQP[1].Q, ctOut.Value[0])
		eval.ringQ.AddLvl(level, ctOut.Value[1], eval.BuffQP[2].Q, ctOut.Value[1])
	}

	ctOut.SetValue(ctOut.Value[:2])
}

// Relinearize relinearizes the ciphertext ctIn of degree > 1 until it is of degree 1, and returns the result in cOut.
//
// It requires a correct evaluation key as additional input:
//
// - it must match the secret-key that was used to create the public key under which the current ct0 is encrypted.
//
// - it must be of degree high enough to relinearize the input ciphertext to degree 1 (e.g., a ciphertext
// of degree 3 will require that the evaluation key stores the keys for both degree 3 and degree 2 ciphertexts).
func (eval *evaluator) Relinearize(ctIn *Ciphertext, ctOut *Ciphertext) {

	if eval.rlk == nil {
		panic("evaluator has no relinearization key")
	}

	if ctIn.Degree()-1 > len(eval.rlk.Keys) {
		panic("input ciphertext degree is too large to allow relinearization with the evluator's relinearization key")
	}

	eval.getElemAndCheckUnary(ctIn, ctOut, 1)

	if ctIn.Degree() < 2 {
		if ctIn != ctOut {
			ctOut.Copy(ctIn.El())
		}
	} else {
		eval.relinearize(ctIn, ctOut)
	}
}

// RelinearizeNew relinearizes the ciphertext ctIn of degree > 1 until it is of degree 1, and creates a new ciphertext to store the result.
//
// Requires a correct evaluation key as additional input:
//
// - it must match the secret-key that was used to create the public key under which the current ct0 is encrypted
//
// - it must be of degree high enough to relinearize the input ciphertext to degree 1 (e.g., a ciphertext
// of degree 3 will require that the evaluation key stores the keys for both degree 3 and degree 2 ciphertexts).
func (eval *evaluator) RelinearizeNew(ctIn *Ciphertext) (ctOut *Ciphertext) {
	ctOut = NewCiphertextLvl(eval.params, 1, ctIn.Level())
	eval.Relinearize(ctIn, ctOut)
	return
}

// SwitchKeys applies the key-switching procedure to the ciphertext ct0 and returns the result in ctOut. It requires as an additional input a valid switching-key:
// it must encrypt the target key under the public key under which ct0 is currently encrypted.
func (eval *evaluator) SwitchKeys(ctIn *Ciphertext, switchKey *rlwe.SwitchingKey, ctOut *Ciphertext) {

	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		panic("cannot SwitchKeys: input and output must be of degree 1 to allow key switching")
	}

	eval.getElemAndCheckUnary(ctIn, ctOut, 1)

	level := utils.MinInt(ctIn.Level(), ctOut.Level())
	eval.SwitchKeysInPlace(level, ctIn.Value[1], switchKey, eval.BuffQP[1].Q, eval.BuffQP[2].Q)
	eval.ringQ.AddLvl(level, ctIn.Value[0], eval.BuffQP[1].Q, ctOut.Value[0])
	ring.CopyValues(eval.BuffQP[2].Q, ctOut.Value[1])
}

// SwitchKeysNew applies the key-switching procedure to the ciphertext ct0 and creates a new ciphertext to store the result. It requires as an additional input a valid switching-key:
// it must encrypt the target key under the public key under which ct0 is currently encrypted.
func (eval *evaluator) SwitchKeysNew(ctIn *Ciphertext, switchkey *rlwe.SwitchingKey) (ctOut *Ciphertext) {
	ctOut = NewCiphertextLvl(eval.params, 1, ctIn.Level())
	eval.SwitchKeys(ctIn, switchkey, ctOut)
	return
}

// RotateColumns rotates the columns of ct0 by k positions to the left and returns the result in ctOut. As an additional input it requires a RotationKeys struct:
//
// - it must either store all the left and right power-of-2 rotations or the specific rotation that is requested.
//
// If only the power-of-two rotations are stored, the numbers k and n/2-k will be decomposed in base-2 and the rotation with the lowest
// hamming weight will be chosen; then the specific rotation will be computed as a sum of powers of two rotations.
func (eval *evaluator) RotateColumns(ctIn *Ciphertext, k int, ctOut *Ciphertext) {

	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		panic("cannot RotateColumns: input and or output must be of degree 1")
	}

	eval.getElemAndCheckUnary(ctIn, ctOut, 1)

	if k == 0 {

		ctOut.Copy(ctIn.El())

	} else {

		galElL := eval.params.GaloisElementForColumnRotationBy(k)
		// Looks in the rotation key if the corresponding rotation has been generated or if the input is a plaintext
		if swk, inSet := eval.rtks.GetRotationKey(galElL); inSet {

			eval.permuteLvl(ctIn.Level(), ctIn, galElL, swk, ctOut)

		} else {
			panic(fmt.Errorf("evaluator has no rotation key for rotation by %d", k))
		}
	}
}

// permuteLvl performs a column rotation on ctIn and returns the result in ctOut
func (eval *evaluator) permuteLvl(level int, ctIn *Ciphertext, generator uint64, switchKey *rlwe.SwitchingKey, ctOut *Ciphertext) {
	eval.SwitchKeysInPlace(level, ctIn.Value[1], switchKey, eval.BuffQP[1].Q, eval.BuffQP[2].Q)
	eval.ringQ.AddLvl(level, eval.BuffQP[1].Q, ctIn.Value[0], eval.BuffQP[1].Q)
	eval.ringQ.PermuteLvl(level, eval.BuffQP[1].Q, generator, ctOut.Value[0])
	eval.ringQ.PermuteLvl(level, eval.BuffQP[2].Q, generator, ctOut.Value[1])
}

// RotateColumnsNew applies RotateColumns and returns the result in a new Ciphertext.
func (eval *evaluator) RotateColumnsNew(ctIn *Ciphertext, k int) (ctOut *Ciphertext) {
	ctOut = NewCiphertextLvl(eval.params, 1, ctIn.Level())
	eval.RotateColumns(ctIn, k, ctOut)
	return
}

// RotateRows rotates the rows of ctIn and returns the result in ctOut.
func (eval *evaluator) RotateRows(ctIn *Ciphertext, ctOut *Ciphertext) {

	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		panic("cannot RotateRows: input and/or output must be of degree 1")
	}

	eval.getElemAndCheckUnary(ctIn, ctOut, 1)

	galEl := eval.params.GaloisElementForRowRotation()

	if key, inSet := eval.rtks.GetRotationKey(galEl); inSet {
		eval.permuteLvl(ctIn.Level(), ctIn, galEl, key, ctOut)
	} else {
		panic("evaluator has no rotation key for row rotation")
	}
}

// RotateRowsNew rotates the rows of ctIn and returns the result a new Ciphertext.
func (eval *evaluator) RotateRowsNew(ctIn *Ciphertext) (ctOut *Ciphertext) {
	ctOut = NewCiphertextLvl(eval.params, 1, ctIn.Level())
	eval.RotateRows(ctIn, ctOut)
	return
}

// InnerSum computes the inner sum of ctIn and returns the result in ctOut. It requires a rotation key that stores all the left powers of two rotations.
// The resulting vector will be of the form [sum, sum, .., sum, sum].
func (eval *evaluator) InnerSum(ctIn *Ciphertext, ctOut *Ciphertext) {

	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		panic("cannot InnerSum: input and output must be of degree 1")
	}

	cTmp := NewCiphertextLvl(eval.params, 1, ctIn.Level())

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
		KeySwitcher:         eval.KeySwitcher.ShallowCopy(),
		evaluatorBuffers:    newEvaluatorBuffer(eval.evaluatorBase),
		basisExtenderQ1toQ2: eval.basisExtenderQ1toQ2.ShallowCopy(),
		rlk:                 eval.rlk,
		rtks:                eval.rtks,
	}
}

// WithKey creates a shallow copy of this evaluator in which the read-only data-structures are
// shared with the receiver but the EvaluationKey is evaluationKey.
func (eval *evaluator) WithKey(evaluationKey rlwe.EvaluationKey) Evaluator {
	return &evaluator{
		evaluatorBase:       eval.evaluatorBase,
		KeySwitcher:         eval.KeySwitcher,
		evaluatorBuffers:    eval.evaluatorBuffers,
		basisExtenderQ1toQ2: eval.basisExtenderQ1toQ2,
		rlk:                 evaluationKey.Rlk,
		rtks:                evaluationKey.Rtks,
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
func (eval *evaluator) BuffPt() *Plaintext {
	return eval.buffPt
}

func (eval *evaluator) getRingQElem(level int, op Operand) *rlwe.Ciphertext {
	switch o := op.(type) {
	case *Ciphertext, *Plaintext:
		return o.El()
	case *PlaintextRingT:
		if eval.tDividesQ {
			ScaleUpTIsQ0VecLvl(level, eval.params.RingQ(), o.Value, eval.buffPt.Value)
		} else {
			ScaleUpTCoprimeWithQVecLvl(level, eval.params.RingQ(), eval.params.RingT(), eval.tInvModQi, eval.BuffQP[0].Q.Coeffs[0], o.Value, eval.buffPt.Value)
		}
		return eval.buffPt.El()
	default:
		panic(fmt.Errorf("cannot getRingQElem: invalid operand type for operation: %T", o))
	}
}

// getElemAndCheckBinary unwraps the elements from the operands and checks that the receiver has sufficiently large degree.
func (eval *evaluator) getElemAndCheckBinary(ctIn *Ciphertext, op1 Operand, ctOut *Ciphertext, opOutMinDegree int, ensureRingQ bool) (el0, el1, elOut *rlwe.Ciphertext) {
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

func (eval *evaluator) getElemAndCheckUnary(ctIn, ctOut *Ciphertext, opOutMinDegree int) (el0, elOut *rlwe.Ciphertext) {

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

	if elOut.Level() > level {
		for i := range elOut.Value {
			elOut.Value[i].Coeffs = elOut.Value[i].Coeffs[:level+1]
		}
	}

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

	if elOut.Level() > level {
		for i := range elOut.Value {
			elOut.Value[i].Coeffs = elOut.Value[i].Coeffs[:level+1]
		}
	}

	for i := range el0.Value {
		evaluate(level, el0.Value[i], elOut.Value[i])
	}
}
