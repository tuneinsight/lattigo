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
	Sub(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)
	SubNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext)
	Neg(ctIn *rlwe.Ciphertext, ctOut *rlwe.Ciphertext)
	NegNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext)
	AddScalar(ctIn *rlwe.Ciphertext, scalar uint64, ctOut *rlwe.Ciphertext)
	MulScalar(ctIn *rlwe.Ciphertext, scalar uint64, ctOut *rlwe.Ciphertext)
	MulScalarThenAdd(ctIn *rlwe.Ciphertext, scalar uint64, ctOut *rlwe.Ciphertext)
	MulScalarNew(ctIn *rlwe.Ciphertext, scalar uint64) (ctOut *rlwe.Ciphertext)
	Rescale(ctIn, ctOut *rlwe.Ciphertext)
	RescaleTo(level int, ctIn, ctOut *rlwe.Ciphertext)
	Mul(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)
	MulNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext)
	MulThenAdd(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)
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
	InnerSum(op0 *rlwe.Ciphertext, batch, n int, ctOut *rlwe.Ciphertext)
	ShallowCopy() Evaluator
	WithKey(rlwe.EvaluationKey) Evaluator

	CheckBinary(op0, op1, opOut rlwe.Operand, opOutMinDegree int) (degree, level int)
	CheckUnary(op0, opOut rlwe.Operand) (degree, level int)
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
	params Parameters

	tInvModQi []uint64
	levelQMul []int      // optimal #QiMul depending on #Qi (variable level)
	qMulHalf  []*big.Int // all prod(QiMul) / 2 depending on #Qi

	tDividesQ bool
}

func newEvaluatorPrecomp(params Parameters) *evaluatorBase {
	ev := new(evaluatorBase)

	ev.params = params

	ringQ := params.RingQ()

	ev.levelQMul = make([]int, params.RingQ().ModuliChainLength())
	for i := range ev.levelQMul {
		ev.levelQMul[i] = int(math.Ceil(float64(ringQ.AtLevel(i).Modulus().BitLen()+params.LogN())/61.0)) - 1
	}

	ringQMul := params.RingQMul()

	ev.qMulHalf = make([]*big.Int, ringQMul.ModuliChainLength())
	for i := range ev.qMulHalf {
		ev.qMulHalf[i] = new(big.Int).Rsh(ringQMul.AtLevel(i).Modulus(), 1)
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
			evb.buffQ[i][j] = eval.params.RingQ().NewPoly()
			evb.buffQMul[i][j] = eval.params.RingQMul().NewPoly()
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

	ringQ := params.RingQ()

	if params.T() != params.Q()[0] {
		ev.tInvModQi = make([]uint64, ringQ.ModuliChainLength())
		for i, s := range ringQ.SubRings {
			ev.tInvModQi[i] = ring.MForm(ring.ModExp(params.T(), s.Modulus-2, s.Modulus), s.Modulus, s.BRedConstant)
		}
	} else {
		ev.tDividesQ = true
	}

	ev.basisExtenderQ1toQ2 = ring.NewBasisExtender(ev.params.RingQ(), ev.params.RingQMul())
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
	_, level := eval.CheckBinary(ctIn, op1, ctOut, utils.MaxInt(ctIn.Degree(), op1.Degree()))
	ctOut.Resize(ctOut.Degree(), level)
	eval.evaluateInPlaceBinary(ctIn, op1.El(), ctOut, eval.params.RingQ().AtLevel(level).Add)
}

// AddNew adds ctIn to op1 and creates a new element ctOut to store the result.
func (eval *evaluator) AddNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, utils.MaxInt(ctIn.Degree(), op1.Degree()), ctIn.Level())
	eval.Add(ctIn, op1, ctOut)
	return
}

// Sub subtracts op1 from ctIn and returns the result in cOut.
func (eval *evaluator) Sub(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext) {
	_, level := eval.CheckBinary(ctIn, op1, ctOut, utils.MaxInt(ctIn.Degree(), op1.Degree()))
	ctOut.Resize(ctOut.Degree(), level)
	eval.evaluateInPlaceBinary(ctIn, op1.El(), ctOut, eval.params.RingQ().AtLevel(level).Sub)

	if ctIn.Degree() < op1.Degree() {
		for i := ctIn.Degree() + 1; i < op1.Degree()+1; i++ {
			eval.params.RingQ().AtLevel(level).Neg(ctOut.Value[i], ctOut.Value[i])
		}
	}
}

// SubNew subtracts op1 from ctIn and creates a new element ctOut to store the result.
func (eval *evaluator) SubNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, utils.MaxInt(ctIn.Degree(), op1.Degree()), ctIn.Level())
	eval.Sub(ctIn, op1, ctOut)
	return
}

// Neg negates ctIn and returns the result in ctOut.
func (eval *evaluator) Neg(ctIn, ctOut *rlwe.Ciphertext) {
	_, level := eval.CheckUnary(ctIn, ctOut)
	ctOut.Resize(ctOut.Degree(), level)
	evaluateInPlaceUnary(ctIn, ctOut, eval.params.RingQ().AtLevel(level).Neg)
}

// NegNew negates ctIn and creates a new element to store the result.
func (eval *evaluator) NegNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, ctIn.Degree(), ctIn.Level())
	eval.Neg(ctIn, ctOut)
	return ctOut
}

// MulScalar multiplies ctIn by a uint64 scalar and returns the result in ctOut.
func (eval *evaluator) MulScalar(ctIn *rlwe.Ciphertext, scalar uint64, ctOut *rlwe.Ciphertext) {
	_, level := eval.CheckUnary(ctIn, ctOut)
	ctOut.Resize(ctOut.Degree(), level)
	evaluateInPlaceUnary(ctIn, ctOut, func(el, elOut *ring.Poly) { eval.params.RingQ().AtLevel(level).MulScalar(el, scalar, elOut) })
}

// MulScalarThenAdd multiplies ctIn by a uint64 scalar and adds the result on ctOut.
func (eval *evaluator) MulScalarThenAdd(ctIn *rlwe.Ciphertext, scalar uint64, ctOut *rlwe.Ciphertext) {
	_, level := eval.CheckUnary(ctIn, ctOut)
	ctOut.Resize(ctOut.Degree(), level)
	evaluateInPlaceUnary(ctIn, ctOut, func(el, elOut *ring.Poly) { eval.params.RingQ().AtLevel(level).MulScalarThenAdd(el, scalar, elOut) })
}

// AddScalar adds the scalar on ctIn and returns the result on ctOut.
func (eval *evaluator) AddScalar(ctIn *rlwe.Ciphertext, scalar uint64, ctOut *rlwe.Ciphertext) {
	_, level := eval.CheckUnary(ctIn, ctOut)
	ctOut.Resize(ctOut.Degree(), level)
	ringQ := eval.params.RingQ().AtLevel(level)
	scalarBigint := new(big.Int).SetUint64(scalar)
	scalarBigint.Mul(scalarBigint, ringQ.Modulus())
	ring.DivRound(scalarBigint, eval.params.RingT().Modulus(), scalarBigint)
	tmp := new(big.Int)

	for i := 0; i < level+1; i++ {
		qi := ringQ.SubRings[i].Modulus
		ctOut.Value[0].Coeffs[i][0] = ring.CRed(ctIn.Value[0].Coeffs[i][0]+tmp.Mod(scalarBigint, new(big.Int).SetUint64(qi)).Uint64(), qi)
	}
}

// MulScalarNew multiplies ctIn by a uint64 scalar and creates a new element ctOut to store the result.
func (eval *evaluator) MulScalarNew(ctIn *rlwe.Ciphertext, scalar uint64) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, ctIn.Degree(), ctIn.Level())
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

	ringQ := eval.params.RingQ().AtLevel(ctIn.Level())

	ringQ.DivRoundByLastModulusMany(ctIn.Level()-level, ctIn.Value[0], eval.buffQ[0][0], ctOut.Value[0])
	ringQ.DivRoundByLastModulusMany(ctIn.Level()-level, ctIn.Value[1], eval.buffQ[0][0], ctOut.Value[1])

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
	eval.tensoreLowDegLvl(level, levelQMul, ct0, ct1)

	eval.quantizeLvl(level, levelQMul, ctOut)
}

func (eval *evaluator) modUpAndNTTLvl(level, levelQMul int, ct *rlwe.Ciphertext, cQ, cQMul []*ring.Poly) {

	ringQ := eval.params.RingQ().AtLevel(level)
	ringQMul := eval.params.RingQMul().AtLevel(levelQMul)

	for i := range ct.Value {
		eval.basisExtenderQ1toQ2.ModUpQtoP(level, levelQMul, ct.Value[i], cQMul[i])
		ringQ.NTTLazy(ct.Value[i], cQ[i])
		ringQMul.NTTLazy(cQMul[i], cQMul[i])
	}
}

func (eval *evaluator) tensoreLowDegLvl(level, levelQMul int, ct0, ct1 *rlwe.Ciphertext) {

	c0Q1 := eval.buffQ[0]    // NTT(ct0) mod Q
	c0Q2 := eval.buffQMul[0] // NTT(ct0) mod QMul

	c1Q1 := eval.buffQ[1]    // NTT(ct1) mod Q
	c1Q2 := eval.buffQMul[1] // NTT(ct1) mod QMul

	c2Q1 := eval.buffQ[2]    //Receiver mod Q
	c2Q2 := eval.buffQMul[2] //Receiver mod QMul

	ringQ := eval.params.RingQ().AtLevel(level)
	ringQMul := eval.params.RingQMul().AtLevel(levelQMul)

	if ct0.Degree() == 1 && ct1.Degree() == 1 {

		c00Q := eval.buffQ[3][0]
		c00Q2 := eval.buffQMul[3][0]

		ringQ.MForm(c0Q1[0], c00Q)
		ringQMul.MForm(c0Q2[0], c00Q2)

		c01Q := eval.buffQ[3][1]
		c01P := eval.buffQMul[3][1]

		ringQ.MForm(c0Q1[1], c01Q)
		ringQMul.MForm(c0Q2[1], c01P)

		// Squaring case
		if ct0 == ct1 {

			// c0 = c0[0]*c0[0]
			ringQ.MulCoeffsMontgomery(c00Q, c0Q1[0], c2Q1[0])
			ringQMul.MulCoeffsMontgomery(c00Q2, c0Q2[0], c2Q2[0])

			// c1 = 2*c0[0]*c0[1]
			ringQ.MulCoeffsMontgomery(c00Q, c0Q1[1], c2Q1[1])
			ringQMul.MulCoeffsMontgomery(c00Q2, c0Q2[1], c2Q2[1])

			ringQ.AddLazy(c2Q1[1], c2Q1[1], c2Q1[1])
			ringQMul.AddLazy(c2Q2[1], c2Q2[1], c2Q2[1])

			// c2 = c0[1]*c0[1]
			ringQ.MulCoeffsMontgomery(c01Q, c0Q1[1], c2Q1[2])
			ringQMul.MulCoeffsMontgomery(c01P, c0Q2[1], c2Q2[2])

			// Normal case
		} else {

			// c0 = c0[0]*c1[0]
			ringQ.MulCoeffsMontgomery(c00Q, c1Q1[0], c2Q1[0])
			ringQMul.MulCoeffsMontgomery(c00Q2, c1Q2[0], c2Q2[0])

			// c1 = c0[0]*c1[1] + c0[1]*c1[0]
			ringQ.MulCoeffsMontgomery(c00Q, c1Q1[1], c2Q1[1])
			ringQMul.MulCoeffsMontgomery(c00Q2, c1Q2[1], c2Q2[1])

			ringQ.MulCoeffsMontgomeryThenAddLazy(c01Q, c1Q1[0], c2Q1[1])
			ringQMul.MulCoeffsMontgomeryThenAddLazy(c01P, c1Q2[0], c2Q2[1])

			// c2 = c0[1]*c1[1]
			ringQ.MulCoeffsMontgomery(c01Q, c1Q1[1], c2Q1[2])
			ringQMul.MulCoeffsMontgomery(c01P, c1Q2[1], c2Q2[2])
		}
	} else {

		c00Q := eval.buffQ[3][0]
		c00Q2 := eval.buffQMul[3][0]

		ringQ.MForm(c1Q1[0], c00Q)
		ringQMul.MForm(c1Q2[0], c00Q2)

		for i := 0; i < ct0.Degree()+1; i++ {
			ringQ.MulCoeffsMontgomery(c00Q, c0Q1[i], c2Q1[i])
			ringQMul.MulCoeffsMontgomery(c00Q2, c0Q2[i], c2Q2[i])
		}
	}
}

func (eval *evaluator) quantizeLvl(level, levelQMul int, ctOut *rlwe.Ciphertext) {

	c2Q1 := eval.buffQ[2]
	c2Q2 := eval.buffQMul[2]

	ringQ := eval.params.RingQ().AtLevel(level)
	ringQMul := eval.params.RingQMul().AtLevel(levelQMul)

	// Applies the inverse NTT to the ciphertext, scales down the ciphertext
	// by t/q and reduces its basis from QP to Q
	for i := range ctOut.Value {
		ringQ.INTTLazy(c2Q1[i], c2Q1[i])
		ringQMul.INTTLazy(c2Q2[i], c2Q2[i])

		// Extends the basis Q of ct(x) to the basis P and Divides (ct(x)Q -> P) by Q
		eval.basisExtenderQ1toQ2.ModDownQPtoP(level, levelQMul, c2Q1[i], c2Q2[i], c2Q2[i]) // QP / Q -> P

		// Centers ct(x)P by (P-1)/2 and extends ct(x)P to the basis Q
		ringQMul.AddScalarBigint(c2Q2[i], eval.qMulHalf[levelQMul], c2Q2[i])
		eval.basisExtenderQ1toQ2.ModUpPtoQ(levelQMul, level, c2Q2[i], ctOut.Value[i])
		ringQ.SubScalarBigint(ctOut.Value[i], eval.qMulHalf[levelQMul], ctOut.Value[i])

		// (ct(x)/Q)*T, doing so only requires that Q*P > Q*Q, faster but adds error ~|T|
		ringQ.MulScalarBigint(ctOut.Value[i], eval.params.RingT().Modulus(), ctOut.Value[i])
	}
}

// Mul multiplies ctIn by op1 and returns the result in ctOut.
func (eval *evaluator) Mul(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext) {
	eval.CheckBinary(ctIn, op1, ctOut, ctIn.Degree()+op1.Degree())
	switch op1 := op1.(type) {
	case *PlaintextMul:
		eval.mulPlaintextMul(ctIn, op1, ctOut)
	case *PlaintextRingT:
		eval.mulPlaintextRingT(ctIn, op1, ctOut)
	case *rlwe.Plaintext, *rlwe.Ciphertext:
		eval.tensorAndRescale(ctIn, op1.El(), ctOut)
	default:
		panic(fmt.Errorf("cannot Mul: invalid rlwe.Operand type for Mul: %T", op1))
	}
}

// MulThenAdd multiplies ctIn with op1 and adds the result on ctOut.
func (eval *evaluator) MulThenAdd(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext) {

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

	ringQ := eval.params.RingQ().AtLevel(utils.MinInt(ctIn.Level(), ctOut.Level()))

	for i := range ctIn.Value {
		ringQ.NTTLazy(ctIn.Value[i], ctOut.Value[i])
		ringQ.MulCoeffsMontgomeryLazy(ctOut.Value[i], ptRt.Value, ctOut.Value[i])
		ringQ.INTT(ctOut.Value[i], ctOut.Value[i])
	}
}

func (eval *evaluator) mulPlaintextRingT(ctIn *rlwe.Ciphertext, ptRt *PlaintextRingT, ctOut *rlwe.Ciphertext) {

	level := utils.MinInt(ctIn.Level(), ctOut.Level())

	ctOut.Resize(ctOut.Degree(), level)

	ringQ := eval.params.RingQ().AtLevel(level)

	coeffs := ptRt.Value.Coeffs[0]
	coeffsNTT := eval.buffQ[0][0].Coeffs[0]

	for i := range ctIn.Value {

		// Copies the inputCT on the outputCT and switches to the NTT domain
		ringQ.NTTLazy(ctIn.Value[i], ctOut.Value[i])

		// Switches the outputCT in the Montgomery domain
		ringQ.MForm(ctOut.Value[i], ctOut.Value[i])

		// For each qi in Q
		for j, s := range ringQ.SubRings[:ctIn.Level()+1] {

			tmp := ctOut.Value[i].Coeffs[j]

			// Transforms the plaintext in the NTT domain of that qi
			s.NTTLazy(coeffs, coeffsNTT)

			// Multiplies NTT_qi(pt) * NTT_qi(ct)
			s.MulCoeffsMontgomery(tmp, coeffsNTT, tmp)

		}

		// Switches the ciphertext out of the NTT domain
		ringQ.INTT(ctOut.Value[i], ctOut.Value[i])
	}
}

// MulNew multiplies ctIn by op1 and creates a new element ctOut to store the result.
func (eval *evaluator) MulNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, ctIn.Degree()+op1.Degree(), ctIn.Level())
	eval.Mul(ctIn, op1, ctOut)
	return
}

// RelinearizeNew relinearizes the ciphertext ctIn of degree > 1 until it is of degree 1, and creates a new ciphertext to store the result.
func (eval *evaluator) RelinearizeNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1, ctIn.Level())
	eval.Relinearize(ctIn, ctOut)
	return
}

// SwitchKeysNew applies the key-switching procedure to the ciphertext ct0 and creates a new ciphertext to store the result. It requires as an additional input a valid switching-key:
// it must encrypt the target key under the public key under which ct0 is currently encrypted.
func (eval *evaluator) SwitchKeysNew(ctIn *rlwe.Ciphertext, switchkey *rlwe.SwitchingKey) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1, ctIn.Level())
	eval.SwitchKeys(ctIn, switchkey, ctOut)
	return
}

// RotateColumns rotates the columns of ct0 by k positions to the left and returns the result in ctOut. As an additional input it requires a RotationKeys struct.
func (eval *evaluator) RotateColumns(ct0 *rlwe.Ciphertext, k int, ctOut *rlwe.Ciphertext) {
	eval.Automorphism(ct0, eval.params.GaloisElementForColumnRotationBy(k), ctOut)
}

// RotateColumnsNew applies RotateColumns and returns the result in a new Ciphertext.
func (eval *evaluator) RotateColumnsNew(ctIn *rlwe.Ciphertext, k int) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1, ctIn.Level())
	eval.RotateColumns(ctIn, k, ctOut)
	return
}

// RotateRows rotates the rows of ct0 and returns the result in ctOut.
func (eval *evaluator) RotateRows(ct0 *rlwe.Ciphertext, ctOut *rlwe.Ciphertext) {
	eval.Automorphism(ct0, eval.params.GaloisElementForRowRotation(), ctOut)
}

// RotateRowsNew rotates the rows of ctIn and returns the result a new Ciphertext.
func (eval *evaluator) RotateRowsNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1, ctIn.Level())
	eval.RotateRows(ctIn, ctOut)
	return
}

func (eval *evaluator) InnerSum(ctIn *rlwe.Ciphertext, batchSize, n int, ctOut *rlwe.Ciphertext) {
	_, level := eval.CheckUnary(ctIn, ctOut)
	eval.params.RingQ().AtLevel(level).NTT(ctIn.Value[0], ctOut.Value[0])
	eval.params.RingQ().AtLevel(level).NTT(ctIn.Value[1], ctOut.Value[1])
	ctOut.IsNTT = true
	eval.Evaluator.InnerSum(ctOut, batchSize, n, ctOut)
	eval.params.RingQ().AtLevel(level).INTT(ctOut.Value[0], ctOut.Value[0])
	eval.params.RingQ().AtLevel(level).INTT(ctOut.Value[1], ctOut.Value[1])
	ctOut.IsNTT = false
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

// evaluateInPlaceBinary applies the provided function in place on el0 and el1 and returns the result in elOut.
func (eval *evaluator) evaluateInPlaceBinary(el0, el1, elOut *rlwe.Ciphertext, evaluate func(*ring.Poly, *ring.Poly, *ring.Poly)) {

	smallest, largest, _ := rlwe.GetSmallestLargest(el0, el1)

	for i := 0; i < smallest.Degree()+1; i++ {
		evaluate(el0.Value[i], el1.Value[i], elOut.Value[i])
	}

	// If the inputs degrees differ, it copies the remaining degree on the receiver.
	if largest != nil && largest != elOut { // checks to avoid unnecessary work.
		for i := smallest.Degree() + 1; i < largest.Degree()+1; i++ {
			elOut.Value[i].Copy(largest.Value[i])
		}
	}
}

// evaluateInPlaceUnary applies the provided function in place on el0 and returns the result in elOut.
func evaluateInPlaceUnary(el0, elOut *rlwe.Ciphertext, evaluate func(*ring.Poly, *ring.Poly)) {
	for i := range el0.Value {
		evaluate(el0.Value[i], elOut.Value[i])
	}
}
