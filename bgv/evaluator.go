package bgv

import (
	"fmt"
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// Evaluator is an interface implementing the public methods of the eval.
type Evaluator interface {

	// Add, Sub, Neg ct-ct & ct-pt
	Add(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)
	AddNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext)
	Sub(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)
	SubNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext)
	Neg(ctIn *rlwe.Ciphertext, ctOut *rlwe.Ciphertext)
	NegNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext)

	// Add, Mul ct-const
	AddScalar(ctIn *rlwe.Ciphertext, scalar uint64, ctOut *rlwe.Ciphertext)
	AddScalarNew(ctIn *rlwe.Ciphertext, scalar uint64) (ctOut *rlwe.Ciphertext)
	MulScalar(ctIn *rlwe.Ciphertext, scalar uint64, ctOut *rlwe.Ciphertext)
	MulScalarNew(ctIn *rlwe.Ciphertext, scalar uint64) (ctOut *rlwe.Ciphertext)
	MulScalarAndAdd(ctIn *rlwe.Ciphertext, scalar uint64, ctOut *rlwe.Ciphertext)

	// Mul ct-ct & ct-pt
	MulNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext)
	Mul(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)
	MulRelinNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext)
	MulRelin(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)

	// MulAndAdd ct-ct & ct-pt
	MulAndAdd(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)
	MulRelinAndAdd(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)

	// Degree Management
	RelinearizeNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext)
	Relinearize(ctIn *rlwe.Ciphertext, ctOut *rlwe.Ciphertext)

	// Error and Level management
	Rescale(ctIn, ctOut *rlwe.Ciphertext) (err error)
	DropLevelNew(ctIn *rlwe.Ciphertext, levels int) (ctOut *rlwe.Ciphertext)
	DropLevel(ctIn *rlwe.Ciphertext, levels int)

	// Column & Rows rotations
	RotateColumnsNew(ctIn *rlwe.Ciphertext, k int) (ctOut *rlwe.Ciphertext)
	RotateColumns(ctIn *rlwe.Ciphertext, k int, ctOut *rlwe.Ciphertext)
	RotateRows(ctIn *rlwe.Ciphertext, ctOut *rlwe.Ciphertext)
	RotateRowsNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext)

	//Polynomial Evaluation
	EvaluatePoly(input interface{}, pol *Polynomial, targetScale rlwe.Scale) (ctOut *rlwe.Ciphertext, err error)
	EvaluatePolyVector(input interface{}, pols []*Polynomial, encoder Encoder, slotIndex map[int][]int, targetScale rlwe.Scale) (ctOut *rlwe.Ciphertext, err error)

	// TODO
	LinearTransformNew(ctIn *rlwe.Ciphertext, linearTransform interface{}) (ctOut []*rlwe.Ciphertext)
	LinearTransform(ctIn *rlwe.Ciphertext, linearTransform interface{}, ctOut []*rlwe.Ciphertext)
	MultiplyByDiagMatrix(ctIn *rlwe.Ciphertext, matrix LinearTransform, c2DecompQP []ringqp.Poly, ctOut *rlwe.Ciphertext)
	MultiplyByDiagMatrixBSGS(ctIn *rlwe.Ciphertext, matrix LinearTransform, c2DecompQP []ringqp.Poly, ctOut *rlwe.Ciphertext)
	InnerSumLog(ctIn *rlwe.Ciphertext, batch, n int, ctOut *rlwe.Ciphertext)
	ReplicateLog(ctIn *rlwe.Ciphertext, batch, n int, ctOut *rlwe.Ciphertext)

	// Key-Switching
	SwitchKeysNew(ctIn *rlwe.Ciphertext, swk *rlwe.SwitchingKey) (ctOut *rlwe.Ciphertext)
	SwitchKeys(ctIn *rlwe.Ciphertext, swk *rlwe.SwitchingKey, ctOut *rlwe.Ciphertext)
	Automorphism(ctIn *rlwe.Ciphertext, galEl uint64, ctOut *rlwe.Ciphertext)
	AutomorphismHoisted(level int, ctIn *rlwe.Ciphertext, c1DecompQP []ringqp.Poly, galEl uint64, ctOut *rlwe.Ciphertext)
	Merge(ctIn map[int]*rlwe.Ciphertext) (ctOut *rlwe.Ciphertext)

	// Others
	GetRLWEEvaluator() *rlwe.Evaluator
	BuffQ() [3]*ring.Poly
	ShallowCopy() Evaluator
	WithKey(rlwe.EvaluationKey) Evaluator
}

// evaluator is a struct that holds the necessary elements to perform the homomorphic operations between ciphertexts and/or plaintexts.
// It also holds a memory buffer used to store intermediate computations.
type evaluator struct {
	*evaluatorBase
	*evaluatorBuffers
	*rlwe.Evaluator
}

type evaluatorBase struct {
	params       Parameters
	qiInvModTNeg []uint64
	qLModqi      [][]uint64
	tInvModQ     []*big.Int
}

func newEvaluatorPrecomp(params Parameters) *evaluatorBase {
	ringQ := params.RingQ()
	ringT := params.RingT()
	t := params.T()

	qiInvModTNeg := make([]uint64, len(ringQ.Modulus))
	tInvModQ := make([]*big.Int, len(ringQ.Modulus))

	for i, qi := range ringQ.Modulus {
		qiInvModTNeg[i] = ring.MForm(t-ring.ModExp(qi, t-2, t), t, ringT.BredParams[0])
		tInvModQ[i] = ring.NewUint(t)
		tInvModQ[i].ModInverse(tInvModQ[i], ringQ.ModulusAtLevel[i])
	}

	qLModqi := make([][]uint64, len(ringQ.Modulus)-1)

	for j := len(ringQ.Modulus) - 1; j > 0; j-- {
		qLModqi[j-1] = make([]uint64, j)
		for i := 0; i < j; i++ {
			qLModqi[j-1][i] = ring.MForm(ringQ.Modulus[j], ringQ.Modulus[i], ringQ.BredParams[i])
		}
	}

	return &evaluatorBase{
		params:       params,
		qiInvModTNeg: qiInvModTNeg,
		qLModqi:      qLModqi,
		tInvModQ:     tInvModQ,
	}
}

type evaluatorBuffers struct {
	buffQ  [3]*ring.Poly
	buffCt *rlwe.Ciphertext
}

// BuffQ returns a pointer to the internal memory buffer buffQ.
func (eval *evaluator) BuffQ() [3]*ring.Poly {
	return eval.buffQ
}

// GetRLWEEvaluator returns the underlying *rlwe.Evaluator.
func (eval *evaluator) GetRLWEEvaluator() *rlwe.Evaluator {
	return eval.Evaluator
}

func newEvaluatorBuffer(eval *evaluatorBase) *evaluatorBuffers {
	ringQ := eval.params.RingQ()
	buffQ := [3]*ring.Poly{ringQ.NewPoly(), ringQ.NewPoly(), ringQ.NewPoly()}
	return &evaluatorBuffers{
		buffQ:  buffQ,
		buffCt: rlwe.NewCiphertext(eval.params.Parameters, 2, eval.params.MaxLevel()),
	}
}

// NewEvaluator creates a new Evaluator, that can be used to do homomorphic
// operations on ciphertexts and/or plaintexts. It stores a memory buffer
// and ciphertexts that will be used for intermediate values.
func NewEvaluator(params Parameters, evaluationKey rlwe.EvaluationKey) Evaluator {
	ev := new(evaluator)
	ev.evaluatorBase = newEvaluatorPrecomp(params)
	ev.evaluatorBuffers = newEvaluatorBuffer(ev.evaluatorBase)
	ev.Evaluator = rlwe.NewEvaluator(params.Parameters, &evaluationKey)

	return ev
}

// ShallowCopy creates a shallow copy of this evaluator in which the read-only data-structures are
// shared with the receiver.
func (eval *evaluator) ShallowCopy() Evaluator {
	return &evaluator{
		evaluatorBase:    eval.evaluatorBase,
		Evaluator:        eval.Evaluator.ShallowCopy(),
		evaluatorBuffers: newEvaluatorBuffer(eval.evaluatorBase),
	}
}

// WithKey creates a shallow copy of this evaluator in which the read-only data-structures are
// shared with the receiver but the EvaluationKey is evaluationKey.
func (eval *evaluator) WithKey(evaluationKey rlwe.EvaluationKey) Evaluator {
	return &evaluator{
		evaluatorBase:    eval.evaluatorBase,
		Evaluator:        eval.Evaluator.WithKey(&evaluationKey),
		evaluatorBuffers: eval.evaluatorBuffers,
	}
}

func (eval *evaluator) checkBinary(op0, op1, opOut rlwe.Operand, opOutMinDegree int) {
	if op0 == nil || op1 == nil || opOut == nil {
		panic("cannot checkBinary: rlwe.Operands cannot be nil")
	}

	if op0.Degree()+op1.Degree() == 0 {
		panic("cannot checkBinary: rlwe.Operands cannot be both plaintext")
	}

	if opOut.Degree() < opOutMinDegree {
		opOut.El().Resize(opOutMinDegree, opOut.Level())
	}

	if op0.Degree() > 2 || op1.Degree() > 2 || opOut.Degree() > 2 {
		panic("cannot checkBinary: rlwe.Operands degree cannot be larger than 2")
	}

	for !op0.El().IsNTT {
		panic("cannot checkBinary: op0 must be in NTT")
	}

	for !op1.El().IsNTT {
		panic("cannot checkBinary: op1 must be in NTT")
	}
}

func (eval *evaluator) evaluateInPlace(el0, el1, elOut *rlwe.Ciphertext, evaluate func(int, *ring.Poly, *ring.Poly, *ring.Poly)) {

	smallest, largest, _ := rlwe.GetSmallestLargest(el0.El(), el1.El())

	level := utils.MinInt(utils.MinInt(el0.Level(), el1.Level()), elOut.Level())

	elOut.Resize(elOut.Degree(), level)

	for i := 0; i < smallest.Degree()+1; i++ {
		evaluate(level, el0.Value[i], el1.Value[i], elOut.Value[i])
	}

	// If the inputs degrees differ, it copies the remaining degree on the receiver.
	if largest != nil && largest != elOut.El() { // checks to avoid unnecessary work.
		for i := smallest.Degree() + 1; i < largest.Degree()+1; i++ {
			elOut.Value[i].Copy(largest.Value[i])
		}
	}

	elOut.Scale = el0.Scale
}

func (eval *evaluator) matchScaleThenEvaluateInPlace(el0, el1, elOut *rlwe.Ciphertext, evaluate func(int, *ring.Poly, uint64, *ring.Poly)) {

	level := utils.MinInt(utils.MinInt(el0.Level(), el1.Level()), elOut.Level())

	r0, r1, _ := eval.matchScalesBinary(el0.Scale.Uint64(), el1.Scale.Uint64())

	for i := range el0.Value {
		eval.params.RingQ().MulScalarLvl(level, el0.Value[i], r0, elOut.Value[i])
	}

	for i := range el1.Value {
		evaluate(level, el1.Value[i], r1, elOut.Value[i])
	}

	elOut.Scale = el0.Scale.Mul(r0, eval.params.T())
}

func (eval *evaluator) newCiphertextBinary(op0, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext) {
	return rlwe.NewCiphertext(eval.params.Parameters, utils.MaxInt(op0.Degree(), op1.Degree()), utils.MinInt(op0.Level(), op1.Level()))
}

// Add adds op1 to ctIn and returns the result in ctOut.
func (eval *evaluator) Add(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext) {
	eval.checkBinary(ctIn, op1, ctOut, utils.MaxInt(ctIn.Degree(), op1.Degree()))

	if ctIn.Scale.Cmp(op1.GetScale()) == 0 {
		eval.evaluateInPlace(ctIn, op1.El(), ctOut, eval.params.RingQ().AddLvl)
	} else {
		eval.matchScaleThenEvaluateInPlace(ctIn, op1.El(), ctOut, eval.params.RingQ().MulScalarAndAddLvl)
	}
}

// AddNew adds op1 to ctIn and returns the result in a new ctOut.
func (eval *evaluator) AddNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext) {
	ctOut = eval.newCiphertextBinary(ctIn, op1)
	eval.Add(ctIn, op1, ctOut)
	return
}

// Sub subtracts op1 to ctIn and returns the result in ctOut.
func (eval *evaluator) Sub(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext) {
	eval.checkBinary(ctIn, op1, ctOut, utils.MaxInt(ctIn.Degree(), op1.Degree()))

	if ctIn.Scale.Cmp(op1.GetScale()) == 0 {
		eval.evaluateInPlace(ctIn, op1.El(), ctOut, eval.params.RingQ().SubLvl)
	} else {
		eval.matchScaleThenEvaluateInPlace(ctIn, op1.El(), ctOut, eval.params.RingQ().MulScalarAndSubLvl)
	}
}

// SubNew subtracts op1 to ctIn and returns the result in a new ctOut.
func (eval *evaluator) SubNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext) {
	ctOut = eval.newCiphertextBinary(ctIn, op1)
	eval.Sub(ctIn, op1, ctOut)
	return
}

// Neg negates ctIn and returns the result in ctOut.
func (eval *evaluator) Neg(ctIn *rlwe.Ciphertext, ctOut *rlwe.Ciphertext) {

	if ctIn.Degree() != ctOut.Degree() {
		panic("cannot Negate: invalid receiver Ciphertext does not match input Ciphertext degree")
	}

	level := utils.MinInt(ctIn.Level(), ctOut.Level())

	for i := range ctIn.Value {
		eval.params.RingQ().NegLvl(level, ctIn.Value[i], ctOut.Value[i])
	}

	ctOut.Scale = ctIn.Scale
}

// NegNew negates ctIn and returns the result in a new ctOut.
func (eval *evaluator) NegNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {
	ctOut = rlwe.NewCiphertext(eval.params.Parameters, ctIn.Degree(), ctIn.Level())
	eval.Neg(ctIn, ctOut)
	return
}

// AddScalar adds a scalar to ctIn and returns the result in ctOut.
func (eval *evaluator) AddScalar(ctIn *rlwe.Ciphertext, scalar uint64, ctOut *rlwe.Ciphertext) {

	ringQ := eval.params.RingQ()
	ringT := eval.params.RingT()

	if ctIn.Scale.Cmp(1) != 0 {
		scalar = ring.BRed(scalar, ctIn.Scale.Uint64(), ringT.Modulus[0], ringT.BredParams[0])
	}

	ringQ.AddScalarLvl(utils.MinInt(ctIn.Level(), ctOut.Level()), ctIn.Value[0], scalar, ctOut.Value[0])

	if ctIn != ctOut {
		for i := 1; i < ctIn.Degree()+1; i++ {
			ring.Copy(ctIn.Value[i], ctOut.Value[i])
		}

		ctOut.MetaData = ctIn.MetaData
	}
}

// AddScalarNew adds a scalar to ctIn and returns the result in a new ctOut.
func (eval *evaluator) AddScalarNew(ctIn *rlwe.Ciphertext, scalar uint64) (ctOut *rlwe.Ciphertext) {
	ctOut = rlwe.NewCiphertext(eval.params.Parameters, ctIn.Degree(), ctIn.Level())
	eval.AddScalar(ctIn, scalar, ctOut)
	return
}

// MulScalar multiplies ctIn with a scalar and returns the result in ctOut.
func (eval *evaluator) MulScalar(ctIn *rlwe.Ciphertext, scalar uint64, ctOut *rlwe.Ciphertext) {
	ringQ := eval.params.RingQ()
	for i := 0; i < ctIn.Degree()+1; i++ {
		ringQ.MulScalarLvl(utils.MinInt(ctIn.Level(), ctOut.Level()), ctIn.Value[i], scalar, ctOut.Value[i])
	}
	ctOut.Scale = ctIn.Scale
}

// MulScalarNew multiplies ctIn with a scalar and returns the result in a new ctOut.
func (eval *evaluator) MulScalarNew(ctIn *rlwe.Ciphertext, scalar uint64) (ctOut *rlwe.Ciphertext) {
	ctOut = rlwe.NewCiphertext(eval.params.Parameters, ctIn.Degree(), ctIn.Level())
	eval.MulScalar(ctIn, scalar, ctOut)
	return
}

// MulScalarAndAdd multiplies ctIn with a scalar adds the result on ctOut.
func (eval *evaluator) MulScalarAndAdd(ctIn *rlwe.Ciphertext, scalar uint64, ctOut *rlwe.Ciphertext) {
	ringQ := eval.params.RingQ()

	// scalar *= (ctOut.scale / ctIn.Scale)
	if ctIn.Scale.Cmp(ctOut.Scale) != 0 {
		ringT := eval.params.RingT()
		ratio := ring.ModExp(ctIn.Scale.Uint64(), ringT.Modulus[0]-2, ringT.Modulus[0])
		ratio = ring.BRed(ratio, ctOut.Scale.Uint64(), ringT.Modulus[0], ringT.BredParams[0])
		scalar = ring.BRed(ratio, scalar, ringT.Modulus[0], ringT.BredParams[0])
	}

	for i := 0; i < ctIn.Degree()+1; i++ {
		ringQ.MulScalarAndAddLvl(utils.MinInt(ctIn.Level(), ctOut.Level()), ctIn.Value[i], scalar, ctOut.Value[i])
	}
}

// DropLevel reduces the level of ctIn by levels and returns the result in ctIn.
// No rescaling is applied during this procedure.
func (eval *evaluator) DropLevel(ctIn *rlwe.Ciphertext, levels int) {
	ctIn.Resize(ctIn.Degree(), ctIn.Level()-levels)
}

// DropLevelNew reduces the level of ctIn by levels and returns the result in a new ctOut.
// No rescaling is applied during this procedure.
func (eval *evaluator) DropLevelNew(ctIn *rlwe.Ciphertext, levels int) (ctOut *rlwe.Ciphertext) {
	ctOut = ctIn.CopyNew()
	eval.DropLevel(ctOut, levels)
	return
}

// Mul multiplies ctIn with op1 without relinearization and returns the result in ctOut.
// The procedure will panic if either ctIn or op1 are have a degree higher than 1.
// The procedure will panic if ctOut.Degree != ctIn.Degree + op1.Degree.
func (eval *evaluator) Mul(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext) {
	eval.mulRelin(ctIn, op1, false, ctOut)
}

// MulNew multiplies ctIn with op1 without relinearization and returns the result in a new ctOut.
// The procedure will panic if either ctIn.Degree or op1.Degree > 1.
func (eval *evaluator) MulNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext) {
	ctOut = rlwe.NewCiphertext(eval.params.Parameters, ctIn.Degree()+op1.Degree(), utils.MinInt(ctIn.Level(), op1.Level()))
	eval.mulRelin(ctIn, op1, false, ctOut)
	return
}

// MulRelinNew multiplies ctIn with op1 with relinearization and returns the result in a new ctOut.
// The procedure will panic if either ctIn.Degree or op1.Degree > 1.
// The procedure will panic if the evaluator was not created with an relinearization key.
func (eval *evaluator) MulRelinNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext) {
	ctOut = rlwe.NewCiphertext(eval.params.Parameters, 1, utils.MinInt(ctIn.Level(), op1.Level()))
	eval.mulRelin(ctIn, op1, true, ctOut)
	return
}

// MulRelin multiplies ctIn with op1 with relinearization and returns the result in ctOut.
// The procedure will panic if either ctIn.Degree or op1.Degree > 1.
// The procedure will panic if ctOut.Degree != ctIn.Degree + op1.Degree.
// The procedure will panic if the evaluator was not created with an relinearization key.
func (eval *evaluator) MulRelin(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext) {
	eval.mulRelin(ctIn, op1, true, ctOut)
}

func (eval *evaluator) mulRelin(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, relin bool, ctOut *rlwe.Ciphertext) {

	eval.checkBinary(ctIn, op1, ctOut, utils.MaxInt(ctIn.Degree(), op1.Degree()))

	level := utils.MinInt(utils.MinInt(ctIn.Level(), op1.Level()), ctOut.Level())

	if ctOut.Level() > level {
		eval.DropLevel(ctOut, ctOut.Level()-level)
	}

	if ctIn.Degree()+op1.Degree() > 2 {
		panic("cannot MulRelin: input elements total degree cannot be larger than 2")
	}

	ctOut.MetaData = ctIn.MetaData
	ctOut.Scale = ctIn.Scale.Mul(op1.GetScale(), eval.params.T())

	ringQ := eval.params.RingQ()

	var c00, c01, c0, c1, c2 *ring.Poly

	// Case Ciphertext (x) Ciphertext
	if ctIn.Degree() == 1 && op1.Degree() == 1 {

		c00 = eval.buffQ[0]
		c01 = eval.buffQ[1]

		c0 = ctOut.Value[0]
		c1 = ctOut.Value[1]

		if !relin {
			if ctOut.Degree() < 2 {
				ctOut.Resize(2, ctOut.Level())
			}
			c2 = ctOut.Value[2]
		} else {
			c2 = eval.buffQ[2]
		}

		// Avoid overwriting if the second input is the output
		var tmp0, tmp1 *rlwe.Ciphertext
		if op1.El() == ctOut.El() {
			tmp0, tmp1 = op1.El(), ctIn.El()
		} else {
			tmp0, tmp1 = ctIn.El(), op1.El()
		}

		ringQ.MFormLvl(level, tmp0.Value[0], c00)
		ringQ.MFormLvl(level, tmp0.Value[1], c01)

		if ctIn.El() == op1.El() { // squaring case
			ringQ.MulCoeffsMontgomeryLvl(level, c00, tmp1.Value[0], c0) // c0 = c[0]*c[0]
			ringQ.MulCoeffsMontgomeryLvl(level, c01, tmp1.Value[1], c2) // c2 = c[1]*c[1]
			ringQ.MulCoeffsMontgomeryLvl(level, c00, tmp1.Value[1], c1) // c1 = 2*c[0]*c[1]
			ringQ.AddLvl(level, c1, c1, c1)

		} else { // regular case
			ringQ.MulCoeffsMontgomeryLvl(level, c00, tmp1.Value[0], c0) // c0 = c0[0]*c0[0]
			ringQ.MulCoeffsMontgomeryLvl(level, c01, tmp1.Value[1], c2) // c2 = c0[1]*c1[1]
			ringQ.MulCoeffsMontgomeryLvl(level, c00, tmp1.Value[1], c1)
			ringQ.MulCoeffsMontgomeryAndAddLvl(level, c01, tmp1.Value[0], c1) // c1 = c0[0]*c1[1] + c0[1]*c1[0]
		}

		if relin {

			// Yup...
			ringQ.MulScalarBigintLvl(level, c2, eval.tInvModQ[level], c2)

			tmpCt := &rlwe.Ciphertext{Value: []*ring.Poly{eval.BuffQP[1].Q, eval.BuffQP[2].Q}}
			tmpCt.IsNTT = true

			eval.GadgetProduct(level, c2, eval.Rlk.Keys[0].GadgetCiphertext, tmpCt)

			// It works >.>
			ringQ.MulScalarAndAddLvl(level, tmpCt.Value[0], eval.params.T(), ctOut.Value[0])
			ringQ.MulScalarAndAddLvl(level, tmpCt.Value[1], eval.params.T(), ctOut.Value[1])
		}

		// Case Plaintext (x) Ciphertext or Ciphertext (x) Plaintext
	} else {

		if ctOut.Degree() < ctIn.Degree() {
			ctOut.Resize(ctIn.Degree(), level)
		}

		c00 := eval.buffQ[0]

		ringQ.MFormLvl(level, op1.El().Value[0], c00)
		for i := range ctOut.Value {
			ringQ.MulCoeffsMontgomeryLvl(level, ctIn.Value[i], c00, ctOut.Value[i])
		}
	}
}

// MulAndAdd multiplies ctIn with op1 (without relinearization)^and adds the result on ctOut.
// The procedure will panic if either ctIn.Degree() or op1.Degree() > 1.
// The procedure will panic if either ctIn == ctOut or op1 == ctOut.
func (eval *evaluator) MulAndAdd(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext) {
	eval.mulRelinAndAdd(ctIn, op1, false, ctOut)
}

// MulRelinAndAdd multiplies ctIn with op1 and adds, relinearize the result on ctOut.
// The procedure will panic if either ctIn.Degree() or op1.Degree() > 1.
// The procedure will panic if either ctIn == ctOut or op1 == ctOut.
func (eval *evaluator) MulRelinAndAdd(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext) {
	eval.mulRelinAndAdd(ctIn, op1, true, ctOut)
}

func (eval *evaluator) mulRelinAndAdd(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, relin bool, ctOut *rlwe.Ciphertext) {

	eval.checkBinary(ctIn, op1, ctOut, utils.MaxInt(ctIn.Degree(), op1.Degree()))

	level := utils.MinInt(utils.MinInt(ctIn.Level(), op1.Level()), ctOut.Level())

	if ctIn.Degree()+op1.Degree() > 2 {
		panic("cannot MulRelinAndAdd: input elements total degree cannot be larger than 2")
	}

	if ctIn.El() == ctOut.El() || op1.El() == ctOut.El() {
		panic("cannot MulRelinAndAdd: ctOut must be different from ctIn and op1")
	}

	ringQ := eval.params.RingQ()
	ringT := eval.params.RingT()

	var c00, c01, c0, c1, c2 *ring.Poly

	// Case Ciphertext (x) Ciphertext
	if ctIn.Degree() == 1 && op1.Degree() == 1 {

		c00 = eval.buffQ[0]
		c01 = eval.buffQ[1]

		c0 = ctOut.Value[0]
		c1 = ctOut.Value[1]

		if !relin {
			ctOut.Resize(2, level)
			c2 = ctOut.Value[2]
		} else {
			ctOut.Resize(1, level)
			c2 = eval.buffQ[2]
		}

		tmp0, tmp1 := ctIn.El(), op1.El()

		var r0 uint64 = 1
		if targetScale := ring.BRed(ctIn.Scale.Uint64(), op1.GetScale().Uint64(), ringT.Modulus[0], ringT.BredParams[0]); ctOut.Scale.Cmp(targetScale) != 0 {
			var r1 uint64
			r0, r1, _ = eval.matchScalesBinary(targetScale, ctOut.Scale.Uint64())

			for i := range ctOut.Value {
				ringQ.MulScalarLvl(level, ctOut.Value[i], r1, ctOut.Value[i])
			}

			ctOut.Scale = ctOut.Scale.Mul(r1, eval.params.T())
		}

		ringQ.MFormLvl(level, tmp0.Value[0], c00)
		ringQ.MFormLvl(level, tmp0.Value[1], c01)

		if r0 != 1 {
			ringQ.MulScalarLvl(level, c00, r0, c00)
			ringQ.MulScalarLvl(level, c01, r0, c01)
		}

		ringQ.MulCoeffsMontgomeryAndAddLvl(level, c00, tmp1.Value[0], c0) // c0 += c[0]*c[0]
		ringQ.MulCoeffsMontgomeryAndAddLvl(level, c00, tmp1.Value[1], c1) // c1 += c[0]*c[1]
		ringQ.MulCoeffsMontgomeryAndAddLvl(level, c01, tmp1.Value[0], c1) // c1 += c[1]*c[0]

		if relin {
			ringQ.MulCoeffsMontgomeryLvl(level, c01, tmp1.Value[1], c2) // c2 += c[1]*c[1]

			ringQ.MulScalarBigintLvl(level, c2, eval.tInvModQ[level], c2)

			tmpCt := &rlwe.Ciphertext{Value: []*ring.Poly{eval.BuffQP[1].Q, eval.BuffQP[2].Q}}
			tmpCt.IsNTT = true

			eval.GadgetProduct(level, c2, eval.Rlk.Keys[0].GadgetCiphertext, tmpCt)

			ringQ.MulScalarAndAddLvl(level, tmpCt.Value[0], eval.params.T(), c0)
			ringQ.MulScalarAndAddLvl(level, tmpCt.Value[1], eval.params.T(), c1)

		} else {
			ringQ.MulCoeffsMontgomeryAndAddLvl(level, c01, tmp1.Value[1], c2) // c2 += c[1]*c[1]
		}

		// Case Plaintext (x) Ciphertext or Ciphertext (x) Plaintext
	} else {

		if ctOut.Degree() < ctIn.Degree() {
			ctOut.Resize(ctIn.Degree(), level)
		}

		c00 := eval.buffQ[0]

		ringQ.MFormLvl(level, op1.El().Value[0], c00)
		for i := range ctIn.Value {
			ringQ.MulCoeffsMontgomeryAndAddLvl(level, ctIn.Value[i], c00, ctOut.Value[i])
		}
	}
}

func (eval *evaluator) rescaleOriginal(ctIn, ctOut *rlwe.Ciphertext) (err error) {

	if ctIn.Level() == 0 {
		return fmt.Errorf("cannot rescale: ctIn already at level 0")
	}

	if ctOut.Level() < ctIn.Level()-1 {
		return fmt.Errorf("cannot rescale: ctOut.Level() < ctIn.Level()-1")
	}

	level := ctIn.Level()
	ringQ := eval.params.RingQ()
	ringT := eval.params.RingT()

	buff0 := eval.BuffQP[0].Q.Coeffs[0]
	buff1 := eval.BuffQP[1].Q.Coeffs[0]
	buff2 := eval.BuffQP[2].Q.Coeffs[0]

	for i, el := range ctIn.Value {

		// buff0 = coeffs[level]
		ringQ.InvNTTSingleLazy(level, el.Coeffs[level], buff0)

		// buff1 = (buff0 * -qL^-1) % t
		ring.MulScalarMontgomeryVec(buff0, buff1, eval.qiInvModTNeg[level], ringT.Modulus[0], ringT.MredParams[0])

		for j := 0; j < level; j++ {

			qi := ringQ.Modulus[j]
			qLModqi := eval.qLModqi[level-1][j]
			mredParams := ringQ.MredParams[j]
			cIn := el.Coeffs[j]
			cOut := ctOut.Value[i].Coeffs[j]

			// buff2 = (buff1 * qL) % qi
			ring.MulScalarMontgomeryVec(buff1, buff2, qLModqi, qi, mredParams)

			// buff2 = buff2 + buff0
			ring.AddVecNoMod(buff2, buff0, buff2)
			ringQ.NTTSingleLazy(j, buff2, buff2)

			// cOut = ((buff2 + 2*qi - cIn) * -qL^-1) % qi
			ring.SubVecAndMulScalarMontgomeryTwoQiVec(buff2, cIn, cOut, ringQ.RescaleParams[level-1][j], qi, mredParams)

		}
	}

	ctOut.Resize(ctOut.Degree(), level-1)
	ctOut.MetaData = ctIn.MetaData
	ctOut.Scale = ctIn.Scale.Div(ringQ.Modulus[level], ringT.Modulus[0])

	return
}

// Rescale divides (rounded) ctIn by the last modulus of the moduli chain and returns the result on ctOut.
// This procedure divides the error by the last modulus of the moduli chain while preserving
// the LSB-plaintext bits.
// The procedure will return an error if:
//  1. ctIn.Level() == 0 (the input ciphertext is already at the last modulus)
//  2. ctOut.Level() < ctIn.Level() - 1 (not enough space to store the result)
func (eval *evaluator) Rescale(ctIn, ctOut *rlwe.Ciphertext) (err error) {

	if ctIn.Level() == 0 {
		return fmt.Errorf("cannot rescale: ctIn already at level 0")
	}

	if ctOut.Level() < ctIn.Level()-1 {
		return fmt.Errorf("cannot rescale: ctOut.Level() < ctIn.Level()-1")
	}

	level := ctIn.Level()
	ringQ := eval.params.RingQ()
	ringT := eval.params.RingT()

	// Center by (qL-1)/2
	qL := ringQ.Modulus[level]
	qLHalf := (qL - 1) >> 1
	T := ringT.Modulus[0]

	buff0 := eval.BuffQP[0].Q
	buff1 := eval.BuffQP[1].Q

	for i := range ctIn.Value {

		ringQ.MulScalarBigintLvl(level, ctIn.Value[i], eval.tInvModQ[level], buff0)
		ringQ.InvNTTSingleLazy(level, buff0.Coeffs[level], buff1.Coeffs[level])
		ring.AddScalarVec(buff1.Coeffs[level], buff1.Coeffs[level], qLHalf, qL)

		for j := 0; j < level; j++ {

			qi := ringQ.Modulus[j]
			bredParams := ringQ.BredParams[j]
			mredParams := ringQ.MredParams[j]
			rescaleParams := ring.BRed(T, ringQ.RescaleParams[level-1][j], qi, bredParams)

			ring.AddScalarNoModVec(buff1.Coeffs[level], buff1.Coeffs[j], qi-ring.BRedAdd(qLHalf, qi, bredParams))
			ringQ.NTTSingleLazy(j, buff1.Coeffs[j], buff1.Coeffs[j])
			ring.SubVecAndMulScalarMontgomeryTwoQiVec(buff1.Coeffs[j], buff0.Coeffs[j], ctOut.Value[i].Coeffs[j], rescaleParams, qi, mredParams)
		}
	}

	ctOut.Resize(ctOut.Degree(), level-1)
	ctOut.MetaData = ctIn.MetaData
	ctOut.Scale = ctIn.Scale.Div(ringQ.Modulus[level], ringT.Modulus[0])
	return
}

// RelinearizeNew applies the relinearization procedure on ctIn and returns the result in a new ctOut.
func (eval *evaluator) RelinearizeNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {
	ctOut = rlwe.NewCiphertext(eval.params.Parameters, 1, ctIn.Level())
	eval.Relinearize(ctIn, ctOut)
	return
}

// Relinearize applies the relinearization procedure on ctIn and returns the result in ctOut.
func (eval *evaluator) Relinearize(ctIn *rlwe.Ciphertext, ctOut *rlwe.Ciphertext) {

	if eval.Rlk == nil || ctIn.Degree()-1 > len(eval.Rlk.Keys) {
		panic("cannot Relinearize: relinearization key missing (or ciphertext degree is too large)")
	}

	level := utils.MinInt(ctIn.Level(), ctOut.Level())

	T := eval.params.T()

	ringQ := eval.params.RingQ()

	tmpCt := &rlwe.Ciphertext{Value: []*ring.Poly{eval.BuffQP[1].Q, eval.BuffQP[2].Q}}
	tmpCt.IsNTT = true

	ringQ.MulScalarBigintLvl(level, ctIn.Value[2], eval.tInvModQ[level], eval.buffQ[0])
	eval.GadgetProduct(level, eval.buffQ[0], eval.Rlk.Keys[0].GadgetCiphertext, tmpCt)

	if ctIn != ctOut {
		ring.Copy(ctIn.Value[0], ctOut.Value[0])
		ring.Copy(ctIn.Value[1], ctOut.Value[1])
		ctOut.MetaData = ctIn.MetaData
	}

	ringQ.MulScalarAndAddLvl(level, eval.BuffQP[1].Q, T, ctOut.Value[0])
	ringQ.MulScalarAndAddLvl(level, eval.BuffQP[2].Q, T, ctOut.Value[1])

	for deg := ctIn.Degree() - 1; deg > 1; deg-- {
		ringQ.MulScalarBigintLvl(level, ctIn.Value[deg], eval.tInvModQ[level], eval.buffQ[0])
		eval.GadgetProduct(level, eval.buffQ[0], eval.Rlk.Keys[deg-2].GadgetCiphertext, tmpCt)
		ringQ.MulScalarAndAddLvl(level, eval.BuffQP[1].Q, T, ctOut.Value[0])
		ringQ.MulScalarAndAddLvl(level, eval.BuffQP[2].Q, T, ctOut.Value[1])
	}

	ctOut.Resize(1, level)

	ctOut.MetaData = ctIn.MetaData
}

// SwitchKeysNew re-encrypts ctIn under a different key and returns the result in a new ctOut.
// It requires a SwitchingKey, which is computed from the key under which the Ciphertext is currently encrypted,
// and the key under which the Ciphertext will be re-encrypted.
// The procedure will panic if either ctIn.Degree() or ctOut.Degree() != 1.
func (eval *evaluator) SwitchKeysNew(ctIn *rlwe.Ciphertext, swk *rlwe.SwitchingKey) (ctOut *rlwe.Ciphertext) {
	ctOut = rlwe.NewCiphertext(eval.params.Parameters, ctIn.Degree(), ctIn.Level())
	eval.SwitchKeys(ctIn, swk, ctOut)
	return
}

// SwitchKeys re-encrypts ctIn under a different key and returns the result in ctOut.
// It requires a SwitchingKey, which is computed from the key under which the Ciphertext is currently encrypted,
// and the key under which the Ciphertext will be re-encrypted.
// The procedure will panic if either ctIn.Degree() or ctOut.Degree() != 1.
func (eval *evaluator) SwitchKeys(ctIn *rlwe.Ciphertext, swk *rlwe.SwitchingKey, ctOut *rlwe.Ciphertext) {

	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		panic("cannot SwitchKeys: input and output Ciphertext must be of degree 1")
	}

	level := utils.MinInt(ctIn.Level(), ctOut.Level())
	ringQ := eval.params.RingQ()

	ringQ.MulScalarBigintLvl(level, ctIn.Value[1], eval.tInvModQ[level], eval.buffQ[0])

	tmpCt := &rlwe.Ciphertext{Value: []*ring.Poly{eval.BuffQP[1].Q, eval.BuffQP[2].Q}}
	tmpCt.IsNTT = true

	eval.GadgetProduct(level, eval.buffQ[0], swk.GadgetCiphertext, tmpCt)

	if ctOut == ctIn {
		ringQ.MulScalarAndAddLvl(level, eval.BuffQP[1].Q, eval.params.T(), ctOut.Value[0])
	} else {
		ringQ.MulScalarLvl(level, eval.BuffQP[1].Q, eval.params.T(), ctOut.Value[0])
		ringQ.AddLvl(level, ctOut.Value[0], ctIn.Value[0], ctOut.Value[0])
	}

	ringQ.MulScalarLvl(level, eval.BuffQP[2].Q, eval.params.T(), ctOut.Value[1])

	ctOut.MetaData = ctIn.MetaData
}

// RotateColumnsNew rotates the columns of ctIn by k positions to the left, and returns the result in a newly created element.
// The procedure will panic if the corresponding Galois key has not been generated and attributed to the evaluator.
// The procedure will panic if ctIn.Degree() != 1.
func (eval *evaluator) RotateColumnsNew(ctIn *rlwe.Ciphertext, k int) (ctOut *rlwe.Ciphertext) {
	ctOut = rlwe.NewCiphertext(eval.params.Parameters, ctIn.Degree(), ctIn.Level())
	eval.RotateColumns(ctIn, k, ctOut)
	return
}

// RotateColumns rotates the columns of ctIn by k positions to the left and returns the result in ctOut.
// The procedure will panic if the corresponding Galois key has not been generated and attributed to the evaluator.
// The procedure will panic if either ctIn.Degree() or ctOut.Degree() != 1.
func (eval *evaluator) RotateColumns(ctIn *rlwe.Ciphertext, k int, ctOut *rlwe.Ciphertext) {
	eval.Automorphism(ctIn, eval.params.GaloisElementForColumnRotationBy(k), ctOut)
}

// RotateRowsNew swaps the rows of ctIn and returns the result in a new ctOut.
// The procedure will panic if the corresponding Galois key has not been generated and attributed to the evaluator.
// The procedure will panic if ctIn.Degree() != 1.
func (eval *evaluator) RotateRowsNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {
	ctOut = rlwe.NewCiphertext(eval.params.Parameters, ctIn.Degree(), ctIn.Level())
	eval.RotateRows(ctIn, ctOut)
	return
}

// RotateRows swaps the rows of ctIn and returns the result in ctOut.
// The procedure will panic if the corresponding Galois key has not been generated and attributed to the evaluator.
// The procedure will panic if either ctIn.Degree() or ctOut.Degree() != 1.
func (eval *evaluator) RotateRows(ctIn *rlwe.Ciphertext, ctOut *rlwe.Ciphertext) {
	eval.Automorphism(ctIn, eval.params.GaloisElementForRowRotation(), ctOut)
}

// Automorphism applies the automorphism X -> X^galEl mod X^N + 1 to ctIn and writes the result on ctOut.
// The procedure will panic if the corresponding Galois key has not been generated and attributed to the evaluator.
// The procedure will panic if ctIn.Degree() != 1.
func (eval *evaluator) Automorphism(ctIn *rlwe.Ciphertext, galEl uint64, ctOut *rlwe.Ciphertext) {

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
		panic(fmt.Sprintf("cannot apply Automorphism: galEl key 5^%d missing", eval.params.InverseGaloisElement(galEl)))
	}

	level := utils.MinInt(ctIn.Level(), ctOut.Level())

	ringQ := eval.params.RingQ()

	ringQ.MulScalarBigintLvl(level, ctIn.Value[1], eval.tInvModQ[level], eval.buffQ[0])

	tmpCt := &rlwe.Ciphertext{Value: []*ring.Poly{eval.BuffQP[1].Q, eval.BuffQP[2].Q}}
	tmpCt.IsNTT = true

	eval.GadgetProduct(level, eval.buffQ[0], rtk.GadgetCiphertext, tmpCt)

	ringQ.MulScalarLvl(level, eval.BuffQP[1].Q, eval.params.T(), eval.BuffQP[1].Q)
	ringQ.MulScalarLvl(level, eval.BuffQP[2].Q, eval.params.T(), eval.BuffQP[2].Q)

	ringQ.AddLvl(level, eval.BuffQP[1].Q, ctIn.Value[0], eval.BuffQP[1].Q)

	ringQ.PermuteNTTWithIndexLvl(level, eval.BuffQP[1].Q, eval.PermuteNTTIndex[galEl], ctOut.Value[0])
	ringQ.PermuteNTTWithIndexLvl(level, eval.BuffQP[2].Q, eval.PermuteNTTIndex[galEl], ctOut.Value[1])

	ctOut.Resize(ctOut.Degree(), level)
	ctOut.MetaData = ctIn.MetaData
}

// Automorphism applies the automorphism X -> X^galEl mod X^N + 1 to ctIn and writes the result on ctOut.
// The procedure also requires the decomposition of ctIn, given as c1DecompQP and which can be computed with rlwe.Evaluator.DecomposeNTT(*).
// The procedure will panic if the corresponding Galois key has not been generated and attributed to the evaluator.
// The procedure will panic if ctIn.Degree() != 1.
func (eval *evaluator) AutomorphismHoisted(level int, ctIn *rlwe.Ciphertext, c1DecompQP []ringqp.Poly, galEl uint64, ctOut *rlwe.Ciphertext) {

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
		panic(fmt.Sprintf("cannot apply AutomorphismHoisted: galEl key 5^%d missing", eval.params.InverseGaloisElement(galEl)))
	}

	ringQ := eval.params.RingQ()

	eval.KeyswitchHoisted(level, c1DecompQP, rtk, eval.BuffQP[0].Q, eval.BuffQP[1].Q, eval.BuffQP[0].P, eval.BuffQP[1].P)

	ringQ.MulScalarLvl(level, eval.BuffQP[0].Q, eval.params.T(), eval.BuffQP[0].Q)
	ringQ.MulScalarLvl(level, eval.BuffQP[1].Q, eval.params.T(), eval.BuffQP[1].Q)

	ringQ.AddLvl(level, eval.BuffQP[0].Q, ctIn.Value[0], eval.BuffQP[0].Q)

	ringQ.PermuteNTTWithIndexLvl(level, eval.BuffQP[0].Q, eval.PermuteNTTIndex[galEl], ctOut.Value[0])
	ringQ.PermuteNTTWithIndexLvl(level, eval.BuffQP[1].Q, eval.PermuteNTTIndex[galEl], ctOut.Value[1])

	ctOut.Resize(ctOut.Degree(), level)

	ctOut.MetaData = ctIn.MetaData
}

func (eval *evaluator) rotateHoistedNoModDownNew(level int, rotations []int, c0 *ring.Poly, c2DecompQP []ringqp.Poly) (cOut map[int]rlwe.CiphertextQP) {
	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()
	cOut = make(map[int]rlwe.CiphertextQP)
	for _, i := range rotations {
		if i != 0 {
			cOut[i] = rlwe.CiphertextQP{Value: [2]ringqp.Poly{{Q: ringQ.NewPolyLvl(level), P: ringP.NewPoly()}, {Q: ringQ.NewPolyLvl(level), P: ringP.NewPoly()}}, MetaData: rlwe.MetaData{IsNTT: true}}
			eval.AutomorphismHoistedNoModDown(level, c0, c2DecompQP, eval.params.GaloisElementForColumnRotationBy(i), cOut[i])
		}
	}

	return
}

// MatchScales updates the both input ciphertexts to ensures that their scale matches.
// To do so it computes t0 * a = ct1 * b such that:
// - ct0.scale * a = ct1.scale: make the scales match.
// - gcd(a, T) == gcd(b, T) == 1: ensure that the new scale is not a zero divisor if T is not prime.
// - |a+b| is minimal: minimize the added noise by the procedure.
func (eval *evaluator) MatchScalesAndLevel(ct0, ct1 *rlwe.Ciphertext) {

	r0, r1, _ := eval.matchScalesBinary(ct0.Scale.Uint64(), ct1.Scale.Uint64())

	level := utils.MinInt(ct0.Level(), ct1.Level())

	ringQ := eval.params.RingQ()

	for _, el := range ct0.Value {
		ringQ.MulScalarLvl(level, el, r0, el)
	}

	ct0.Resize(ct0.Degree(), level)
	ct0.Scale = ct0.Scale.Mul(r0, eval.params.T())

	for _, el := range ct1.Value {
		ringQ.MulScalarLvl(level, el, r1, el)
	}

	ct1.Resize(ct1.Degree(), level)
	ct1.Scale = ct1.Scale.Mul(r1, eval.params.T())
}

func (eval *evaluator) matchScalesBinary(scale0, scale1 uint64) (r0, r1, e uint64) {

	ringT := eval.params.RingT()

	t := ringT.Modulus[0]
	tHalf := t >> 1
	bredParams := ringT.BredParams[0]

	if utils.GCD(scale0, t) != 1 {
		panic("cannot matchScalesBinary: invalid ciphertext scale: gcd(scale, t) != 1")
	}

	var a = ringT.Modulus[0]
	var b uint64 = 0
	var A = ring.BRed(ring.ModExp(scale0, t-2, t), scale1, t, bredParams)
	var B uint64 = 1

	e = center(A, tHalf, t) + 1

	for A != 0 {
		q := a / A
		a, A = A, a%A
		b, B = B, ring.CRed(t+b-ring.BRed(B, q, t, bredParams), t)

		if A != 0 && utils.GCD(A, t) == 1 {
			tmp := center(A, tHalf, t) + center(B, tHalf, t)
			if tmp < e {
				e = tmp
				r0, r1 = A, B
			}
		}
	}

	return
}

func center(x, thalf, t uint64) uint64 {
	if x >= thalf {
		return t - x
	}
	return x
}

// MergeE merges a batch of Ciphertexts, packing the first coefficient of each Ciphertext into a single Ciphertext.
// The operation will require N/gap + log(gap) key-switches, where gap is the minimum gap between
// two non-zero coefficients of the final ciphertext.
// The method takes as input a map of Ciphertexts, indexing in which coefficient of the final
// Ciphertext the first coefficient of each Ciphertext of the map must be packed.
// The method will update the scale of all the Ciphertexts in the map to match the scale
// of the Ciphertext with the smallest index
func (eval *evaluator) Merge(ctIn map[int]*rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {

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

	ciphertextslist := make([]*rlwe.Ciphertext, ringQ.N)

	for i := range ctIn {
		ciphertextslist[i] = ctIn[i]
	}

	var scale0 uint64

	// Fetch the scale of the first ct
	for _, ct := range ciphertextslist {
		if ct != nil {
			scale0 = ct.Scale.Uint64()
			break
		}
	}

	t := params.T()
	bredParams := params.RingT().BredParams[0]
	for _, ct := range ciphertextslist {

		if ct != nil {
			if scale1 := ct.Scale.Uint64(); scale1 != 0 {

				update := ring.BRed(ring.ModExp(scale0, t-2, t), scale1, t, bredParams)

				ringQ.MulScalarLvl(ct.Level(), ct.Value[0], update, ct.Value[0])
				ringQ.MulScalarLvl(ct.Level(), ct.Value[1], update, ct.Value[1])

				ct.Scale = rlwe.NewScale(scale0)

			}
		}
	}

	if ciphertextslist[0] == nil {
		ciphertextslist[0] = rlwe.NewCiphertext(params.Parameters, 1, levelQ)
		ciphertextslist[0].Scale = rlwe.NewScale(scale0)
	}

	return eval.mergeRLWERecurse(ciphertextslist, xPow2)
}

func (eval *evaluator) mergeRLWERecurse(ciphertexts []*rlwe.Ciphertext, xPow []*ring.Poly) *rlwe.Ciphertext {

	ringQ := eval.params.RingQ()

	L := bits.Len64(uint64(len(ciphertexts))) - 1

	if L == 0 {
		return ciphertexts[0]
	}

	odd := make([]*rlwe.Ciphertext, len(ciphertexts)>>1)
	even := make([]*rlwe.Ciphertext, len(ciphertexts)>>1)

	for i := 0; i < len(ciphertexts)>>1; i++ {
		odd[i] = ciphertexts[2*i]
		even[i] = ciphertexts[2*i+1]
	}

	ctEven := eval.mergeRLWERecurse(odd, xPow)
	ctOdd := eval.mergeRLWERecurse(even, xPow)

	if ctEven == nil && ctOdd == nil {
		return nil
	}

	var tmpEven *rlwe.Ciphertext
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
