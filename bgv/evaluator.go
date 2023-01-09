package bgv

import (
	"fmt"
	"math/big"

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
	MulScalarThenAdd(ctIn *rlwe.Ciphertext, scalar uint64, ctOut *rlwe.Ciphertext)

	// Mul ct-ct & ct-pt
	MulNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext)
	Mul(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)
	MulRelinNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext)
	MulRelin(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)

	// MulThenAdd ct-ct & ct-pt
	MulThenAdd(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)
	MulRelinThenAdd(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext)

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
	InnerSum(ctIn *rlwe.Ciphertext, batch, n int, ctOut *rlwe.Ciphertext)
	Replicate(ctIn *rlwe.Ciphertext, batch, n int, ctOut *rlwe.Ciphertext)

	// Key-Switching
	SwitchKeysNew(ctIn *rlwe.Ciphertext, swk *rlwe.SwitchingKey) (ctOut *rlwe.Ciphertext)
	SwitchKeys(ctIn *rlwe.Ciphertext, swk *rlwe.SwitchingKey, ctOut *rlwe.Ciphertext)
	Automorphism(ctIn *rlwe.Ciphertext, galEl uint64, ctOut *rlwe.Ciphertext)
	AutomorphismHoisted(level int, ctIn *rlwe.Ciphertext, c1DecompQP []ringqp.Poly, galEl uint64, ctOut *rlwe.Ciphertext)
	RotateHoistedLazyNew(level int, rotations []int, c0 *ring.Poly, c2DecompQP []ringqp.Poly) (cOut map[int]rlwe.CiphertextQP)
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
	params   Parameters
	tInvModQ []*big.Int
}

func newEvaluatorPrecomp(params Parameters) *evaluatorBase {
	ringQ := params.RingQ()
	t := params.T()

	tInvModQ := make([]*big.Int, ringQ.ModuliChainLength())

	for i := range tInvModQ {
		tInvModQ[i] = ring.NewUint(t)
		tInvModQ[i].ModInverse(tInvModQ[i], ringQ.ModulusAtLevel[i])
	}

	return &evaluatorBase{
		params:   params,
		tInvModQ: tInvModQ,
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
		buffCt: NewCiphertext(eval.params, 2, eval.params.MaxLevel()),
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

func (eval *evaluator) evaluateInPlace(level int, el0, el1, elOut *rlwe.Ciphertext, evaluate func(*ring.Poly, *ring.Poly, *ring.Poly)) {

	smallest, largest, _ := rlwe.GetSmallestLargest(el0.El(), el1.El())

	elOut.Resize(elOut.Degree(), level)

	for i := 0; i < smallest.Degree()+1; i++ {
		evaluate(el0.Value[i], el1.Value[i], elOut.Value[i])
	}

	// If the inputs degrees differ, it copies the remaining degree on the receiver.
	if largest != nil && largest != elOut.El() { // checks to avoid unnecessary work.
		for i := smallest.Degree() + 1; i < largest.Degree()+1; i++ {
			elOut.Value[i].Copy(largest.Value[i])
		}
	}

	elOut.MetaData = el0.MetaData
}

func (eval *evaluator) matchScaleThenEvaluateInPlace(level int, el0, el1, elOut *rlwe.Ciphertext, evaluate func(*ring.Poly, uint64, *ring.Poly)) {

	r0, r1, _ := eval.matchScalesBinary(el0.Scale.Uint64(), el1.Scale.Uint64())

	for i := range el0.Value {
		eval.params.RingQ().AtLevel(level).MulScalar(el0.Value[i], r0, elOut.Value[i])
	}

	for i := range el1.Value {
		evaluate(el1.Value[i], r1, elOut.Value[i])
	}

	elOut.MetaData = el0.MetaData
	elOut.Scale = el0.Scale.Mul(eval.params.NewScale(r0))
}

func (eval *evaluator) newCiphertextBinary(op0, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext) {
	return NewCiphertext(eval.params, utils.MaxInt(op0.Degree(), op1.Degree()), utils.MinInt(op0.Level(), op1.Level()))
}

// Add adds op1 to ctIn and returns the result in ctOut.
func (eval *evaluator) Add(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext) {
	_, level := eval.params.CheckBinary(ctIn, op1, ctOut, utils.MaxInt(ctIn.Degree(), op1.Degree()))

	if ctIn.Scale.Cmp(op1.GetScale()) == 0 {
		eval.evaluateInPlace(level, ctIn, op1.El(), ctOut, eval.params.RingQ().AtLevel(level).Add)
	} else {
		eval.matchScaleThenEvaluateInPlace(level, ctIn, op1.El(), ctOut, eval.params.RingQ().AtLevel(level).MulScalarThenAdd)
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
	_, level := eval.params.CheckBinary(ctIn, op1, ctOut, utils.MaxInt(ctIn.Degree(), op1.Degree()))

	if ctIn.Scale.Cmp(op1.GetScale()) == 0 {
		eval.evaluateInPlace(level, ctIn, op1.El(), ctOut, eval.params.RingQ().AtLevel(level).Sub)
	} else {
		eval.matchScaleThenEvaluateInPlace(level, ctIn, op1.El(), ctOut, eval.params.RingQ().AtLevel(level).MulScalarThenSub)
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
		eval.params.RingQ().AtLevel(level).Neg(ctIn.Value[i], ctOut.Value[i])
	}

	ctOut.MetaData = ctIn.MetaData
}

// NegNew negates ctIn and returns the result in a new ctOut.
func (eval *evaluator) NegNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, ctIn.Degree(), ctIn.Level())
	eval.Neg(ctIn, ctOut)
	return
}

// AddScalar adds a scalar to ctIn and returns the result in ctOut.
func (eval *evaluator) AddScalar(ctIn *rlwe.Ciphertext, scalar uint64, ctOut *rlwe.Ciphertext) {

	ringT := eval.params.RingT()

	level := utils.MinInt(ctIn.Level(), ctOut.Level())

	if ctIn.Scale.Cmp(eval.params.NewScale(1)) != 0 {
		scalar = ring.BRed(scalar, ctIn.Scale.Uint64(), ringT.SubRings[0].Modulus, ringT.SubRings[0].BRedConstant)
	}

	scalarBig := new(big.Int).SetUint64(scalar)

	scalarBig.Mul(scalarBig, eval.tInvModQ[level])

	eval.params.RingQ().AtLevel(level).AddScalarBigint(ctIn.Value[0], scalarBig, ctOut.Value[0])

	if ctIn != ctOut {
		for i := 1; i < ctIn.Degree()+1; i++ {
			ring.Copy(ctIn.Value[i], ctOut.Value[i])
		}

		ctOut.MetaData = ctIn.MetaData
	}
}

// AddScalarNew adds a scalar to ctIn and returns the result in a new ctOut.
func (eval *evaluator) AddScalarNew(ctIn *rlwe.Ciphertext, scalar uint64) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, ctIn.Degree(), ctIn.Level())
	eval.AddScalar(ctIn, scalar, ctOut)
	return
}

// MulScalar multiplies ctIn with a scalar and returns the result in ctOut.
func (eval *evaluator) MulScalar(ctIn *rlwe.Ciphertext, scalar uint64, ctOut *rlwe.Ciphertext) {
	ringQ := eval.params.RingQ().AtLevel(utils.MinInt(ctIn.Level(), ctOut.Level()))
	for i := 0; i < ctIn.Degree()+1; i++ {
		ringQ.MulScalar(ctIn.Value[i], scalar, ctOut.Value[i])
	}
	ctOut.MetaData = ctIn.MetaData
}

// MulScalarNew multiplies ctIn with a scalar and returns the result in a new ctOut.
func (eval *evaluator) MulScalarNew(ctIn *rlwe.Ciphertext, scalar uint64) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, ctIn.Degree(), ctIn.Level())
	eval.MulScalar(ctIn, scalar, ctOut)
	return
}

// MulScalarThenAdd multiplies ctIn with a scalar adds the result on ctOut.
func (eval *evaluator) MulScalarThenAdd(ctIn *rlwe.Ciphertext, scalar uint64, ctOut *rlwe.Ciphertext) {
	ringQ := eval.params.RingQ().AtLevel(utils.MinInt(ctIn.Level(), ctOut.Level()))

	// scalar *= (ctOut.scale / ctIn.Scale)
	if ctIn.Scale.Cmp(ctOut.Scale) != 0 {
		ringT := eval.params.RingT()
		ratio := ring.ModExp(ctIn.Scale.Uint64(), ringT.SubRings[0].Modulus-2, ringT.SubRings[0].Modulus)
		ratio = ring.BRed(ratio, ctOut.Scale.Uint64(), ringT.SubRings[0].Modulus, ringT.SubRings[0].BRedConstant)
		scalar = ring.BRed(ratio, scalar, ringT.SubRings[0].Modulus, ringT.SubRings[0].BRedConstant)
	}

	for i := 0; i < ctIn.Degree()+1; i++ {
		ringQ.MulScalarThenAdd(ctIn.Value[i], scalar, ctOut.Value[i])
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
	ctOut = NewCiphertext(eval.params, ctIn.Degree()+op1.Degree(), utils.MinInt(ctIn.Level(), op1.Level()))
	eval.mulRelin(ctIn, op1, false, ctOut)
	return
}

// MulRelinNew multiplies ctIn with op1 with relinearization and returns the result in a new ctOut.
// The procedure will panic if either ctIn.Degree or op1.Degree > 1.
// The procedure will panic if the evaluator was not created with an relinearization key.
func (eval *evaluator) MulRelinNew(ctIn *rlwe.Ciphertext, op1 rlwe.Operand) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1, utils.MinInt(ctIn.Level(), op1.Level()))
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
	ctOut.Scale = ctIn.Scale.Mul(op1.GetScale())

	ringQ := eval.params.RingQ().AtLevel(level)

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

		ringQ.MForm(tmp0.Value[0], c00)
		ringQ.MForm(tmp0.Value[1], c01)

		ringQ.MulScalar(c00, eval.params.T(), c00)
		ringQ.MulScalar(c01, eval.params.T(), c01)

		if ctIn.El() == op1.El() { // squaring case
			ringQ.MulCoeffsMontgomery(c00, tmp1.Value[0], c0) // c0 = c[0]*c[0]
			ringQ.MulCoeffsMontgomery(c01, tmp1.Value[1], c2) // c2 = c[1]*c[1]
			ringQ.MulCoeffsMontgomery(c00, tmp1.Value[1], c1) // c1 = 2*c[0]*c[1]
			ringQ.Add(c1, c1, c1)

		} else { // regular case
			ringQ.MulCoeffsMontgomery(c00, tmp1.Value[0], c0) // c0 = c0[0]*c0[0]
			ringQ.MulCoeffsMontgomery(c01, tmp1.Value[1], c2) // c2 = c0[1]*c1[1]
			ringQ.MulCoeffsMontgomery(c00, tmp1.Value[1], c1)
			ringQ.MulCoeffsMontgomeryThenAdd(c01, tmp1.Value[0], c1) // c1 = c0[0]*c1[1] + c0[1]*c1[0]
		}

		if relin {

			if eval.Rlk == nil {
				panic("cannot MulRelin: relinerization key is missing")
			}

			tmpCt := &rlwe.Ciphertext{Value: []*ring.Poly{eval.BuffQP[1].Q, eval.BuffQP[2].Q}}
			tmpCt.IsNTT = true

			eval.GadgetProduct(level, c2, eval.Rlk.Keys[0].GadgetCiphertext, tmpCt)

			ringQ.Add(ctOut.Value[0], tmpCt.Value[0], ctOut.Value[0])
			ringQ.Add(ctOut.Value[1], tmpCt.Value[1], ctOut.Value[1])
		}

		// Case Plaintext (x) Ciphertext or Ciphertext (x) Plaintext
	} else {

		if ctOut.Degree() < ctIn.Degree() {
			ctOut.Resize(ctIn.Degree(), level)
		}

		c00 := eval.buffQ[0]

		ringQ.MForm(op1.El().Value[0], c00)
		ringQ.MulScalar(c00, eval.params.T(), c00)
		for i := range ctOut.Value {
			ringQ.MulCoeffsMontgomery(ctIn.Value[i], c00, ctOut.Value[i])
		}
	}
}

// MulThenAdd multiplies ctIn with op1 (without relinearization)^and adds the result on ctOut.
// The procedure will panic if either ctIn.Degree() or op1.Degree() > 1.
// The procedure will panic if either ctIn == ctOut or op1 == ctOut.
func (eval *evaluator) MulThenAdd(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext) {
	eval.mulRelinThenAdd(ctIn, op1, false, ctOut)
}

// MulRelinThenAdd multiplies ctIn with op1 and adds, relinearize the result on ctOut.
// The procedure will panic if either ctIn.Degree() or op1.Degree() > 1.
// The procedure will panic if either ctIn == ctOut or op1 == ctOut.
func (eval *evaluator) MulRelinThenAdd(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, ctOut *rlwe.Ciphertext) {
	eval.mulRelinThenAdd(ctIn, op1, true, ctOut)
}

func (eval *evaluator) mulRelinThenAdd(ctIn *rlwe.Ciphertext, op1 rlwe.Operand, relin bool, ctOut *rlwe.Ciphertext) {

	eval.checkBinary(ctIn, op1, ctOut, utils.MaxInt(ctIn.Degree(), op1.Degree()))

	level := utils.MinInt(utils.MinInt(ctIn.Level(), op1.Level()), ctOut.Level())

	if ctIn.Degree()+op1.Degree() > 2 {
		panic("cannot MulRelinThenAdd: input elements total degree cannot be larger than 2")
	}

	if ctIn.El() == ctOut.El() || op1.El() == ctOut.El() {
		panic("cannot MulRelinThenAdd: ctOut must be different from ctIn and op1")
	}

	ringQ := eval.params.RingQ().AtLevel(level)
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
		if targetScale := ring.BRed(ctIn.Scale.Uint64(), op1.GetScale().Uint64(), ringT.SubRings[0].Modulus, ringT.SubRings[0].BRedConstant); ctOut.Scale.Cmp(eval.params.NewScale(targetScale)) != 0 {
			var r1 uint64
			r0, r1, _ = eval.matchScalesBinary(targetScale, ctOut.Scale.Uint64())

			for i := range ctOut.Value {
				ringQ.MulScalar(ctOut.Value[i], r1, ctOut.Value[i])
			}

			ctOut.Scale = ctOut.Scale.Mul(eval.params.NewScale(r1))
		}

		ringQ.MForm(tmp0.Value[0], c00)
		ringQ.MForm(tmp0.Value[1], c01)

		ringQ.MulScalar(c00, eval.params.T(), c00)
		ringQ.MulScalar(c01, eval.params.T(), c01)

		if r0 != 1 {
			ringQ.MulScalar(c00, r0, c00)
			ringQ.MulScalar(c01, r0, c01)
		}

		ringQ.MulCoeffsMontgomeryThenAdd(c00, tmp1.Value[0], c0) // c0 += c[0]*c[0]
		ringQ.MulCoeffsMontgomeryThenAdd(c00, tmp1.Value[1], c1) // c1 += c[0]*c[1]
		ringQ.MulCoeffsMontgomeryThenAdd(c01, tmp1.Value[0], c1) // c1 += c[1]*c[0]

		if relin {

			if eval.Rlk == nil {
				panic("cannot MulRelinThenAdd: relinerization key is missing")
			}

			ringQ.MulCoeffsMontgomery(c01, tmp1.Value[1], c2) // c2 += c[1]*c[1]

			tmpCt := &rlwe.Ciphertext{Value: []*ring.Poly{eval.BuffQP[1].Q, eval.BuffQP[2].Q}}
			tmpCt.IsNTT = true

			eval.GadgetProduct(level, c2, eval.Rlk.Keys[0].GadgetCiphertext, tmpCt)

			ringQ.Add(ctOut.Value[0], tmpCt.Value[0], ctOut.Value[0])
			ringQ.Add(ctOut.Value[1], tmpCt.Value[1], ctOut.Value[1])

		} else {
			ringQ.MulCoeffsMontgomeryThenAdd(c01, tmp1.Value[1], c2) // c2 += c[1]*c[1]
		}

		// Case Plaintext (x) Ciphertext or Ciphertext (x) Plaintext
	} else {

		if ctOut.Degree() < ctIn.Degree() {
			ctOut.Resize(ctIn.Degree(), level)
		}

		c00 := eval.buffQ[0]

		ringQ.MForm(op1.El().Value[0], c00)
		ringQ.MulScalar(c00, eval.params.T(), c00)

		var r0 = uint64(1)
		if targetScale := ring.BRed(ctIn.Scale.Uint64(), op1.GetScale().Uint64(), ringT.SubRings[0].Modulus, ringT.SubRings[0].BRedConstant); ctOut.Scale.Cmp(eval.params.NewScale(targetScale)) != 0 {
			var r1 uint64
			r0, r1, _ = eval.matchScalesBinary(targetScale, ctOut.Scale.Uint64())

			for i := range ctOut.Value {
				ringQ.MulScalar(ctOut.Value[i], r1, ctOut.Value[i])
			}

			ctOut.Scale = ctOut.Scale.Mul(eval.params.NewScale(r1))
		}

		if r0 != 1 {
			ringQ.MulScalar(c00, r0, c00)
		}

		for i := range ctIn.Value {
			ringQ.MulCoeffsMontgomeryThenAdd(ctIn.Value[i], c00, ctOut.Value[i])
		}
	}
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
	ringQ := eval.params.RingQ().AtLevel(level)

	for i := range ctOut.Value {
		ringQ.DivRoundByLastModulusNTT(ctIn.Value[i], eval.buffQ[0], ctOut.Value[i])
	}

	ctOut.Resize(ctOut.Degree(), level-1)
	ctOut.MetaData = ctIn.MetaData
	ctOut.Scale = ctIn.Scale.Div(eval.params.NewScale(ringQ.SubRings[level].Modulus))
	return
}

// RelinearizeNew applies the relinearization procedure on ctIn and returns the result in a new ctOut.
func (eval *evaluator) RelinearizeNew(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1, ctIn.Level())
	eval.Relinearize(ctIn, ctOut)
	return
}

// SwitchKeysNew re-encrypts ctIn under a different key and returns the result in a new ctOut.
// It requires a SwitchingKey, which is computed from the key under which the Ciphertext is currently encrypted,
// and the key under which the Ciphertext will be re-encrypted.
// The procedure will panic if either ctIn.Degree() or ctOut.Degree() != 1.
func (eval *evaluator) SwitchKeysNew(ctIn *rlwe.Ciphertext, swk *rlwe.SwitchingKey) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, ctIn.Degree(), ctIn.Level())
	eval.SwitchKeys(ctIn, swk, ctOut)
	return
}

// RotateColumnsNew rotates the columns of ctIn by k positions to the left, and returns the result in a newly created element.
// The procedure will panic if the corresponding Galois key has not been generated and attributed to the evaluator.
// The procedure will panic if ctIn.Degree() != 1.
func (eval *evaluator) RotateColumnsNew(ctIn *rlwe.Ciphertext, k int) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, ctIn.Degree(), ctIn.Level())
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
	ctOut = NewCiphertext(eval.params, ctIn.Degree(), ctIn.Level())
	eval.RotateRows(ctIn, ctOut)
	return
}

// RotateRows swaps the rows of ctIn and returns the result in ctOut.
// The procedure will panic if the corresponding Galois key has not been generated and attributed to the evaluator.
// The procedure will panic if either ctIn.Degree() or ctOut.Degree() != 1.
func (eval *evaluator) RotateRows(ctIn *rlwe.Ciphertext, ctOut *rlwe.Ciphertext) {
	eval.Automorphism(ctIn, eval.params.GaloisElementForRowRotation(), ctOut)
}

func (eval *evaluator) RotateHoistedLazyNew(level int, rotations []int, c0 *ring.Poly, c2DecompQP []ringqp.Poly) (cOut map[int]rlwe.CiphertextQP) {
	cOut = make(map[int]rlwe.CiphertextQP)
	for _, i := range rotations {
		if i != 0 {
			cOut[i] = rlwe.NewCiphertextQP(eval.params.Parameters, level, eval.params.PCount()-1)
			eval.AutomorphismHoistedLazy(level, c0, c2DecompQP, eval.params.GaloisElementForColumnRotationBy(i), cOut[i])
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

	ringQ := eval.params.RingQ().AtLevel(level)

	for _, el := range ct0.Value {
		ringQ.MulScalar(el, r0, el)
	}

	ct0.Resize(ct0.Degree(), level)
	ct0.Scale = ct0.Scale.Mul(eval.params.NewScale(r0))

	for _, el := range ct1.Value {
		ringQ.MulScalar(el, r1, el)
	}

	ct1.Resize(ct1.Degree(), level)
	ct1.Scale = ct1.Scale.Mul(eval.params.NewScale(r1))
}

func (eval *evaluator) matchScalesBinary(scale0, scale1 uint64) (r0, r1, e uint64) {

	ringT := eval.params.RingT()

	t := ringT.SubRings[0].Modulus
	tHalf := t >> 1
	BRedConstant := ringT.SubRings[0].BRedConstant

	if utils.GCD(scale0, t) != 1 {
		panic("cannot matchScalesBinary: invalid ciphertext scale: gcd(scale, t) != 1")
	}

	var a = ringT.SubRings[0].Modulus
	var b uint64 = 0
	var A = ring.BRed(ring.ModExp(scale0, t-2, t), scale1, t, BRedConstant)
	var B uint64 = 1

	e = center(A, tHalf, t) + 1

	for A != 0 {
		q := a / A
		a, A = A, a%A
		b, B = B, ring.CRed(t+b-ring.BRed(B, q, t, BRedConstant), t)

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
