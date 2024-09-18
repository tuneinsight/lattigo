package bgv

import (
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
	"github.com/tuneinsight/lattigo/v6/utils"
)

// Evaluator is a struct that holds the necessary elements to perform the homomorphic operations between ciphertexts and/or plaintexts.
// It also holds a memory buffer used to store intermediate computations.
// The [Evaluator.ScaleInvariant] flag needs to be set in order to use a BFV-style
// version of the evaluator.
type Evaluator struct {
	*evaluatorBase
	*evaluatorBuffers
	*rlwe.Evaluator
	*Encoder

	// ScaleInvariant is a flag indicating whether the evaluator executes
	// scale-invariant multiplications (transforming the BGV evaluator into
	// BFV evaluator).
	ScaleInvariant bool
}

type evaluatorBase struct {
	tMontgomery         ring.RNSScalar
	levelQMul           []int      // optimal #QiMul depending on #Qi (variable level)
	pHalf               []*big.Int // all prod(QiMul) / 2 depending on #Qi
	basisExtenderQ1toQ2 *ring.BasisExtender
}

func (eval evaluatorBase) ShallowCopy() *evaluatorBase {
	return &evaluatorBase{
		tMontgomery:         eval.tMontgomery,
		levelQMul:           eval.levelQMul,
		pHalf:               eval.pHalf,
		basisExtenderQ1toQ2: eval.basisExtenderQ1toQ2.ShallowCopy(),
	}
}

func newEvaluatorPrecomp(parameters Parameters) *evaluatorBase {
	ringQ := parameters.RingQ()
	ringQMul := parameters.RingQMul()
	t := parameters.PlaintextModulus()

	levelQMul := make([]int, ringQ.ModuliChainLength())
	Q := new(big.Int).SetUint64(1)
	for i := range levelQMul {
		Q.Mul(Q, new(big.Int).SetUint64(ringQ.SubRings[i].Modulus))
		levelQMul[i] = int(math.Ceil(float64(Q.BitLen()+parameters.LogN())/61.0)) - 1
	}

	pHalf := make([]*big.Int, ringQMul.ModuliChainLength())

	QMul := new(big.Int).SetUint64(1)
	for i := range pHalf {
		QMul.Mul(QMul, new(big.Int).SetUint64(ringQMul.SubRings[i].Modulus))
		pHalf[i] = new(big.Int).Rsh(QMul, 1)
	}

	basisExtenderQ1toQ2 := ring.NewBasisExtender(ringQ, ringQMul)

	// PlaintextModulus * 2^{64} mod Q
	tMontgomery := ringQ.NewRNSScalarFromBigint(new(big.Int).Lsh(new(big.Int).SetUint64(t), 64))
	ringQ.MFormRNSScalar(tMontgomery, tMontgomery)

	return &evaluatorBase{
		tMontgomery:         tMontgomery,
		levelQMul:           levelQMul,
		pHalf:               pHalf,
		basisExtenderQ1toQ2: basisExtenderQ1toQ2,
	}
}

type evaluatorBuffers struct {
	buffQ    [3]ring.Poly
	buffQMul [9]ring.Poly
}

// BuffQ returns a pointer to the internal memory buffer buffQ.
func (eval Evaluator) BuffQ() [3]ring.Poly {
	return eval.buffQ
}

func newEvaluatorBuffer(params Parameters) *evaluatorBuffers {

	ringQ := params.RingQ()
	buffQ := [3]ring.Poly{
		ringQ.NewPoly(),
		ringQ.NewPoly(),
		ringQ.NewPoly(),
	}

	ringQMul := params.RingQMul()

	buffQMul := [9]ring.Poly{
		ringQMul.NewPoly(),
		ringQMul.NewPoly(),
		ringQMul.NewPoly(),
		ringQMul.NewPoly(),
		ringQMul.NewPoly(),
		ringQMul.NewPoly(),
		ringQMul.NewPoly(),
		ringQMul.NewPoly(),
		ringQMul.NewPoly(),
	}

	return &evaluatorBuffers{
		buffQ:    buffQ,
		buffQMul: buffQMul,
	}
}

// NewEvaluator creates a new [Evaluator], that can be used to do homomorphic
// operations on ciphertexts and/or plaintexts. It stores a memory buffer
// and ciphertexts that will be used for intermediate values.
// The evaluator can optionally be initialized as scale-invariant, which
// transforms it into a BFV evaluator. See `schemes/bfv/README.md` for more
// information.
func NewEvaluator(parameters Parameters, evk rlwe.EvaluationKeySet, scaleInvariant ...bool) *Evaluator {
	ev := new(Evaluator)
	ev.evaluatorBase = newEvaluatorPrecomp(parameters)
	ev.evaluatorBuffers = newEvaluatorBuffer(parameters)
	ev.Evaluator = rlwe.NewEvaluator(parameters.Parameters, evk)
	ev.Encoder = NewEncoder(parameters)
	ev.ScaleInvariant = len(scaleInvariant) > 0 && scaleInvariant[0]

	return ev
}

// GetParameters returns a pointer to the underlying [bgv.Parameters].
func (eval Evaluator) GetParameters() *Parameters {
	return &eval.Encoder.parameters
}

// ShallowCopy creates a shallow copy of this [Evaluator] in which the read-only data-structures are
// shared with the receiver.
func (eval Evaluator) ShallowCopy() *Evaluator {
	return &Evaluator{
		evaluatorBase:    eval.evaluatorBase.ShallowCopy(),
		Evaluator:        eval.Evaluator.ShallowCopy(),
		evaluatorBuffers: newEvaluatorBuffer(*eval.GetParameters()),
		Encoder:          eval.Encoder.ShallowCopy(),
		ScaleInvariant:   eval.ScaleInvariant,
	}
}

// WithKey creates a shallow copy of this [Evaluator] in which the read-only data-structures are
// shared with the receiver but the evaluation key is set to the provided [rlwe.EvaluationKeySet].
func (eval Evaluator) WithKey(evk rlwe.EvaluationKeySet) *Evaluator {
	return &Evaluator{
		evaluatorBase:    eval.evaluatorBase,
		Evaluator:        eval.Evaluator.WithKey(evk),
		evaluatorBuffers: eval.evaluatorBuffers,
		Encoder:          eval.Encoder,
	}
}

// Add adds op1 to op0 and returns the result in opOut.
// inputs:
//   - op0: an *[rlwe.Ciphertext]
//   - op1:
//   - [rlwe.ElementInterface][[ring.Poly]]
//   - *big.Int, uint64, int64, int
//   - []uint64 or []int64 (of size at most N where N is the smallest integer satisfying PlaintextModulus = 1 mod 2N)
//   - opOut: an *[rlwe.Ciphertext]
//
// If op1 is an [rlwe.ElementInterface][[ring.Poly]] and the scales of op0, op1 and opOut do not match, then a scale matching operation will
// be automatically carried out to ensure that addition is performed between operands of the same scale.
// This scale matching operation will increase the noise by a small factor.
// For this reason it is preferable to ensure that all operands are already at the same scale when calling this method.
func (eval Evaluator) Add(op0 *rlwe.Ciphertext, op1 rlwe.Operand, opOut *rlwe.Ciphertext) (err error) {

	ringQ := eval.parameters.RingQ()

	switch op1 := op1.(type) {
	case rlwe.ElementInterface[ring.Poly]:

		degree, level, err := eval.InitOutputBinaryOp(op0.El(), op1.El(), op0.Degree()+op1.Degree(), opOut.El())
		if err != nil {
			return fmt.Errorf("cannot Add: %w", err)
		}

		opOut.Resize(degree, level)

		if op0.Scale.Cmp(op1.El().Scale) == 0 {
			eval.evaluateInPlace(level, op0, op1.El(), opOut, ringQ.AtLevel(level).Add)
		} else {
			eval.matchScaleThenEvaluateInPlace(level, op0, op1.El(), opOut, ringQ.AtLevel(level).MulScalarThenAdd)
		}

	case *big.Int:

		_, level, err := eval.InitOutputUnaryOp(op0.El(), opOut.El())
		if err != nil {
			return fmt.Errorf("cannot Add: %w", err)
		}

		opOut.Resize(op0.Degree(), level)

		TBig := eval.parameters.RingT().ModulusAtLevel[0]

		// Sets op1 to the scale of op0
		op1.Mul(op1, new(big.Int).SetUint64(op0.Scale.Uint64()))

		op1.Mod(op1, TBig)

		// If op1 > T/2 -> op1 -= T
		if op1.Cmp(new(big.Int).Rsh(TBig, 1)) == 1 {
			op1.Sub(op1, TBig)
		}

		// Scales op0 by T^{-1} mod Q
		op1.Mul(op1, eval.tInvModQ[level])

		ringQ.AtLevel(level).AddScalarBigint(op0.Value[0], op1, opOut.Value[0])

		if op0 != opOut {
			for i := 1; i < op0.Degree()+1; i++ {
				opOut.Value[i].CopyLvl(level, op0.Value[i])
			}
		}

	case uint64:
		return eval.Add(op0, new(big.Int).SetUint64(op1), opOut)
	case int64:
		return eval.Add(op0, new(big.Int).SetInt64(op1), opOut)
	case int:
		return eval.Add(op0, new(big.Int).SetInt64(int64(op1)), opOut)
	case []uint64, []int64:

		_, level, err := eval.InitOutputUnaryOp(op0.El(), opOut.El())
		if err != nil {
			return fmt.Errorf("cannot Add: %w", err)
		}

		opOut.Resize(op0.Degree(), level)

		// Instantiates new plaintext from buffer
		pt, err := rlwe.NewPlaintextAtLevelFromPoly(level, eval.buffQ[0])

		// This error should not happen, unless the evaluator's buffer were
		// improperly tempered with. If it does happen, there is no way to
		// recover from it.
		if err != nil {
			panic(err)
		}

		pt.MetaData = op0.MetaData // Sets the metadata, notably matches scalses

		// Encodes the vector on the plaintext
		if err = eval.Encoder.Encode(op1, pt); err != nil {
			return err
		}

		// Generic in place evaluation
		eval.evaluateInPlace(level, op0, pt.El(), opOut, eval.parameters.RingQ().AtLevel(level).Add)
	default:
		return fmt.Errorf("invalid op1.(Type), expected rlwe.ElementInterface[ring.Poly], []uint64, []int64, *big.Int, uint64, int64 or int, but got %T", op1)
	}

	return
}

func (eval Evaluator) evaluateInPlace(level int, el0 *rlwe.Ciphertext, el1 *rlwe.Element[ring.Poly], elOut *rlwe.Ciphertext, evaluate func(ring.Poly, ring.Poly, ring.Poly)) {

	smallest, largest, _ := rlwe.GetSmallestLargest(el0.El(), el1.El())

	for i := 0; i < smallest.Degree()+1; i++ {
		evaluate(el0.Value[i], el1.Value[i], elOut.Value[i])
	}

	elOut.Scale = el0.Scale.Max(el1.Scale)

	// If the inputs degrees differ, it copies the remaining degree on the receiver.
	if largest != nil && largest != elOut.El() { // checks to avoid unnecessary work.
		for i := smallest.Degree() + 1; i < largest.Degree()+1; i++ {
			elOut.Value[i].CopyLvl(level, largest.Value[i])
		}
	}
}

func (eval Evaluator) matchScaleThenEvaluateInPlace(level int, el0 *rlwe.Ciphertext, el1 *rlwe.Element[ring.Poly], elOut *rlwe.Ciphertext, evaluate func(ring.Poly, uint64, ring.Poly)) {

	r0, r1, _ := eval.matchScalesBinary(el0.Scale.Uint64(), el1.Scale.Uint64())

	for i := range el0.Value {
		eval.parameters.RingQ().AtLevel(level).MulScalar(el0.Value[i], r0, elOut.Value[i])
	}

	for i := el0.Degree() + 1; i < elOut.Degree()+1; i++ {
		elOut.Value[i].Zero()
	}

	for i := range el1.Value {
		evaluate(el1.Value[i], r1, elOut.Value[i])
	}

	elOut.Scale = el0.Scale.Mul(eval.parameters.NewScale(r0))
}

func (eval Evaluator) newCiphertextBinary(op0, op1 rlwe.ElementInterface[ring.Poly]) (opOut *rlwe.Ciphertext) {
	return NewCiphertext(*eval.GetParameters(), utils.Max(op0.Degree(), op1.Degree()), utils.Min(op0.Level(), op1.Level()))
}

// AddNew adds op1 to op0 and returns the result on a new *[rlwe.Ciphertext] opOut.
// inputs:
//   - op0: an *[rlwe.Ciphertext]
//   - op1:
//   - [rlwe.ElementInterface][[ring.Poly]]
//   - *big.Int, uint64, int64, int
//   - []uint64 or []int64 (of size at most N where N is the smallest integer satisfying PlaintextModulus = 1 mod 2N)
//
// If op1 is an [rlwe.ElementInterface][[ring.Poly]] and the scales of op0 and op1 not match, then a scale matching operation will
// be automatically carried out to ensure that addition is performed between operands of the same scale.
// This scale matching operation will increase the noise by a small factor.
// For this reason it is preferable to ensure that all operands are already at the same scale when calling this method.
func (eval Evaluator) AddNew(op0 *rlwe.Ciphertext, op1 rlwe.Operand) (opOut *rlwe.Ciphertext, err error) {

	switch op1 := op1.(type) {
	case rlwe.ElementInterface[ring.Poly]:
		opOut = eval.newCiphertextBinary(op0, op1)
	default:
		opOut = NewCiphertext(eval.parameters, op0.Degree(), op0.Level())
	}

	return opOut, eval.Add(op0, op1, opOut)
}

// Sub subtracts op1 to op0 and returns the result in opOut.
// inputs:
//   - op0: an *[rlwe.Ciphertext]
//   - op1:
//   - [rlwe.ElementInterface][[ring.Poly]]
//   - *big.Int, uint64, int64, int
//   - []uint64 or []int64 (of size at most N where N is the smallest integer satisfying PlaintextModulus = 1 mod 2N)
//   - opOut: an *[rlwe.Ciphertext]
//
// If op1 is an [rlwe.ElementInterface][[ring.Poly]] and the scales of op0, op1 and opOut do not match, then a scale matching operation will
// be automatically carried out to ensure that the subtraction is performed between operands of the same scale.
// This scale matching operation will increase the noise by a small factor.
// For this reason it is preferable to ensure that all operands are already at the same scale when calling this method.
func (eval Evaluator) Sub(op0 *rlwe.Ciphertext, op1 rlwe.Operand, opOut *rlwe.Ciphertext) (err error) {

	switch op1 := op1.(type) {
	case rlwe.ElementInterface[ring.Poly]:

		degree, level, err := eval.InitOutputBinaryOp(op0.El(), op1.El(), op0.Degree()+op1.Degree(), opOut.El())
		if err != nil {
			return fmt.Errorf("cannot Sub: %w", err)
		}

		opOut.Resize(degree, level)

		ringQ := eval.parameters.RingQ()

		if op0.Scale.Cmp(op1.El().Scale) == 0 {
			eval.evaluateInPlace(level, op0, op1.El(), opOut, ringQ.AtLevel(level).Sub)
		} else {
			eval.matchScaleThenEvaluateInPlace(level, op0, op1.El(), opOut, ringQ.AtLevel(level).MulScalarThenSub)
		}

	case *big.Int:
		return eval.Add(op0, new(big.Int).Neg(op1), opOut)
	case uint64:
		return eval.Sub(op0, new(big.Int).SetUint64(op1), opOut)
	case int64:
		return eval.Sub(op0, new(big.Int).SetInt64(op1), opOut)
	case int:
		return eval.Sub(op0, new(big.Int).SetInt64(int64(op1)), opOut)
	case []uint64, []int64:

		_, level, err := eval.InitOutputUnaryOp(op0.El(), opOut.El())
		if err != nil {
			return fmt.Errorf("cannot Sub: %w", err)
		}

		opOut.Resize(op0.Degree(), level)

		// Instantiates new plaintext from buffer
		pt, err := rlwe.NewPlaintextAtLevelFromPoly(level, eval.buffQ[0])

		// This error should not happen, unless the evaluator's buffer were
		// improperly tempered with. If it does happen, there is no way to
		// recover from it.
		if err != nil {
			panic(err)
		}

		pt.MetaData = op0.MetaData // Sets the metadata, notably matches scales

		// Encodes the vector on the plaintext
		if err := eval.Encoder.Encode(op1, pt); err != nil {
			return fmt.Errorf("cannot Sub: %w", err)
		}

		// Generic in place evaluation
		eval.evaluateInPlace(level, op0, pt.El(), opOut, eval.parameters.RingQ().AtLevel(level).Sub)
	default:
		return fmt.Errorf("invalid op1.(Type), expected rlwe.ElementInterface[ring.Poly], []uint64, []int64, *big.Int, uint64, int64 or int, but got %T", op1)
	}

	return
}

// SubNew subtracts op1 to op0 and returns the result in a new *[rlwe.Ciphertext] opOut.
// inputs:
//   - op0: an *[rlwe.Ciphertext]
//   - op1:
//   - [rlwe.ElementInterface][[ring.Poly]]
//   - *big.Int, uint64, int64, int
//   - []uint64 or []int64 (of size at most N where N is the smallest integer satisfying PlaintextModulus = 1 mod 2N)
//
// If op1 is an [rlwe.ElementInterface][[ring.Poly]] and the scales of op0, op1 and opOut do not match, then a scale matching operation will
// be automatically carried out to ensure that the subtraction is performed between operands of the same scale.
// This scale matching operation will increase the noise by a small factor.
// For this reason it is preferable to ensure that all operands are already at the same scale when calling this method.
func (eval Evaluator) SubNew(op0 *rlwe.Ciphertext, op1 rlwe.Operand) (opOut *rlwe.Ciphertext, err error) {
	switch op1 := op1.(type) {
	case rlwe.ElementInterface[ring.Poly]:
		opOut = eval.newCiphertextBinary(op0, op1)
	default:
		opOut = NewCiphertext(eval.parameters, op0.Degree(), op0.Level())
	}

	return opOut, eval.Sub(op0, op1, opOut)
}

// DropLevel reduces the level of op0 by levels.
// No rescaling is applied during this procedure.
func (eval Evaluator) DropLevel(op0 *rlwe.Ciphertext, levels int) {
	op0.Resize(op0.Degree(), op0.Level()-levels)
}

// Mul multiplies op0 with op1 without relinearization using either standard tensoring (BGV/CKKS-style) when [Evaluator.ScaleInvariant]
// is set to false or scale-invariant tensoring (BFV-style) otherwise, i.e., [Evaluator.MulScaleInvariant], and returns the result in opOut.
// This tensoring increases the noise by a multiplicative factor of the plaintext and noise norms of the operands and will usually
// require to be followed by a rescaling operation to avoid an exponential growth of the noise from subsequent multiplications.
// The procedure will return an error if either op0 or op1 are have a degree higher than 1.
// The procedure will return an error if opOut.Degree != op0.Degree + op1.Degree.
//
// inputs:
//   - op0: an *[rlwe.Ciphertext]
//   - op1:
//   - [rlwe.ElementInterface][[ring.Poly]]
//   - *big.Int, uint64, int64, int
//   - []uint64 or []int64 (of size at most N where N is the smallest integer satisfying PlaintextModulus = 1 mod 2N)
//   - opOut: an *[rlwe.Ciphertext]
//
// If op1 is an [rlwe.ElementInterface][[ring.Poly]]:
//   - the level of opOut will be updated to min(op0.Level(), op1.Level())
//   - the scale of opOut will be updated to op0.Scale * op1.Scale
func (eval Evaluator) Mul(op0 *rlwe.Ciphertext, op1 rlwe.Operand, opOut *rlwe.Ciphertext) (err error) {

	if eval.ScaleInvariant {
		switch op1 := op1.(type) {
		case rlwe.ElementInterface[ring.Poly], []uint64, []int64:
			return eval.MulScaleInvariant(op0, op1, opOut)
		}
	}

	switch op1 := op1.(type) {
	case rlwe.ElementInterface[ring.Poly]:

		_, level, err := eval.InitOutputBinaryOp(op0.El(), op1.El(), 2, opOut.El())
		if err != nil {
			return fmt.Errorf("cannot Mul: %w", err)
		}

		opOut.Resize(opOut.Degree(), level)

		if err = eval.tensorStandard(op0, op1.El(), false, opOut); err != nil {
			return fmt.Errorf("cannot Mul: %w", err)
		}

	case *big.Int:

		_, level, err := eval.InitOutputUnaryOp(op0.El(), opOut.El())
		if err != nil {
			return fmt.Errorf("cannot Mul: %w", err)
		}

		opOut.Resize(op0.Degree(), level)

		ringQ := eval.parameters.RingQ().AtLevel(level)

		TBig := eval.parameters.RingT().ModulusAtLevel[0]

		op1.Mod(op1, TBig)

		// If op1 > T/2 then subtract T to minimize the noise
		if op1.Cmp(new(big.Int).Rsh(TBig, 1)) == 1 {
			op1.Sub(op1, TBig)
		}

		for i := 0; i < op0.Degree()+1; i++ {
			ringQ.MulScalarBigint(op0.Value[i], op1, opOut.Value[i])
		}

	case uint64:
		return eval.Mul(op0, new(big.Int).SetUint64(op1), opOut)
	case int:
		return eval.Mul(op0, new(big.Int).SetInt64(int64(op1)), opOut)
	case int64:
		return eval.Mul(op0, new(big.Int).SetInt64(op1), opOut)
	case []uint64, []int64:

		_, level, err := eval.InitOutputUnaryOp(op0.El(), opOut.El())
		if err != nil {
			return fmt.Errorf("cannot Mul: %w", err)
		}

		// Instantiates new plaintext from buffer
		pt, err := rlwe.NewPlaintextAtLevelFromPoly(level, eval.buffQ[0])

		// This error should not happen, unless the evaluator's buffer were
		// improperly tempered with. If it does happen, there is no way to
		// recover from it.
		if err != nil {
			panic(err)
		}

		pt.MetaData = op0.MetaData.CopyNew() // Sets the metadata, notably matches scales
		pt.Scale = rlwe.NewScale(1)

		// Encodes the vector on the plaintext
		if err := eval.Encoder.Encode(op1, pt); err != nil {
			return err
		}

		if err = eval.Mul(op0, pt, opOut); err != nil {
			return fmt.Errorf("cannot Mul: %w", err)
		}
	default:
		return fmt.Errorf("invalid op1.(Type), expected rlwe.ElementInterface[ring.Poly], []uint64, []int64, *big.Int, uint64, int64 or int, but got %T", op1)
	}

	return
}

// MulNew multiplies op0 with op1 without relinearization using standard tensoring (BGV/CKKS-style) when [Evaluator.ScaleInvariant]
// is set to false or scale-invariant tensoring (BFV-style) otherwise, i.e., [Evaluator.MulScaleInvariantNew], and returns
// the result in a new *[rlwe.Ciphertext] opOut. This tensoring increases the noise by a multiplicative factor of the plaintext
// and noise norms of the operands and will usually require to be followed by a rescaling operation to avoid an exponential
// growth of the noise from subsequent multiplications. The procedure will return an error if either op0 or op1 are have a
// degree higher than 1.
//
// inputs:
//   - op0: an *[rlwe.Ciphertext]
//   - op1:
//   - [rlwe.ElementInterface][[ring.Poly]]
//   - *big.Int, uint64, int64, int
//   - []uint64 or []int64 (of size at most N where N is the smallest integer satisfying PlaintextModulus = 1 mod 2N)
//
// If op1 is an [rlwe.ElementInterface][[ring.Poly]]:
//   - the degree of opOut will be op0.Degree() + op1.Degree()
//   - the level of opOut will be to min(op0.Level(), op1.Level())
//   - the scale of opOut will be to op0.Scale * op1.Scale
func (eval Evaluator) MulNew(op0 *rlwe.Ciphertext, op1 rlwe.Operand) (opOut *rlwe.Ciphertext, err error) {

	if eval.ScaleInvariant {
		switch op1 := op1.(type) {
		case rlwe.ElementInterface[ring.Poly], []uint64, []int64:
			return eval.MulScaleInvariantNew(op0, op1)
		}
	}

	switch op1 := op1.(type) {
	case rlwe.ElementInterface[ring.Poly]:
		opOut = NewCiphertext(eval.parameters, op0.Degree()+op1.Degree(), utils.Min(op0.Level(), op1.Level()))
	default:
		opOut = NewCiphertext(eval.parameters, op0.Degree(), op0.Level())
	}

	return opOut, eval.Mul(op0, op1, opOut)
}

// MulRelin multiplies op0 with op1 with relinearization using standard tensoring (BGV/CKKS-style) when [Evaluator.ScaleInvariant]
// is set to false or scale-invariant tensoring (BFV-style) otherwise, i.e., [Evaluator.MulRelinScaleInvariant], and returns the result in
// opOut. This tensoring increases the noise by a multiplicative factor of the plaintext and noise norms of the operands and will usually
// require to be followed by a rescaling operation to avoid an exponential growth of the noise from subsequent multiplications.
// The procedure will return an error if either op0.Degree or op1.Degree > 1.
// The procedure will return an error if opOut.Degree != op0.Degree + op1.Degree.
// The procedure will return an error if the evaluator was not created with an relinearization key.
//
// inputs:
//   - op0: an *[rlwe.Ciphertext]
//   - op1:
//   - [rlwe.ElementInterface][[ring.Poly]]
//   - *big.Int, uint64, int64, int
//   - []uint64 or []int64 (of size at most N where N is the smallest integer satisfying PlaintextModulus = 1 mod 2N)
//   - opOut: an *[rlwe.Ciphertext]
//
// If op1 is an [rlwe.ElementInterface][[ring.Poly]]:
//   - the level of opOut will be updated to min(op0.Level(), op1.Level())
//   - the scale of opOut will be updated to op0.Scale * op1.Scale
func (eval Evaluator) MulRelin(op0 *rlwe.Ciphertext, op1 rlwe.Operand, opOut *rlwe.Ciphertext) (err error) {

	if eval.ScaleInvariant {
		return eval.MulRelinScaleInvariant(op0, op1, opOut)
	}

	switch op1 := op1.(type) {
	case rlwe.ElementInterface[ring.Poly]:

		_, level, err := eval.InitOutputBinaryOp(op0.El(), op1.El(), 2, opOut.El())
		if err != nil {
			return fmt.Errorf("cannot MulRelin: %w", err)
		}

		opOut.Resize(opOut.Degree(), level)

		if err = eval.tensorStandard(op0, op1.El(), true, opOut); err != nil {
			return fmt.Errorf("cannot MulRelin: %w", err)
		}

	default:
		if err = eval.Mul(op0, op1, opOut); err != nil {
			return fmt.Errorf("cannot MulRelin: %w", err)
		}
	}

	return
}

// MulRelinNew multiplies op0 with op1 with relinearization using standard tensoring (BGV/CKKS-style) when [Evaluator.ScaleInvariant]
// is set to false or scale-invariant tensoring (BFV-style) otherwise, i.e., [Evaluator.MulRelinScaleInvariantNew], and returns the result
// in a new *[rlwe.Ciphertext] opOut. This tensoring increases the noise by a multiplicative factor of the plaintext and noise norms
// of the operands and will usually require to be followed by a rescaling operation to avoid an exponential growth of the noise from
// subsequent multiplications.
// The procedure will return an error if either op0.Degree or op1.Degree > 1.
// The procedure will return an error if the evaluator was not created with an relinearization key.
//
// inputs:
//   - op0: an *[rlwe.Ciphertext]
//   - op1:
//   - [rlwe.ElementInterface][[ring.Poly]]
//   - *big.Int, uint64, int64, int
//   - []uint64 or []int64 (of size at most N where N is the smallest integer satisfying PlaintextModulus = 1 mod 2N)
//
// If op1 is an [rlwe.ElementInterface][[ring.Poly]]:
//   - the level of opOut will be to min(op0.Level(), op1.Level())
//   - the scale of opOut will be to op0.Scale * op1.Scale
func (eval Evaluator) MulRelinNew(op0 *rlwe.Ciphertext, op1 rlwe.Operand) (opOut *rlwe.Ciphertext, err error) {

	if eval.ScaleInvariant {
		return eval.MulRelinScaleInvariantNew(op0, op1)
	}

	switch op1 := op1.(type) {
	case rlwe.ElementInterface[ring.Poly]:
		opOut = NewCiphertext(eval.parameters, 1, utils.Min(op0.Level(), op1.Level()))
	default:
		opOut = NewCiphertext(eval.parameters, 1, op0.Level())
	}

	return opOut, eval.MulRelin(op0, op1, opOut)
}

func (eval Evaluator) tensorStandard(op0 *rlwe.Ciphertext, op1 *rlwe.Element[ring.Poly], relin bool, opOut *rlwe.Ciphertext) (err error) {

	level := opOut.Level()

	opOut.Scale = op0.Scale.Mul(op1.Scale)

	ringQ := eval.parameters.RingQ().AtLevel(level)

	var c00, c01, c0, c1, c2 ring.Poly

	// Case Ciphertext (x) Ciphertext
	if op0.Degree() == 1 && op1.Degree() == 1 {

		c00 = eval.buffQ[0]
		c01 = eval.buffQ[1]

		c0 = opOut.Value[0]
		c1 = opOut.Value[1]

		if !relin {
			opOut.Resize(2, opOut.Level())
			c2 = opOut.Value[2]
		} else {
			opOut.Resize(1, opOut.Level())
			c2 = eval.buffQ[2]
		}

		// Avoid overwriting if the second input is the output
		var tmp0, tmp1 *rlwe.Element[ring.Poly]
		if op1.El() == opOut.El() {
			tmp0, tmp1 = op1.El(), op0.El()
		} else {
			tmp0, tmp1 = op0.El(), op1.El()
		}

		// Multiply by T * 2^{64} * 2^{64} -> result multipled by T and switched in the Montgomery domain
		ringQ.MulRNSScalarMontgomery(tmp0.Value[0], eval.tMontgomery, c00)
		ringQ.MulRNSScalarMontgomery(tmp0.Value[1], eval.tMontgomery, c01)

		if op0.El() == op1.El() { // squaring case
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

			var rlk *rlwe.RelinearizationKey
			var err error
			if rlk, err = eval.CheckAndGetRelinearizationKey(); err != nil {
				return fmt.Errorf("cannot Tensor: cannot Relinearize: %w", err)
			}

			tmpCt := &rlwe.Ciphertext{}
			tmpCt.Value = []ring.Poly{eval.BuffQP[1].Q, eval.BuffQP[2].Q}
			tmpCt.MetaData = &rlwe.MetaData{}
			tmpCt.IsNTT = true

			eval.GadgetProduct(level, c2, &rlk.GadgetCiphertext, tmpCt)

			ringQ.Add(opOut.Value[0], tmpCt.Value[0], opOut.Value[0])
			ringQ.Add(opOut.Value[1], tmpCt.Value[1], opOut.Value[1])
		}

		// Case Plaintext (x) Ciphertext or Ciphertext (x) Plaintext
	} else {

		opOut.Resize(op0.Degree(), level)

		c00 := eval.buffQ[0]

		// Multiply by T * 2^{64} * 2^{64} -> result multipled by T and switched in the Montgomery domain
		ringQ.MulRNSScalarMontgomery(op1.El().Value[0], eval.tMontgomery, c00)
		for i := range opOut.Value {
			ringQ.MulCoeffsMontgomery(op0.Value[i], c00, opOut.Value[i])
		}
	}

	return
}

// MulScaleInvariant multiplies op0 with op1 without relinearization and using scale invariant tensoring (BFV-style), and returns the result in opOut.
// This tensoring increases the noise by a constant factor regardless of the current noise, thus no rescaling is required with subsequent multiplications if they are
// performed with the invariant tensoring procedure. Rescaling can still be useful to reduce the size of the ciphertext, once the noise is higher than the prime
// that will be used for the rescaling or to ensure that the noise is minimal before using the regular tensoring.
// The procedure will return an error if either op0.Degree or op1.Degree > 1.
// The procedure will return an error if the evaluator was not created with an relinearization key.
//
// inputs:
//   - op0: an *[rlwe.Ciphertext]
//   - op1:
//   - [rlwe.ElementInterface][[ring.Poly]]
//   - *big.Int, uint64, int64, int
//   - []uint64 or []int64 (of size at most N where N is the smallest integer satisfying PlaintextModulus = 1 mod 2N)
//   - opOut: an *[rlwe.Ciphertext]
//
// If op1 is an [rlwe.ElementInterface][[ring.Poly]]:
//   - the level of opOut will be updated to min(op0.Level(), op1.Level())
//   - the scale of opOut will be to op0.Scale * op1.Scale * (-Q mod T)^{-1} mod T
func (eval Evaluator) MulScaleInvariant(op0 *rlwe.Ciphertext, op1 rlwe.Operand, opOut *rlwe.Ciphertext) (err error) {
	switch op1 := op1.(type) {
	case rlwe.ElementInterface[ring.Poly]:

		_, level, err := eval.InitOutputBinaryOp(op0.El(), op1.El(), 2, opOut.El())
		if err != nil {
			return fmt.Errorf("cannot MulInvariant: %w", err)
		}

		opOut.Resize(opOut.Degree(), level)

		if op1.Degree() == 0 {

			if err = eval.tensorStandard(op0, op1.El(), false, opOut); err != nil {
				return fmt.Errorf("cannot MulInvariant: %w", err)
			}

		} else {

			if err = eval.tensorScaleInvariant(op0, op1.El(), false, opOut); err != nil {
				return fmt.Errorf("cannot MulInvariant: %w", err)
			}
		}
	case []uint64, []int64:

		_, level, err := eval.InitOutputUnaryOp(op0.El(), opOut.El())
		if err != nil {
			return fmt.Errorf("cannot MulInvariant: %w", err)
		}

		opOut.Resize(op0.Degree(), level)

		// Instantiates new plaintext from buffer
		pt, err := rlwe.NewPlaintextAtLevelFromPoly(level, eval.buffQ[0])

		// This error should not happen, unless the evaluator's buffer were
		// improperly tempered with. If it does happen, there is no way to
		// recover from it.
		if err != nil {
			panic(err)
		}
		pt.MetaData = op0.MetaData.CopyNew() // Sets the metadata, notably matches scales
		pt.Scale = rlwe.NewScale(1)

		// Encodes the vector on the plaintext
		if err := eval.Encoder.Encode(op1, pt); err != nil {
			return err
		}

		if err = eval.tensorStandard(op0, pt.El(), false, opOut); err != nil {
			return fmt.Errorf("cannot MulInvariant: %w", err)
		}

	default:
		if err = eval.Mul(op0, op1, opOut); err != nil {
			return fmt.Errorf("cannot MulInvariant: %w", err)
		}
	}
	return
}

// MulScaleInvariantNew multiplies op0 with op1 without relinearization and using scale invariant tensoring (BFV-style), and returns the result in a new *[rlwe.Ciphertext] opOut.
// This tensoring increases the noise by a constant factor regardless of the current noise, thus no rescaling is required with subsequent multiplications if they are
// performed with the invariant tensoring procedure. Rescaling can still be useful to reduce the size of the ciphertext, once the noise is higher than the prime
// that will be used for the rescaling or to ensure that the noise is minimal before using the regular tensoring.
// The procedure will return an error if either op0.Degree or op1.Degree > 1.
// The procedure will return an error if the evaluator was not created with an relinearization key.
//
// inputs:
//   - op0: an *[rlwe.Ciphertext]
//   - op1:
//   - [rlwe.ElementInterface][[ring.Poly]]
//   - *big.Int, uint64, int64, int
//   - []uint64 or []int64 (of size at most N where N is the smallest integer satisfying PlaintextModulus = 1 mod 2N)
//
// If op1 is an [rlwe.ElementInterface][[ring.Poly]]:
//   - the level of opOut will be to min(op0.Level(), op1.Level())
//   - the scale of opOut will be to op0.Scale * op1.Scale * (-Q mod PlaintextModulus)^{-1} mod PlaintextModulus
func (eval Evaluator) MulScaleInvariantNew(op0 *rlwe.Ciphertext, op1 rlwe.Operand) (opOut *rlwe.Ciphertext, err error) {
	switch op1 := op1.(type) {
	case rlwe.ElementInterface[ring.Poly]:
		opOut = NewCiphertext(eval.parameters, op0.Degree()+op1.Degree(), utils.Min(op0.Level(), op1.Level()))
	default:
		opOut = NewCiphertext(eval.parameters, op0.Degree(), op0.Level())
	}
	return opOut, eval.MulScaleInvariant(op0, op1, opOut)
}

// MulRelinScaleInvariant multiplies op0 with op1 with relinearization and using scale invariant tensoring (BFV-style), and returns the result in opOut.
// This tensoring increases the noise by a constant factor regardless of the current noise, thus no rescaling is required with subsequent multiplications if they are
// performed with the invariant tensoring procedure. Rescaling can still be useful to reduce the size of the ciphertext, once the noise is higher than the prime
// that will be used for the rescaling or to ensure that the noise is minimal before using the regular tensoring.
// The procedure will return an error if either op0.Degree or op1.Degree > 1.
// The procedure will return an error if the evaluator was not created with an relinearization key.
//
// inputs:
//   - op0: an *[rlwe.Ciphertext]
//   - op1:
//   - [rlwe.ElementInterface][[ring.Poly]]
//   - *big.Int, uint64, int64, int
//   - []uint64 or []int64 (of size at most N where N is the smallest integer satisfying PlaintextModulus = 1 mod 2N)
//   - opOut: an *[rlwe.Ciphertext]
//
// If op1 is an [rlwe.ElementInterface][[ring.Poly]]:
//   - the level of opOut will be updated to min(op0.Level(), op1.Level())
//   - the scale of opOut will be to op0.Scale * op1.Scale * (-Q mod PlaintextModulus)^{-1} mod PlaintextModulus
func (eval Evaluator) MulRelinScaleInvariant(op0 *rlwe.Ciphertext, op1 rlwe.Operand, opOut *rlwe.Ciphertext) (err error) {
	switch op1 := op1.(type) {
	case rlwe.ElementInterface[ring.Poly]:

		_, level, err := eval.InitOutputBinaryOp(op0.El(), op1.El(), 2, opOut.El())
		if err != nil {
			return fmt.Errorf("cannot MulRelinInvariant: %w", err)
		}

		opOut.Resize(opOut.Degree(), level)

		if op1.Degree() == 0 {

			if err = eval.tensorStandard(op0, op1.El(), true, opOut); err != nil {
				return fmt.Errorf("cannot MulRelinInvariant: %w", err)
			}

		} else {

			if err = eval.tensorScaleInvariant(op0, op1.El(), true, opOut); err != nil {
				return fmt.Errorf("cannot MulRelinInvariant: %w", err)
			}
		}

	case []uint64, []int64:

		_, level, err := eval.InitOutputUnaryOp(op0.El(), opOut.El())

		if err != nil {
			return fmt.Errorf("cannot MulRelinInvariant: %w", err)
		}

		opOut.Resize(op0.Degree(), level)

		// Instantiates new plaintext from buffer
		pt, err := rlwe.NewPlaintextAtLevelFromPoly(level, eval.buffQ[0])

		// This error should not happen, unless the evaluator's buffer were
		// improperly tempered with. If it does happen, there is no way to
		// recover from it.
		if err != nil {
			panic(err)
		}

		pt.MetaData = op0.MetaData.CopyNew() // Sets the metadata, notably matches scales
		pt.Scale = rlwe.NewScale(1)

		// Encodes the vector on the plaintext
		if err := eval.Encoder.Encode(op1, pt); err != nil {
			return fmt.Errorf("cannot MulRelinInvariant: %w", err)
		}

		if err = eval.tensorStandard(op0, pt.El(), true, opOut); err != nil {
			return fmt.Errorf("cannot MulRelinInvariant: %w", err)
		}

	case uint64, int64, int, *big.Int:
		if err = eval.Mul(op0, op1, opOut); err != nil {
			return fmt.Errorf("cannot MulRelinInvariant: %w", err)
		}
	default:
		return fmt.Errorf("cannot MulRelinInvariant: invalid op1.(Type), expected rlwe.ElementInterface[ring.Poly], []uint64, []int64, uint64, int64 or int, but got %T", op1)
	}
	return
}

// MulRelinScaleInvariantNew multiplies op0 with op1 with relinearization and using scale invariant tensoring (BFV-style), and returns the result in a new *[rlwe.Ciphertext] opOut.
// This tensoring increases the noise by a constant factor regardless of the current noise, thus no rescaling is required with subsequent multiplications if they are
// performed with the invariant tensoring procedure. Rescaling can still be useful to reduce the size of the ciphertext, once the noise is higher than the prime
// that will be used for the rescaling or to ensure that the noise is minimal before using the regular tensoring.
// The procedure will return an error if either op0.Degree or op1.Degree > 1.
// The procedure will return an error if the evaluator was not created with an relinearization key.
//
// inputs:
//   - op0: an *[rlwe.Ciphertext]
//   - op1:
//   - [rlwe.ElementInterface][[ring.Poly]]
//   - *big.Int, uint64, int64, int
//   - []uint64 or []int64 (of size at most N where N is the smallest integer satisfying PlaintextModulus = 1 mod 2N)
//
// If op1 is an [rlwe.ElementInterface][[ring.Poly]]:
//   - the level of opOut will be to min(op0.Level(), op1.Level())
//   - the scale of opOut will be to op0.Scale * op1.Scale * (-Q mod PlaintextModulus)^{-1} mod PlaintextModulus
func (eval Evaluator) MulRelinScaleInvariantNew(op0 *rlwe.Ciphertext, op1 rlwe.Operand) (opOut *rlwe.Ciphertext, err error) {
	switch op1 := op1.(type) {
	case rlwe.ElementInterface[ring.Poly]:
		opOut = NewCiphertext(eval.parameters, 1, utils.Min(op0.Level(), op1.Level()))
	default:
		opOut = NewCiphertext(eval.parameters, op0.Degree(), op0.Level())
	}

	if err = eval.MulRelinScaleInvariant(op0, op1, opOut); err != nil {
		return nil, fmt.Errorf("cannot MulRelinInvariantNew: %w", err)
	}
	return
}

// tensorScaleInvariant computes (ct0 x ct1) * (t/Q) and stores the result in opOut.
func (eval Evaluator) tensorScaleInvariant(ct0 *rlwe.Ciphertext, ct1 *rlwe.Element[ring.Poly], relin bool, opOut *rlwe.Ciphertext) (err error) {

	level := opOut.Level()

	levelQMul := eval.levelQMul[level]

	// Avoid overwriting if the second input is the output
	var tmp0Q0, tmp1Q0 *rlwe.Element[ring.Poly]
	if ct1 == opOut.El() {
		tmp0Q0, tmp1Q0 = ct1, ct0.El()
	} else {
		tmp0Q0, tmp1Q0 = ct0.El(), ct1
	}

	tmp0Q1 := &rlwe.Element[ring.Poly]{Value: eval.buffQMul[0:3]}
	tmp1Q1 := &rlwe.Element[ring.Poly]{Value: eval.buffQMul[3:5]}
	tmp2Q1 := tmp0Q1

	eval.modUpAndNTT(level, levelQMul, tmp0Q0, tmp0Q1)

	if tmp0Q0 != tmp1Q0 {
		eval.modUpAndNTT(level, levelQMul, tmp1Q0, tmp1Q1)
	}

	var c2 ring.Poly
	if !relin {
		opOut.Resize(2, opOut.Level())
		c2 = opOut.Value[2]
	} else {
		opOut.Resize(1, opOut.Level())
		c2 = eval.buffQ[2]
	}

	tmp2Q0 := &rlwe.Element[ring.Poly]{Value: []ring.Poly{opOut.Value[0], opOut.Value[1], c2}}

	eval.tensorLowDeg(level, levelQMul, tmp0Q0, tmp1Q0, tmp2Q0, tmp0Q1, tmp1Q1, tmp2Q1)

	eval.quantize(level, levelQMul, tmp2Q0.Value[0], tmp2Q1.Value[0])
	eval.quantize(level, levelQMul, tmp2Q0.Value[1], tmp2Q1.Value[1])
	eval.quantize(level, levelQMul, tmp2Q0.Value[2], tmp2Q1.Value[2])

	if relin {

		var rlk *rlwe.RelinearizationKey

		if rlk, err = eval.GetRelinearizationKey(); err != nil {
			return fmt.Errorf("cannot TensorInvariant: %w", err)
		}

		tmpCt := &rlwe.Ciphertext{}
		tmpCt.Value = []ring.Poly{eval.BuffQP[1].Q, eval.BuffQP[2].Q}
		tmpCt.MetaData = &rlwe.MetaData{}
		tmpCt.IsNTT = true

		eval.GadgetProduct(level, c2, &rlk.GadgetCiphertext, tmpCt)

		ringQ := eval.parameters.RingQ().AtLevel(level)

		ringQ.Add(opOut.Value[0], tmpCt.Value[0], opOut.Value[0])
		ringQ.Add(opOut.Value[1], tmpCt.Value[1], opOut.Value[1])
	}

	opOut.Scale = MulScaleInvariant(eval.parameters, ct0.Scale, tmp1Q0.Scale, level)

	return
}

// MulScaleInvariant returns c = a * b / (-Q[level] mod PlaintextModulus), where a, b are the input scale,
// level the level at which the operation is carried out and and c is the new scale after performing the
// invariant tensoring (BFV-style).
func MulScaleInvariant(params Parameters, a, b rlwe.Scale, level int) (c rlwe.Scale) {
	c = a.Mul(b)
	qModTNeg := new(big.Int).Mod(params.RingQ().ModulusAtLevel[level], new(big.Int).SetUint64(params.PlaintextModulus())).Uint64()
	qModTNeg = params.PlaintextModulus() - qModTNeg
	c = c.Div(params.NewScale(qModTNeg))
	return
}

func (eval Evaluator) modUpAndNTT(level, levelQMul int, ctQ0, ctQ1 *rlwe.Element[ring.Poly]) {
	ringQ, ringQMul := eval.parameters.RingQ().AtLevel(level), eval.parameters.RingQMul().AtLevel(levelQMul)
	for i := range ctQ0.Value {
		ringQ.INTT(ctQ0.Value[i], eval.buffQ[0])
		eval.basisExtenderQ1toQ2.ModUpQtoP(level, levelQMul, eval.buffQ[0], ctQ1.Value[i])
		ringQMul.NTTLazy(ctQ1.Value[i], ctQ1.Value[i])
	}
}

func (eval Evaluator) tensorLowDeg(level, levelQMul int, ct0Q0, ct1Q0, ct2Q0, ct0Q1, ct1Q1, ct2Q1 *rlwe.Element[ring.Poly]) {

	ringQ, ringQMul := eval.parameters.RingQ().AtLevel(level), eval.parameters.RingQMul().AtLevel(levelQMul)

	c00 := eval.buffQ[0]
	c01 := eval.buffQ[1]

	ringQ.MForm(ct0Q0.Value[0], c00)
	ringQ.MForm(ct0Q0.Value[1], c01)

	c00M := eval.buffQMul[5]
	c01M := eval.buffQMul[6]

	ringQMul.MForm(ct0Q1.Value[0], c00M)
	ringQMul.MForm(ct0Q1.Value[1], c01M)

	// Squaring case
	if ct0Q0 == ct1Q0 {
		ringQ.MulCoeffsMontgomery(c00, ct0Q0.Value[0], ct2Q0.Value[0]) // c0 = c0[0]*c0[0]
		ringQ.MulCoeffsMontgomery(c01, ct0Q0.Value[1], ct2Q0.Value[2]) // c2 = c0[1]*c0[1]
		ringQ.MulCoeffsMontgomery(c00, ct0Q0.Value[1], ct2Q0.Value[1]) // c1 = 2*c0[0]*c0[1]
		ringQ.AddLazy(ct2Q0.Value[1], ct2Q0.Value[1], ct2Q0.Value[1])

		ringQMul.MulCoeffsMontgomery(c00M, ct0Q1.Value[0], ct2Q1.Value[0])
		ringQMul.MulCoeffsMontgomery(c01M, ct0Q1.Value[1], ct2Q1.Value[2])
		ringQMul.MulCoeffsMontgomery(c00M, ct0Q1.Value[1], ct2Q1.Value[1])
		ringQMul.AddLazy(ct2Q1.Value[1], ct2Q1.Value[1], ct2Q1.Value[1])

		// Normal case
	} else {
		ringQ.MulCoeffsMontgomery(c00, ct1Q0.Value[0], ct2Q0.Value[0]) // c0 = c0[0]*c1[0]
		ringQ.MulCoeffsMontgomery(c01, ct1Q0.Value[1], ct2Q0.Value[2]) // c2 = c0[1]*c1[1]
		ringQ.MulCoeffsMontgomery(c00, ct1Q0.Value[1], ct2Q0.Value[1]) // c1 = c0[0]*c1[1] + c0[1]*c1[0]
		ringQ.MulCoeffsMontgomeryThenAddLazy(c01, ct1Q0.Value[0], ct2Q0.Value[1])

		ringQMul.MulCoeffsMontgomery(c00M, ct1Q1.Value[0], ct2Q1.Value[0])
		ringQMul.MulCoeffsMontgomery(c01M, ct1Q1.Value[1], ct2Q1.Value[2])
		ringQMul.MulCoeffsMontgomery(c00M, ct1Q1.Value[1], ct2Q1.Value[1])
		ringQMul.MulCoeffsMontgomeryThenAddLazy(c01M, ct1Q1.Value[0], ct2Q1.Value[1])
	}
}

func (eval Evaluator) quantize(level, levelQMul int, c2Q1, c2Q2 ring.Poly) {

	ringQ, ringQMul := eval.parameters.RingQ().AtLevel(level), eval.parameters.RingQMul().AtLevel(levelQMul)

	// Applies the inverse NTT to the ciphertext, scales down the ciphertext
	// by t/q and reduces its basis from QP to Q

	ringQ.INTTLazy(c2Q1, c2Q1)
	ringQMul.INTTLazy(c2Q2, c2Q2)

	// Extends the basis Q of ct(x) to the basis P and Divides (ct(x)Q -> P) by Q
	eval.basisExtenderQ1toQ2.ModDownQPtoP(level, levelQMul, c2Q1, c2Q2, c2Q2) // QP / Q -> P

	// Centers ct(x)P by (P-1)/2 and extends ct(x)P to the basis Q
	eval.basisExtenderQ1toQ2.ModUpPtoQ(levelQMul, level, c2Q2, c2Q1)

	// (ct(x)/Q)*T, doing so only requires that Q*P > Q*Q, faster but adds error ~|T|
	ringQ.MulScalar(c2Q1, eval.parameters.PlaintextModulus(), c2Q1)

	ringQ.NTT(c2Q1, c2Q1)
}

// MulThenAdd multiplies op0 with op1 using standard tensoring and without relinearization, and adds the result on opOut.
// The procedure will return an error if either op0.Degree() or op1.Degree() > 1.
// The procedure will return an error if either op0 == opOut or op1 == opOut.
//
// inputs:
//   - op0: an *[rlwe.Ciphertext]
//   - op1:
//   - [rlwe.ElementInterface][[ring.Poly]]
//   - *big.Int, uint64, int64, int
//   - []uint64 or []int64 (of size at most N where N is the smallest integer satisfying PlaintextModulus = 1 mod 2N)
//   - opOut: an *[rlwe.Ciphertext]
//
// If op1 is an [rlwe.ElementInterface][[ring.Poly]] and opOut.Scale != op1.Scale * op0.Scale, then a scale matching operation will
// be automatically carried out to ensure that addition is performed between operands of the same scale.
// This scale matching operation will increase the noise by a small factor.
// For this reason it is preferable to ensure that opOut.Scale == op1.Scale * op0.Scale when calling this method.
func (eval Evaluator) MulThenAdd(op0 *rlwe.Ciphertext, op1 rlwe.Operand, opOut *rlwe.Ciphertext) (err error) {

	switch op1 := op1.(type) {
	case rlwe.ElementInterface[ring.Poly]:

		_, level, err := eval.InitOutputBinaryOp(op0.El(), op1.El(), 2, opOut.El())
		if err != nil {
			return fmt.Errorf("cannot MulThenAdd: %w", err)
		}

		if op0.El() == opOut.El() || op1.El() == opOut.El() {
			return fmt.Errorf("cannot MulThenAdd: opOut must be different from op0 and op1")
		}

		opOut.Resize(opOut.Degree(), level)

		if err = eval.mulRelinThenAdd(op0, op1.El(), false, opOut); err != nil {
			return fmt.Errorf("cannot MulThenAdd: %w", err)
		}

	case *big.Int:

		_, level, err := eval.InitOutputUnaryOp(op0.El(), opOut.El())

		if err != nil {
			return fmt.Errorf("cannot MulThenAdd: %w", err)
		}

		opOut.Resize(op0.Degree(), opOut.Level())

		ringQ := eval.parameters.RingQ().AtLevel(level)

		s := eval.parameters.RingT().SubRings[0]

		// op1 *= (op1.Scale / opOut.Scale)
		if op0.Scale.Cmp(opOut.Scale) != 0 {
			ratio := ring.ModExp(op0.Scale.Uint64(), s.Modulus-2, s.Modulus)
			ratio = ring.BRed(ratio, opOut.Scale.Uint64(), s.Modulus, s.BRedConstant)
			op1.Mul(op1, new(big.Int).SetUint64(ratio))
		}

		TBig := eval.parameters.RingT().ModulusAtLevel[0]

		op1.Mod(op1, TBig)

		// If op1 > T/2 then subtract T to minimize the noise
		if op1.Cmp(new(big.Int).Rsh(TBig, 1)) == 1 {
			op1.Sub(op1, TBig)
		}

		for i := 0; i < op0.Degree()+1; i++ {
			ringQ.MulScalarBigintThenAdd(op0.Value[i], op1, opOut.Value[i])
		}

	case int:
		return eval.MulThenAdd(op0, new(big.Int).SetInt64(int64(op1)), opOut)
	case int64:
		return eval.MulThenAdd(op0, new(big.Int).SetInt64(op1), opOut)
	case uint64:
		return eval.MulThenAdd(op0, new(big.Int).SetUint64(op1), opOut)
	case []uint64, []int64:

		_, level, err := eval.InitOutputUnaryOp(op0.El(), opOut.El())

		if err != nil {
			return fmt.Errorf("cannot MulThenAdd: %w", err)
		}

		opOut.Resize(op0.Degree(), opOut.Level())

		// Instantiates new plaintext from buffer
		pt, err := rlwe.NewPlaintextAtLevelFromPoly(level, eval.buffQ[0])

		// This error should not happen, unless the evaluator's buffer were
		// improperly tempered with. If it does happen, there is no way to
		// recover from it.
		if err != nil {
			panic(err)
		}
		pt.MetaData = op0.MetaData.CopyNew() // Sets the metadata, notably matches scales

		// op1 *= (op1.Scale / opOut.Scale)
		if op0.Scale.Cmp(opOut.Scale) != 0 {
			s := eval.parameters.RingT().SubRings[0]
			ratio := ring.ModExp(op0.Scale.Uint64(), s.Modulus-2, s.Modulus)
			pt.Scale = rlwe.NewScale(ring.BRed(ratio, opOut.Scale.Uint64(), s.Modulus, s.BRedConstant))
		} else {
			pt.Scale = rlwe.NewScale(1)
		}

		// Encodes the vector on the plaintext
		if err := eval.Encoder.Encode(op1, pt); err != nil {
			return fmt.Errorf("cannot MulThenAdd: %w", err)
		}

		if err = eval.MulThenAdd(op0, pt, opOut); err != nil {
			return fmt.Errorf("cannot MulThenAdd: %w", err)
		}

	default:
		return fmt.Errorf("cannot MulThenAdd: invalid op1.(Type), expected rlwe.ElementInterface[ring.Poly], []uint64, []int64, *big.Int, uint64, int64 or int, but got %T", op1)
	}

	return
}

// MulRelinThenAdd multiplies op0 with op1 using standard tensoring and with relinearization, and adds the result on opOut.
// The procedure will return an error if either op0.Degree() or op1.Degree() > 1.
// The procedure will return an error if either op0 == opOut or op1 == opOut.
//
// inputs:
//   - op0: an *[rlwe.Ciphertext]
//   - op1:
//   - [rlwe.ElementInterface][[ring.Poly]]
//   - *big.Int, uint64, int64, int
//   - []uint64 or []int64 (of size at most N where N is the smallest integer satisfying PlaintextModulus = 1 mod 2N)
//   - opOut: an *[rlwe.Ciphertext]
//
// If op1 is an [rlwe.ElementInterface][[ring.Poly]] and opOut.Scale != op1.Scale * op0.Scale, then a scale matching operation will
// be automatically carried out to ensure that addition is performed between operands of the same scale.
// This scale matching operation will increase the noise by a small factor.
// For this reason it is preferable to ensure that opOut.Scale == op1.Scale * op0.Scale when calling this method.
func (eval Evaluator) MulRelinThenAdd(op0 *rlwe.Ciphertext, op1 rlwe.Operand, opOut *rlwe.Ciphertext) (err error) {
	switch op1 := op1.(type) {
	case rlwe.ElementInterface[ring.Poly]:
		if op1.Degree() == 0 {
			return eval.MulThenAdd(op0, op1, opOut)
		} else {

			_, level, err := eval.InitOutputBinaryOp(op0.El(), op1.El(), 2, opOut.El())
			if err != nil {
				return fmt.Errorf("cannot MulThenAdd: %w", err)
			}

			if op0.El() == opOut.El() || op1.El() == opOut.El() {
				return fmt.Errorf("cannot MulThenAdd: opOut must be different from op0 and op1")
			}

			opOut.Resize(opOut.Degree(), level)

			return eval.mulRelinThenAdd(op0, op1.El(), true, opOut)
		}
	default:
		return eval.MulThenAdd(op0, op1, opOut)
	}
}

func (eval Evaluator) mulRelinThenAdd(op0 *rlwe.Ciphertext, op1 *rlwe.Element[ring.Poly], relin bool, opOut *rlwe.Ciphertext) (err error) {

	level := opOut.Level()

	ringQ := eval.parameters.RingQ().AtLevel(level)
	sT := eval.parameters.RingT().SubRings[0]

	var c00, c01, c0, c1, c2 ring.Poly

	// Case Ciphertext (x) Ciphertext
	if op0.Degree() == 1 && op1.Degree() == 1 {

		c00 = eval.buffQ[0]
		c01 = eval.buffQ[1]

		c0 = opOut.Value[0]
		c1 = opOut.Value[1]

		if !relin {
			opOut.Resize(2, level)
			c2 = opOut.Value[2]
		} else {
			opOut.Resize(utils.Max(1, opOut.Degree()), level)
			c2 = eval.buffQ[2]
		}

		tmp0, tmp1 := op0.El(), op1.El()

		// If op0.Scale * op1.Scale != opOut.Scale then
		// updates op1.Scale and opOut.Scale
		var r0 uint64 = 1
		if targetScale := ring.BRed(op0.Scale.Uint64(), op1.Scale.Uint64(), sT.Modulus, sT.BRedConstant); opOut.Scale.Cmp(eval.parameters.NewScale(targetScale)) != 0 {
			var r1 uint64
			r0, r1, _ = eval.matchScalesBinary(targetScale, opOut.Scale.Uint64())

			for i := range opOut.Value {
				ringQ.MulScalar(opOut.Value[i], r1, opOut.Value[i])
			}

			opOut.Scale = opOut.Scale.Mul(eval.parameters.NewScale(r1))
		}

		// Multiply by T * 2^{64} * 2^{64} -> result multipled by T and switched in the Montgomery domain
		ringQ.MulRNSScalarMontgomery(tmp0.Value[0], eval.tMontgomery, c00)
		ringQ.MulRNSScalarMontgomery(tmp0.Value[1], eval.tMontgomery, c01)

		// Scales the input to the output scale
		if r0 != 1 {
			ringQ.MulScalar(c00, r0, c00)
			ringQ.MulScalar(c01, r0, c01)
		}

		ringQ.MulCoeffsMontgomeryThenAdd(c00, tmp1.Value[0], c0) // c0 += c[0]*c[0]
		ringQ.MulCoeffsMontgomeryThenAdd(c00, tmp1.Value[1], c1) // c1 += c[0]*c[1]
		ringQ.MulCoeffsMontgomeryThenAdd(c01, tmp1.Value[0], c1) // c1 += c[1]*c[0]

		if relin {

			var rlk *rlwe.RelinearizationKey
			var err error
			if rlk, err = eval.CheckAndGetRelinearizationKey(); err != nil {
				return fmt.Errorf("cannot Relinearize: %w", err)
			}

			ringQ.MulCoeffsMontgomery(c01, tmp1.Value[1], c2) // c2 += c[1]*c[1]

			tmpCt := &rlwe.Ciphertext{}
			tmpCt.Value = []ring.Poly{eval.BuffQP[1].Q, eval.BuffQP[2].Q}
			tmpCt.MetaData = &rlwe.MetaData{}
			tmpCt.IsNTT = true

			eval.GadgetProduct(level, c2, &rlk.GadgetCiphertext, tmpCt)

			ringQ.Add(opOut.Value[0], tmpCt.Value[0], opOut.Value[0])
			ringQ.Add(opOut.Value[1], tmpCt.Value[1], opOut.Value[1])

		} else {
			ringQ.MulCoeffsMontgomeryThenAdd(c01, tmp1.Value[1], c2) // c2 += c[1]*c[1]
		}

		// Case Plaintext (x) Ciphertext or Ciphertext (x) Plaintext
	} else {

		opOut.Resize(utils.Max(op0.Degree(), opOut.Degree()), level)

		c00 := eval.buffQ[0]

		// Multiply by T * 2^{64} * 2^{64} -> result multipled by T and switched in the Montgomery domain
		ringQ.MulRNSScalarMontgomery(op1.El().Value[0], eval.tMontgomery, c00)

		// If op0.Scale * op1.Scale != opOut.Scale then
		// updates op1.Scale and opOut.Scale
		var r0 = uint64(1)
		if targetScale := ring.BRed(op0.Scale.Uint64(), op1.Scale.Uint64(), sT.Modulus, sT.BRedConstant); opOut.Scale.Cmp(eval.parameters.NewScale(targetScale)) != 0 {
			var r1 uint64
			r0, r1, _ = eval.matchScalesBinary(targetScale, opOut.Scale.Uint64())

			for i := range opOut.Value {
				ringQ.MulScalar(opOut.Value[i], r1, opOut.Value[i])
			}

			opOut.Scale = opOut.Scale.Mul(eval.parameters.NewScale(r1))
		}

		if r0 != 1 {
			ringQ.MulScalar(c00, r0, c00)
		}

		for i := range op0.Value {
			ringQ.MulCoeffsMontgomeryThenAdd(op0.Value[i], c00, opOut.Value[i])
		}
	}

	return
}

// Rescale divides (rounded) op0 by the last prime of the moduli chain and returns the result on opOut.
// This procedure divides the noise by the last prime of the moduli chain while preserving
// the MSB-plaintext bits.
// The procedure will return an error if:
//   - op0.Level() == 0 (the input ciphertext is already at the last prime)
//   - opOut.Level() < op0.Level() - 1 (not enough space to store the result)
//
// The scale of opOut will be updated to op0.Scale * qi^{-1} mod PlaintextModulus where qi is the prime consumed by
// the rescaling operation.
// Note that if the evaluator has been instantiated as scale-invariant (BFV-style), then Rescale is a nop.
func (eval Evaluator) Rescale(op0, opOut *rlwe.Ciphertext) (err error) {

	if eval.ScaleInvariant {
		return nil
	}

	if op0.MetaData == nil || opOut.MetaData == nil {
		return fmt.Errorf("cannot Rescale: op0.MetaData or opOut.MetaData is nil")
	}

	if op0.Level() == 0 {
		return fmt.Errorf("cannot rescale: op0 already at level 0")
	}

	if opOut.Level() < op0.Level()-1 {
		return fmt.Errorf("cannot rescale: opOut.Level() < op0.Level()-1")
	}

	level := op0.Level()
	ringQ := eval.parameters.RingQ().AtLevel(level)

	for i := range opOut.Value {
		ringQ.DivRoundByLastModulusNTT(op0.Value[i], eval.buffQ[0], opOut.Value[i])
	}

	opOut.Resize(opOut.Degree(), level-1)

	*opOut.MetaData = *op0.MetaData
	opOut.Scale = op0.Scale.Div(eval.parameters.NewScale(ringQ.SubRings[level].Modulus))
	return
}

// RelinearizeNew applies the relinearization procedure on op0 and returns the result in a new opOut.
func (eval Evaluator) RelinearizeNew(op0 *rlwe.Ciphertext) (opOut *rlwe.Ciphertext, err error) {
	opOut = NewCiphertext(eval.parameters, 1, op0.Level())
	return opOut, eval.Relinearize(op0, opOut)
}

// ApplyEvaluationKeyNew re-encrypts op0 under a different key and returns the result in a new opOut.
// It requires a [rlwe.EvaluationKey], which is computed from the key under which the Ciphertext is currently encrypted,
// and the key under which the Ciphertext will be re-encrypted.
// The procedure will return an error if either op0.Degree() or opOut.Degree() != 1.
func (eval Evaluator) ApplyEvaluationKeyNew(op0 *rlwe.Ciphertext, evk *rlwe.EvaluationKey) (opOut *rlwe.Ciphertext, err error) {
	opOut = NewCiphertext(eval.parameters, op0.Degree(), op0.Level())
	return opOut, eval.ApplyEvaluationKey(op0, evk, opOut)
}

// RotateColumnsNew rotates the columns of op0 by k positions to the left, and returns the result in a newly created element.
// The procedure will return an error if the corresponding Galois key has not been generated and attributed to the evaluator.
// The procedure will return an error if op0.Degree() != 1.
func (eval Evaluator) RotateColumnsNew(op0 *rlwe.Ciphertext, k int) (opOut *rlwe.Ciphertext, err error) {
	opOut = NewCiphertext(eval.parameters, op0.Degree(), op0.Level())
	return opOut, eval.RotateColumns(op0, k, opOut)
}

// RotateColumns rotates the columns of op0 by k positions to the left and returns the result in opOut.
// The procedure will return an error if the corresponding Galois key has not been generated and attributed to the evaluator.
// The procedure will return an error if either op0.Degree() or opOut.Degree() != 1.
func (eval Evaluator) RotateColumns(op0 *rlwe.Ciphertext, k int, opOut *rlwe.Ciphertext) (err error) {
	return eval.Automorphism(op0, eval.parameters.GaloisElement(k), opOut)
}

// RotateRowsNew swaps the rows of op0 and returns the result in a new opOut.
// The procedure will return an error if the corresponding Galois key has not been generated and attributed to the evaluator.
// The procedure will return an error if op0.Degree() != 1.
func (eval Evaluator) RotateRowsNew(op0 *rlwe.Ciphertext) (opOut *rlwe.Ciphertext, err error) {
	opOut = NewCiphertext(eval.parameters, op0.Degree(), op0.Level())
	return opOut, eval.RotateRows(op0, opOut)
}

// RotateRows swaps the rows of op0 and returns the result in op1.
// The procedure will return an error if the corresponding Galois key has not been generated and attributed to the evaluator.
// The procedure will return an error if either op0.Degree() or op1.Degree() != 1.
func (eval Evaluator) RotateRows(op0, opOut *rlwe.Ciphertext) (err error) {
	return eval.Automorphism(op0, eval.parameters.GaloisElementForRowRotation(), opOut)
}

// RotateHoistedLazyNew applies a series of rotations on the same ciphertext and returns each different rotation in a map indexed by the rotation.
// Results are not rescaled by P.
func (eval Evaluator) RotateHoistedLazyNew(level int, rotations []int, op0 *rlwe.Ciphertext, c2DecompQP []ringqp.Poly) (opOut map[int]*rlwe.Element[ringqp.Poly], err error) {
	opOut = make(map[int]*rlwe.Element[ringqp.Poly])
	for _, i := range rotations {
		if i != 0 {
			opOut[i] = rlwe.NewElementExtended(eval.parameters, 1, level, eval.parameters.MaxLevelP())
			if err = eval.AutomorphismHoistedLazy(level, op0, c2DecompQP, eval.parameters.GaloisElement(i), opOut[i]); err != nil {
				return
			}
		}
	}

	return
}

// MatchScalesAndLevel updates the both input ciphertexts to ensures that their scale matches.
// To do so it computes t0 * a = opOut * b such that:
//   - ct0.Scale * a = opOut.Scale: make the scales match.
//   - gcd(a, PlaintextModulus) == gcd(b, PlaintextModulus) == 1: ensure that the new scale is not a zero divisor if PlaintextModulus is not prime.
//   - |a+b| is minimal: minimize the added noise by the procedure.
func (eval Evaluator) MatchScalesAndLevel(ct0, opOut *rlwe.Ciphertext) {

	r0, r1, _ := eval.matchScalesBinary(ct0.Scale.Uint64(), opOut.Scale.Uint64())

	level := utils.Min(ct0.Level(), opOut.Level())

	ringQ := eval.parameters.RingQ().AtLevel(level)

	for _, el := range ct0.Value {
		ringQ.MulScalar(el, r0, el)
	}

	ct0.Resize(ct0.Degree(), level)
	ct0.Scale = ct0.Scale.Mul(eval.parameters.NewScale(r0))

	for _, el := range opOut.Value {
		ringQ.MulScalar(el, r1, el)
	}

	opOut.Resize(opOut.Degree(), level)
	opOut.Scale = opOut.Scale.Mul(eval.parameters.NewScale(r1))
}

func (eval Evaluator) GetRLWEParameters() *rlwe.Parameters {
	return eval.Evaluator.GetRLWEParameters()
}

func (eval Evaluator) matchScalesBinary(scale0, scale1 uint64) (r0, r1, e uint64) {

	ringT := eval.parameters.RingT()

	t := ringT.SubRings[0].Modulus
	tHalf := t >> 1
	GenBRedConstant := ringT.SubRings[0].BRedConstant

	// This should never happen and if it were to happen,
	// there is no way to recover from it.
	if utils.GCD(scale0, t) != 1 {
		panic("cannot matchScalesBinary: invalid ciphertext scale: gcd(scale, t) != 1")
	}

	var a = ringT.SubRings[0].Modulus
	var b uint64 = 0
	var A = ring.BRed(ring.ModExp(scale0, t-2, t), scale1, t, GenBRedConstant)
	var B uint64 = 1

	r0, r1 = A, B

	e = center(A, tHalf, t) + 1

	for A != 0 {

		q := a / A
		a, A = A, a%A
		b, B = B, ring.CRed(t+b-ring.BRed(B, q, t, GenBRedConstant), t)

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
