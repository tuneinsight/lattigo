// Package bfv is a depreciated placeholder package wrapping the bgv package for backward compatibility. This package will be removed in the next major version.
package bfv

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/bgv"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
)

func NewPlaintext(params Parameters, level int) (pt *rlwe.Plaintext) {
	return rlwe.NewPlaintext(params.Parameters, level)
}

func NewCiphertext(params Parameters, degree, level int) (ct *rlwe.Ciphertext) {
	return rlwe.NewCiphertext(params.Parameters, degree, level)
}

func NewEncryptor(params Parameters, key interface{}) rlwe.Encryptor {
	return rlwe.NewEncryptor(params.Parameters, key)
}

func NewDecryptor(params Parameters, key *rlwe.SecretKey) rlwe.Decryptor {
	return rlwe.NewDecryptor(params.Parameters, key)
}

func NewKeyGenerator(params Parameters) *rlwe.KeyGenerator {
	return rlwe.NewKeyGenerator(params.Parameters)
}

func NewPRNGEncryptor(params Parameters, key *rlwe.SecretKey) rlwe.PRNGEncryptor {
	return rlwe.NewPRNGEncryptor(params.Parameters, key)
}

type Encoder bgv.Encoder

func NewEncoder(params Parameters) Encoder {
	return bgv.NewEncoder(bgv.Parameters(params))
}

// Evaluator is an interface implementing the public methods of the eval.
type Evaluator interface {

	// Add: ct-ct & ct-pt & ct-scalar
	Add(op0 *rlwe.Ciphertext, op1 interface{}, op2 *rlwe.Ciphertext)
	AddNew(op0 *rlwe.Ciphertext, op1 interface{}) (op2 *rlwe.Ciphertext)

	// Sub: ct-ct & ct-pt & ct-scalar
	Sub(op0 *rlwe.Ciphertext, op1 interface{}, op2 *rlwe.Ciphertext)
	SubNew(op0 *rlwe.Ciphertext, op1 interface{}) (op2 *rlwe.Ciphertext)

	// Neg
	Neg(op0 *rlwe.Ciphertext, op1 *rlwe.Ciphertext)
	NegNew(op0 *rlwe.Ciphertext) (op1 *rlwe.Ciphertext)

	// Mul ct-ct & ct-pt & ct-scalar
	Mul(op0 *rlwe.Ciphertext, op1 interface{}, op2 *rlwe.Ciphertext)
	MulNew(op0 *rlwe.Ciphertext, op1 interface{}) (op2 *rlwe.Ciphertext)
	MulRelin(op0 *rlwe.Ciphertext, op1 interface{}, op2 *rlwe.Ciphertext)
	MulRelinNew(op0 *rlwe.Ciphertext, op1 interface{}) (op2 *rlwe.Ciphertext)

	// MulThenAdd ct-ct & ct-pt
	MulThenAdd(op0 *rlwe.Ciphertext, op1 interface{}, op2 *rlwe.Ciphertext)
	MulRelinThenAdd(op0 *rlwe.Ciphertext, op1 interface{}, op2 *rlwe.Ciphertext)

	// Degree Management
	RelinearizeNew(op0 *rlwe.Ciphertext) (op1 *rlwe.Ciphertext)
	Relinearize(op0 *rlwe.Ciphertext, op1 *rlwe.Ciphertext)

	// Error and Level management
	Rescale(op0, op1 *rlwe.Ciphertext) (err error)
	DropLevelNew(op0 *rlwe.Ciphertext, levels int) (op1 *rlwe.Ciphertext)
	DropLevel(op0 *rlwe.Ciphertext, levels int)

	// Column & Rows rotations
	RotateColumnsNew(op0 *rlwe.Ciphertext, k int) (op1 *rlwe.Ciphertext)
	RotateColumns(op0 *rlwe.Ciphertext, k int, op1 *rlwe.Ciphertext)
	RotateRows(op0 *rlwe.Ciphertext, op1 *rlwe.Ciphertext)
	RotateRowsNew(op0 *rlwe.Ciphertext) (op1 *rlwe.Ciphertext)

	//Polynomial Evaluation
	EvaluatePoly(op0 interface{}, pol *Polynomial, targetScale rlwe.Scale) (op1 *rlwe.Ciphertext, err error)
	EvaluatePolyVector(op0 interface{}, pols []*Polynomial, encoder Encoder, slotIndex map[int][]int, targetScale rlwe.Scale) (op1 *rlwe.Ciphertext, err error)

	// TODO
	LinearTransformNew(op0 *rlwe.Ciphertext, linearTransform interface{}) (op1 []*rlwe.Ciphertext)
	LinearTransform(op0 *rlwe.Ciphertext, linearTransform interface{}, op1 []*rlwe.Ciphertext)
	MultiplyByDiagMatrix(op0 *rlwe.Ciphertext, matrix LinearTransform, c2DecompQP []ringqp.Poly, op1 *rlwe.Ciphertext)
	MultiplyByDiagMatrixBSGS(op0 *rlwe.Ciphertext, matrix LinearTransform, c2DecompQP []ringqp.Poly, op1 *rlwe.Ciphertext)
	InnerSum(op0 *rlwe.Ciphertext, batch, n int, ctOut *rlwe.Ciphertext)
	Replicate(op0 *rlwe.Ciphertext, batch, n int, ctOut *rlwe.Ciphertext)

	// Key-Switching
	ApplyEvaluationKeyNew(op0 *rlwe.Ciphertext, evk *rlwe.EvaluationKey) (op1 *rlwe.Ciphertext)
	ApplyEvaluationKey(op0 *rlwe.Ciphertext, evk *rlwe.EvaluationKey, op1 *rlwe.Ciphertext)
	Automorphism(op0 *rlwe.Ciphertext, galEl uint64, op1 *rlwe.Ciphertext)
	AutomorphismHoisted(level int, op0 *rlwe.Ciphertext, c1DecompQP []ringqp.Poly, galEl uint64, op1 *rlwe.Ciphertext)
	RotateHoistedLazyNew(level int, rotations []int, op0 *rlwe.Ciphertext, c2DecompQP []ringqp.Poly) (op1 map[int]*rlwe.OperandQP)

	// Others
	GetRLWEEvaluator() *rlwe.Evaluator
	BuffQ() [3]*ring.Poly
	ShallowCopy() Evaluator
	WithKey(evk rlwe.EvaluationKeySetInterface) (eval Evaluator)
}

type evaluator struct {
	bgv.Evaluator
}

func NewEvaluator(params Parameters, evk rlwe.EvaluationKeySetInterface) Evaluator {
	return &evaluator{bgv.NewEvaluator(bgv.Parameters(params), evk)}
}

func (eval *evaluator) WithKey(evk rlwe.EvaluationKeySetInterface) Evaluator {
	return &evaluator{eval.Evaluator.WithKey(evk)}
}

func (eval *evaluator) ShallowCopy() Evaluator {
	return &evaluator{eval.Evaluator.ShallowCopy()}
}

// Mul multiplies op0 with op1 without relinearization and returns the result in op2.
// The procedure will panic if either op0 or op1 are have a degree higher than 1.
// The procedure will panic if op2.Degree != op0.Degree + op1.Degree.
func (eval *evaluator) Mul(op0 *rlwe.Ciphertext, op1 interface{}, op2 *rlwe.Ciphertext) {
	switch op1 := op1.(type) {
	case rlwe.Operand:
		eval.Evaluator.MulInvariant(op0, op1, op2)
	case uint64:
		eval.Evaluator.Mul(op0, op1, op0)
	default:
		panic(fmt.Sprintf("invalid op1.(Type), expected rlwe.Operand or uint64, but got %T", op1))
	}

}

func (eval *evaluator) MulNew(op0 *rlwe.Ciphertext, op1 interface{}) (op2 *rlwe.Ciphertext) {
	switch op1 := op1.(type) {
	case rlwe.Operand:
		return eval.Evaluator.MulInvariantNew(op0, op1)
	case uint64:
		return eval.Evaluator.MulNew(op0, op1)
	default:
		panic(fmt.Sprintf("invalid op1.(Type), expected rlwe.Operand or uint64, but got %T", op1))
	}
}

func (eval *evaluator) MulRelinNew(op0 *rlwe.Ciphertext, op1 interface{}) (op2 *rlwe.Ciphertext) {
	return eval.Evaluator.MulRelinInvariantNew(op0, op1)
}

func (eval *evaluator) MulRelin(op0 *rlwe.Ciphertext, op1 interface{}, op2 *rlwe.Ciphertext) {
	eval.Evaluator.MulRelinInvariant(op0, op1, op2)
}

type Polynomial bgv.Polynomial

func NewPoly(coeffs []uint64) (p *Polynomial) {
	poly := Polynomial(*bgv.NewPoly(coeffs))
	return &poly
}

type PowerBasis *bgv.PowerBasis

func NewPowerBasis(ct *rlwe.Ciphertext) (p *PowerBasis) {
	pb := PowerBasis(bgv.NewPowerBasis(ct))
	return &pb
}

func (eval *evaluator) EvaluatePoly(input interface{}, pol *Polynomial, targetScale rlwe.Scale) (opOut *rlwe.Ciphertext, err error) {
	poly := bgv.Polynomial(*pol)
	return eval.Evaluator.EvaluatePolyInvariant(input, &poly, targetScale)
}

func (eval *evaluator) EvaluatePolyVector(input interface{}, pols []*Polynomial, encoder Encoder, slotsIndex map[int][]int, targetScale rlwe.Scale) (opOut *rlwe.Ciphertext, err error) {

	polys := make([]*bgv.Polynomial, len(pols))

	for i := range polys {
		p := bgv.Polynomial(*pols[i])
		polys[i] = &p
	}

	return eval.Evaluator.EvaluatePolyVectorInvariant(input, polys, encoder, slotsIndex, targetScale)
}

type LinearTransform bgv.LinearTransform

func (lt *LinearTransform) Rotations() (rotations []int) {
	ll := bgv.LinearTransform(*lt)
	return ll.Rotations()
}

func NewLinearTransform(params Parameters, nonZeroDiags []int, level int, BSGSRatio float64) LinearTransform {
	return LinearTransform(bgv.NewLinearTransform(bgv.Parameters(params), nonZeroDiags, level, BSGSRatio))
}

func GenLinearTransform(ecd Encoder, dMat map[int][]uint64, level int, scale rlwe.Scale) LinearTransform {
	return LinearTransform(bgv.GenLinearTransform(bgv.Encoder(ecd), dMat, level, scale))
}

func GenLinearTransformBSGS(ecd Encoder, dMat map[int][]uint64, level int, scale rlwe.Scale, BSGSRatio float64) LinearTransform {
	return LinearTransform(bgv.GenLinearTransformBSGS(bgv.Encoder(ecd), dMat, level, scale, BSGSRatio))
}

func (eval *evaluator) LinearTransformNew(ctIn *rlwe.Ciphertext, linearTransform interface{}) (ctOut []*rlwe.Ciphertext) {

	var LTs []bgv.LinearTransform

	switch linearTransform := linearTransform.(type) {
	case LinearTransform:
		LTs = []bgv.LinearTransform{bgv.LinearTransform(linearTransform)}
	case []LinearTransform:
		LTs := make([]bgv.LinearTransform, len(linearTransform))
		for i := range LTs {
			LTs[i] = bgv.LinearTransform(linearTransform[i])
		}
	}

	return eval.Evaluator.LinearTransformNew(ctIn, LTs)
}

func (eval *evaluator) LinearTransform(ctIn *rlwe.Ciphertext, linearTransform interface{}, ctOut []*rlwe.Ciphertext) {
	var LTs []bgv.LinearTransform

	switch linearTransform := linearTransform.(type) {
	case LinearTransform:
		LTs = []bgv.LinearTransform{bgv.LinearTransform(linearTransform)}
	case []LinearTransform:
		LTs := make([]bgv.LinearTransform, len(linearTransform))
		for i := range LTs {
			LTs[i] = bgv.LinearTransform(linearTransform[i])
		}
	}

	eval.Evaluator.LinearTransform(ctIn, LTs, ctOut)
}

func (eval *evaluator) MultiplyByDiagMatrix(ctIn *rlwe.Ciphertext, matrix LinearTransform, BuffDecompQP []ringqp.Poly, ctOut *rlwe.Ciphertext) {
	eval.Evaluator.MultiplyByDiagMatrix(ctIn, bgv.LinearTransform(matrix), BuffDecompQP, ctOut)
}

func (eval *evaluator) MultiplyByDiagMatrixBSGS(ctIn *rlwe.Ciphertext, matrix LinearTransform, BuffDecompQP []ringqp.Poly, ctOut *rlwe.Ciphertext) {
	eval.Evaluator.MultiplyByDiagMatrixBSGS(ctIn, bgv.LinearTransform(matrix), BuffDecompQP, ctOut)
}
