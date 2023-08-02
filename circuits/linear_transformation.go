package circuits

import (
	"fmt"
	"sort"

	"github.com/tuneinsight/lattigo/v4/bgv"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

type Numeric interface {
	ckks.Float | bgv.Integer
}

type EvaluatorForLinearTransform interface {
	rlwe.ParameterProvider
	// TODO: separated int
	DecomposeNTT(levelQ, levelP, nbPi int, c2 ring.Poly, c2IsNTT bool, decompQP []ringqp.Poly)
	CheckAndGetGaloisKey(galEl uint64) (evk *rlwe.GaloisKey, err error)
	GadgetProductLazy(levelQ int, cx ring.Poly, gadgetCt *rlwe.GadgetCiphertext, ct *rlwe.Operand[ringqp.Poly])
	GadgetProductHoistedLazy(levelQ int, BuffQPDecompQP []ringqp.Poly, gadgetCt *rlwe.GadgetCiphertext, ct *rlwe.Operand[ringqp.Poly])
	AutomorphismHoistedLazy(levelQ int, ctIn *rlwe.Ciphertext, c1DecompQP []ringqp.Poly, galEl uint64, ctQP *rlwe.Operand[ringqp.Poly]) (err error)
	ModDownQPtoQNTT(levelQ, levelP int, p1Q, p1P, p2Q ring.Poly)
	AutomorphismIndex(uint64) []uint64

	GetEvaluatorBuffer() *rlwe.EvaluatorBuffers // TODO extract
}

type LinearTransformEvaluator struct {
	EvaluatorForLinearTransform
	*rlwe.EvaluatorBuffers
}

// EncoderInterface defines a set of common and scheme agnostic method provided by an Encoder struct.
type EncoderInterface[T Numeric, U *ring.Poly | ringqp.Poly | *rlwe.Plaintext] interface {
	Encode(values []T, metaData *rlwe.MetaData, output U) (err error)
}

// NewEvaluator instantiates a new LinearTransformEvaluator from an EvaluatorForLinearTransform.
// The method is allocation free if the underlying EvaluatorForLinearTransform returns a non-nil
// *rlwe.EvaluatorBuffers.
func NewEvaluator(eval EvaluatorForLinearTransform) (linTransEval *LinearTransformEvaluator) {
	linTransEval = new(LinearTransformEvaluator)
	linTransEval.EvaluatorForLinearTransform = eval
	linTransEval.EvaluatorBuffers = eval.GetEvaluatorBuffer()
	if linTransEval.EvaluatorBuffers == nil {
		linTransEval.EvaluatorBuffers = rlwe.NewEvaluatorBuffers(*eval.GetRLWEParameters())
	}
	return
}

// LinearTransformationParameters is a struct storing the parameterization of a
// linear transformation.
//
// # A homomorphic linear transformations on a ciphertext acts as evaluating
//
// Ciphertext([1 x n] vector) <- Ciphertext([1 x n] vector) x Plaintext([n x n] matrix)
//
// where n is the number of plaintext slots.
//
// The diagonal representation of a linear transformations is defined by first expressing
// the linear transformation through its nxn matrix and then reading the matrix diagonally.
//
// For example, the following nxn for n=4 matrix:
//
// 0 1 2 3 (diagonal index)
// | 1 2 3 0 |
// | 0 1 2 3 |
// | 3 0 1 2 |
// | 2 3 0 1 |
//
// its diagonal representation is comprised of 3 non-zero diagonals at indexes [0, 1, 2]:
// 0: [1, 1, 1, 1]
// 1: [2, 2, 2, 2]
// 2: [3, 3, 3, 3]
//
// Note that negative indexes can be used and will be interpreted modulo the matrix dimension.
//
// The diagonal representation is well suited for two reasons:
//  1. It is the effective format used during the homomorphic evaluation.
//  2. It enables on average a more compact and efficient representation of linear transformations
//     than their matrix representation by being able to only store the non-zero diagonals.
//
// Finally, some metrics about the time and storage complexity of homomorphic linear transformations:
//   - Storage: #diagonals polynomials mod Q_level * P
//   - Evaluation: #diagonals multiplications and 2sqrt(#diagonals) ciphertexts rotations.
type LinearTransformationParameters struct {
	// DiagonalsIndexList is the list of the non-zero diagonals of the square matrix.
	// A non zero diagonals is a diagonal with a least one non-zero element.
	DiagonalsIndexList []int

	// Level is the level at which to encode the linear transformation.
	Level int

	// Scale is the plaintext scale at which to encode the linear transformation.
	Scale rlwe.Scale

	// LogDimensions is the log2 dimensions of the matrix that can be SIMD packed
	// in a single plaintext polynomial.
	// This method is equivalent to params.PlaintextDimensions().
	// Note that the linear transformation is evaluated independently on each rows of
	// the SIMD packed matrix.
	LogDimensions ring.Dimensions

	// LogBabyStepGianStepRatio is the log2 of the ratio n1/n2 for n = n1 * n2 and
	// n is the dimension of the linear transformation. The number of Galois keys required
	// is minimized when this value is 0 but the overall complexity of the homomorphic evaluation
	// can be reduced by increasing the ratio (at the expanse of increasing the number of keys required).
	// If the value returned is negative, then the baby-step giant-step algorithm is not used
	// and the evaluation complexity (as well as the number of keys) becomes O(n) instead of O(sqrt(n)).
	LogBabyStepGianStepRatio int
}

type Diagonals[T Numeric] map[int][]T

// DiagonalsIndexList returns the list of the non-zero diagonals of the square matrix.
// A non zero diagonals is a diagonal with a least one non-zero element.
func (m Diagonals[T]) DiagonalsIndexList() (indexes []int) {
	indexes = make([]int, 0, len(m))
	for k := range m {
		indexes = append(indexes, k)
	}
	return indexes
}

// At returns the i-th non-zero diagonal.
// Method accepts negative values with the equivalency -i = n - i.
func (m Diagonals[T]) At(i, slots int) ([]T, error) {

	v, ok := m[i]

	if !ok {

		var j int
		if i > 0 {
			j = i - slots
		} else if j < 0 {
			j = i + slots
		} else {
			return nil, fmt.Errorf("cannot At[0]: diagonal does not exist")
		}

		v, ok := m[j]

		if !ok {
			return nil, fmt.Errorf("cannot At[%d or %d]: diagonal does not exist", i, j)
		}

		return v, nil
	}

	return v, nil
}

// LinearTransformation is a type for linear transformations on ciphertexts.
// It stores a plaintext matrix in diagonal form and
// can be evaluated on a ciphertext by using the evaluator.LinearTransformation method.
type LinearTransformation struct {
	*rlwe.MetaData
	LogBSGSRatio int
	N1           int                 // N1 is the number of inner loops of the baby-step giant-step algorithm used in the evaluation (if N1 == 0, BSGS is not used).
	Level        int                 // Level is the level at which the matrix is encoded (can be circuit dependent)
	Vec          map[int]ringqp.Poly // Vec is the matrix, in diagonal form, where each entry of vec is an indexed non-zero diagonal.
}

// GaloisElements returns the list of Galois elements needed for the evaluation of the linear transformation.
func (LT LinearTransformation) GaloisElements(params rlwe.ParameterProvider) (galEls []uint64) {
	return galoisElementsForLinearTransformation(params, utils.GetKeys(LT.Vec), LT.LogDimensions.Cols, LT.LogBSGSRatio)
}

// GaloisElementsForLinearTransformation returns the list of Galois elements required to evaluate the linear transformation.
func GaloisElementsForLinearTransformation(params rlwe.ParameterProvider, lt LinearTransformationParameters) (galEls []uint64) {
	return galoisElementsForLinearTransformation(params, lt.DiagonalsIndexList, 1<<lt.LogDimensions.Cols, lt.LogBabyStepGianStepRatio)
}

func galoisElementsForLinearTransformation(params rlwe.ParameterProvider, diags []int, slots, logbsgs int) (galEls []uint64) {

	p := params.GetRLWEParameters()

	if logbsgs < 0 {

		_, _, rotN2 := BSGSIndex(diags, slots, slots)

		galEls = make([]uint64, len(rotN2))
		for i := range rotN2 {
			galEls[i] = p.GaloisElement(rotN2[i])
		}

		return
	}

	N1 := FindBestBSGSRatio(diags, slots, logbsgs)

	_, rotN1, rotN2 := BSGSIndex(diags, slots, N1)

	return p.GaloisElements(utils.GetDistincts(append(rotN1, rotN2...)))
}

// NewLinearTransformation allocates a new LinearTransformation with zero values according to the parameters specified by the LinearTranfromationParameters.
func NewLinearTransformation(params rlwe.ParameterProvider, lt LinearTransformationParameters) LinearTransformation {

	p := params.GetRLWEParameters()

	vec := make(map[int]ringqp.Poly)
	cols := 1 << lt.LogDimensions.Cols
	logBSGS := lt.LogBabyStepGianStepRatio
	levelQ := lt.Level
	levelP := p.MaxLevelP()
	ringQP := p.RingQP().AtLevel(levelQ, levelP)

	diagslislt := lt.DiagonalsIndexList

	var N1 int
	if logBSGS < 0 {
		N1 = 0
		for _, i := range diagslislt {
			idx := i
			if idx < 0 {
				idx += cols
			}
			vec[idx] = ringQP.NewPoly()
		}
	} else {
		N1 = FindBestBSGSRatio(diagslislt, cols, logBSGS)
		index, _, _ := BSGSIndex(diagslislt, cols, N1)
		for j := range index {
			for _, i := range index[j] {
				vec[j+i] = ringQP.NewPoly()
			}
		}
	}

	metadata := &rlwe.MetaData{
		PlaintextMetaData: rlwe.PlaintextMetaData{
			LogDimensions: lt.LogDimensions,
			Scale:         lt.Scale,
			IsBatched:     true,
		},
		CiphertextMetaData: rlwe.CiphertextMetaData{
			IsNTT:        true,
			IsMontgomery: true,
		},
	}

	return LinearTransformation{MetaData: metadata, LogBSGSRatio: logBSGS, N1: N1, Level: levelQ, Vec: vec}
}

// EncodeLinearTransformation encodes on a pre-allocated LinearTransformation a set of non-zero diagonaes of a matrix representing a linear transformation.
//
// inputs:
//   - allocated: a pre-allocated LinearTransformation using `NewLinearTransformation`
//   - diagonals: linear transformation parameters
//   - encoder: an struct complying to the EncoderInterface
func EncodeLinearTransformation[T Numeric](params LinearTransformationParameters, encoder EncoderInterface[T, ringqp.Poly], diagonals Diagonals[T], allocated LinearTransformation) (err error) {

	if allocated.LogDimensions != params.LogDimensions {
		return fmt.Errorf("cannot EncodeLinearTransformation: LogDimensions between allocated and parameters do not match (%v != %v)", allocated.LogDimensions, params.LogDimensions)
	}

	rows := 1 << params.LogDimensions.Rows
	cols := 1 << params.LogDimensions.Cols
	N1 := allocated.N1

	diags := params.DiagonalsIndexList

	buf := make([]T, rows*cols)

	metaData := allocated.MetaData

	metaData.Scale = params.Scale

	var v []T

	if N1 == 0 {
		for _, i := range diags {

			idx := i
			if idx < 0 {
				idx += cols
			}

			if vec, ok := allocated.Vec[idx]; !ok {
				return fmt.Errorf("cannot EncodeLinearTransformation: error encoding on LinearTransformation: plaintext diagonal [%d] does not exist", idx)
			} else {

				if v, err = diagonals.At(i, cols); err != nil {
					return fmt.Errorf("cannot EncodeLinearTransformation: %w", err)
				}

				if err = rotateAndEncodeDiagonal(v, encoder, 0, metaData, buf, vec); err != nil {
					return
				}
			}
		}
	} else {

		index, _, _ := BSGSIndex(diags, cols, N1)

		for j := range index {

			rot := -j & (cols - 1)

			for _, i := range index[j] {

				if vec, ok := allocated.Vec[i+j]; !ok {
					return fmt.Errorf("cannot Encode: error encoding on LinearTransformation BSGS: input does not match the same non-zero diagonals")
				} else {

					if v, err = diagonals.At(i+j, cols); err != nil {
						return fmt.Errorf("cannot EncodeLinearTransformation: %w", err)
					}

					if err = rotateAndEncodeDiagonal(v, encoder, rot, metaData, buf, vec); err != nil {
						return
					}
				}
			}
		}
	}

	return
}

func rotateAndEncodeDiagonal[T Numeric](v []T, encoder EncoderInterface[T, ringqp.Poly], rot int, metaData *rlwe.MetaData, buf []T, poly ringqp.Poly) (err error) {

	rows := 1 << metaData.LogDimensions.Rows
	cols := 1 << metaData.LogDimensions.Cols

	rot &= (cols - 1)

	var values []T
	if rot != 0 {

		values = buf

		for i := 0; i < rows; i++ {
			utils.RotateSliceAllocFree(v[i*cols:(i+1)*cols], rot, values[i*cols:(i+1)*cols])
		}

	} else {
		values = v
	}

	return encoder.Encode(values, metaData, poly)
}

// LinearTransformationsNew takes as input a ciphertext ctIn and a list of linear transformations [M0, M1, M2, ...] and returns opOut:[M0(ctIn), M1(ctIn), M2(ctInt), ...].
func (eval LinearTransformEvaluator) LinearTransformationsNew(ctIn *rlwe.Ciphertext, linearTransformations []LinearTransformation) (opOut []*rlwe.Ciphertext, err error) {

	params := eval.GetRLWEParameters()
	opOut = make([]*rlwe.Ciphertext, len(linearTransformations))
	for i := range opOut {
		opOut[i] = rlwe.NewCiphertext(params, 1, linearTransformations[i].Level)
	}

	return opOut, eval.LinearTransformations(ctIn, linearTransformations, opOut)
}

// LinearTransformationNew takes as input a ciphertext ctIn and a linear transformation M and evaluate and returns opOut: M(ctIn).
func (eval LinearTransformEvaluator) LinearTransformationNew(ctIn *rlwe.Ciphertext, linearTransformation LinearTransformation) (opOut *rlwe.Ciphertext, err error) {
	cts, err := eval.LinearTransformationsNew(ctIn, []LinearTransformation{linearTransformation})
	return cts[0], err
}

// LinearTransformation takes as input a ciphertext ctIn, a linear transformation M and evaluates opOut: M(ctIn).
func (eval LinearTransformEvaluator) LinearTransformation(ctIn *rlwe.Ciphertext, linearTransformation LinearTransformation, opOut *rlwe.Ciphertext) (err error) {
	return eval.LinearTransformations(ctIn, []LinearTransformation{linearTransformation}, []*rlwe.Ciphertext{opOut})
}

// LinearTransformations takes as input a ciphertext ctIn, a list of linear transformations [M0, M1, M2, ...] and a list of pre-allocated receiver opOut
// and evaluates opOut: [M0(ctIn), M1(ctIn), M2(ctIn), ...]
func (eval LinearTransformEvaluator) LinearTransformations(ctIn *rlwe.Ciphertext, linearTransformations []LinearTransformation, opOut []*rlwe.Ciphertext) (err error) {

	params := eval.GetRLWEParameters()

	if len(opOut) < len(linearTransformations) {
		return fmt.Errorf("output *rlwe.Ciphertext slice is too small")
	}
	for i := range linearTransformations {
		if opOut[i] == nil {
			return fmt.Errorf("output slice contains unallocated ciphertext")
		}
	}

	var level int
	for _, lt := range linearTransformations {
		level = utils.Max(level, lt.Level)
	}
	level = utils.Min(level, ctIn.Level())

	eval.DecomposeNTT(level, params.MaxLevelP(), params.PCount(), ctIn.Value[1], ctIn.IsNTT, eval.BuffDecompQP)
	for i, lt := range linearTransformations {
		if lt.N1 == 0 {
			if err = eval.MultiplyByDiagMatrix(ctIn, lt, eval.BuffDecompQP, opOut[i]); err != nil {
				return
			}
		} else {
			if err = eval.MultiplyByDiagMatrixBSGS(ctIn, lt, eval.BuffDecompQP, opOut[i]); err != nil {
				return
			}
		}
	}
	return
}

// MultiplyByDiagMatrix multiplies the Ciphertext "ctIn" by the plaintext matrix "matrix" and returns the result on the Ciphertext
// "opOut". Memory buffers for the decomposed ciphertext BuffDecompQP, BuffDecompQP must be provided, those are list of poly of ringQ and ringP
// respectively, each of size params.Beta().
// The naive approach is used (single hoisting and no baby-step giant-step), which is faster than MultiplyByDiagMatrixBSGS
// for matrix of only a few non-zero diagonals but uses more keys.
func (eval LinearTransformEvaluator) MultiplyByDiagMatrix(ctIn *rlwe.Ciphertext, matrix LinearTransformation, BuffDecompQP []ringqp.Poly, opOut *rlwe.Ciphertext) (err error) {

	*opOut.MetaData = *ctIn.MetaData
	opOut.Scale = opOut.Scale.Mul(matrix.Scale)

	params := eval.GetRLWEParameters()

	levelQ := utils.Min(opOut.Level(), utils.Min(ctIn.Level(), matrix.Level))
	levelP := params.RingP().MaxLevel()

	ringQP := params.RingQP().AtLevel(levelQ, levelP)
	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	opOut.Resize(opOut.Degree(), levelQ)

	QiOverF := params.QiOverflowMargin(levelQ)
	PiOverF := params.PiOverflowMargin(levelP)

	c0OutQP := ringqp.Poly{Q: opOut.Value[0], P: eval.BuffQP[5].Q}
	c1OutQP := ringqp.Poly{Q: opOut.Value[1], P: eval.BuffQP[5].P}

	ct0TimesP := eval.BuffQP[0].Q // ct0 * P mod Q
	tmp0QP := eval.BuffQP[1]
	tmp1QP := eval.BuffQP[2]

	cQP := &rlwe.Operand[ringqp.Poly]{}
	cQP.Value = []ringqp.Poly{eval.BuffQP[3], eval.BuffQP[4]}
	cQP.MetaData = &rlwe.MetaData{}
	cQP.MetaData.IsNTT = true

	ring.Copy(ctIn.Value[0], eval.BuffCt.Value[0])
	ring.Copy(ctIn.Value[1], eval.BuffCt.Value[1])
	ctInTmp0, ctInTmp1 := eval.BuffCt.Value[0], eval.BuffCt.Value[1]

	ringQ.MulScalarBigint(ctInTmp0, ringP.ModulusAtLevel[levelP], ct0TimesP) // P*c0

	slots := 1 << matrix.LogDimensions.Cols

	keys := utils.GetSortedKeys(matrix.Vec)

	var state bool
	if keys[0] == 0 {
		state = true
		keys = keys[1:]
	}

	for i, k := range keys {

		k &= (slots - 1)

		galEl := params.GaloisElement(k)

		var evk *rlwe.GaloisKey
		var err error
		if evk, err = eval.CheckAndGetGaloisKey(galEl); err != nil {
			return fmt.Errorf("cannot MultiplyByDiagMatrix: Automorphism: CheckAndGetGaloisKey: %w", err)
		}

		index := eval.AutomorphismIndex(galEl)

		eval.GadgetProductHoistedLazy(levelQ, BuffDecompQP, &evk.GadgetCiphertext, cQP)
		ringQ.Add(cQP.Value[0].Q, ct0TimesP, cQP.Value[0].Q)
		ringQP.AutomorphismNTTWithIndex(cQP.Value[0], index, tmp0QP)
		ringQP.AutomorphismNTTWithIndex(cQP.Value[1], index, tmp1QP)

		pt := matrix.Vec[k]

		if i == 0 {
			// keyswitch(c1_Q) = (d0_QP, d1_QP)
			ringQP.MulCoeffsMontgomery(pt, tmp0QP, c0OutQP)
			ringQP.MulCoeffsMontgomery(pt, tmp1QP, c1OutQP)
		} else {
			// keyswitch(c1_Q) = (d0_QP, d1_QP)
			ringQP.MulCoeffsMontgomeryThenAdd(pt, tmp0QP, c0OutQP)
			ringQP.MulCoeffsMontgomeryThenAdd(pt, tmp1QP, c1OutQP)
		}

		if i%QiOverF == QiOverF-1 {
			ringQ.Reduce(c0OutQP.Q, c0OutQP.Q)
			ringQ.Reduce(c1OutQP.Q, c1OutQP.Q)
		}

		if i%PiOverF == PiOverF-1 {
			ringP.Reduce(c0OutQP.P, c0OutQP.P)
			ringP.Reduce(c1OutQP.P, c1OutQP.P)
		}
	}

	if len(keys)%QiOverF == 0 {
		ringQ.Reduce(c0OutQP.Q, c0OutQP.Q)
		ringQ.Reduce(c1OutQP.Q, c1OutQP.Q)
	}

	if len(keys)%PiOverF == 0 {
		ringP.Reduce(c0OutQP.P, c0OutQP.P)
		ringP.Reduce(c1OutQP.P, c1OutQP.P)
	}

	eval.ModDownQPtoQNTT(levelQ, levelP, c0OutQP.Q, c0OutQP.P, c0OutQP.Q) // sum(phi(c0 * P + d0_QP))/P
	eval.ModDownQPtoQNTT(levelQ, levelP, c1OutQP.Q, c1OutQP.P, c1OutQP.Q) // sum(phi(d1_QP))/P

	if state { // Rotation by zero
		ringQ.MulCoeffsMontgomeryThenAdd(matrix.Vec[0].Q, ctInTmp0, c0OutQP.Q) // opOut += c0_Q * plaintext
		ringQ.MulCoeffsMontgomeryThenAdd(matrix.Vec[0].Q, ctInTmp1, c1OutQP.Q) // opOut += c1_Q * plaintext
	}

	return
}

// MultiplyByDiagMatrixBSGS multiplies the Ciphertext "ctIn" by the plaintext matrix "matrix" and returns the result on the Ciphertext
// "opOut". Memory buffers for the decomposed Ciphertext BuffDecompQP, BuffDecompQP must be provided, those are list of poly of ringQ and ringP
// respectively, each of size params.Beta().
// The BSGS approach is used (double hoisting with baby-step giant-step), which is faster than MultiplyByDiagMatrix
// for matrix with more than a few non-zero diagonals and uses significantly less keys.
func (eval LinearTransformEvaluator) MultiplyByDiagMatrixBSGS(ctIn *rlwe.Ciphertext, matrix LinearTransformation, BuffDecompQP []ringqp.Poly, opOut *rlwe.Ciphertext) (err error) {

	params := eval.GetRLWEParameters()

	*opOut.MetaData = *ctIn.MetaData
	opOut.Scale = opOut.Scale.Mul(matrix.Scale)

	levelQ := utils.Min(opOut.Level(), utils.Min(ctIn.Level(), matrix.Level))
	levelP := params.MaxLevelP()

	ringQP := params.RingQP().AtLevel(levelQ, levelP)
	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	opOut.Resize(opOut.Degree(), levelQ)

	QiOverF := params.QiOverflowMargin(levelQ) >> 1
	PiOverF := params.PiOverflowMargin(levelP) >> 1

	// Computes the N2 rotations indexes of the non-zero rows of the diagonalized DFT matrix for the baby-step giant-step algorithm
	index, _, rotN2 := BSGSIndex(utils.GetKeys(matrix.Vec), 1<<matrix.LogDimensions.Cols, matrix.N1)

	ring.Copy(ctIn.Value[0], eval.BuffCt.Value[0])
	ring.Copy(ctIn.Value[1], eval.BuffCt.Value[1])

	ctInTmp0, ctInTmp1 := eval.BuffCt.Value[0], eval.BuffCt.Value[1]

	// Pre-rotates ciphertext for the baby-step giant-step algorithm, does not divide by P yet
	ctInRotQP := map[int]*rlwe.Operand[ringqp.Poly]{}
	for _, i := range rotN2 {
		if i != 0 {
			ctInRotQP[i] = rlwe.NewOperandQP(params, 1, levelQ, levelP)
			if err = eval.AutomorphismHoistedLazy(levelQ, ctIn, BuffDecompQP, params.GaloisElement(i), ctInRotQP[i]); err != nil {
				return
			}
		}
	}

	// Accumulator inner loop
	tmp0QP := eval.BuffQP[1]
	tmp1QP := eval.BuffQP[2]

	// Accumulator outer loop
	cQP := &rlwe.Operand[ringqp.Poly]{}
	cQP.Value = []ringqp.Poly{eval.BuffQP[3], eval.BuffQP[4]}
	cQP.MetaData = &rlwe.MetaData{}
	cQP.IsNTT = true

	// Result in QP
	c0OutQP := ringqp.Poly{Q: opOut.Value[0], P: eval.BuffQP[5].Q}
	c1OutQP := ringqp.Poly{Q: opOut.Value[1], P: eval.BuffQP[5].P}

	ringQ.MulScalarBigint(ctInTmp0, ringP.ModulusAtLevel[levelP], ctInTmp0) // P*c0
	ringQ.MulScalarBigint(ctInTmp1, ringP.ModulusAtLevel[levelP], ctInTmp1) // P*c1

	keys := utils.GetSortedKeys(index)

	// OUTER LOOP
	var cnt0 int
	for _, j := range keys {

		// INNER LOOP
		var cnt1 int
		for _, i := range index[j] {

			pt := matrix.Vec[j+i]
			ct := ctInRotQP[i]

			if i == 0 {
				if cnt1 == 0 {
					ringQ.MulCoeffsMontgomeryLazy(pt.Q, ctInTmp0, tmp0QP.Q)
					ringQ.MulCoeffsMontgomeryLazy(pt.Q, ctInTmp1, tmp1QP.Q)
					tmp0QP.P.Zero()
					tmp1QP.P.Zero()
				} else {
					ringQ.MulCoeffsMontgomeryLazyThenAddLazy(pt.Q, ctInTmp0, tmp0QP.Q)
					ringQ.MulCoeffsMontgomeryLazyThenAddLazy(pt.Q, ctInTmp1, tmp1QP.Q)
				}
			} else {
				if cnt1 == 0 {
					ringQP.MulCoeffsMontgomeryLazy(pt, ct.Value[0], tmp0QP)
					ringQP.MulCoeffsMontgomeryLazy(pt, ct.Value[1], tmp1QP)
				} else {
					ringQP.MulCoeffsMontgomeryLazyThenAddLazy(pt, ct.Value[0], tmp0QP)
					ringQP.MulCoeffsMontgomeryLazyThenAddLazy(pt, ct.Value[1], tmp1QP)
				}
			}

			if cnt1%QiOverF == QiOverF-1 {
				ringQ.Reduce(tmp0QP.Q, tmp0QP.Q)
				ringQ.Reduce(tmp1QP.Q, tmp1QP.Q)
			}

			if cnt1%PiOverF == PiOverF-1 {
				ringP.Reduce(tmp0QP.P, tmp0QP.P)
				ringP.Reduce(tmp1QP.P, tmp1QP.P)
			}

			cnt1++
		}

		if cnt1%QiOverF != 0 {
			ringQ.Reduce(tmp0QP.Q, tmp0QP.Q)
			ringQ.Reduce(tmp1QP.Q, tmp1QP.Q)
		}

		if cnt1%PiOverF != 0 {
			ringP.Reduce(tmp0QP.P, tmp0QP.P)
			ringP.Reduce(tmp1QP.P, tmp1QP.P)
		}

		// If j != 0, then rotates ((tmp0QP.Q, tmp0QP.P), (tmp1QP.Q, tmp1QP.P)) by N1*j and adds the result on ((cQP.Value[0].Q, cQP.Value[0].P), (cQP.Value[1].Q, cQP.Value[1].P))
		if j != 0 {

			// Hoisting of the ModDown of sum(sum(phi(d1) * plaintext))
			eval.ModDownQPtoQNTT(levelQ, levelP, tmp1QP.Q, tmp1QP.P, tmp1QP.Q) // c1 * plaintext + sum(phi(d1) * plaintext) + phi(c1) * plaintext mod Q

			galEl := params.GaloisElement(j)

			var evk *rlwe.GaloisKey
			var err error
			if evk, err = eval.CheckAndGetGaloisKey(galEl); err != nil {
				return fmt.Errorf("cannot MultiplyByDiagMatrix: Automorphism: CheckAndGetGaloisKey: %w", err)
			}

			rotIndex := eval.AutomorphismIndex(galEl)

			eval.GadgetProductLazy(levelQ, tmp1QP.Q, &evk.GadgetCiphertext, cQP) // EvaluationKey(P*phi(tmpRes_1)) = (d0, d1) in base QP
			ringQP.Add(cQP.Value[0], tmp0QP, cQP.Value[0])

			// Outer loop rotations
			if cnt0 == 0 {
				ringQP.AutomorphismNTTWithIndex(cQP.Value[0], rotIndex, c0OutQP)
				ringQP.AutomorphismNTTWithIndex(cQP.Value[1], rotIndex, c1OutQP)
			} else {
				ringQP.AutomorphismNTTWithIndexThenAddLazy(cQP.Value[0], rotIndex, c0OutQP)
				ringQP.AutomorphismNTTWithIndexThenAddLazy(cQP.Value[1], rotIndex, c1OutQP)
			}

			// Else directly adds on ((cQP.Value[0].Q, cQP.Value[0].P), (cQP.Value[1].Q, cQP.Value[1].P))
		} else {
			if cnt0 == 0 {
				ringqp.CopyLvl(levelQ, levelP, tmp0QP, c0OutQP)
				ringqp.CopyLvl(levelQ, levelP, tmp1QP, c1OutQP)
			} else {
				ringQP.AddLazy(c0OutQP, tmp0QP, c0OutQP)
				ringQP.AddLazy(c1OutQP, tmp1QP, c1OutQP)
			}
		}

		if cnt0%QiOverF == QiOverF-1 {
			ringQ.Reduce(opOut.Value[0], opOut.Value[0])
			ringQ.Reduce(opOut.Value[1], opOut.Value[1])
		}

		if cnt0%PiOverF == PiOverF-1 {
			ringP.Reduce(c0OutQP.P, c0OutQP.P)
			ringP.Reduce(c1OutQP.P, c1OutQP.P)
		}

		cnt0++
	}

	if cnt0%QiOverF != 0 {
		ringQ.Reduce(opOut.Value[0], opOut.Value[0])
		ringQ.Reduce(opOut.Value[1], opOut.Value[1])
	}

	if cnt0%PiOverF != 0 {
		ringP.Reduce(c0OutQP.P, c0OutQP.P)
		ringP.Reduce(c1OutQP.P, c1OutQP.P)
	}

	eval.ModDownQPtoQNTT(levelQ, levelP, opOut.Value[0], c0OutQP.P, opOut.Value[0]) // sum(phi(c0 * P + d0_QP))/P
	eval.ModDownQPtoQNTT(levelQ, levelP, opOut.Value[1], c1OutQP.P, opOut.Value[1]) // sum(phi(d1_QP))/P

	return
}

// FindBestBSGSRatio finds the best N1*N2 = N for the baby-step giant-step algorithm for matrix multiplication.
func FindBestBSGSRatio(nonZeroDiags []int, maxN int, logMaxRatio int) (minN int) {

	maxRatio := float64(int(1 << logMaxRatio))

	for N1 := 1; N1 < maxN; N1 <<= 1 {

		_, rotN1, rotN2 := BSGSIndex(nonZeroDiags, maxN, N1)

		nbN1, nbN2 := len(rotN1)-1, len(rotN2)-1

		if float64(nbN2)/float64(nbN1) == maxRatio {
			return N1
		}

		if float64(nbN2)/float64(nbN1) > maxRatio {
			return N1 / 2
		}
	}

	return 1
}

// BSGSIndex returns the index map and needed rotation for the BSGS matrix-vector multiplication algorithm.
func BSGSIndex(nonZeroDiags []int, slots, N1 int) (index map[int][]int, rotN1, rotN2 []int) {
	index = make(map[int][]int)
	rotN1Map := make(map[int]bool)
	rotN2Map := make(map[int]bool)

	for _, rot := range nonZeroDiags {
		rot &= (slots - 1)
		idxN1 := ((rot / N1) * N1) & (slots - 1)
		idxN2 := rot & (N1 - 1)
		if index[idxN1] == nil {
			index[idxN1] = []int{idxN2}
		} else {
			index[idxN1] = append(index[idxN1], idxN2)
		}
		rotN1Map[idxN1] = true
		rotN2Map[idxN2] = true
	}

	for k := range index {
		sort.Ints(index[k])
	}

	return index, utils.GetSortedKeys(rotN1Map), utils.GetSortedKeys(rotN2Map)
}
