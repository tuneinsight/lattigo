package rlwe

import (
	"fmt"
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// LinearTranfromationParameters is an interface defining a set of methods
// for structs representing and parameterizing a linear transformation.
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
// - Storage: #diagonals polynomials mod Q_level * P
// - Evaluation: #diagonals multiplications and 2sqrt(#diagonals) ciphertexts rotations.
type LinearTranfromationParameters[T any] interface {

	// DiagonalsList returns the list of the non-zero diagonals of the square matrix.
	// A non zero diagonals is a diagonal with a least one non-zero element.
	GetDiagonalsList() []int

	// Diagonals returns all non-zero diagonals of the square matrix in a map indexed
	// by their position.
	GetDiagonals() map[int][]T

	// At returns the i-th non-zero diagonal.
	// Method must accept negative values with the equivalency -i = n - i.
	At(i int) ([]T, error)

	// Level returns level at which to encode the linear transformation.
	GetLevel() int

	// PlaintextScale returns the plaintext scale at which to encode the linear transformation.
	GetPlaintextScale() Scale

	// PlaintextLogDimensions returns log2 dimensions of the matrix that can be SIMD packed
	// in a single plaintext polynomial.
	// This method is equivalent to params.PlaintextDimensions().
	// Note that the linear transformation is evaluated independently on each rows of
	// the SIMD packed matrix.
	GetPlaintextLogDimensions() [2]int

	// LogBabyStepGianStepRatio return the log2 of the ratio n1/n2 for n = n1 * n2 and
	// n is the dimension of the linear transformation. The number of Galois keys required
	// is minimized when this value is 0 but the overall complexity of the homomorphic evaluation
	// can be reduced by increasing the ratio (at the expanse of increasing the number of keys required).
	// If the value returned is negative, then the baby-step giant-step algorithm is not used
	// and the evaluation complexity (as well as the number of keys) becomes O(n) instead of O(sqrt(n)).
	GetLogBabyStepGianStepRatio() int
}

type MemLinearTransformationParameters[T any] struct {
	Diagonals                map[int][]T
	Level                    int
	PlaintextScale           Scale
	PlaintextLogDimensions   [2]int
	LogBabyStepGianStepRatio int
}

func (m MemLinearTransformationParameters[T]) GetDiagonalsList() []int {
	return utils.GetKeys(m.Diagonals)
}

func (m MemLinearTransformationParameters[T]) GetDiagonals() map[int][]T {
	return m.Diagonals
}

func (m MemLinearTransformationParameters[T]) At(i int) ([]T, error) {

	slots := 1 << m.PlaintextLogDimensions[1]

	v, ok := m.Diagonals[i]

	if !ok {

		var j int
		if i > 0 {
			j = i - slots
		} else if j < 0 {
			j = i + slots
		} else {
			return nil, fmt.Errorf("cannot At[0]: diagonal does not exist")
		}

		v, ok := m.Diagonals[j]

		if !ok {
			return nil, fmt.Errorf("cannot At[%d or %d]: diagonal does not exist", i, j)
		}

		return v, nil
	}

	return v, nil
}

func (m MemLinearTransformationParameters[T]) GetLevel() int {
	return m.Level
}

func (m MemLinearTransformationParameters[T]) GetPlaintextScale() Scale {
	return m.PlaintextScale
}

func (m MemLinearTransformationParameters[T]) GetPlaintextLogDimensions() [2]int {
	return m.PlaintextLogDimensions
}

func (m MemLinearTransformationParameters[T]) GetLogBabyStepGianStepRatio() int {
	return m.LogBabyStepGianStepRatio
}

// LinearTransformation is a type for linear transformations on ciphertexts.
// It stores a plaintext matrix in diagonal form and
// can be evaluated on a ciphertext by using the evaluator.LinearTransformation method.
type LinearTransformation struct {
	*MetaData
	LogBSGSRatio int
	N1           int                 // N1 is the number of inner loops of the baby-step giant-step algorithm used in the evaluation (if N1 == 0, BSGS is not used).
	Level        int                 // Level is the level at which the matrix is encoded (can be circuit dependent)
	Vec          map[int]ringqp.Poly // Vec is the matrix, in diagonal form, where each entry of vec is an indexed non-zero diagonal.
}

// GaloisElements returns the list of Galois elements needed for the evaluation of the linear transformation.
func (LT LinearTransformation) GaloisElements(params ParametersInterface) (galEls []uint64) {
	return galoisElementsForLinearTransformation(params, utils.GetKeys(LT.Vec), LT.PlaintextLogDimensions[1], LT.LogBSGSRatio)
}

// GaloisElementsForLinearTransformation returns the list of Galois elements required to perform a linear transform
// with the provided non-zero diagonals.
func GaloisElementsForLinearTransformation[T any](params ParametersInterface, lt LinearTranfromationParameters[T]) (galEls []uint64) {
	return galoisElementsForLinearTransformation(params, lt.GetDiagonalsList(), 1<<lt.GetPlaintextLogDimensions()[1], lt.GetLogBabyStepGianStepRatio())
}

func galoisElementsForLinearTransformation(params ParametersInterface, diags []int, slots, logbsgs int) (galEls []uint64) {

	if logbsgs < 0 {

		_, _, rotN2 := BSGSIndex(diags, slots, slots)

		galEls = make([]uint64, len(rotN2))

		for i := range rotN2 {
			galEls[i] = params.GaloisElement(rotN2[i])
		}

		return
	}

	N1 := FindBestBSGSRatio(diags, slots, logbsgs)

	_, rotN1, rotN2 := BSGSIndex(diags, slots, N1)

	return params.GaloisElements(utils.GetDistincts(append(rotN1, rotN2...)))
}

// NewLinearTransformation allocates a new LinearTransformation with zero values according to the parameters specified by the LinearTranfromationParameters.
func NewLinearTransformation[T any](params ParametersInterface, lt LinearTranfromationParameters[T]) LinearTransformation {
	vec := make(map[int]ringqp.Poly)
	cols := 1 << lt.GetPlaintextLogDimensions()[1]
	logBSGS := lt.GetLogBabyStepGianStepRatio()
	levelQ := lt.GetLevel()
	levelP := params.MaxLevelP()
	ringQP := params.RingQP().AtLevel(levelQ, levelP)

	diagslislt := lt.GetDiagonalsList()

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

	metadata := &MetaData{
		PlaintextLogDimensions: lt.GetPlaintextLogDimensions(),
		PlaintextScale:         lt.GetPlaintextScale(),
		EncodingDomain:         FrequencyDomain,
		IsNTT:                  true,
		IsMontgomery:           true,
	}

	return LinearTransformation{MetaData: metadata, LogBSGSRatio: logBSGS, N1: N1, Level: levelQ, Vec: vec}
}

// EncodeLinearTransformation encodes on a pre-allocated LinearTransformation a set of non-zero diagonaes of a matrix representing a linear transformation.
//
// inputs:
// - allocated: a pre-allocated LinearTransformation using `NewLinearTransformation`
// - diagonals: linear transformation parameters
// - encoder: an struct complying to the EncoderInterface
func EncodeLinearTransformation[T any](allocated LinearTransformation, params LinearTranfromationParameters[T], encoder EncoderInterface[T, ringqp.Poly]) (err error) {

	if allocated.PlaintextLogDimensions != params.GetPlaintextLogDimensions() {
		return fmt.Errorf("cannot EncodeLinearTransformation: PlaintextLogDimensions between allocated and parameters do not match (%v != %v)", allocated.PlaintextLogDimensions, params.GetPlaintextLogDimensions())
	}

	rows := 1 << params.GetPlaintextLogDimensions()[0]
	cols := 1 << params.GetPlaintextLogDimensions()[1]
	N1 := allocated.N1

	diags := params.GetDiagonalsList()

	buf := make([]T, rows*cols)

	metaData := allocated.MetaData

	metaData.PlaintextScale = params.GetPlaintextScale()

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

				if v, err = params.At(i); err != nil {
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

					if v, err = params.At(i + j); err != nil {
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

func rotateAndEncodeDiagonal[T any](v []T, encoder EncoderInterface[T, ringqp.Poly], rot int, metaData *MetaData, buf []T, poly ringqp.Poly) (err error) {

	rows := 1 << metaData.PlaintextLogDimensions[0]
	cols := 1 << metaData.PlaintextLogDimensions[1]

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

// LinearTransformationNew evaluates a linear transform on the pre-allocated Ciphertexts.
// The LinearTransformation can either be an (ordered) list of LinearTransformation or a single LinearTransformation.
// In either case a list of Ciphertext is returned (the second case returning a list containing a single Ciphertext).
func (eval Evaluator) LinearTransformationNew(ctIn *Ciphertext, linearTransformation interface{}) (opOut []*Ciphertext, err error) {

	switch LTs := linearTransformation.(type) {
	case []LinearTransformation:
		opOut = make([]*Ciphertext, len(LTs))

		var maxLevel int
		for _, LT := range LTs {
			maxLevel = utils.Max(maxLevel, LT.Level)
		}

		minLevel := utils.Min(maxLevel, ctIn.Level())
		eval.DecomposeNTT(minLevel, eval.params.MaxLevelP(), eval.params.PCount(), ctIn.Value[1], ctIn.IsNTT, eval.BuffDecompQP)

		for i, LT := range LTs {
			opOut[i] = NewCiphertext(eval.params, 1, minLevel)

			if LT.N1 == 0 {
				if err = eval.MultiplyByDiagMatrix(ctIn, LT, eval.BuffDecompQP, opOut[i]); err != nil {
					return
				}
			} else {
				if err = eval.MultiplyByDiagMatrixBSGS(ctIn, LT, eval.BuffDecompQP, opOut[i]); err != nil {
					return
				}
			}
		}

	case LinearTransformation:

		minLevel := utils.Min(LTs.Level, ctIn.Level())
		eval.DecomposeNTT(minLevel, eval.params.MaxLevelP(), eval.params.PCount(), ctIn.Value[1], ctIn.IsNTT, eval.BuffDecompQP)

		opOut = []*Ciphertext{NewCiphertext(eval.params, 1, minLevel)}

		if LTs.N1 == 0 {
			if err = eval.MultiplyByDiagMatrix(ctIn, LTs, eval.BuffDecompQP, opOut[0]); err != nil {
				return
			}
		} else {
			if err = eval.MultiplyByDiagMatrixBSGS(ctIn, LTs, eval.BuffDecompQP, opOut[0]); err != nil {
				return
			}
		}
	}
	return
}

// LinearTransformation evaluates a linear transform on the pre-allocated Ciphertexts.
// The LinearTransformation can either be an (ordered) list of LinearTransformation or a single LinearTransformation.
// In either case a list of Ciphertext is returned (the second case returning a list containing a single Ciphertext).
func (eval Evaluator) LinearTransformation(ctIn *Ciphertext, linearTransformation interface{}, opOut []*Ciphertext) (err error) {

	switch LTs := linearTransformation.(type) {
	case []LinearTransformation:
		var maxLevel int
		for _, LT := range LTs {
			maxLevel = utils.Max(maxLevel, LT.Level)
		}

		minLevel := utils.Min(maxLevel, ctIn.Level())
		eval.DecomposeNTT(minLevel, eval.params.MaxLevelP(), eval.params.PCount(), ctIn.Value[1], true, eval.BuffDecompQP)

		for i, LT := range LTs {
			if LT.N1 == 0 {
				if err = eval.MultiplyByDiagMatrix(ctIn, LT, eval.BuffDecompQP, opOut[i]); err != nil {
					return
				}
			} else {
				if err = eval.MultiplyByDiagMatrixBSGS(ctIn, LT, eval.BuffDecompQP, opOut[i]); err != nil {
					return
				}
			}
		}

	case LinearTransformation:
		minLevel := utils.Min(LTs.Level, ctIn.Level())
		eval.DecomposeNTT(minLevel, eval.params.MaxLevelP(), eval.params.PCount(), ctIn.Value[1], true, eval.BuffDecompQP)
		if LTs.N1 == 0 {
			if err = eval.MultiplyByDiagMatrix(ctIn, LTs, eval.BuffDecompQP, opOut[0]); err != nil {
				return
			}
		} else {
			if err = eval.MultiplyByDiagMatrixBSGS(ctIn, LTs, eval.BuffDecompQP, opOut[0]); err != nil {
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
func (eval Evaluator) MultiplyByDiagMatrix(ctIn *Ciphertext, matrix LinearTransformation, BuffDecompQP []ringqp.Poly, opOut *Ciphertext) (err error) {

	*opOut.MetaData = *ctIn.MetaData
	opOut.PlaintextScale = opOut.PlaintextScale.Mul(matrix.PlaintextScale)

	levelQ := utils.Min(opOut.Level(), utils.Min(ctIn.Level(), matrix.Level))
	levelP := eval.params.RingP().MaxLevel()

	ringQP := eval.params.RingQP().AtLevel(levelQ, levelP)
	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	opOut.Resize(opOut.Degree(), levelQ)

	QiOverF := eval.params.QiOverflowMargin(levelQ)
	PiOverF := eval.params.PiOverflowMargin(levelP)

	c0OutQP := ringqp.Poly{Q: opOut.Value[0], P: eval.BuffQP[5].Q}
	c1OutQP := ringqp.Poly{Q: opOut.Value[1], P: eval.BuffQP[5].P}

	ct0TimesP := eval.BuffQP[0].Q // ct0 * P mod Q
	tmp0QP := eval.BuffQP[1]
	tmp1QP := eval.BuffQP[2]

	cQP := &Operand[ringqp.Poly]{}
	cQP.Value = []ringqp.Poly{eval.BuffQP[3], eval.BuffQP[4]}
	cQP.MetaData = &MetaData{IsNTT: true}

	ring.Copy(ctIn.Value[0], eval.BuffCt.Value[0])
	ring.Copy(ctIn.Value[1], eval.BuffCt.Value[1])
	ctInTmp0, ctInTmp1 := eval.BuffCt.Value[0], eval.BuffCt.Value[1]

	ringQ.MulScalarBigint(ctInTmp0, ringP.ModulusAtLevel[levelP], ct0TimesP) // P*c0

	slots := 1 << matrix.PlaintextLogDimensions[1]

	keys := utils.GetSortedKeys(matrix.Vec)

	var state bool
	if keys[0] == 0 {
		state = true
		keys = keys[1:]
	}

	for i, k := range keys {

		k &= (slots - 1)

		galEl := eval.params.GaloisElement(k)

		var evk *GaloisKey
		var err error
		if evk, err = eval.CheckAndGetGaloisKey(galEl); err != nil {
			return fmt.Errorf("cannot MultiplyByDiagMatrix: Automorphism: CheckAndGetGaloisKey: %w", err)
		}

		index := eval.AutomorphismIndex[galEl]

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

	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c0OutQP.Q, c0OutQP.P, c0OutQP.Q) // sum(phi(c0 * P + d0_QP))/P
	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c1OutQP.Q, c1OutQP.P, c1OutQP.Q) // sum(phi(d1_QP))/P

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
func (eval Evaluator) MultiplyByDiagMatrixBSGS(ctIn *Ciphertext, matrix LinearTransformation, BuffDecompQP []ringqp.Poly, opOut *Ciphertext) (err error) {

	*opOut.MetaData = *ctIn.MetaData
	opOut.PlaintextScale = opOut.PlaintextScale.Mul(matrix.PlaintextScale)

	levelQ := utils.Min(opOut.Level(), utils.Min(ctIn.Level(), matrix.Level))
	levelP := eval.Parameters().MaxLevelP()

	ringQP := eval.params.RingQP().AtLevel(levelQ, levelP)
	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	opOut.Resize(opOut.Degree(), levelQ)

	QiOverF := eval.params.QiOverflowMargin(levelQ) >> 1
	PiOverF := eval.params.PiOverflowMargin(levelP) >> 1

	// Computes the N2 rotations indexes of the non-zero rows of the diagonalized DFT matrix for the baby-step giant-step algorithm
	index, _, rotN2 := BSGSIndex(utils.GetKeys(matrix.Vec), 1<<matrix.PlaintextLogDimensions[1], matrix.N1)

	ring.Copy(ctIn.Value[0], eval.BuffCt.Value[0])
	ring.Copy(ctIn.Value[1], eval.BuffCt.Value[1])

	ctInTmp0, ctInTmp1 := eval.BuffCt.Value[0], eval.BuffCt.Value[1]

	// Pre-rotates ciphertext for the baby-step giant-step algorithm, does not divide by P yet
	ctInRotQP := map[int]*Operand[ringqp.Poly]{}
	for _, i := range rotN2 {
		if i != 0 {
			ctInRotQP[i] = NewOperandQP(eval.Parameters(), 1, levelQ, levelP)
			if err = eval.AutomorphismHoistedLazy(levelQ, ctIn, BuffDecompQP, eval.Parameters().GaloisElement(i), ctInRotQP[i]); err != nil {
				return
			}
		}
	}

	// Accumulator inner loop
	tmp0QP := eval.BuffQP[1]
	tmp1QP := eval.BuffQP[2]

	// Accumulator outer loop
	cQP := &Operand[ringqp.Poly]{}
	cQP.Value = []ringqp.Poly{eval.BuffQP[3], eval.BuffQP[4]}
	cQP.MetaData = &MetaData{IsNTT: true}

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
			eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, tmp1QP.Q, tmp1QP.P, tmp1QP.Q) // c1 * plaintext + sum(phi(d1) * plaintext) + phi(c1) * plaintext mod Q

			galEl := eval.params.GaloisElement(j)

			var evk *GaloisKey
			var err error
			if evk, err = eval.CheckAndGetGaloisKey(galEl); err != nil {
				return fmt.Errorf("cannot MultiplyByDiagMatrix: Automorphism: CheckAndGetGaloisKey: %w", err)
			}

			rotIndex := eval.AutomorphismIndex[galEl]

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

	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, opOut.Value[0], c0OutQP.P, opOut.Value[0]) // sum(phi(c0 * P + d0_QP))/P
	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, opOut.Value[1], c1OutQP.P, opOut.Value[1]) // sum(phi(d1_QP))/P

	return
}

// Trace maps X -> sum((-1)^i * X^{i*n+1}) for n <= i < N
// Monomial X^k vanishes if k is not divisible by (N/n), otherwise it is multiplied by (N/n).
// Ciphertext is pre-multiplied by (N/n)^-1 to remove the (N/n) factor.
// Examples of full Trace for [0 + 1X + 2X^2 + 3X^3 + 4X^4 + 5X^5 + 6X^6 + 7X^7]
//
// 1.
//
//	  [1 + 2X + 3X^2 + 4X^3 + 5X^4 + 6X^5 + 7X^6 + 8X^7]
//	+ [1 - 6X - 3X^2 + 8X^3 + 5X^4 + 2X^5 - 7X^6 - 4X^7]  {X-> X^(i * 5^1)}
//	= [2 - 4X + 0X^2 +12X^3 +10X^4 + 8X^5 - 0X^6 + 4X^7]
//
// 2.
//
//	  [2 - 4X + 0X^2 +12X^3 +10X^4 + 8X^5 - 0X^6 + 4X^7]
//	+ [2 + 4X + 0X^2 -12X^3 +10X^4 - 8X^5 + 0X^6 - 4X^7]  {X-> X^(i * 5^2)}
//	= [4 + 0X + 0X^2 - 0X^3 +20X^4 + 0X^5 + 0X^6 - 0X^7]
//
// 3.
//
//	  [4 + 0X + 0X^2 - 0X^3 +20X^4 + 0X^5 + 0X^6 - 0X^7]
//	+ [4 + 0X + 0X^2 - 0X^3 -20X^4 + 0X^5 + 0X^6 - 0X^7]  {X-> X^(i * -1)}
//	= [8 + 0X + 0X^2 - 0X^3 + 0X^4 + 0X^5 + 0X^6 - 0X^7]
//
// The method will return an error if the input and output ciphertexts degree is not one.
func (eval Evaluator) Trace(ctIn *Ciphertext, logN int, opOut *Ciphertext) (err error) {

	if ctIn.Degree() != 1 || opOut.Degree() != 1 {
		return fmt.Errorf("ctIn.Degree() != 1 or opOut.Degree() != 1")
	}

	level := utils.Min(ctIn.Level(), opOut.Level())

	opOut.Resize(opOut.Degree(), level)

	*opOut.MetaData = *ctIn.MetaData

	gap := 1 << (eval.params.LogN() - logN - 1)

	if logN == 0 {
		gap <<= 1
	}

	if gap > 1 {

		ringQ := eval.params.RingQ().AtLevel(level)

		if ringQ.Type() == ring.ConjugateInvariant {
			gap >>= 1 // We skip the last step that applies phi(5^{-1})
		}

		NInv := new(big.Int).SetUint64(uint64(gap))
		NInv.ModInverse(NInv, ringQ.ModulusAtLevel[level])

		// pre-multiplication by (N/n)^-1
		ringQ.MulScalarBigint(ctIn.Value[0], NInv, opOut.Value[0])
		ringQ.MulScalarBigint(ctIn.Value[1], NInv, opOut.Value[1])

		if !ctIn.IsNTT {
			ringQ.NTT(opOut.Value[0], opOut.Value[0])
			ringQ.NTT(opOut.Value[1], opOut.Value[1])
			opOut.IsNTT = true
		}

		buff, err := NewCiphertextAtLevelFromPoly(level, []ring.Poly{eval.BuffQP[3].Q, eval.BuffQP[4].Q})

		if err != nil {
			panic(err)
		}

		buff.IsNTT = true

		for i := logN; i < eval.params.LogN()-1; i++ {

			if err = eval.Automorphism(opOut, eval.params.GaloisElement(1<<i), buff); err != nil {
				return err
			}

			ringQ.Add(opOut.Value[0], buff.Value[0], opOut.Value[0])
			ringQ.Add(opOut.Value[1], buff.Value[1], opOut.Value[1])
		}

		if logN == 0 && ringQ.Type() == ring.Standard {

			if err = eval.Automorphism(opOut, ringQ.NthRoot()-1, buff); err != nil {
				return err
			}

			ringQ.Add(opOut.Value[0], buff.Value[0], opOut.Value[0])
			ringQ.Add(opOut.Value[1], buff.Value[1], opOut.Value[1])
		}

		if !ctIn.IsNTT {
			ringQ.INTT(opOut.Value[0], opOut.Value[0])
			ringQ.INTT(opOut.Value[1], opOut.Value[1])
			opOut.IsNTT = false
		}

	} else {
		if ctIn != opOut {
			opOut.Copy(ctIn)
		}
	}

	return
}

// Expand expands a RLWE Ciphertext encrypting sum ai * X^i to 2^logN ciphertexts,
// each encrypting ai * X^0 for 0 <= i < 2^LogN. That is, it extracts the first 2^logN
// coefficients, whose degree is a multiple of 2^logGap, of ctIn and returns an RLWE
// Ciphertext for each coefficient extracted.
//
// The method will return an error if:
//   - The input ciphertext degree is not one
//   - The ring type is not ring.Standard
func (eval Evaluator) Expand(ctIn *Ciphertext, logN, logGap int) (opOut []*Ciphertext, err error) {

	if ctIn.Degree() != 1 {
		return nil, fmt.Errorf("cannot Expand: ctIn.Degree() != 1")
	}

	if eval.params.RingType() != ring.Standard {
		return nil, fmt.Errorf("cannot Expand: method is only supported for ring.Type = ring.Standard (X^{-2^{i}} does not exist in the sub-ring Z[X + X^{-1}])")
	}

	params := eval.params

	level := ctIn.Level()

	ringQ := params.RingQ().AtLevel(level)

	// Compute X^{-2^{i}} from 1 to LogN
	xPow2 := genXPow2(ringQ, logN, true)

	opOut = make([]*Ciphertext, 1<<(logN-logGap))
	opOut[0] = ctIn.CopyNew()
	opOut[0].PlaintextLogDimensions = [2]int{0, 0}

	if ct := opOut[0]; !ctIn.IsNTT {
		ringQ.NTT(ct.Value[0], ct.Value[0])
		ringQ.NTT(ct.Value[1], ct.Value[1])
		ct.IsNTT = true
	}

	// Multiplies by 2^{-logN} mod Q
	NInv := new(big.Int).SetUint64(1 << logN)
	NInv.ModInverse(NInv, ringQ.ModulusAtLevel[level])

	ringQ.MulScalarBigint(opOut[0].Value[0], NInv, opOut[0].Value[0])
	ringQ.MulScalarBigint(opOut[0].Value[1], NInv, opOut[0].Value[1])

	gap := 1 << logGap

	tmp, err := NewCiphertextAtLevelFromPoly(level, []ring.Poly{eval.BuffCt.Value[0], eval.BuffCt.Value[1]})

	if err != nil {
		panic(err)
	}

	tmp.MetaData = ctIn.MetaData

	for i := 0; i < logN; i++ {

		n := 1 << i

		galEl := uint64(ringQ.N()/n + 1)

		half := n / gap

		for j := 0; j < (n+gap-1)/gap; j++ {

			c0 := opOut[j]

			// X -> X^{N/n + 1}
			//[a, b, c, d] -> [a, -b, c, -d]
			if err = eval.Automorphism(c0, galEl, tmp); err != nil {
				return
			}

			if j+half > 0 {

				c1 := opOut[j].CopyNew()

				// Zeroes odd coeffs: [a, b, c, d] + [a, -b, c, -d] -> [2a, 0, 2b, 0]
				ringQ.Add(c0.Value[0], tmp.Value[0], c0.Value[0])
				ringQ.Add(c0.Value[1], tmp.Value[1], c0.Value[1])

				// Zeroes even coeffs: [a, b, c, d] - [a, -b, c, -d] -> [0, 2b, 0, 2d]
				ringQ.Sub(c1.Value[0], tmp.Value[0], c1.Value[0])
				ringQ.Sub(c1.Value[1], tmp.Value[1], c1.Value[1])

				// c1 * X^{-2^{i}}: [0, 2b, 0, 2d] * X^{-n} -> [2b, 0, 2d, 0]
				ringQ.MulCoeffsMontgomery(c1.Value[0], xPow2[i], c1.Value[0])
				ringQ.MulCoeffsMontgomery(c1.Value[1], xPow2[i], c1.Value[1])

				opOut[j+half] = c1

			} else {

				// Zeroes odd coeffs: [a, b, c, d] + [a, -b, c, -d] -> [2a, 0, 2b, 0]
				ringQ.Add(c0.Value[0], tmp.Value[0], c0.Value[0])
				ringQ.Add(c0.Value[1], tmp.Value[1], c0.Value[1])
			}
		}
	}

	for _, ct := range opOut {
		if ct != nil && !ctIn.IsNTT {
			ringQ.INTT(ct.Value[0], ct.Value[0])
			ringQ.INTT(ct.Value[1], ct.Value[1])
			ct.IsNTT = false
		}
	}
	return
}

// Pack packs a batch of RLWE ciphertexts, packing the batch of ciphertexts into a single ciphertext.
// The number of key-switching operations is inputLogGap - log2(gap) + len(cts), where log2(gap) is the
// minimum distance between two keys of the map cts[int]*Ciphertext.
//
// Input:
//
//		cts: a map of Ciphertext, where the index in the map is the future position of the first coefficient
//		      of the indexed ciphertext in the final ciphertext (see example). Ciphertexts can be in or out of the NTT domain.
//		logGap: all coefficients of the input ciphertexts that are not a multiple of X^{2^{logGap}} will be zeroed
//		        during the merging (see example). This is equivalent to skipping the first 2^{logGap} steps of the
//		        algorithm, i.e. having as input ciphertexts that are already partially packed.
//	 zeroGarbageSlots: if set to true, slots which are not multiples of X^{2^{logGap}} will be zeroed during the procedure.
//	                   this will greatly increase the noise and increase the number of key-switching operations to inputLogGap + len(cts).
//
// Output: a ciphertext packing all input ciphertexts
//
// Example: we want to pack 4 ciphertexts into one, and keep only coefficients which are a multiple of X^{4}.
//
//	To do so, we must set logGap = 2.
//	Here the `X` slots are treated as garbage slots that we want to discard during the procedure.
//
//	input: map[int]{
//	   0: [x00, X, X, X, x01, X, X, X],   with logGap = 2
//	   1: [x10, X, X, X, x11, X, X, X],
//	   2: [x20, X, X, X, x21, X, X, X],
//	   3: [x30, X, X, X, x31, X, X, X],
//		}
//
//	 Step 1:
//	         map[0]: 2^{-1} * (map[0] + X^2 * map[2] + phi_{5^2}(map[0] - X^2 * map[2]) = [x00, X, x20, X, x01, X, x21, X]
//	         map[1]: 2^{-1} * (map[1] + X^2 * map[3] + phi_{5^2}(map[1] - X^2 * map[3]) = [x10, X, x30, X, x11, X, x31, X]
//	 Step 2:
//	         map[0]: 2^{-1} * (map[0] + X^1 * map[1] + phi_{5^4}(map[0] - X^1 * map[1]) = [x00, x10, x20, x30, x01, x11, x21, x22]
func (eval Evaluator) Pack(cts map[int]*Ciphertext, inputLogGap int, zeroGarbageSlots bool) (ct *Ciphertext, err error) {

	params := eval.Parameters()

	if params.RingType() != ring.Standard {
		return nil, fmt.Errorf("cannot Pack: procedure is only supported for ring.Type = ring.Standard (X^{2^{i}} does not exist in the sub-ring Z[X + X^{-1}])")
	}

	if len(cts) < 2 {
		return nil, fmt.Errorf("cannot Pack: #cts must be at least 2")
	}

	keys := utils.GetSortedKeys(cts)

	gap := keys[1] - keys[0]
	level := cts[keys[0]].Level()

	for i, key := range keys[1:] {
		level = utils.Min(level, cts[key].Level())

		if i < len(keys)-1 {
			gap = utils.Min(gap, keys[i+1]-keys[i])
		}
	}

	logN := params.LogN()
	ringQ := params.RingQ().AtLevel(level)

	logStart := logN - inputLogGap
	logEnd := logN

	if !zeroGarbageSlots {
		if gap > 0 {
			logEnd -= bits.Len64(uint64(gap - 1))
		}
	}

	if logStart >= logEnd {
		return nil, fmt.Errorf("cannot Pack: gaps between ciphertexts is smaller than inputLogGap > N")
	}

	xPow2 := genXPow2(ringQ.AtLevel(level), params.LogN(), false) // log(N) polynomial to generate, quick

	NInv := new(big.Int).SetUint64(uint64(1 << (logEnd - logStart)))
	NInv.ModInverse(NInv, ringQ.ModulusAtLevel[level])

	for _, key := range keys {

		ct := cts[key]

		if ct.Degree() != 1 {
			return nil, fmt.Errorf("cannot Pack: cts[%d].Degree() != 1", key)
		}

		if !ct.IsNTT {
			ringQ.NTT(ct.Value[0], ct.Value[0])
			ringQ.NTT(ct.Value[1], ct.Value[1])
			ct.IsNTT = true
		}

		ringQ.MulScalarBigint(ct.Value[0], NInv, ct.Value[0])
		ringQ.MulScalarBigint(ct.Value[1], NInv, ct.Value[1])
	}

	tmpa := &Ciphertext{}
	tmpa.Value = []ring.Poly{ringQ.NewPoly(), ringQ.NewPoly()}
	tmpa.MetaData = &MetaData{IsNTT: true}

	for i := logStart; i < logEnd; i++ {

		t := 1 << (logN - 1 - i)

		for jx, jy := 0, t; jx < t; jx, jy = jx+1, jy+1 {

			a := cts[jx]
			b := cts[jy]

			if b != nil {

				//X^(N/2^L)
				ringQ.MulCoeffsMontgomery(b.Value[0], xPow2[len(xPow2)-i-1], b.Value[0])
				ringQ.MulCoeffsMontgomery(b.Value[1], xPow2[len(xPow2)-i-1], b.Value[1])

				if a != nil {

					// tmpa = phi(a - b * X^{N/2^{i}}, 2^{i-1})
					ringQ.Sub(a.Value[0], b.Value[0], tmpa.Value[0])
					ringQ.Sub(a.Value[1], b.Value[1], tmpa.Value[1])

					// a = a + b * X^{N/2^{i}}
					ringQ.Add(a.Value[0], b.Value[0], a.Value[0])
					ringQ.Add(a.Value[1], b.Value[1], a.Value[1])

				} else {
					// if ct[jx] == nil, then simply re-assigns
					cts[jx] = cts[jy]
				}
			}

			if a != nil {

				var galEl uint64

				if i == 0 {
					galEl = ringQ.NthRoot() - 1
				} else {
					galEl = eval.Parameters().GaloisElement(1 << (i - 1))
				}

				if b != nil {
					if err = eval.Automorphism(tmpa, galEl, tmpa); err != nil {
						return
					}
				} else {
					if err = eval.Automorphism(a, galEl, tmpa); err != nil {
						return
					}
				}

				// a + b * X^{N/2^{i}} + phi(a - b * X^{N/2^{i}}, 2^{i-1})
				ringQ.Add(a.Value[0], tmpa.Value[0], a.Value[0])
				ringQ.Add(a.Value[1], tmpa.Value[1], a.Value[1])
			}
		}
	}

	return cts[0], nil
}

func genXPow2(r *ring.Ring, logN int, div bool) (xPow []ring.Poly) {

	// Compute X^{-n} from 0 to LogN
	xPow = make([]ring.Poly, logN)

	moduli := r.ModuliChain()[:r.Level()+1]
	BRC := r.BRedConstants()

	var idx int
	for i := 0; i < logN; i++ {

		idx = 1 << i

		if div {
			idx = r.N() - idx
		}

		xPow[i] = r.NewPoly()

		if i == 0 {

			for j := range moduli {
				xPow[i].Coeffs[j][idx] = ring.MForm(1, moduli[j], BRC[j])
			}

			r.NTT(xPow[i], xPow[i])

		} else {
			r.MulCoeffsMontgomery(xPow[i-1], xPow[i-1], xPow[i]) // X^{n} = X^{1} * X^{n-1}
		}
	}

	if div {
		r.Neg(xPow[0], xPow[0])
	}

	return
}

// InnerSum applies an optimized inner sum on the Ciphertext (log2(n) + HW(n) rotations with double hoisting).
// The operation assumes that `ctIn` encrypts SlotCount/`batchSize` sub-vectors of size `batchSize` which it adds together (in parallel) in groups of `n`.
// It outputs in opOut a Ciphertext for which the "leftmost" sub-vector of each group is equal to the sum of the group.
func (eval Evaluator) InnerSum(ctIn *Ciphertext, batchSize, n int, opOut *Ciphertext) (err error) {

	levelQ := ctIn.Level()
	levelP := eval.params.PCount() - 1

	ringQP := eval.params.RingQP().AtLevel(ctIn.Level(), levelP)

	ringQ := ringQP.RingQ

	opOut.Resize(opOut.Degree(), levelQ)
	*opOut.MetaData = *ctIn.MetaData

	ctInNTT, err := NewCiphertextAtLevelFromPoly(levelQ, eval.BuffCt.Value[:2])

	if err != nil {
		panic(err)
	}

	ctInNTT.MetaData = &MetaData{IsNTT: true}

	if !ctIn.IsNTT {
		ringQ.NTT(ctIn.Value[0], ctInNTT.Value[0])
		ringQ.NTT(ctIn.Value[1], ctInNTT.Value[1])
	} else {
		ring.CopyLvl(levelQ, ctIn.Value[0], ctInNTT.Value[0])
		ring.CopyLvl(levelQ, ctIn.Value[1], ctInNTT.Value[1])
	}

	if n == 1 {
		if ctIn != opOut {
			ring.CopyLvl(levelQ, ctIn.Value[0], opOut.Value[0])
			ring.CopyLvl(levelQ, ctIn.Value[1], opOut.Value[1])
		}
	} else {

		// BuffQP[0:2] are used by AutomorphismHoistedLazy

		// Accumulator mod QP (i.e. opOut Mod QP)
		accQP := &Operand[ringqp.Poly]{Value: []ringqp.Poly{eval.BuffQP[2], eval.BuffQP[3]}}
		accQP.MetaData = ctInNTT.MetaData

		// Buffer mod QP (i.e. to store the result of lazy gadget products)
		cQP := &Operand[ringqp.Poly]{Value: []ringqp.Poly{eval.BuffQP[4], eval.BuffQP[5]}}
		cQP.MetaData = ctInNTT.MetaData

		// Buffer mod Q (i.e. to store the result of gadget products)
		cQ, err := NewCiphertextAtLevelFromPoly(levelQ, []ring.Poly{cQP.Value[0].Q, cQP.Value[1].Q})

		if err != nil {
			panic(err)
		}

		cQ.MetaData = ctInNTT.MetaData

		state := false
		copy := true
		// Binary reading of the input n
		for i, j := 0, n; j > 0; i, j = i+1, j>>1 {

			// Starts by decomposing the input ciphertext
			eval.DecomposeNTT(levelQ, levelP, levelP+1, ctInNTT.Value[1], true, eval.BuffDecompQP)

			// If the binary reading scans a 1 (j is odd)
			if j&1 == 1 {

				k := n - (n & ((2 << i) - 1))
				k *= batchSize

				// If the rotation is not zero
				if k != 0 {

					rot := eval.params.GaloisElement(k)

					// opOutQP = opOutQP + Rotate(ctInNTT, k)
					if copy {
						if err = eval.AutomorphismHoistedLazy(levelQ, ctInNTT, eval.BuffDecompQP, rot, accQP); err != nil {
							return err
						}
						copy = false
					} else {
						if err = eval.AutomorphismHoistedLazy(levelQ, ctInNTT, eval.BuffDecompQP, rot, cQP); err != nil {
							return err
						}
						ringQP.Add(accQP.Value[0], cQP.Value[0], accQP.Value[0])
						ringQP.Add(accQP.Value[1], cQP.Value[1], accQP.Value[1])
					}

					// j is even
				} else {

					state = true

					// if n is not a power of two, then at least one j was odd, and thus the buffer opOutQP is not empty
					if n&(n-1) != 0 {

						// opOut = opOutQP/P + ctInNTT
						eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, accQP.Value[0].Q, accQP.Value[0].P, opOut.Value[0]) // Division by P
						eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, accQP.Value[1].Q, accQP.Value[1].P, opOut.Value[1]) // Division by P

						ringQ.Add(opOut.Value[0], ctInNTT.Value[0], opOut.Value[0])
						ringQ.Add(opOut.Value[1], ctInNTT.Value[1], opOut.Value[1])

					} else {
						ring.CopyLvl(levelQ, ctInNTT.Value[0], opOut.Value[0])
						ring.CopyLvl(levelQ, ctInNTT.Value[1], opOut.Value[1])
					}
				}
			}

			if !state {

				rot := eval.params.GaloisElement((1 << i) * batchSize)

				// ctInNTT = ctInNTT + Rotate(ctInNTT, 2^i)
				if err = eval.AutomorphismHoisted(levelQ, ctInNTT, eval.BuffDecompQP, rot, cQ); err != nil {
					return err
				}
				ringQ.Add(ctInNTT.Value[0], cQ.Value[0], ctInNTT.Value[0])
				ringQ.Add(ctInNTT.Value[1], cQ.Value[1], ctInNTT.Value[1])
			}
		}
	}

	if !ctIn.IsNTT {
		ringQ.INTT(opOut.Value[0], opOut.Value[0])
		ringQ.INTT(opOut.Value[1], opOut.Value[1])
	}

	return
}

// Replicate applies an optimized replication on the Ciphertext (log2(n) + HW(n) rotations with double hoisting).
// It acts as the inverse of a inner sum (summing elements from left to right).
// The replication is parameterized by the size of the sub-vectors to replicate "batchSize" and
// the number of times 'n' they need to be replicated.
// To ensure correctness, a gap of zero values of size batchSize * (n-1) must exist between
// two consecutive sub-vectors to replicate.
// This method is faster than Replicate when the number of rotations is large and it uses log2(n) + HW(n) instead of 'n'.
func (eval Evaluator) Replicate(ctIn *Ciphertext, batchSize, n int, opOut *Ciphertext) (err error) {
	return eval.InnerSum(ctIn, -batchSize, n, opOut)
}
