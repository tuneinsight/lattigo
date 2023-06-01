package rlwe

import (
	"fmt"
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"

	"runtime"
)

// LinearTransform is a type for linear transformations on ciphertexts.
// It stores a plaintext matrix in diagonal form and
// can be evaluated on a ciphertext by using the evaluator.LinearTransform method.
type LinearTransform struct {
	MetaData
	LogBSGSRatio int
	N1           int                 // N1 is the number of inner loops of the baby-step giant-step algorithm used in the evaluation (if N1 == 0, BSGS is not used).
	Level        int                 // Level is the level at which the matrix is encoded (can be circuit dependent)
	Vec          map[int]ringqp.Poly // Vec is the matrix, in diagonal form, where each entry of vec is an indexed non-zero diagonal.
}

// NewLinearTransform allocates a new LinearTransform with zero plaintexts at the specified level.
//
// inputs:
// - params: a struct compliant to the ParametersInterface
// - nonZeroDiags: the list of the indexes of the non-zero diagonals
// - level: the level of the encoded diagonals
// - plaintextScale: the scaling factor of the encoded diagonals
// - plaintextLogDimensions: the log2 dimension of the plaintext matrix (e.g. [1, x] for BFV/BGV and [0, x] for CKKS)
// - logBSGSRatio: the log2 ratio outer/inner loops of the BSGS linear transform evaluation algorithm. Set to -1 to not use the BSGS algorithm.
func NewLinearTransform(params ParametersInterface, nonZeroDiags []int, level int, plaintextScale Scale, plaintextLogDimensions [2]int, LogBSGSRatio int) LinearTransform {
	vec := make(map[int]ringqp.Poly)
	cols := 1 << plaintextLogDimensions[1]
	levelQ := level
	levelP := params.MaxLevelP()
	ringQP := params.RingQP().AtLevel(levelQ, levelP)
	var N1 int
	if LogBSGSRatio < 0 {
		N1 = 0
		for _, i := range nonZeroDiags {
			idx := i
			if idx < 0 {
				idx += cols
			}
			vec[idx] = *ringQP.NewPoly()
		}
	} else {
		N1 = FindBestBSGSRatio(nonZeroDiags, cols, LogBSGSRatio)
		index, _, _ := BSGSIndex(nonZeroDiags, cols, N1)
		for j := range index {
			for _, i := range index[j] {
				vec[j+i] = *ringQP.NewPoly()
			}
		}
	}

	metadata := MetaData{
		PlaintextLogDimensions: plaintextLogDimensions,
		PlaintextScale:         plaintextScale,
		EncodingDomain:         FrequencyDomain,
		IsNTT:                  true,
		IsMontgomery:           true,
	}

	return LinearTransform{MetaData: metadata, LogBSGSRatio: LogBSGSRatio, N1: N1, Level: level, Vec: vec}
}

// GaloisElements returns the list of Galois elements needed for the evaluation of the linear transformation.
func (LT *LinearTransform) GaloisElements(params ParametersInterface) (galEls []uint64) {
	return params.GaloisElementsForLinearTransform(utils.GetKeys(LT.Vec), LT.PlaintextLogDimensions[1], LT.LogBSGSRatio)
}

// EncodeLinearTransform encodes on a pre-allocated LinearTransform a set of non-zero diagonales of a matrix representing a linear transformation.
//
// inputs:
// - LT: a pre-allocated LinearTransform using `NewLinearTransform`
// - diagonals: the set of non-zero diagonals
// - encoder: an struct complying to the EncoderInterface
func EncodeLinearTransform[T any](LT LinearTransform, diagonals map[int][]T, encoder EncoderInterface[T, ringqp.Poly]) (err error) {

	scale := LT.PlaintextScale
	PlaintextLogDimensions := LT.PlaintextLogDimensions
	rows := 1 << PlaintextLogDimensions[0]
	cols := 1 << PlaintextLogDimensions[1]
	N1 := LT.N1

	keys := utils.GetKeys(diagonals)

	buf := make([]T, rows*cols)

	metaData := MetaData{
		PlaintextLogDimensions: PlaintextLogDimensions,
		IsNTT:                  true,
		IsMontgomery:           true,
		PlaintextScale:         scale,
	}

	if N1 == 0 {
		for _, i := range keys {

			idx := i
			if idx < 0 {
				idx += cols
			}

			if vec, ok := LT.Vec[idx]; !ok {
				return fmt.Errorf("cannot Encode: error encoding on LinearTransform: plaintext diagonal [%d] does not exist", idx)
			} else {
				if err = rotateAndEncodeDiagonal(diagonals, encoder, i, 0, metaData, buf, vec); err != nil {
					return
				}
			}
		}
	} else {

		index, _, _ := BSGSIndex(keys, cols, N1)

		for j := range index {

			rot := -j & (cols - 1)

			for _, i := range index[j] {

				if vec, ok := LT.Vec[i+j]; !ok {
					return fmt.Errorf("cannot Encode: error encoding on LinearTransform BSGS: input does not match the same non-zero diagonals")
				} else {
					if err = rotateAndEncodeDiagonal(diagonals, encoder, i+j, rot, metaData, buf, vec); err != nil {
						return
					}
				}
			}
		}
	}

	return
}

func rotateAndEncodeDiagonal[T any](diagonals map[int][]T, encoder EncoderInterface[T, ringqp.Poly], i, rot int, metaData MetaData, buf []T, poly ringqp.Poly) error {

	rows := 1 << metaData.PlaintextLogDimensions[0]
	cols := 1 << metaData.PlaintextLogDimensions[1]

	// manages inputs that have rotation between 0 and cols-1 or between -cols/2 and cols/2-1
	v, ok := diagonals[i]
	if !ok {
		if v, ok = diagonals[i-cols]; !ok {
			return fmt.Errorf("cannot EncodeLinearTransformDiagonalNaive: diagonal [%d] doesn't exist", i)
		}
	}

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

// GenLinearTransform allocates a new LinearTransform encoding the provided set of non-zero diagonals of a matrix representing a linear transformation.
//
// inputs:
// - diagonals: the set of non-zero diagonals
// - encoder: an struct complying to the EncoderInterface
// - level: the level of the encoded diagonals
// - plaintextScale: the scaling factor of the encoded diagonals
// - plaintextLogDimensions: the log2 dimension of the plaintext matrix (e.g. [1, x] for BFV/BGV and [0, x] for CKKS)
// - logBSGSRatio: the log2 ratio outer/inner loops of the BSGS linear transform evaluation algorithm. Set to -1 to not use the BSGS algorithm.
func GenLinearTransform[T any](diagonals map[int][]T, encoder EncoderInterface[T, ringqp.Poly], level int, plaintextScale Scale, plaintextLogDimensions [2]int, logBSGSRatio int) (LT LinearTransform, err error) {

	params := encoder.Parameters()

	ringQP := params.RingQP().AtLevel(level, params.MaxLevelP())

	rows := 1 << plaintextLogDimensions[0]
	cols := 1 << plaintextLogDimensions[1]

	keys := utils.GetKeys(diagonals)

	buf := make([]T, cols*rows)

	vec := make(map[int]ringqp.Poly)

	metaData := MetaData{
		PlaintextLogDimensions: plaintextLogDimensions,
		EncodingDomain:         FrequencyDomain,
		IsNTT:                  true,
		IsMontgomery:           true,
		PlaintextScale:         plaintextScale,
	}

	var N1 int

	if logBSGSRatio < 0 {

		for _, i := range keys {

			idx := i
			if idx < 0 {
				idx += cols
			}

			pt := *ringQP.NewPoly()

			if err = rotateAndEncodeDiagonal(diagonals, encoder, i, 0, metaData, buf, pt); err != nil {
				return
			}

			vec[idx] = pt
		}

	} else {

		// N1*N2 = N
		N1 = FindBestBSGSRatio(keys, cols, logBSGSRatio)
		index, _, _ := BSGSIndex(keys, cols, N1)

		for j := range index {

			rot := -j & (cols - 1)

			for _, i := range index[j] {

				pt := *ringQP.NewPoly()

				if err = rotateAndEncodeDiagonal(diagonals, encoder, i+j, rot, metaData, buf, pt); err != nil {
					return
				}

				vec[i+j] = pt

			}
		}
	}

	return LinearTransform{MetaData: metaData, LogBSGSRatio: logBSGSRatio, N1: N1, Vec: vec, Level: level}, nil
}

// LinearTransformNew evaluates a linear transform on the pre-allocated Ciphertexts.
// The linearTransform can either be an (ordered) list of LinearTransform or a single LinearTransform.
// In either case a list of Ciphertext is returned (the second case returning a list containing a single Ciphertext).
func (eval *Evaluator) LinearTransformNew(ctIn *Ciphertext, linearTransform interface{}) (ctOut []*Ciphertext) {

	switch LTs := linearTransform.(type) {
	case []LinearTransform:
		ctOut = make([]*Ciphertext, len(LTs))

		var maxLevel int
		for _, LT := range LTs {
			maxLevel = utils.Max(maxLevel, LT.Level)
		}

		minLevel := utils.Min(maxLevel, ctIn.Level())
		eval.DecomposeNTT(minLevel, eval.params.MaxLevelP(), eval.params.PCount(), ctIn.Value[1], ctIn.IsNTT, eval.BuffDecompQP)

		for i, LT := range LTs {
			ctOut[i] = NewCiphertext(eval.params, 1, minLevel)

			if LT.N1 == 0 {
				eval.MultiplyByDiagMatrix(ctIn, LT, eval.BuffDecompQP, ctOut[i])
			} else {
				eval.MultiplyByDiagMatrixBSGS(ctIn, LT, eval.BuffDecompQP, ctOut[i])
			}
		}

	case LinearTransform:

		minLevel := utils.Min(LTs.Level, ctIn.Level())
		eval.DecomposeNTT(minLevel, eval.params.MaxLevelP(), eval.params.PCount(), ctIn.Value[1], ctIn.IsNTT, eval.BuffDecompQP)

		ctOut = []*Ciphertext{NewCiphertext(eval.params, 1, minLevel)}

		if LTs.N1 == 0 {
			eval.MultiplyByDiagMatrix(ctIn, LTs, eval.BuffDecompQP, ctOut[0])
		} else {
			eval.MultiplyByDiagMatrixBSGS(ctIn, LTs, eval.BuffDecompQP, ctOut[0])
		}
	}
	return
}

// LinearTransform evaluates a linear transform on the pre-allocated Ciphertexts.
// The linearTransform can either be an (ordered) list of LinearTransform or a single LinearTransform.
// In either case a list of Ciphertext is returned (the second case returning a list containing a single Ciphertext).
func (eval *Evaluator) LinearTransform(ctIn *Ciphertext, linearTransform interface{}, ctOut []*Ciphertext) {

	switch LTs := linearTransform.(type) {
	case []LinearTransform:
		var maxLevel int
		for _, LT := range LTs {
			maxLevel = utils.Max(maxLevel, LT.Level)
		}

		minLevel := utils.Min(maxLevel, ctIn.Level())
		eval.DecomposeNTT(minLevel, eval.params.MaxLevelP(), eval.params.PCount(), ctIn.Value[1], true, eval.BuffDecompQP)

		for i, LT := range LTs {
			if LT.N1 == 0 {
				eval.MultiplyByDiagMatrix(ctIn, LT, eval.BuffDecompQP, ctOut[i])
			} else {
				eval.MultiplyByDiagMatrixBSGS(ctIn, LT, eval.BuffDecompQP, ctOut[i])
			}
		}

	case LinearTransform:
		minLevel := utils.Min(LTs.Level, ctIn.Level())
		eval.DecomposeNTT(minLevel, eval.params.MaxLevelP(), eval.params.PCount(), ctIn.Value[1], true, eval.BuffDecompQP)
		if LTs.N1 == 0 {
			eval.MultiplyByDiagMatrix(ctIn, LTs, eval.BuffDecompQP, ctOut[0])
		} else {
			eval.MultiplyByDiagMatrixBSGS(ctIn, LTs, eval.BuffDecompQP, ctOut[0])
		}
	}
}

// MultiplyByDiagMatrix multiplies the Ciphertext "ctIn" by the plaintext matrix "matrix" and returns the result on the Ciphertext
// "ctOut". Memory buffers for the decomposed ciphertext BuffDecompQP, BuffDecompQP must be provided, those are list of poly of ringQ and ringP
// respectively, each of size params.Beta().
// The naive approach is used (single hoisting and no baby-step giant-step), which is faster than MultiplyByDiagMatrixBSGS
// for matrix of only a few non-zero diagonals but uses more keys.
func (eval *Evaluator) MultiplyByDiagMatrix(ctIn *Ciphertext, matrix LinearTransform, BuffDecompQP []ringqp.Poly, ctOut *Ciphertext) {

	ctOut.MetaData = ctIn.MetaData
	ctOut.PlaintextScale = ctOut.PlaintextScale.Mul(matrix.PlaintextScale)

	levelQ := utils.Min(ctOut.Level(), utils.Min(ctIn.Level(), matrix.Level))
	levelP := eval.params.RingP().MaxLevel()

	ringQP := eval.params.RingQP().AtLevel(levelQ, levelP)
	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	ctOut.Resize(ctOut.Degree(), levelQ)

	QiOverF := eval.params.QiOverflowMargin(levelQ)
	PiOverF := eval.params.PiOverflowMargin(levelP)

	c0OutQP := ringqp.Poly{Q: ctOut.Value[0], P: eval.BuffQP[5].Q}
	c1OutQP := ringqp.Poly{Q: ctOut.Value[1], P: eval.BuffQP[5].P}

	ct0TimesP := eval.BuffQP[0].Q // ct0 * P mod Q
	tmp0QP := eval.BuffQP[1]
	tmp1QP := eval.BuffQP[2]

	cQP := &OperandQP{}
	cQP.Value = []*ringqp.Poly{&eval.BuffQP[3], &eval.BuffQP[4]}
	cQP.IsNTT = true

	ring.Copy(ctIn.Value[0], eval.BuffCt.Value[0])
	ring.Copy(ctIn.Value[1], eval.BuffCt.Value[1])
	ctInTmp0, ctInTmp1 := eval.BuffCt.Value[0], eval.BuffCt.Value[1]

	ringQ.MulScalarBigint(ctInTmp0, ringP.ModulusAtLevel[levelP], ct0TimesP) // P*c0

	slots := 1 << matrix.PlaintextLogDimensions[1]

	var state bool
	var cnt int
	for k := range matrix.Vec {

		k &= (slots - 1)

		if k == 0 {
			state = true
		} else {

			galEl := eval.params.GaloisElement(k)

			var evk *GaloisKey
			var err error
			if evk, err = eval.CheckAndGetGaloisKey(galEl); err != nil {
				panic(fmt.Errorf("cannot apply Automorphism: %w", err))
			}

			index := eval.AutomorphismIndex[galEl]

			eval.GadgetProductHoistedLazy(levelQ, BuffDecompQP, &evk.GadgetCiphertext, cQP)
			ringQ.Add(cQP.Value[0].Q, ct0TimesP, cQP.Value[0].Q)
			ringQP.AutomorphismNTTWithIndex(cQP.Value[0], index, &tmp0QP)
			ringQP.AutomorphismNTTWithIndex(cQP.Value[1], index, &tmp1QP)

			pt := matrix.Vec[k]

			if cnt == 0 {
				// keyswitch(c1_Q) = (d0_QP, d1_QP)
				ringQP.MulCoeffsMontgomery(&pt, &tmp0QP, &c0OutQP)
				ringQP.MulCoeffsMontgomery(&pt, &tmp1QP, &c1OutQP)
			} else {
				// keyswitch(c1_Q) = (d0_QP, d1_QP)
				ringQP.MulCoeffsMontgomeryThenAdd(&pt, &tmp0QP, &c0OutQP)
				ringQP.MulCoeffsMontgomeryThenAdd(&pt, &tmp1QP, &c1OutQP)
			}

			if cnt%QiOverF == QiOverF-1 {
				ringQ.Reduce(c0OutQP.Q, c0OutQP.Q)
				ringQ.Reduce(c1OutQP.Q, c1OutQP.Q)
			}

			if cnt%PiOverF == PiOverF-1 {
				ringP.Reduce(c0OutQP.P, c0OutQP.P)
				ringP.Reduce(c1OutQP.P, c1OutQP.P)
			}

			cnt++
		}
	}

	if cnt%QiOverF == 0 {
		ringQ.Reduce(c0OutQP.Q, c0OutQP.Q)
		ringQ.Reduce(c1OutQP.Q, c1OutQP.Q)
	}

	if cnt%PiOverF == 0 {
		ringP.Reduce(c0OutQP.P, c0OutQP.P)
		ringP.Reduce(c1OutQP.P, c1OutQP.P)
	}

	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c0OutQP.Q, c0OutQP.P, c0OutQP.Q) // sum(phi(c0 * P + d0_QP))/P
	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c1OutQP.Q, c1OutQP.P, c1OutQP.Q) // sum(phi(d1_QP))/P

	if state { // Rotation by zero
		ringQ.MulCoeffsMontgomeryThenAdd(matrix.Vec[0].Q, ctInTmp0, c0OutQP.Q) // ctOut += c0_Q * plaintext
		ringQ.MulCoeffsMontgomeryThenAdd(matrix.Vec[0].Q, ctInTmp1, c1OutQP.Q) // ctOut += c1_Q * plaintext
	}
}

// MultiplyByDiagMatrixBSGS multiplies the Ciphertext "ctIn" by the plaintext matrix "matrix" and returns the result on the Ciphertext
// "ctOut". Memory buffers for the decomposed Ciphertext BuffDecompQP, BuffDecompQP must be provided, those are list of poly of ringQ and ringP
// respectively, each of size params.Beta().
// The BSGS approach is used (double hoisting with baby-step giant-step), which is faster than MultiplyByDiagMatrix
// for matrix with more than a few non-zero diagonals and uses significantly less keys.
func (eval *Evaluator) MultiplyByDiagMatrixBSGS(ctIn *Ciphertext, matrix LinearTransform, BuffDecompQP []ringqp.Poly, ctOut *Ciphertext) {

	ctOut.MetaData = ctIn.MetaData
	ctOut.PlaintextScale = ctOut.PlaintextScale.Mul(matrix.PlaintextScale)

	levelQ := utils.Min(ctOut.Level(), utils.Min(ctIn.Level(), matrix.Level))
	levelP := eval.Parameters().MaxLevelP()

	ringQP := eval.params.RingQP().AtLevel(levelQ, levelP)
	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	ctOut.Resize(ctOut.Degree(), levelQ)

	QiOverF := eval.params.QiOverflowMargin(levelQ) >> 1
	PiOverF := eval.params.PiOverflowMargin(levelP) >> 1

	// Computes the N2 rotations indexes of the non-zero rows of the diagonalized DFT matrix for the baby-step giant-step algorithm
	index, _, rotN2 := BSGSIndex(utils.GetKeys(matrix.Vec), 1<<matrix.PlaintextLogDimensions[1], matrix.N1)

	ring.Copy(ctIn.Value[0], eval.BuffCt.Value[0])
	ring.Copy(ctIn.Value[1], eval.BuffCt.Value[1])

	ctInTmp0, ctInTmp1 := eval.BuffCt.Value[0], eval.BuffCt.Value[1]

	// Pre-rotates ciphertext for the baby-step giant-step algorithm, does not divide by P yet
	ctInRotQP := map[int]*OperandQP{}
	for _, i := range rotN2 {
		if i != 0 {
			ctInRotQP[i] = NewOperandQP(eval.Parameters(), 1, levelQ, levelP)
			eval.AutomorphismHoistedLazy(levelQ, ctIn, BuffDecompQP, eval.Parameters().GaloisElement(i), ctInRotQP[i])
		}
	}

	// Accumulator inner loop
	tmp0QP := eval.BuffQP[1]
	tmp1QP := eval.BuffQP[2]

	// Accumulator outer loop
	cQP := &OperandQP{}
	cQP.Value = []*ringqp.Poly{&eval.BuffQP[3], &eval.BuffQP[4]}
	cQP.IsNTT = true

	// Result in QP
	c0OutQP := ringqp.Poly{Q: ctOut.Value[0], P: eval.BuffQP[5].Q}
	c1OutQP := ringqp.Poly{Q: ctOut.Value[1], P: eval.BuffQP[5].P}

	ringQ.MulScalarBigint(ctInTmp0, ringP.ModulusAtLevel[levelP], ctInTmp0) // P*c0
	ringQ.MulScalarBigint(ctInTmp1, ringP.ModulusAtLevel[levelP], ctInTmp1) // P*c1

	// OUTER LOOP
	var cnt0 int
	for j := range index {

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
					ringQP.MulCoeffsMontgomeryLazy(&pt, ct.Value[0], &tmp0QP)
					ringQP.MulCoeffsMontgomeryLazy(&pt, ct.Value[1], &tmp1QP)
				} else {
					ringQP.MulCoeffsMontgomeryLazyThenAddLazy(&pt, ct.Value[0], &tmp0QP)
					ringQP.MulCoeffsMontgomeryLazyThenAddLazy(&pt, ct.Value[1], &tmp1QP)
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
				panic(fmt.Errorf("cannot apply Automorphism: %w", err))
			}

			rotIndex := eval.AutomorphismIndex[galEl]

			eval.GadgetProductLazy(levelQ, tmp1QP.Q, &evk.GadgetCiphertext, cQP) // EvaluationKey(P*phi(tmpRes_1)) = (d0, d1) in base QP
			ringQP.Add(cQP.Value[0], &tmp0QP, cQP.Value[0])

			// Outer loop rotations
			if cnt0 == 0 {
				ringQP.AutomorphismNTTWithIndex(cQP.Value[0], rotIndex, &c0OutQP)
				ringQP.AutomorphismNTTWithIndex(cQP.Value[1], rotIndex, &c1OutQP)
			} else {
				ringQP.AutomorphismNTTWithIndexThenAddLazy(cQP.Value[0], rotIndex, &c0OutQP)
				ringQP.AutomorphismNTTWithIndexThenAddLazy(cQP.Value[1], rotIndex, &c1OutQP)
			}

			// Else directly adds on ((cQP.Value[0].Q, cQP.Value[0].P), (cQP.Value[1].Q, cQP.Value[1].P))
		} else {
			if cnt0 == 0 {
				ringqp.CopyLvl(levelQ, levelP, &tmp0QP, &c0OutQP)
				ringqp.CopyLvl(levelQ, levelP, &tmp1QP, &c1OutQP)
			} else {
				ringQP.AddLazy(&c0OutQP, &tmp0QP, &c0OutQP)
				ringQP.AddLazy(&c1OutQP, &tmp1QP, &c1OutQP)
			}
		}

		if cnt0%QiOverF == QiOverF-1 {
			ringQ.Reduce(ctOut.Value[0], ctOut.Value[0])
			ringQ.Reduce(ctOut.Value[1], ctOut.Value[1])
		}

		if cnt0%PiOverF == PiOverF-1 {
			ringP.Reduce(c0OutQP.P, c0OutQP.P)
			ringP.Reduce(c1OutQP.P, c1OutQP.P)
		}

		cnt0++
	}

	if cnt0%QiOverF != 0 {
		ringQ.Reduce(ctOut.Value[0], ctOut.Value[0])
		ringQ.Reduce(ctOut.Value[1], ctOut.Value[1])
	}

	if cnt0%PiOverF != 0 {
		ringP.Reduce(c0OutQP.P, c0OutQP.P)
		ringP.Reduce(c1OutQP.P, c1OutQP.P)
	}

	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, ctOut.Value[0], c0OutQP.P, ctOut.Value[0]) // sum(phi(c0 * P + d0_QP))/P
	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, ctOut.Value[1], c1OutQP.P, ctOut.Value[1]) // sum(phi(d1_QP))/P

	ctInRotQP = nil
	runtime.GC()
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
func (eval *Evaluator) Trace(ctIn *Ciphertext, logN int, ctOut *Ciphertext) {

	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		panic("ctIn.Degree() != 1 or ctOut.Degree() != 1")
	}

	level := utils.Min(ctIn.Level(), ctOut.Level())

	ctOut.Resize(ctOut.Degree(), level)

	ctOut.MetaData = ctIn.MetaData

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
		ringQ.MulScalarBigint(ctIn.Value[0], NInv, ctOut.Value[0])
		ringQ.MulScalarBigint(ctIn.Value[1], NInv, ctOut.Value[1])

		if !ctIn.IsNTT {
			ringQ.NTT(ctOut.Value[0], ctOut.Value[0])
			ringQ.NTT(ctOut.Value[1], ctOut.Value[1])
			ctOut.IsNTT = true
		}

		buff := NewCiphertextAtLevelFromPoly(level, []*ring.Poly{eval.BuffQP[3].Q, eval.BuffQP[4].Q})
		buff.IsNTT = true

		for i := logN; i < eval.params.LogN()-1; i++ {
			eval.Automorphism(ctOut, eval.params.GaloisElement(1<<i), buff)
			ringQ.Add(ctOut.Value[0], buff.Value[0], ctOut.Value[0])
			ringQ.Add(ctOut.Value[1], buff.Value[1], ctOut.Value[1])
		}

		if logN == 0 && ringQ.Type() == ring.Standard {
			eval.Automorphism(ctOut, ringQ.NthRoot()-1, buff)
			ringQ.Add(ctOut.Value[0], buff.Value[0], ctOut.Value[0])
			ringQ.Add(ctOut.Value[1], buff.Value[1], ctOut.Value[1])
		}

		if !ctIn.IsNTT {
			ringQ.INTT(ctOut.Value[0], ctOut.Value[0])
			ringQ.INTT(ctOut.Value[1], ctOut.Value[1])
			ctOut.IsNTT = false
		}

	} else {
		if ctIn != ctOut {
			ctOut.Copy(ctIn)
		}
	}
}

// Expand expands a RLWE Ciphertext encrypting sum ai * X^i to 2^logN ciphertexts,
// each encrypting ai * X^0 for 0 <= i < 2^LogN. That is, it extracts the first 2^logN
// coefficients, whose degree is a multiple of 2^logGap, of ctIn and returns an RLWE
// Ciphertext for each coefficient extracted.
func (eval *Evaluator) Expand(ctIn *Ciphertext, logN, logGap int) (ctOut []*Ciphertext) {

	if ctIn.Degree() != 1 {
		panic("ctIn.Degree() != 1")
	}

	if eval.params.RingType() != ring.Standard {
		panic("Expand is only supported for ring.Type = ring.Standard (X^{-2^{i}} does not exist in the sub-ring Z[X + X^{-1}])")
	}

	params := eval.params

	level := ctIn.Level()

	ringQ := params.RingQ().AtLevel(level)

	// Compute X^{-2^{i}} from 1 to LogN
	xPow2 := genXPow2(ringQ, logN, true)

	ctOut = make([]*Ciphertext, 1<<(logN-logGap))
	ctOut[0] = ctIn.CopyNew()
	ctOut[0].PlaintextLogDimensions = [2]int{0, 0}

	if ct := ctOut[0]; !ctIn.IsNTT {
		ringQ.NTT(ct.Value[0], ct.Value[0])
		ringQ.NTT(ct.Value[1], ct.Value[1])
		ct.IsNTT = true
	}

	// Multiplies by 2^{-logN} mod Q
	NInv := new(big.Int).SetUint64(1 << logN)
	NInv.ModInverse(NInv, ringQ.ModulusAtLevel[level])

	ringQ.MulScalarBigint(ctOut[0].Value[0], NInv, ctOut[0].Value[0])
	ringQ.MulScalarBigint(ctOut[0].Value[1], NInv, ctOut[0].Value[1])

	gap := 1 << logGap

	tmp := NewCiphertextAtLevelFromPoly(level, []*ring.Poly{eval.BuffCt.Value[0], eval.BuffCt.Value[1]})
	tmp.MetaData = ctIn.MetaData

	for i := 0; i < logN; i++ {

		n := 1 << i

		galEl := uint64(ringQ.N()/n + 1)

		half := n / gap

		for j := 0; j < (n+gap-1)/gap; j++ {

			c0 := ctOut[j]

			// X -> X^{N/n + 1}
			//[a, b, c, d] -> [a, -b, c, -d]
			eval.Automorphism(c0, galEl, tmp)

			if j+half > 0 {

				c1 := ctOut[j].CopyNew()

				// Zeroes odd coeffs: [a, b, c, d] + [a, -b, c, -d] -> [2a, 0, 2b, 0]
				ringQ.Add(c0.Value[0], tmp.Value[0], c0.Value[0])
				ringQ.Add(c0.Value[1], tmp.Value[1], c0.Value[1])

				// Zeroes even coeffs: [a, b, c, d] - [a, -b, c, -d] -> [0, 2b, 0, 2d]
				ringQ.Sub(c1.Value[0], tmp.Value[0], c1.Value[0])
				ringQ.Sub(c1.Value[1], tmp.Value[1], c1.Value[1])

				// c1 * X^{-2^{i}}: [0, 2b, 0, 2d] * X^{-n} -> [2b, 0, 2d, 0]
				ringQ.MulCoeffsMontgomery(c1.Value[0], xPow2[i], c1.Value[0])
				ringQ.MulCoeffsMontgomery(c1.Value[1], xPow2[i], c1.Value[1])

				ctOut[j+half] = c1

			} else {

				// Zeroes odd coeffs: [a, b, c, d] + [a, -b, c, -d] -> [2a, 0, 2b, 0]
				ringQ.Add(c0.Value[0], tmp.Value[0], c0.Value[0])
				ringQ.Add(c0.Value[1], tmp.Value[1], c0.Value[1])
			}
		}
	}

	for _, ct := range ctOut {
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
func (eval *Evaluator) Pack(cts map[int]*Ciphertext, inputLogGap int, zeroGarbageSlots bool) (ct *Ciphertext) {

	params := eval.Parameters()

	if params.RingType() != ring.Standard {
		panic(fmt.Errorf("cannot Pack: procedure is only supported for ring.Type = ring.Standard (X^{2^{i}} does not exist in the sub-ring Z[X + X^{-1}])"))
	}

	if len(cts) < 2 {
		panic(fmt.Errorf("cannot Pack: #cts must be at least 2"))
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
		panic(fmt.Errorf("cannot PackRLWE: gaps between ciphertexts is smaller than inputLogGap > N"))
	}

	xPow2 := genXPow2(ringQ.AtLevel(level), params.LogN(), false) // log(N) polynomial to generate, quick

	NInv := new(big.Int).SetUint64(uint64(1 << (logEnd - logStart)))
	NInv.ModInverse(NInv, ringQ.ModulusAtLevel[level])

	for _, key := range keys {

		ct := cts[key]

		if ct.Degree() != 1 {
			panic(fmt.Errorf("cannot PackRLWE: cts[%d].Degree() != 1", key))
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
	tmpa.Value = []*ring.Poly{ringQ.NewPoly(), ringQ.NewPoly()}
	tmpa.IsNTT = true

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
					eval.Automorphism(tmpa, galEl, tmpa)
				} else {
					eval.Automorphism(a, galEl, tmpa)
				}

				// a + b * X^{N/2^{i}} + phi(a - b * X^{N/2^{i}}, 2^{i-1})
				ringQ.Add(a.Value[0], tmpa.Value[0], a.Value[0])
				ringQ.Add(a.Value[1], tmpa.Value[1], a.Value[1])
			}
		}
	}

	return cts[0]
}

func genXPow2(r *ring.Ring, logN int, div bool) (xPow []*ring.Poly) {

	// Compute X^{-n} from 0 to LogN
	xPow = make([]*ring.Poly, logN)

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
// It outputs in ctOut a Ciphertext for which the "leftmost" sub-vector of each group is equal to the sum of the group.
func (eval *Evaluator) InnerSum(ctIn *Ciphertext, batchSize, n int, ctOut *Ciphertext) {

	levelQ := ctIn.Level()
	levelP := eval.params.PCount() - 1

	ringQP := eval.params.RingQP().AtLevel(ctIn.Level(), levelP)

	ringQ := ringQP.RingQ

	ctOut.Resize(ctOut.Degree(), levelQ)
	ctOut.MetaData = ctIn.MetaData

	ctInNTT := NewCiphertextAtLevelFromPoly(levelQ, eval.BuffCt.Value[:2])
	ctInNTT.IsNTT = true

	if !ctIn.IsNTT {
		ringQ.NTT(ctIn.Value[0], ctInNTT.Value[0])
		ringQ.NTT(ctIn.Value[1], ctInNTT.Value[1])
	} else {
		ring.CopyLvl(levelQ, ctIn.Value[0], ctInNTT.Value[0])
		ring.CopyLvl(levelQ, ctIn.Value[1], ctInNTT.Value[1])
	}

	if n == 1 {
		if ctIn != ctOut {
			ring.CopyLvl(levelQ, ctIn.Value[0], ctOut.Value[0])
			ring.CopyLvl(levelQ, ctIn.Value[1], ctOut.Value[1])
		}
	} else {

		// BuffQP[0:2] are used by AutomorphismHoistedLazy

		// Accumulator mod QP (i.e. ctOut Mod QP)
		accQP := &OperandQP{Value: []*ringqp.Poly{&eval.BuffQP[2], &eval.BuffQP[3]}}
		accQP.IsNTT = true

		// Buffer mod QP (i.e. to store the result of lazy gadget products)
		cQP := &OperandQP{Value: []*ringqp.Poly{&eval.BuffQP[4], &eval.BuffQP[5]}}
		cQP.IsNTT = true

		// Buffer mod Q (i.e. to store the result of gadget products)
		cQ := NewCiphertextAtLevelFromPoly(levelQ, []*ring.Poly{cQP.Value[0].Q, cQP.Value[1].Q})
		cQ.IsNTT = true

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

					// ctOutQP = ctOutQP + Rotate(ctInNTT, k)
					if copy {
						eval.AutomorphismHoistedLazy(levelQ, ctInNTT, eval.BuffDecompQP, rot, accQP)
						copy = false
					} else {
						eval.AutomorphismHoistedLazy(levelQ, ctInNTT, eval.BuffDecompQP, rot, cQP)
						ringQP.Add(accQP.Value[0], cQP.Value[0], accQP.Value[0])
						ringQP.Add(accQP.Value[1], cQP.Value[1], accQP.Value[1])
					}

					// j is even
				} else {

					state = true

					// if n is not a power of two, then at least one j was odd, and thus the buffer ctOutQP is not empty
					if n&(n-1) != 0 {

						// ctOut = ctOutQP/P + ctInNTT
						eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, accQP.Value[0].Q, accQP.Value[0].P, ctOut.Value[0]) // Division by P
						eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, accQP.Value[1].Q, accQP.Value[1].P, ctOut.Value[1]) // Division by P

						ringQ.Add(ctOut.Value[0], ctInNTT.Value[0], ctOut.Value[0])
						ringQ.Add(ctOut.Value[1], ctInNTT.Value[1], ctOut.Value[1])

					} else {
						ring.CopyLvl(levelQ, ctInNTT.Value[0], ctOut.Value[0])
						ring.CopyLvl(levelQ, ctInNTT.Value[1], ctOut.Value[1])
					}
				}
			}

			if !state {

				rot := eval.params.GaloisElement((1 << i) * batchSize)

				// ctInNTT = ctInNTT + Rotate(ctInNTT, 2^i)
				eval.AutomorphismHoisted(levelQ, ctInNTT, eval.BuffDecompQP, rot, cQ)
				ringQ.Add(ctInNTT.Value[0], cQ.Value[0], ctInNTT.Value[0])
				ringQ.Add(ctInNTT.Value[1], cQ.Value[1], ctInNTT.Value[1])
			}
		}
	}

	if !ctIn.IsNTT {
		ringQ.INTT(ctOut.Value[0], ctOut.Value[0])
		ringQ.INTT(ctOut.Value[1], ctOut.Value[1])
	}
}

// Replicate applies an optimized replication on the Ciphertext (log2(n) + HW(n) rotations with double hoisting).
// It acts as the inverse of a inner sum (summing elements from left to right).
// The replication is parameterized by the size of the sub-vectors to replicate "batchSize" and
// the number of times 'n' they need to be replicated.
// To ensure correctness, a gap of zero values of size batchSize * (n-1) must exist between
// two consecutive sub-vectors to replicate.
// This method is faster than Replicate when the number of rotations is large and it uses log2(n) + HW(n) instead of 'n'.
func (eval *Evaluator) Replicate(ctIn *Ciphertext, batchSize, n int, ctOut *Ciphertext) {
	eval.InnerSum(ctIn, -batchSize, n, ctOut)
}
