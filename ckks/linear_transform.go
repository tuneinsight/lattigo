package ckks

import (
	"fmt"
	"runtime"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// TraceNew maps X -> sum((-1)^i * X^{i*n+1}) for 0 <= i < N and returns the result on a new ciphertext.
// For log(n) = logSlots.
func (eval *evaluator) TraceNew(ctIn *rlwe.Ciphertext, logSlots int) (ctOut *rlwe.Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1, ctIn.Level())
	eval.Trace(ctIn, logSlots, ctOut)
	return
}

// Average returns the average of vectors of batchSize elements.
// The operation assumes that ctIn encrypts SlotCount/'batchSize' sub-vectors of size 'batchSize'.
// It then replaces all values of those sub-vectors by the component-wise average between all the sub-vectors.
// Example for batchSize=4 and slots=8: [{a, b, c, d}, {e, f, g, h}] -> [0.5*{a+e, b+f, c+g, d+h}, 0.5*{a+e, b+f, c+g, d+h}]
// Operation requires log2(SlotCout/'batchSize') rotations.
// Required rotation keys can be generated with 'RotationsForInnerSumLog(batchSize, SlotCount/batchSize)”
func (eval *evaluator) Average(ctIn *rlwe.Ciphertext, logBatchSize int, ctOut *rlwe.Ciphertext) {

	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		panic("ctIn.Degree() != 1 or ctOut.Degree() != 1")
	}

	if logBatchSize > eval.params.LogSlots() {
		panic("cannot Average: batchSize must be smaller or equal to the number of slots")
	}

	ringQ := eval.params.RingQ()

	level := utils.Min(ctIn.Level(), ctOut.Level())

	n := eval.params.Slots() / (1 << logBatchSize)

	// pre-multiplication by n^-1
	for i, s := range ringQ.SubRings[:level+1] {

		invN := ring.ModExp(uint64(n), s.Modulus-2, s.Modulus)
		invN = ring.MForm(invN, s.Modulus, s.BRedConstant)

		s.MulScalarMontgomery(ctIn.Value[0].Coeffs[i], invN, ctOut.Value[0].Coeffs[i])
		s.MulScalarMontgomery(ctIn.Value[1].Coeffs[i], invN, ctOut.Value[1].Coeffs[i])
	}

	eval.InnerSum(ctOut, 1<<logBatchSize, n, ctOut)
}

// RotateHoistedNew takes an input Ciphertext and a list of rotations and returns a map of Ciphertext, where each element of the map is the input Ciphertext
// rotation by one element of the list. It is much faster than sequential calls to Rotate.
func (eval *evaluator) RotateHoistedNew(ctIn *rlwe.Ciphertext, rotations []int) (ctOut map[int]*rlwe.Ciphertext) {
	ctOut = make(map[int]*rlwe.Ciphertext)
	for _, i := range rotations {
		ctOut[i] = NewCiphertext(eval.params, 1, ctIn.Level())
	}
	eval.RotateHoisted(ctIn, rotations, ctOut)
	return
}

// RotateHoisted takes an input Ciphertext and a list of rotations and populates a map of pre-allocated Ciphertexts,
// where each element of the map is the input Ciphertext rotation by one element of the list.
// It is much faster than sequential calls to Rotate.
func (eval *evaluator) RotateHoisted(ctIn *rlwe.Ciphertext, rotations []int, ctOut map[int]*rlwe.Ciphertext) {
	levelQ := ctIn.Level()
	eval.DecomposeNTT(levelQ, eval.params.MaxLevelP(), eval.params.PCount(), ctIn.Value[1], ctIn.IsNTT, eval.BuffDecompQP)
	for _, i := range rotations {
		eval.AutomorphismHoisted(levelQ, ctIn, eval.BuffDecompQP, eval.params.GaloisElementForColumnRotationBy(i), ctOut[i])
	}
}

// LinearTransform is a type for linear transformations on ciphertexts.
// It stores a plaintext matrix diagonalized in diagonal form and
// can be evaluated on a ciphertext by using the evaluator.LinearTransform method.
type LinearTransform struct {
	LogSlots int                 // Log of the number of slots of the plaintext (needed to compute the appropriate rotation keys)
	N1       int                 // N1 is the number of inner loops of the baby-step giant-step algorithm used in the evaluation (if N1 == 0, BSGS is not used).
	Level    int                 // Level is the level at which the matrix is encoded (can be circuit dependent)
	Scale    rlwe.Scale          // Scale is the scale at which the matrix is encoded (can be circuit dependent)
	Vec      map[int]ringqp.Poly // Vec is the matrix, in diagonal form, where each entry of vec is an indexed non-zero diagonal.
}

// NewLinearTransform allocates a new LinearTransform with zero plaintexts at the specified level.
// If LogBSGSRatio < 0, the LinearTransform is set to not use the BSGS approach.
func NewLinearTransform(params Parameters, nonZeroDiags []int, level, logSlots int, LogBSGSRatio int) LinearTransform {
	vec := make(map[int]ringqp.Poly)
	slots := 1 << logSlots
	levelQ := level
	levelP := params.PCount() - 1
	ringQP := params.RingQP().AtLevel(levelQ, levelP)
	var N1 int
	if LogBSGSRatio < 0 {
		N1 = 0
		for _, i := range nonZeroDiags {
			idx := i
			if idx < 0 {
				idx += slots
			}
			vec[idx] = *ringQP.NewPoly()
		}
	} else {
		N1 = FindBestBSGSRatio(nonZeroDiags, slots, LogBSGSRatio)
		index, _, _ := BSGSIndex(nonZeroDiags, slots, N1)
		for j := range index {
			for _, i := range index[j] {
				vec[j+i] = *ringQP.NewPoly()
			}
		}
	}

	return LinearTransform{LogSlots: logSlots, N1: N1, Level: level, Vec: vec}
}

// GaloisElements returns the list of Galois elements needed for the evaluation
// of the linear transform.
func (LT *LinearTransform) GaloisElements(params Parameters) (galEls []uint64) {
	slots := 1 << LT.LogSlots

	rotIndex := make(map[int]bool)

	var index int

	N1 := LT.N1

	if LT.N1 == 0 {

		for j := range LT.Vec {
			rotIndex[j] = true
		}

	} else {

		for j := range LT.Vec {

			index = ((j / N1) * N1) & (slots - 1)
			rotIndex[index] = true

			index = j & (N1 - 1)
			rotIndex[index] = true
		}
	}

	galEls = make([]uint64, len(rotIndex))
	var i int
	for j := range rotIndex {
		galEls[i] = params.GaloisElementForColumnRotationBy(j)
		i++
	}

	return
}

// Encode encodes on a pre-allocated LinearTransform the linear transforms' matrix in diagonal form `value`.
// values.(type) can be either map[int][]complex128 or map[int][]float64.
// User must ensure that 1 <= len([]complex128\[]float64) <= 2^logSlots < 2^logN.
// It can then be evaluated on a ciphertext using evaluator.LinearTransform.
// Evaluation will use the naive approach (single hoisting and no baby-step giant-step).
// Faster if there is only a few non-zero diagonals but uses more keys.
func (LT *LinearTransform) Encode(encoder Encoder, value interface{}, scale rlwe.Scale) {

	enc, ok := encoder.(*encoderComplex128)
	if !ok {
		panic("cannot Encode: encoder should be an encoderComplex128")
	}

	dMat := interfaceMapToMapOfInterface(value)
	slots := 1 << LT.LogSlots
	N1 := LT.N1

	if N1 == 0 {
		for i := range dMat {
			idx := i
			if idx < 0 {
				idx += slots
			}

			if _, ok := LT.Vec[idx]; !ok {
				panic("cannot Encode: error encoding on LinearTransform: input does not match the same non-zero diagonals")
			}

			enc.Embed(dMat[i], LT.LogSlots, scale, true, LT.Vec[idx])
		}
	} else {

		index, _, _ := BSGSIndex(value, slots, N1)

		var values interface{}
		switch value.(type) {
		case map[int][]complex128:
			values = make([]complex128, slots)
		case map[int][]float64:
			values = make([]float64, slots)
		}

		for j := range index {

			rot := -j & (slots - 1)

			for _, i := range index[j] {
				// manages inputs that have rotation between 0 and slots-1 or between -slots/2 and slots/2-1
				v, ok := dMat[j+i]
				if !ok {
					v = dMat[j+i-slots]
				}

				if _, ok := LT.Vec[j+i]; !ok {
					panic("cannot Encode: error encoding on LinearTransform BSGS: input does not match the same non-zero diagonals")
				}

				copyRotInterface(values, v, rot)

				enc.Embed(values, LT.LogSlots, scale, true, LT.Vec[j+i])
			}
		}
	}

	LT.Scale = scale
}

// GenLinearTransform allocates and encode a new LinearTransform struct from the linear transforms' matrix in diagonal form `value`.
// values.(type) can be either map[int][]complex128 or map[int][]float64.
// User must ensure that 1 <= len([]complex128\[]float64) <= 2^logSlots < 2^logN.
// It can then be evaluated on a ciphertext using evaluator.LinearTransform.
// Evaluation will use the naive approach (single hoisting and no baby-step giant-step).
// Faster if there is only a few non-zero diagonals but uses more keys.
func GenLinearTransform(encoder Encoder, value interface{}, level int, scale rlwe.Scale, logslots int) LinearTransform {

	enc, ok := encoder.(*encoderComplex128)
	if !ok {
		panic("cannot GenLinearTransform: encoder should be an encoderComplex128")
	}

	params := enc.params
	dMat := interfaceMapToMapOfInterface(value)
	vec := make(map[int]ringqp.Poly)
	slots := 1 << logslots
	levelQ := level
	levelP := params.PCount() - 1
	ringQP := params.RingQP().AtLevel(levelQ, levelP)
	for i := range dMat {

		idx := i
		if idx < 0 {
			idx += slots
		}
		vec[idx] = *ringQP.NewPoly()
		enc.Embed(dMat[i], logslots, scale, true, vec[idx])
	}

	return LinearTransform{LogSlots: logslots, N1: 0, Vec: vec, Level: level, Scale: scale}
}

// GenLinearTransformBSGS allocates and encodes a new LinearTransform struct from the linear transforms' matrix in diagonal form `value` for evaluation with a baby-step giant-step approach.
// values.(type) can be either map[int][]complex128 or map[int][]float64.
// User must ensure that 1 <= len([]complex128\[]float64) <= 2^logSlots < 2^logN.
// LinearTransform types can be be evaluated on a ciphertext using evaluator.LinearTransform.
// Evaluation will use the optimized approach (double hoisting and baby-step giant-step).
// Faster if there is more than a few non-zero diagonals.
// LogBSGSRatio is the log of the maximum ratio between the inner and outer loop of the baby-step giant-step algorithm used in evaluator.LinearTransform.
// Optimal LogBSGSRatio value is between 0 and 4 depending on the sparsity of the matrix.
func GenLinearTransformBSGS(encoder Encoder, value interface{}, level int, scale rlwe.Scale, LogBSGSRatio int, logSlots int) (LT LinearTransform) {

	enc, ok := encoder.(*encoderComplex128)
	if !ok {
		panic("cannot GenLinearTransformBSGS: encoder should be an encoderComplex128")
	}

	params := enc.params

	slots := 1 << logSlots

	// N1*N2 = N
	N1 := FindBestBSGSRatio(value, slots, LogBSGSRatio)

	index, _, _ := BSGSIndex(value, slots, N1)

	vec := make(map[int]ringqp.Poly)

	dMat := interfaceMapToMapOfInterface(value)
	levelQ := level
	levelP := params.PCount() - 1

	ringQP := params.RingQP().AtLevel(levelQ, levelP)

	var values interface{}
	switch value.(type) {
	case map[int][]complex128:
		values = make([]complex128, slots)
	case map[int][]float64:
		values = make([]float64, slots)
	}

	for j := range index {

		rot := -j & (slots - 1)

		for _, i := range index[j] {

			// manages inputs that have rotation between 0 and slots-1 or between -slots/2 and slots/2-1
			v, ok := dMat[j+i]
			if !ok {
				v = dMat[j+i-slots]
			}
			vec[j+i] = *ringQP.NewPoly()

			copyRotInterface(values, v, rot)

			enc.Embed(values, logSlots, scale, true, vec[j+i])
		}
	}

	return LinearTransform{LogSlots: logSlots, N1: N1, Vec: vec, Level: level, Scale: scale}
}

func copyRotInterface(a, b interface{}, rot int) {
	switch a.(type) {
	case []complex128:

		ac128 := a.([]complex128)
		bc128 := b.([]complex128)

		n := len(ac128)

		if len(bc128) >= rot {
			copy(ac128[:n-rot], bc128[rot:])
			copy(ac128[n-rot:], bc128[:rot])
		} else {
			copy(ac128[n-rot:], bc128)
		}
	case []float64:

		af64 := a.([]float64)
		bf64 := b.([]float64)

		n := len(af64)

		if len(bf64) >= rot {
			copy(af64[:n-rot], bf64[rot:])
			copy(af64[n-rot:], bf64[:rot])
		} else {
			copy(af64[n-rot:], bf64)
		}
	}
}

// BSGSIndex returns the index map and needed rotation for the BSGS matrix-vector multiplication algorithm.
func BSGSIndex(el interface{}, slots, N1 int) (index map[int][]int, rotN1, rotN2 []int) {
	index = make(map[int][]int)
	rotN1Map := make(map[int]bool)
	rotN2Map := make(map[int]bool)
	var nonZeroDiags []int
	switch element := el.(type) {
	case map[int][]complex128:
		nonZeroDiags = make([]int, len(element))
		var i int
		for key := range element {
			nonZeroDiags[i] = key
			i++
		}
	case map[int][]float64:
		nonZeroDiags = make([]int, len(element))
		var i int
		for key := range element {
			nonZeroDiags[i] = key
			i++
		}
	case map[int]bool:
		nonZeroDiags = make([]int, len(element))
		var i int
		for key := range element {
			nonZeroDiags[i] = key
			i++
		}
	case map[int]ringqp.Poly:
		nonZeroDiags = make([]int, len(element))
		var i int
		for key := range element {
			nonZeroDiags[i] = key
			i++
		}
	case []int:
		nonZeroDiags = element
	}

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

	rotN1 = []int{}
	for i := range rotN1Map {
		rotN1 = append(rotN1, i)
	}

	rotN2 = []int{}
	for i := range rotN2Map {
		rotN2 = append(rotN2, i)
	}

	return
}

func interfaceMapToMapOfInterface(m interface{}) map[int]interface{} {
	d := make(map[int]interface{})
	switch el := m.(type) {
	case map[int][]complex128:
		for i := range el {
			d[i] = el[i]
		}
	case map[int][]float64:
		for i := range el {
			d[i] = el[i]
		}
	default:
		panic("cannot interfaceMapToMapOfInterface: invalid input, must be map[int][]complex128 or map[int][]float64")
	}
	return d
}

// FindBestBSGSRatio finds the best N1*N2 = N for the baby-step giant-step algorithm for matrix multiplication.
func FindBestBSGSRatio(diagMatrix interface{}, maxN int, logMaxRatio int) (minN int) {

	maxRatio := float64(int(1 << logMaxRatio))

	for N1 := 1; N1 < maxN; N1 <<= 1 {

		_, rotN1, rotN2 := BSGSIndex(diagMatrix, maxN, N1)

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

// LinearTransformNew evaluates a linear transform on the Ciphertext "ctIn" and returns the result on a new Ciphertext.
// The linearTransform can either be an (ordered) list of PtDiagMatrix or a single PtDiagMatrix.
// In either case, a list of Ciphertext is returned (the second case returning a list
// containing a single Ciphertext). A PtDiagMatrix is a diagonalized plaintext matrix constructed with an Encoder using
// the method encoder.EncodeDiagMatrixAtLvl(*).
func (eval *evaluator) LinearTransformNew(ctIn *rlwe.Ciphertext, linearTransform interface{}) (ctOut []*rlwe.Ciphertext) {

	switch LTs := linearTransform.(type) {
	case []LinearTransform:
		ctOut = make([]*rlwe.Ciphertext, len(LTs))

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

			ctOut[i].MetaData = ctIn.MetaData
			ctOut[i].Scale = ctIn.Scale.Mul(LT.Scale)
		}

	case LinearTransform:

		minLevel := utils.Min(LTs.Level, ctIn.Level())
		eval.DecomposeNTT(minLevel, eval.params.MaxLevelP(), eval.params.PCount(), ctIn.Value[1], ctIn.IsNTT, eval.BuffDecompQP)

		ctOut = []*rlwe.Ciphertext{NewCiphertext(eval.params, 1, minLevel)}

		if LTs.N1 == 0 {
			eval.MultiplyByDiagMatrix(ctIn, LTs, eval.BuffDecompQP, ctOut[0])
		} else {
			eval.MultiplyByDiagMatrixBSGS(ctIn, LTs, eval.BuffDecompQP, ctOut[0])
		}

		ctOut[0].MetaData = ctIn.MetaData
		ctOut[0].Scale = ctIn.Scale.Mul(LTs.Scale)
	}
	return
}

// LinearTransform evaluates a linear transform on the pre-allocated Ciphertexts.
// The linearTransform can either be an (ordered) list of PtDiagMatrix or a single PtDiagMatrix.
// In either case a list of Ciphertext is returned (the second case returning a list
// containing a single Ciphertext). A PtDiagMatrix is a diagonalized plaintext matrix constructed with an Encoder using
// the method encoder.EncodeDiagMatrixAtLvl(*).
func (eval *evaluator) LinearTransform(ctIn *rlwe.Ciphertext, linearTransform interface{}, ctOut []*rlwe.Ciphertext) {

	switch LTs := linearTransform.(type) {
	case []LinearTransform:
		var maxLevel int
		for _, LT := range LTs {
			maxLevel = utils.Max(maxLevel, LT.Level)
		}

		minLevel := utils.Min(maxLevel, ctIn.Level())
		eval.DecomposeNTT(minLevel, eval.params.MaxLevelP(), eval.params.PCount(), ctIn.Value[1], ctIn.IsNTT, eval.BuffDecompQP)

		for i, LT := range LTs {
			if LT.N1 == 0 {
				eval.MultiplyByDiagMatrix(ctIn, LT, eval.BuffDecompQP, ctOut[i])
			} else {
				eval.MultiplyByDiagMatrixBSGS(ctIn, LT, eval.BuffDecompQP, ctOut[i])
			}

			ctOut[i].MetaData = ctIn.MetaData
			ctOut[i].Scale = ctIn.Scale.Mul(LT.Scale)
		}

	case LinearTransform:
		minLevel := utils.Min(LTs.Level, ctIn.Level())
		eval.DecomposeNTT(minLevel, eval.params.MaxLevelP(), eval.params.PCount(), ctIn.Value[1], ctIn.IsNTT, eval.BuffDecompQP)
		if LTs.N1 == 0 {
			eval.MultiplyByDiagMatrix(ctIn, LTs, eval.BuffDecompQP, ctOut[0])
		} else {
			eval.MultiplyByDiagMatrixBSGS(ctIn, LTs, eval.BuffDecompQP, ctOut[0])
		}

		ctOut[0].MetaData = ctIn.MetaData
		ctOut[0].Scale = ctIn.Scale.Mul(LTs.Scale)
	}
}

// MultiplyByDiagMatrix multiplies the Ciphertext "ctIn" by the plaintext matrix "matrix" and returns the result on the Ciphertext
// "ctOut". Memory buffers for the decomposed ciphertext BuffDecompQP, BuffDecompQP must be provided, those are list of poly of ringQ and ringP
// respectively, each of size params.Beta().
// The naive approach is used (single hoisting and no baby-step giant-step), which is faster than MultiplyByDiagMatrixBSGS
// for matrix of only a few non-zero diagonals but uses more keys.
func (eval *evaluator) MultiplyByDiagMatrix(ctIn *rlwe.Ciphertext, matrix LinearTransform, BuffDecompQP []ringqp.Poly, ctOut *rlwe.Ciphertext) {

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
	ksRes0QP := eval.BuffQP[3]
	ksRes1QP := eval.BuffQP[4]

	ksRes := &rlwe.OperandQP{}
	ksRes.Value = []*ringqp.Poly{
		&eval.BuffQP[3],
		&eval.BuffQP[4],
	}
	ksRes.MetaData.IsNTT = true

	ring.Copy(ctIn.Value[0], eval.buffCt.Value[0])
	ring.Copy(ctIn.Value[1], eval.buffCt.Value[1])
	ctInTmp0, ctInTmp1 := eval.buffCt.Value[0], eval.buffCt.Value[1]

	ringQ.MulScalarBigint(ctInTmp0, ringP.Modulus(), ct0TimesP) // P*c0

	var state bool
	var cnt int
	for k := range matrix.Vec {

		k &= int((ringQ.NthRoot() >> 2) - 1)

		if k == 0 {
			state = true
		} else {

			galEl := eval.params.GaloisElementForColumnRotationBy(k)

			var evk *rlwe.GaloisKey
			var err error
			if evk, err = eval.CheckAndGetGaloisKey(galEl); err != nil {
				panic(fmt.Errorf("cannot apply Automorphism: %w", err))
			}

			index := eval.AutomorphismIndex[galEl]

			eval.GadgetProductHoistedLazy(levelQ, BuffDecompQP, &evk.GadgetCiphertext, ksRes)
			ringQ.Add(ksRes0QP.Q, ct0TimesP, ksRes0QP.Q)

			ringQP.AutomorphismNTTWithIndex(&ksRes0QP, index, &tmp0QP)
			ringQP.AutomorphismNTTWithIndex(&ksRes1QP, index, &tmp1QP)

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
func (eval *evaluator) MultiplyByDiagMatrixBSGS(ctIn *rlwe.Ciphertext, matrix LinearTransform, PoolDecompQP []ringqp.Poly, ctOut *rlwe.Ciphertext) {

	levelQ := utils.Min(ctOut.Level(), utils.Min(ctIn.Level(), matrix.Level))
	levelP := eval.params.RingP().MaxLevel()

	ringQP := eval.params.RingQP().AtLevel(levelQ, levelP)

	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	ctOut.Resize(ctOut.Degree(), levelQ)

	QiOverF := eval.params.QiOverflowMargin(levelQ) >> 1
	PiOverF := eval.params.PiOverflowMargin(levelP) >> 1

	// Computes the N2 rotations indexes of the non-zero rows of the diagonalized DFT matrix for the baby-step giant-step algorithm
	index, _, rotN2 := BSGSIndex(matrix.Vec, 1<<matrix.LogSlots, matrix.N1)

	ring.Copy(ctIn.Value[0], eval.buffCt.Value[0])
	ring.Copy(ctIn.Value[1], eval.buffCt.Value[1])

	ctInTmp0, ctInTmp1 := eval.buffCt.Value[0], eval.buffCt.Value[1]

	// Pre-rotates ciphertext for the baby-step giant-step algorithm, does not divide by P yet
	ctInRotQP := eval.RotateHoistedLazyNew(levelQ, rotN2, eval.buffCt, eval.BuffDecompQP)

	// Accumulator inner loop
	tmp0QP := eval.BuffQP[1]
	tmp1QP := eval.BuffQP[2]

	// Accumulator outer loop
	cQP := &rlwe.OperandQP{}
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

			// Hoisting of the Lazy of sum(sum(phi(d1) * plaintext))
			eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, tmp1QP.Q, tmp1QP.P, tmp1QP.Q) // c1 * plaintext + sum(phi(d1) * plaintext) + phi(c1) * plaintext mod Q

			galEl := eval.params.GaloisElementForColumnRotationBy(j)

			var evk *rlwe.GaloisKey
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
