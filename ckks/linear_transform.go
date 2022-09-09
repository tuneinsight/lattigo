package ckks

import (
	"encoding/binary"
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// Trace maps X -> sum((-1)^i * X^{i*n+1}) for 0 <= i < N
// For log(n) = logSlots.
// Monomial X^k vanishes if k is not divisible by (N/n), otherwise it is multiplied by (N/n).
// Ciphertext is pre-multiplied by (N/n)^-1 to remove the (N/n) factor.
// Examples of full Trace for [0 + 1X + 2X^2 + 3X^3 + 4X^4 + 5X^5 + 6X^6 + 7X^7]
//
// 1)	   [1 + 2X + 3X^2 + 4X^3 + 5X^4 + 6X^5 + 7X^6 + 8X^7]
//       + [1 - 6X - 3X^2 + 8X^3 + 5X^4 + 2X^5 - 7X^6 - 4X^7]  {X-> X^(i * 5^1)}
//   	 = [2 - 4X + 0X^2 +12X^3 +10X^4 + 8X^5 - 0X^6 + 4X^7]
//
// 2)      [2 - 4X + 0X^2 +12X^3 +10X^4 + 8X^5 - 0X^6 + 4X^7]
//       + [2 + 4X + 0X^2 -12X^3 +10X^4 - 8X^5 + 0X^6 - 4X^7]  {X-> X^(i * 5^2)}
//       = [4 + 0X + 0X^2 - 0X^3 +20X^4 + 0X^5 + 0X^6 - 0X^7]
//
// 3)      [4 + 0X + 0X^2 - 0X^3 +20X^4 + 0X^5 + 0X^6 - 0X^7]
//       + [4 + 4X + 0X^2 - 0X^3 -20X^4 + 0X^5 + 0X^6 - 0X^7]  {X-> X^(i * -1)}
//       = [8 + 0X + 0X^2 - 0X^3 + 0X^4 + 0X^5 + 0X^6 - 0X^7]
func (eval *evaluator) Trace(ctIn *Ciphertext, logSlots int, ctOut *Ciphertext) {
	eval.Evaluator.Trace(ctIn.Ciphertext, logSlots, ctOut.Ciphertext)
	ctOut.scale = ctIn.scale
}

// TraceNew maps X -> sum((-1)^i * X^{i*n+1}) for 0 <= i < N and returns the result on a new ciphertext.
// For log(n) = logSlots.
func (eval *evaluator) TraceNew(ctIn *Ciphertext, logSlots int) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1, ctIn.Level(), ctIn.scale)
	eval.Trace(ctIn, logSlots, ctOut)
	return
}

// RotateHoistedNew takes an input Ciphertext and a list of rotations and returns a map of Ciphertext, where each element of the map is the input Ciphertext
// rotation by one element of the list. It is much faster than sequential calls to Rotate.
func (eval *evaluator) RotateHoistedNew(ctIn *Ciphertext, rotations []int) (ctOut map[int]*Ciphertext) {
	ctOut = make(map[int]*Ciphertext)
	for _, i := range rotations {
		ctOut[i] = NewCiphertext(eval.params, 1, ctIn.Level(), ctIn.scale)
	}
	eval.RotateHoisted(ctIn, rotations, ctOut)
	return
}

// RotateHoisted takes an input Ciphertext and a list of rotations and populates a map of pre-allocated Ciphertexts,
// where each element of the map is the input Ciphertext rotation by one element of the list.
// It is much faster than sequential calls to Rotate.
func (eval *evaluator) RotateHoisted(ctIn *Ciphertext, rotations []int, ctOut map[int]*Ciphertext) {
	levelQ := ctIn.Level()
	eval.DecomposeNTT(levelQ, eval.params.PCount()-1, eval.params.PCount(), ctIn.Value[1], eval.BuffDecompQP)
	for _, i := range rotations {
		eval.AutomorphismHoisted(levelQ, ctIn.Ciphertext, eval.BuffDecompQP, eval.params.GaloisElementForColumnRotationBy(i), ctOut[i].Ciphertext)
		ctOut[i].scale = ctIn.scale
	}
}

// LinearTransform is a type for linear transformations on ciphertexts.
// It stores a plaintext matrix diagonalized in diagonal form and
// can be evaluated on a ciphertext by using the evaluator.LinearTransform method.
type LinearTransform struct {
	rlwe.LinearTransform
	Scale float64 // Scale is the scale at which the matrix is encoded (can be circuit dependent)
}

func (LT *LinearTransform) IsPlaintext() bool {
	for _, el := range LT.Vec {
		switch el.(type) {
		case *rlwe.PlaintextQP:
			return true
		default:
			return false
		}
	}
	return false
}

// GetDataLen64 returns the size in bytes of the target LinearTransform when
// encoded on a slice of bytes using MarshalBinary.
func (LT *LinearTransform) GetDataLen64(WithMetaData bool) (dataLen int) {
	dataLen += 4*8 + 1

	var opQPLen int
	for i := range LT.Vec {
		opQPLen = LT.Vec[i].GetDataLen64(WithMetaData)
		break
	}

	return dataLen + len(LT.Vec)*(opQPLen+8)
}

// MarshalBinary encodes the target LinearTransform on a slice of bytes.
func (LT *LinearTransform) MarshalBinary() (data []byte, err error) {

	data = make([]byte, LT.GetDataLen64(true))

	var ptr int

	binary.LittleEndian.PutUint64(data[ptr:], uint64(LT.LogSlots))
	ptr += 8
	binary.LittleEndian.PutUint64(data[ptr:], uint64(LT.N1))
	ptr += 8
	binary.LittleEndian.PutUint64(data[ptr:], uint64(LT.Level))
	ptr += 8
	binary.LittleEndian.PutUint64(data[ptr:], math.Float64bits(LT.Scale))
	ptr += 8

	for _, diag := range LT.Vec {
		if diag.Degree() > 0 {
			data[ptr] = 1
		}
		break
	}
	ptr++

	var inc int
	for rot, diag := range LT.Vec {

		binary.BigEndian.PutUint64(data[ptr:], uint64(rot))
		ptr += 8

		if inc, err = diag.WriteTo64(data[ptr:]); err != nil {
			return nil, err
		}

		ptr += inc
	}

	return
}

// UnmarshalBinary decodes the input slice of bytes on the target LinearTransform.
func (LT *LinearTransform) UnmarshalBinary(data []byte) (err error) {

	var ptr int

	LT.LogSlots = int(binary.LittleEndian.Uint64(data[ptr:]))
	ptr += 8
	LT.N1 = int(binary.LittleEndian.Uint64(data[ptr:]))
	ptr += 8
	LT.Level = int(binary.LittleEndian.Uint64(data[ptr:]))
	ptr += 8
	LT.Scale = math.Float64frombits(binary.LittleEndian.Uint64(data[ptr:]))
	ptr += 8

	LT.Vec = make(map[int]rlwe.OperandQP)

	encrypted := data[ptr] == 1
	ptr++

	var inc int
	for ptr < len(data) {

		rot := int(binary.BigEndian.Uint64(data[ptr:]))
		ptr += 8

		var diag rlwe.OperandQP
		if encrypted {
			diag = &rlwe.CiphertextQP{}
		} else {
			diag = &rlwe.PlaintextQP{}
		}

		if inc, err = diag.Decode64(data[ptr:]); err != nil {
			return
		}
		ptr += inc

		LT.Vec[rot] = diag
	}

	if ptr != len(data) {
		return fmt.Errorf("remaining unparsed data")
	}

	return
}

// NewLinearTransform allocates a new LinearTransform with zero plaintexts at the specified level.
// If BSGSRatio == 0, the LinearTransform is set to not use the BSGS approach.
// Method will panic if BSGSRatio < 0.
func NewLinearTransform(params Parameters, nonZeroDiags []int, level, logSlots int, BSGSRatio float64) LinearTransform {
	vec := make(map[int]rlwe.OperandQP)
	slots := 1 << logSlots
	levelQ := level
	levelP := params.PCount() - 1
	var N1 int
	if BSGSRatio == 0 {
		N1 = 0
		for _, i := range nonZeroDiags {
			idx := i
			if idx < 0 {
				idx += slots
			}
			vec[idx] = &rlwe.PlaintextQP{Value: params.RingQP().NewPolyLvl(levelQ, levelP)}
		}
	} else if BSGSRatio > 0 {
		N1 = rlwe.FindBestBSGSSplit(nonZeroDiags, slots, BSGSRatio)
		index, _, _ := rlwe.BsgsIndex(nonZeroDiags, slots, N1)
		for j := range index {
			for _, i := range index[j] {
				vec[j+i] = &rlwe.PlaintextQP{Value: params.RingQP().NewPolyLvl(levelQ, levelP)}
			}
		}
	} else {
		panic("BSGS ratio cannot be negative")
	}

	return LinearTransform{
		LinearTransform: rlwe.LinearTransform{
			LogSlots: logSlots,
			N1:       N1,
			Vec:      vec,
			Level:    level,
		},
	}
}

// Rotations returns the list of rotations needed for the evaluation
// of the linear transform.
func (LT *LinearTransform) Rotations() (rotations []int) {
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

	rotations = make([]int, len(rotIndex))
	var i int
	for j := range rotIndex {
		rotations[i] = j
		i++
	}

	return rotations
}

// Encode encodes on a pre-allocated LinearTransform the linear transforms' matrix in diagonal form `value`.
// values.(type) can be either map[int][]complex128 or map[int][]float64.
// User must ensure that 1 <= len([]complex128\[]float64) <= 2^logSlots < 2^logN.
// It can then be evaluated on a ciphertext using evaluator.LinearTransform.
// Evaluation will use the naive approach (single hoisting and no baby-step giant-step).
// Faster if there is only a few non-zero diagonals but uses more keys.
func (LT *LinearTransform) Encode(encoder Encoder, value interface{}, scale float64) {

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

			enc.Embed(dMat[i], LT.LogSlots, scale, true, LT.Vec[idx].El().Value[0])
		}
	} else {
		index, _, _ := rlwe.BsgsIndex(value, slots, N1)

		for j := range index {
			for _, i := range index[j] {
				// manages inputs that have rotation between 0 and slots-1 or between -slots/2 and slots/2-1
				v, ok := dMat[j+i]
				if !ok {
					v = dMat[j+i-slots]
				}

				if _, ok := LT.Vec[j+i]; !ok {
					panic("cannot Encode: error encoding on LinearTransform BSGS: input does not match the same non-zero diagonals")
				}

				enc.Embed(utils.RotateSlice(v, -j), LT.LogSlots, scale, true, LT.Vec[j+i].El().Value[0])
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
func GenLinearTransform(encoder Encoder, value interface{}, level int, scale float64, logslots int) LinearTransform {

	enc, ok := encoder.(*encoderComplex128)
	if !ok {
		panic("cannot GenLinearTransform: encoder should be an encoderComplex128")
	}

	params := enc.params
	dMat := interfaceMapToMapOfInterface(value)
	vec := make(map[int]rlwe.OperandQP)
	slots := 1 << logslots
	levelQ := level
	levelP := params.PCount() - 1
	for i := range dMat {

		idx := i
		if idx < 0 {
			idx += slots
		}
		vec[idx] = &rlwe.PlaintextQP{Value: params.RingQP().NewPolyLvl(levelQ, levelP)}
		enc.Embed(dMat[i], logslots, scale, true, vec[idx].El().Value[0])
	}

	return LinearTransform{
		LinearTransform: rlwe.LinearTransform{
			N1:    0,
			Vec:   vec,
			Level: level,
		},
		Scale: scale,
	}
}

// GenLinearTransformEncrypted allocates and encrypts a new LinearTransform struct from the linear transforms' matrix in diagonal form `value`.
// values.(type) can be either map[int][]complex128 or map[int][]float64.
// User must ensure that 1 <= len([]complex128\[]float64) <= 2^logSlots < 2^logN.
// It can then be evaluated on a ciphertext using evaluator.LinearTransform.
// Evaluation will use the naive approach (single hoisting and no baby-step giant-step).
// Faster if there is only a few non-zero diagonals but uses more keys.
func GenLinearTransformEncrypted(encoder Encoder, enc Encryptor, value interface{}, level int, scale float64, logslots int) LinearTransform {

	ecd, ok := encoder.(*encoderComplex128)
	if !ok {
		panic("encoder should be an encoderComplex128")
	}

	params := ecd.params
	dMat := interfaceMapToMapOfInterface(value)
	vec := make(map[int]rlwe.OperandQP)
	slots := 1 << logslots
	levelQ := level
	levelP := params.PCount() - 1
	pt := params.RingQP().NewPolyLvl(level, params.PCount()-1)
	for i := range dMat {

		idx := i
		if idx < 0 {
			idx += slots
		}

		ct := &rlwe.CiphertextQP{
			Value: []ringqp.Poly{
				params.RingQP().NewPolyLvl(levelQ, levelP),
				params.RingQP().NewPolyLvl(levelQ, levelP),
			},
		}

		enc.GetRLWEEncryptor().EncryptZero(ct)

		ecd.Embed(dMat[i], logslots, scale, true, pt)
		params.RingQP().AddLvl(levelQ, levelP, ct.Value[0], pt, ct.Value[0])

		vec[idx] = ct
	}

	return LinearTransform{
		LinearTransform: rlwe.LinearTransform{
			N1:    0,
			Vec:   vec,
			Level: level,
		},
		Scale: scale,
	}
}

// GenLinearTransformBSGS allocates and encodes a new LinearTransform struct from the linear transforms' matrix in diagonal form `value` for evaluation with a baby-step giant-step approach.
// values.(type) can be either map[int][]complex128 or map[int][]float64.
// User must ensure that 1 <= len([]complex128\[]float64) <= 2^logSlots < 2^logN.
// LinearTransform types can be be evaluated on a ciphertext using evaluator.LinearTransform.
// Evaluation will use the optimized approach (double hoisting and baby-step giant-step).
// Faster if there is more than a few non-zero diagonals.
// BSGSRatio is the maximum ratio between the inner and outer loop of the baby-step giant-step algorithm used in evaluator.LinearTransform.
// Optimal BSGSRatio value is between 4 and 16 depending on the sparsity of the matrix.
func GenLinearTransformBSGS(encoder Encoder, value interface{}, level int, scale, BSGSRatio float64, logSlots int) (LT LinearTransform) {

	enc, ok := encoder.(*encoderComplex128)
	if !ok {
		panic("cannot GenLinearTransformBSGS: encoder should be an encoderComplex128")
	}

	params := enc.params

	slots := 1 << logSlots

	// N1*N2 = N
	N1 := rlwe.FindBestBSGSSplit(value, slots, BSGSRatio)

	index, _, _ := rlwe.BsgsIndex(value, slots, N1)

	vec := make(map[int]rlwe.OperandQP)

	dMat := interfaceMapToMapOfInterface(value)
	levelQ := level
	levelP := params.PCount() - 1
	for j := range index {

		for _, i := range index[j] {

			// manages inputs that have rotation between 0 and slots-1 or between -slots/2 and slots/2-1
			v, ok := dMat[j+i]
			if !ok {
				v = dMat[j+i-slots]
			}
			vec[j+i] = &rlwe.PlaintextQP{Value: params.RingQP().NewPolyLvl(levelQ, levelP)}
			enc.Embed(utils.RotateSlice(v, -j), logSlots, scale, true, vec[j+i].El().Value[0])
		}
	}

	return LinearTransform{
		LinearTransform: rlwe.LinearTransform{
			LogSlots: logSlots,
			N1:       N1,
			Vec:      vec,
			Level:    level,
		},
		Scale: scale,
	}
}

// GenLinearTransformBSGSEncrypted allocates and encrypts a new LinearTransform struct from the linear transforms' matrix in diagonal form `value` for evaluation with a baby-step giant-step approach.
// values.(type) can be either map[int][]complex128 or map[int][]float64.
// User must ensure that 1 <= len([]complex128\[]float64) <= 2^logSlots < 2^logN.
// LinearTransform types can be be evaluated on a ciphertext using evaluator.LinearTransform.
// Evaluation will use the optimized approach (double hoisting and baby-step giant-step).
// Faster if there is more than a few non-zero diagonals.
// BSGSRatio is the maximum ratio between the inner and outer loop of the baby-step giant-step algorithm used in evaluator.LinearTransform.
// Optimal BSGSRatio value is between 4 and 16 depending on the sparsity of the matrix.
func GenLinearTransformBSGSEncrypted(encoder Encoder, enc Encryptor, value interface{}, level int, scale, BSGSRatio float64, logSlots int) (LT LinearTransform) {

	ecd, ok := encoder.(*encoderComplex128)
	if !ok {
		panic("cannot GenLinearTransformBSGS: encoder should be an encoderComplex128")
	}

	params := ecd.params

	slots := 1 << logSlots

	// N1*N2 = N
	N1 := rlwe.FindBestBSGSSplit(value, slots, BSGSRatio)

	index, _, _ := rlwe.BsgsIndex(value, slots, N1)

	vec := make(map[int]rlwe.OperandQP)

	dMat := interfaceMapToMapOfInterface(value)
	levelQ := level
	levelP := params.PCount() - 1
	pt := params.RingQP().NewPolyLvl(level, params.PCount()-1)
	for j := range index {

		for _, i := range index[j] {

			// manages inputs that have rotation between 0 and slots-1 or between -slots/2 and slots/2-1
			v, ok := dMat[j+i]
			if !ok {
				v = dMat[j+i-slots]
			}

			ct := &rlwe.CiphertextQP{
				Value: []ringqp.Poly{
					params.RingQP().NewPolyLvl(levelQ, levelP),
					params.RingQP().NewPolyLvl(levelQ, levelP),
				},
			}

			ct.Value[0].Q.IsNTT = true
			ct.Value[0].P.IsNTT = true
			ct.Value[1].Q.IsNTT = true
			ct.Value[1].P.IsNTT = true

			enc.GetRLWEEncryptor().EncryptZero(ct)

			ecd.Embed(utils.RotateSlice(v, -j), logSlots, scale, true, pt)
			params.RingQP().AddLvl(levelQ, levelP, ct.Value[0], pt, ct.Value[0])

			vec[j+i] = ct
		}
	}

	return LinearTransform{
		LinearTransform: rlwe.LinearTransform{
			LogSlots: logSlots,
			N1:       N1,
			Vec:      vec,
			Level:    level,
		},
		Scale: scale,
	}
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

// LinearTransformNew evaluates a linear transform on the ciphertext and returns the result on a new ciphertext.
// The linearTransform can either be an (ordered) list of PtDiagMatrix or a single PtDiagMatrix.
// In either case a list of ciphertext is returned (the second case returning a list of
// containing a single ciphertext. A PtDiagMatrix is a diagonalized plaintext matrix constructed with an Encoder using
// the method encoder.EncodeDiagMatrixAtLvl(*).
func (eval *evaluator) LinearTransformNew(ctIn *Ciphertext, linearTransform interface{}) (ctOut []*Ciphertext) {

	switch LTs := linearTransform.(type) {
	case []LinearTransform:
		ctOut = make([]*Ciphertext, len(LTs))

		var maxLevel int
		for _, LT := range LTs {
			maxLevel = utils.MaxInt(maxLevel, LT.Level)
		}

		minLevel := utils.MinInt(maxLevel, ctIn.Level())

		for i := range LTs {
			ctOut[i] = NewCiphertext(eval.params, 1, minLevel, ctIn.scale)
		}

	case LinearTransform:
		ctOut = []*Ciphertext{NewCiphertext(eval.params, 1, utils.MinInt(LTs.Level, ctIn.Level()), ctIn.scale)}
	}

	eval.LinearTransform(ctIn, linearTransform, ctOut)

	return
}

// LinearTransformNew evaluates a linear transform on the pre-allocated ciphertexts.
// The linearTransform can either be an (ordered) list of PtDiagMatrix or a single PtDiagMatrix.
// In either case a list of ciphertext is returned (the second case returning a list of
// containing a single ciphertext. A PtDiagMatrix is a diagonalized plaintext matrix constructed with an Encoder using
// the method encoder.EncodeDiagMatrixAtLvl(*).
func (eval *evaluator) LinearTransform(ctIn *Ciphertext, linearTransform interface{}, ctOut []*Ciphertext) {

	switch LTs := linearTransform.(type) {
	case []LinearTransform:
		var maxLevel int
		for _, LT := range LTs {
			maxLevel = utils.MaxInt(maxLevel, LT.Level)
		}

		minLevel := utils.MinInt(maxLevel, ctIn.Level())

		eval.DecomposeNTT(minLevel, eval.params.PCount()-1, eval.params.PCount(), ctIn.Value[1], eval.BuffDecompQP)

		for i, LT := range LTs {
			if LT.N1 == 0 {
				eval.MultiplyByDiagMatrix(ctIn.Ciphertext, LT.LinearTransform, eval.BuffDecompQP, ctOut[i].Ciphertext)
			} else {
				eval.MultiplyByDiagMatrixBSGS(ctIn.Ciphertext, LT.LinearTransform, eval.BuffDecompQP, ctOut[i].Ciphertext)
			}

			ctOut[i].SetScale(ctOut[i].Scale() * LT.Scale)
		}

	case LinearTransform:
		minLevel := utils.MinInt(LTs.Level, ctIn.Level())
		eval.DecomposeNTT(minLevel, eval.params.PCount()-1, eval.params.PCount(), ctIn.Value[1], eval.BuffDecompQP)
		if LTs.N1 == 0 {
			eval.MultiplyByDiagMatrix(ctIn.Ciphertext, LTs.LinearTransform, eval.BuffDecompQP, ctOut[0].Ciphertext)
		} else {
			eval.MultiplyByDiagMatrixBSGS(ctIn.Ciphertext, LTs.LinearTransform, eval.BuffDecompQP, ctOut[0].Ciphertext)
		}

		ctOut[0].SetScale(ctOut[0].Scale() * LTs.Scale)
	}
}

// Average returns the average of vectors of batchSize elements.
// The operation assumes that ctIn encrypts SlotCount/'batchSize' sub-vectors of size 'batchSize'.
// It then replaces all values of those sub-vectors by the component-wise average between all the sub-vectors.
// Example for batchSize=4 and slots=8: [{a, b, c, d}, {e, f, g, h}] -> [0.5*{a+e, b+f, c+g, d+h}, 0.5*{a+e, b+f, c+g, d+h}]
// Operation requires log2(SlotCout/'batchSize') rotations.
// Required rotation keys can be generated with 'RotationsForInnerSumLog(batchSize, SlotCount/batchSize)''
func (eval *evaluator) Average(ctIn *Ciphertext, logBatchSize int, ctOut *Ciphertext) {

	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		panic("ctIn.Degree() != 1 or ctOut.Degree() != 1")
	}

	if logBatchSize > eval.params.LogSlots() {
		panic("cannot Average: batchSize must be smaller or equal to the number of slots")
	}

	ringQ := eval.params.RingQ()

	level := utils.MinInt(ctIn.Level(), ctOut.Level())

	n := eval.params.Slots() / (1 << logBatchSize)

	// pre-multiplication by n^-1
	for i := 0; i < level+1; i++ {
		Q := ringQ.Modulus[i]
		bredParams := ringQ.BredParams[i]
		mredparams := ringQ.MredParams[i]
		invN := ring.ModExp(uint64(n), Q-2, Q)
		invN = ring.MForm(invN, Q, bredParams)

		ring.MulScalarMontgomeryVec(ctIn.Value[0].Coeffs[i], ctOut.Value[0].Coeffs[i], invN, Q, mredparams)
		ring.MulScalarMontgomeryVec(ctIn.Value[1].Coeffs[i], ctOut.Value[1].Coeffs[i], invN, Q, mredparams)
	}

	eval.InnerSumLog(ctOut, 1<<logBatchSize, n, ctOut)
}

// InnerSumLog applies an optimized inner sum on the ciphertext (log2(n) + HW(n) rotations with double hoisting).
// The operation assumes that `ctIn` encrypts SlotCount/`batchSize` sub-vectors of size `batchSize` which it adds together (in parallel) by groups of `n`.
// It outputs in ctOut a ciphertext for which the "leftmost" sub-vector of each group is equal to the sum of the group.
// This method is faster than InnerSum when the number of rotations is large and uses log2(n) + HW(n) instead of 'n' keys.
func (eval *evaluator) InnerSumLog(ctIn *Ciphertext, batchSize, n int, ctOut *Ciphertext) {

	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		panic("ctIn.Degree() != 1 or ctOut.Degree() != 1")
	}

	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()
	ringQP := eval.params.RingQP()

	levelQ := ctIn.Level()
	levelP := len(ringP.Modulus) - 1

	ctOut.Resize(ctOut.Degree(), levelQ)
	ctOut.scale = ctIn.scale

	if n == 1 {
		if ctIn != ctOut {
			ring.CopyValuesLvl(levelQ, ctIn.Value[0], ctOut.Value[0])
			ring.CopyValuesLvl(levelQ, ctIn.Value[1], ctOut.Value[1])
		}
	} else {

		// Memory buffer for ctIn = ctIn + rot(ctIn, 2^i) in Q
		tmpc0 := eval.buffQ[0] // unused memory buffer from evaluator
		tmpc1 := eval.buffQ[1] // unused memory buffer from evaluator
		tmpc0.IsNTT = true
		tmpc1.IsNTT = true
		tmpct := NewCiphertextAtLevelFromPoly(levelQ, [2]*ring.Poly{tmpc0, tmpc1})
		c0OutQP := eval.BuffQP[2]
		c1OutQP := eval.BuffQP[3]

		accQP := rlwe.CiphertextQP{Value: []ringqp.Poly{eval.BuffQP[4], eval.BuffQP[5]}}
		accQP.Value[0].Q.IsNTT = true
		accQP.Value[1].Q.IsNTT = true
		ctqp := NewCiphertextAtLevelFromPoly(levelQ, [2]*ring.Poly{accQP.Value[0].Q, accQP.Value[1].Q})

		state := false
		copy := true
		// Binary reading of the input n
		for i, j := 0, n; j > 0; i, j = i+1, j>>1 {

			// Starts by decomposing the input ciphertext
			if i == 0 {
				// If first iteration, then copies directly from the input ciphertext that hasn't been rotated
				eval.DecomposeNTT(levelQ, levelP, levelP+1, ctIn.Value[1], eval.BuffDecompQP)
			} else {
				// Else copies from the rotated input ciphertext
				eval.DecomposeNTT(levelQ, levelP, levelP+1, tmpc1, eval.BuffDecompQP)
			}

			// If the binary reading scans a 1
			if j&1 == 1 {

				k := n - (n & ((2 << i) - 1))
				k *= batchSize

				// If the rotation is not zero
				if k != 0 {

					// Rotate((tmpc0, tmpc1), k)
					if i == 0 {
						eval.AutomorphismHoistedNoModDown(levelQ, ctIn.Value[0], eval.BuffDecompQP, eval.params.GaloisElementForColumnRotationBy(k), accQP)
					} else {
						eval.AutomorphismHoistedNoModDown(levelQ, tmpc0, eval.BuffDecompQP, eval.params.GaloisElementForColumnRotationBy(k), accQP)
					}

					// ctOut += Rotate((tmpc0, tmpc1), k)
					if copy {
						ringQP.CopyValuesLvl(levelQ, levelP, accQP.Value[0], c0OutQP)
						ringQP.CopyValuesLvl(levelQ, levelP, accQP.Value[1], c1OutQP)
						copy = false
					} else {
						ringQP.AddLvl(levelQ, levelP, c0OutQP, accQP.Value[0], c0OutQP)
						ringQP.AddLvl(levelQ, levelP, c1OutQP, accQP.Value[1], c1OutQP)
					}
				} else {

					state = true

					// if n is not a power of two
					if n&(n-1) != 0 {
						eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c0OutQP.Q, c0OutQP.P, c0OutQP.Q) // Division by P
						eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c1OutQP.Q, c1OutQP.P, c1OutQP.Q) // Division by P

						// ctOut += (tmpc0, tmpc1)
						ringQ.AddLvl(levelQ, c0OutQP.Q, tmpc0, ctOut.Value[0])
						ringQ.AddLvl(levelQ, c1OutQP.Q, tmpc1, ctOut.Value[1])

					} else {
						ring.CopyValuesLvl(levelQ, tmpc0, ctOut.Value[0])
						ring.CopyValuesLvl(levelQ, tmpc1, ctOut.Value[1])
						ctOut.Value[0].IsNTT = true
						ctOut.Value[1].IsNTT = true
					}
				}
			}

			if !state {

				rot := eval.params.GaloisElementForColumnRotationBy((1 << i) * batchSize)
				if i == 0 {
					eval.AutomorphismHoisted(levelQ, ctIn.Ciphertext, eval.BuffDecompQP, rot, tmpct.Ciphertext)
					ringQ.AddLvl(levelQ, tmpc0, ctIn.Value[0], tmpc0)
					ringQ.AddLvl(levelQ, tmpc1, ctIn.Value[1], tmpc1)
				} else {

					// (tmpc0, tmpc1) = Rotate((tmpc0, tmpc1), 2^i)
					eval.AutomorphismHoisted(levelQ, tmpct.Ciphertext, eval.BuffDecompQP, rot, ctqp.Ciphertext)
					ringQ.AddLvl(levelQ, tmpc0, accQP.Value[0].Q, tmpc0)
					ringQ.AddLvl(levelQ, tmpc1, accQP.Value[1].Q, tmpc1)
				}
			}
		}
	}
}

// InnerSum applies an naive inner sum on the ciphertext (n rotations with single hoisting).
// The operation assumes that `ctIn` encrypts SlotCount/`batchSize` sub-vectors of size `batchSize` which it adds together (in parallel) by groups of `n`.
// It outputs in ctOut a ciphertext for which the "leftmost" sub-vector of each group is equal to the sum of the group.
// This method is faster than InnerSumLog when the number of rotations is small but uses 'n' keys instead of log(n) + HW(n).
func (eval *evaluator) InnerSum(ctIn *Ciphertext, batchSize, n int, ctOut *Ciphertext) {

	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		panic("ctIn.Degree() != 1 or ctOut.Degree() != 1")
	}

	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()
	ringQP := eval.params.RingQP()

	levelQ := ctIn.Level()
	levelP := len(ringP.Modulus) - 1

	QiOverF := eval.params.QiOverflowMargin(levelQ) >> 1
	PiOverF := eval.params.PiOverflowMargin(levelP) >> 1

	ctOut.Resize(ctOut.Degree(), levelQ)
	ctOut.scale = ctIn.scale

	// If sum with only the first element, then returns the input
	if n == 1 {

		// If the input-output points differ, copies on the output
		if ctIn != ctOut {
			ring.CopyValuesLvl(levelQ, ctIn.Value[0], ctOut.Value[0])
			ring.CopyValuesLvl(levelQ, ctIn.Value[1], ctOut.Value[1])
		}
		// If sum on at least two elements
	} else {

		// List of n-2 rotations
		rotations := []int{}
		for i := 1; i < n; i++ {
			rotations = append(rotations, i*batchSize)
		}

		// Memory buffer
		tmp0QP := eval.BuffQP[1]
		tmp1QP := eval.BuffQP[2]

		// Basis decomposition
		eval.DecomposeNTT(levelQ, levelP, levelP+1, ctIn.Value[1], eval.BuffDecompQP)

		// Pre-rotates all [1, ..., n-1] rotations
		// Hoisted rotation without division by P
		ctInRotQP := eval.RotateHoistedNoModDownNew(levelQ, rotations, ctIn.Value[0], eval.BuffDecompQP)

		// P*c0 -> tmp0QP.Q
		ringQ.MulScalarBigintLvl(levelQ, ctIn.Value[0], ringP.ModulusAtLevel[levelP], tmp0QP.Q)

		// Adds phi_k(P*c0) on each of the vecRotQ
		// Does not need to add on the vecRotP because mod P === 0
		var reduce int
		// Sums elements [2, ..., n-1]
		for i := 1; i < n; i++ {

			j := i * batchSize

			if i == 1 {
				ringQP.CopyValuesLvl(levelQ, levelP, ctInRotQP[j].Value[0], tmp0QP)
				ringQP.CopyValuesLvl(levelQ, levelP, ctInRotQP[j].Value[1], tmp1QP)
			} else {
				ringQP.AddNoModLvl(levelQ, levelP, tmp0QP, ctInRotQP[j].Value[0], tmp0QP)
				ringQP.AddNoModLvl(levelQ, levelP, tmp1QP, ctInRotQP[j].Value[1], tmp1QP)
			}

			if reduce%QiOverF == QiOverF-1 {
				ringQ.ReduceLvl(levelQ, tmp0QP.Q, tmp0QP.Q)
				ringQ.ReduceLvl(levelQ, tmp1QP.Q, tmp1QP.Q)
			}

			if reduce%PiOverF == PiOverF-1 {
				ringP.ReduceLvl(levelP, tmp0QP.P, tmp0QP.P)
				ringP.ReduceLvl(levelP, tmp1QP.P, tmp1QP.P)
			}

			reduce++
		}

		if reduce%QiOverF != 0 {
			ringQ.ReduceLvl(levelQ, tmp0QP.Q, tmp0QP.Q)
			ringQ.ReduceLvl(levelQ, tmp1QP.Q, tmp1QP.Q)
		}

		if reduce%PiOverF != 0 {
			ringP.ReduceLvl(levelP, tmp0QP.P, tmp0QP.P)
			ringP.ReduceLvl(levelP, tmp1QP.P, tmp1QP.P)
		}

		// Division by P of sum(elements [2, ..., n-1] )
		eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, tmp0QP.Q, tmp0QP.P, tmp0QP.Q) // sum_{i=1, n-1}(phi(d0))/P
		eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, tmp1QP.Q, tmp1QP.P, tmp1QP.Q) // sum_{i=1, n-1}(phi(d1))/P

		// Adds element[1] (which did not require rotation)
		ringQ.AddLvl(levelQ, ctIn.Value[0], tmp0QP.Q, ctOut.Value[0]) // sum_{i=1, n-1}(phi(d0))/P + ct0
		ringQ.AddLvl(levelQ, ctIn.Value[1], tmp1QP.Q, ctOut.Value[1]) // sum_{i=1, n-1}(phi(d1))/P + ct1
	}
}

// ReplicateLog applies an optimized replication on the ciphertext (log2(n) + HW(n) rotations with double hoisting).
// It acts as the inverse of a inner sum (summing elements from left to right).
// The replication is parameterized by the size of the sub-vectors to replicate "batchSize" and
// the number of time "n" they need to be replicated.
// To ensure correctness, a gap of zero values of size batchSize * (n-1) must exist between
// two consecutive sub-vectors to replicate.
// This method is faster than Replicate when the number of rotations is large and uses log2(n) + HW(n) instead of 'n'.
func (eval *evaluator) ReplicateLog(ctIn *Ciphertext, batchSize, n int, ctOut *Ciphertext) {
	eval.InnerSumLog(ctIn, -batchSize, n, ctOut)
}

// Replicate applies naive replication on the ciphertext (n rotations with single hoisting).
// It acts as the inverse of a inner sum (summing elements from left to right).
// The replication is parameterized by the size of the sub-vectors to replicate "batchSize" and
// the number of time "n" they need to be replicated.
// To ensure correctness, a gap of zero values of size batchSize * (n-1) must exist between
// two consecutive sub-vectors to replicate.
// This method is faster than ReplicateLog when the number of rotations is small but uses 'n' keys instead of log2(n) + HW(n).
func (eval *evaluator) Replicate(ctIn *Ciphertext, batchSize, n int, ctOut *Ciphertext) {
	eval.InnerSum(ctIn, -batchSize, n, ctOut)
}
