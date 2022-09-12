package ckks

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// Trace maps X -> sum((-1)^i * X^{i*n+1}) for 0 <= i < N and writes the result on the output ciphertext.
// For log(n) = logSlots.
func (eval *evaluator) Trace(ctIn *Ciphertext, logSlots int, ctOut *Ciphertext) {
	eval.Evaluator.Trace(ctIn.Ciphertext, logSlots, ctOut.Ciphertext)
}

// TraceNew maps X -> sum((-1)^i * X^{i*n+1}) for 0 <= i < N and returns the result on a new ciphertext.
// For log(n) = logSlots.
func (eval *evaluator) TraceNew(ctIn *Ciphertext, logSlots int) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1, ctIn.Level())
	eval.Trace(ctIn, logSlots, ctOut)
	return
}

// RotateHoistedNew takes an input Ciphertext and a list of rotations and returns a map of Ciphertext, where each element of the map is the input Ciphertext
// rotation by one element of the list. It is much faster than sequential calls to Rotate.
func (eval *evaluator) RotateHoistedNew(ctIn *Ciphertext, rotations []int) (ctOut map[int]*Ciphertext) {
	ctOut = make(map[int]*Ciphertext)
	for _, i := range rotations {
		ctOut[i] = NewCiphertext(eval.params, 1, ctIn.Level())
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
		ctOut[i].Scale().Set(ctIn.Scale())
	}
}

// LinearTransform is a type for linear transformations on ciphertexts.
// It stores a plaintext matrix diagonalized in diagonal form and
// can be evaluated on a ciphertext by using the evaluator.LinearTransform method.
type LinearTransform struct {
	rlwe.LinearTransform
}

func NewLinearTransform(params Parameters, nonZeroDiags []int, level, logSlots int, BSGSRatio float64) LinearTransform {
	LT := rlwe.NewLinearTransform(params.Parameters, nonZeroDiags, level, logSlots, BSGSRatio)
	LT.Scale = NewScale(0)
	return LinearTransform{LinearTransform: LT}
}

// Encode encodes on a pre-allocated LinearTransform the linear transforms' matrix in diagonal form `value`.
// values.(type) can be either map[int][]complex128 or map[int][]float64.
// User must ensure that 1 <= len([]complex128\[]float64) <= 2^logSlots < 2^logN.
// It can then be evaluated on a ciphertext using evaluator.LinearTransform.
// Evaluation will use the naive approach (single hoisting and no baby-step giant-step).
// Faster if there is only a few non-zero diagonals but uses more keys.
func (LT *LinearTransform) Encode(ecd Encoder, value interface{}, scale rlwe.Scale) {

	encode := func(value interface{}, opQP rlwe.OperandQP) {
		ecd.Embed(value, LT.LogDimension, scale, true, opQP.El().Value[0])
	}

	LT.LinearTransform.Encode(encode, interfaceMapToMapOfInterface(value))
	LT.Scale.Set(scale)
}

func (LT *LinearTransform) MarshalBinary() (data []byte, err error) {
	return LT.LinearTransform.MarshalBinary()
}

func (LT *LinearTransform) UnmarshalBinary(data []byte) (err error) {
	LT.LinearTransform = rlwe.LinearTransform{}
	LT.Scale = NewScale(0)
	return LT.LinearTransform.UnmarshalBinary(data)
}

// GenLinearTransform allocates and encrypts a new LinearTransform struct from the linear transforms' matrix in diagonal form `value` for evaluation with a baby-step giant-step approach.
// values.(type) can be either map[int][]complex128 or map[int][]float64.
// User must ensure that 1 <= len([]complex128\[]float64) <= 2^logSlots < 2^logN.
// LinearTransform types can be be evaluated on a ciphertext using evaluator.LinearTransform.
// Evaluation will use the optimized approach (double hoisting and baby-step giant-step).
// Faster if there is more than a few non-zero diagonals.
// BSGSRatio is the maximum ratio between the inner and outer loop of the baby-step giant-step algorithm used in evaluator.LinearTransform.
// Optimal BSGSRatio value is between 4 and 16 depending on the sparsity of the matrix.
func GenLinearTransform(ecd Encoder, enc Encryptor, value interface{}, level int, scale rlwe.Scale, BSGSRatio float64, logSlots int) LinearTransform {

	if ecd == nil {
		panic("GenLinearTransformBSGS: ecd cannot be nil")
	}

	params := ecd.Parameters()
	levelQ := level
	levelP := params.PCount() - 1

	var embed func(value interface{}) (opQP rlwe.OperandQP)

	if enc != nil {
		pt := params.RingQP().NewPolyLvl(level, params.PCount()-1)

		embed = func(value interface{}) (opQP rlwe.OperandQP) {
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

			ecd.Embed(value, logSlots, scale, true, pt)
			params.RingQP().AddLvl(levelQ, levelP, ct.Value[0], pt, ct.Value[0])

			return ct
		}

	} else {

		embed = func(value interface{}) (opQP rlwe.OperandQP) {
			opQP = &rlwe.PlaintextQP{Value: params.RingQP().NewPolyLvl(levelQ, levelP)}
			ecd.Embed(value, logSlots, scale, true, opQP.El().Value[0])
			return
		}
	}

	LT := rlwe.GenLinearTransform(embed, interfaceMapToMapOfInterface(value), BSGSRatio, logSlots)
	LT.Scale = scale.CopyNew()

	return LinearTransform{LinearTransform: LT}
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
			ctOut[i] = NewCiphertext(eval.params, 1, minLevel)
		}

	case LinearTransform:
		ctOut = []*Ciphertext{NewCiphertext(eval.params, 1, utils.MinInt(LTs.Level, ctIn.Level()))}
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
		}

	case LinearTransform:
		minLevel := utils.MinInt(LTs.Level, ctIn.Level())
		eval.DecomposeNTT(minLevel, eval.params.PCount()-1, eval.params.PCount(), ctIn.Value[1], eval.BuffDecompQP)
		if LTs.N1 == 0 {
			eval.MultiplyByDiagMatrix(ctIn.Ciphertext, LTs.LinearTransform, eval.BuffDecompQP, ctOut[0].Ciphertext)
		} else {
			eval.MultiplyByDiagMatrixBSGS(ctIn.Ciphertext, LTs.LinearTransform, eval.BuffDecompQP, ctOut[0].Ciphertext)
		}
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

	eval.InnerSum(ctOut, 1<<logBatchSize, n, ctOut)
}

// InnerSumLog applies an optimized inner sum on the ciphertext (log2(n) + HW(n) rotations with double hoisting).
// The operation assumes that `ctIn` encrypts SlotCount/`batchSize` sub-vectors of size `batchSize` which it adds together (in parallel) by groups of `n`.
// It outputs in ctOut a ciphertext for which the "leftmost" sub-vector of each group is equal to the sum of the group.
// This method is faster than InnerSum when the number of rotations is large and uses log2(n) + HW(n) instead of 'n' keys.
func (eval *evaluator) InnerSum(ctIn *Ciphertext, batchSize, n int, ctOut *Ciphertext) {
	eval.Evaluator.InnerSum(ctIn.Ciphertext, batchSize, n, ctOut.Ciphertext)
}

// ReplicateLog applies an optimized replication on the ciphertext (log2(n) + HW(n) rotations with double hoisting).
// It acts as the inverse of a inner sum (summing elements from left to right).
// The replication is parameterized by the size of the sub-vectors to replicate "batchSize" and
// the number of time "n" they need to be replicated.
// To ensure correctness, a gap of zero values of size batchSize * (n-1) must exist between
// two consecutive sub-vectors to replicate.
func (eval *evaluator) Replicate(ctIn *Ciphertext, batchSize, n int, ctOut *Ciphertext) {
	eval.InnerSum(ctIn, -batchSize, n, ctOut)
}
