package bgv

import (
	"encoding/binary"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// InnerSumLog applies an optimized inner sum on the ciphertext (log2(n) + HW(n) rotations with double hoisting).
// The operation assumes that `ctIn` encrypts SlotCount/`batchSize` sub-vectors of size `batchSize` which it adds together (in parallel) by groups of `n`.
// It outputs in ctOut a ciphertext for which the "leftmost" sub-vector of each group is equal to the sum of the group.
// This method is faster than InnerSum when the number of rotations is large and uses log2(n) + HW(n) instead of 'n' keys.
func (eval *evaluator) InnerSumLog(ctIn *Ciphertext, batchSize, n int, ctOut *Ciphertext) {

	ringQ := eval.params.RingQ()
	levelQ := utils.MinInt(ctIn.Level(), ctOut.Level())

	ct1 := rlwe.NewCiphertextAtLevelFromPoly(levelQ, [2]*ring.Poly{eval.buffCt.Value[0], eval.buffCt.Value[1]})
	ct1.Value[0].IsNTT = true
	ct1.Value[1].IsNTT = true

	ringQ.MulScalarBigintLvl(levelQ, ctIn.Value[0], eval.tInvModQ[levelQ], ct1.Value[0])
	ringQ.MulScalarBigintLvl(levelQ, ctIn.Value[1], eval.tInvModQ[levelQ], ct1.Value[1])

	eval.Evaluator.InnerSumLog(ct1, batchSize, n, ctOut.Ciphertext)

	ringQ.MulScalarLvl(levelQ, ctOut.Value[0], eval.params.T(), ctOut.Value[0])
	ringQ.MulScalarLvl(levelQ, ctOut.Value[1], eval.params.T(), ctOut.Value[1])
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

type LinearTransform struct {
	rlwe.LinearTransform
}

// MarshalBinary encodes the target LinearTransform on a slice of bytes.
func (LT *LinearTransform) MarshalBinary() (data []byte, err error) {
	return LT.LinearTransform.MarshalBinary()
}

// UnmarshalBinary decodes the input slice of bytes on the target LinearTransform.
func (LT *LinearTransform) UnmarshalBinary(data []byte) (err error) {
	LT.LinearTransform = rlwe.LinearTransform{}
	LT.Scale = NewScale(0)
	return LT.LinearTransform.UnmarshalBinary(data)
}

// NewLinearTransform allocates a new LinearTransform with zero plaintexts at the specified level.
// If BSGSRatio == 0, the LinearTransform is set to not use the BSGS approach.
// Method will panic if BSGSRatio < 0.
func NewLinearTransform(params Parameters, nonZeroDiags []int, level int, BSGSRatio float64) LinearTransform {
	return LinearTransform{
		LinearTransform: rlwe.NewLinearTransform(params.Parameters, nonZeroDiags, level, params.LogN()-1, BSGSRatio),
	}
}

// Encode encodes on a pre-allocated LinearTransform the linear transforms' matrix in diagonal form `value`.
// values.(type) can be either map[int][]complex128 or map[int][]float64.
// User must ensure that 1 <= len([]complex128\[]float64) <= 2^logSlots < 2^logN.
// It can then be evaluated on a ciphertext using evaluator.LinearTransform.
// Evaluation will use the naive approach (single hoisting and no baby-step giant-step).
// Faster if there is only a few non-zero diagonals but uses more keys.
func (LT *LinearTransform) Encode(ecd Encoder, dMat map[int][]uint64, scale rlwe.Scale) {

	params := ecd.Parameters()
	ringQP := params.RingQP()

	levelQ := LT.Level
	levelP := params.PCount() - 1

	buffT := params.RingT().NewPoly()

	encode := func(value interface{}, opQP rlwe.OperandQP) {

		ecd.EncodeRingT(value, scale, buffT)
		ecd.RingT2Q(levelQ, buffT, opQP.El().Value[0].Q)
		ecd.RingT2Q(levelP, buffT, opQP.El().Value[0].P)

		ringQP.NTTLvl(levelQ, levelP, opQP.El().Value[0], opQP.El().Value[0])
		ringQP.MFormLvl(levelQ, levelP, opQP.El().Value[0], opQP.El().Value[0])
	}

	LT.LinearTransform.Encode(encode, mapUint64ToMapOfInterface(dMat))
	LT.Scale = scale.CopyNew()
}

func mapUint64ToMapOfInterface(m map[int][]uint64) map[int]interface{} {
	d := make(map[int]interface{})
	for i := range m {
		d[i] = m[i]
	}
	return d
}

// GenLinearTransformBSGS allocates and encodes a new LinearTransform struct from the linear transforms' matrix in diagonal form `value` for evaluation with a baby-step giant-step approach.
// values.(type) can be either map[int][]complex128 or map[int][]float64.
// User must ensure that 1 <= len([]complex128\[]float64) <= 2^logSlots < 2^logN.
// LinearTransform types can be be evaluated on a ciphertext using evaluator.LinearTransform.
// Evaluation will use the optimized approach (double hoisting and baby-step giant-step).
// Faster if there is more than a few non-zero diagonals.
// BSGSRatio is the maximum ratio between the inner and outer loop of the baby-step giant-step algorithm used in evaluator.LinearTransform.
// Optimal BSGSRatio value is between 4 and 16 depending on the sparsity of the matrix.
func GenLinearTransform(ecd Encoder, enc Encryptor, dMat map[int][]uint64, level int, scale rlwe.Scale, BSGSRatio float64) LinearTransform {

	if ecd == nil {
		panic("GenLinearTransformBSGS: ecd cannot be nil")
	}

	params := ecd.Parameters()
	ringQP := params.RingQP()
	levelQ := level
	levelP := params.PCount() - 1

	var embed func(value interface{}) (opQP rlwe.OperandQP)

	buffT := params.RingT().NewPoly()

	if enc != nil {

		pt := &rlwe.PlaintextQP{Value: params.RingQP().NewPolyLvl(levelQ, levelP)}

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

			ecd.EncodeRingT(value, scale, buffT)
			ecd.RingT2Q(levelQ, buffT, pt.Value.Q)
			ecd.RingT2Q(levelP, buffT, pt.Value.P)

			ringQP.NTTLvl(levelQ, levelP, pt.Value, pt.Value)
			ringQP.MFormLvl(levelQ, levelP, pt.Value, pt.Value)

			params.RingQP().AddLvl(levelQ, levelP, ct.Value[0], pt.Value, ct.Value[0])

			return ct
		}

	} else {
		embed = func(value interface{}) (opQP rlwe.OperandQP) {

			pt := &rlwe.PlaintextQP{Value: params.RingQP().NewPolyLvl(levelQ, levelP)}

			ecd.EncodeRingT(value, scale, buffT)
			ecd.RingT2Q(levelQ, buffT, pt.Value.Q)
			ecd.RingT2Q(levelP, buffT, pt.Value.P)

			ringQP.NTTLvl(levelQ, levelP, pt.Value, pt.Value)
			ringQP.MFormLvl(levelQ, levelP, pt.Value, pt.Value)

			return pt
		}
	}

	LT := rlwe.GenLinearTransform(embed, mapUint64ToMapOfInterface(dMat), BSGSRatio, params.LogN()-1)
	LT.Scale = scale.CopyNew()

	return LinearTransform{
		LinearTransform: LT,
	}
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
			ctOut[i] = NewCiphertext(eval.params, 1, minLevel, NewScale(0))
		}

	case LinearTransform:
		ctOut = []*Ciphertext{NewCiphertext(eval.params, 1, utils.MinInt(LTs.Level, ctIn.Level()), NewScale(0))}

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

	ringQ := eval.params.RingQ()

	switch LTs := linearTransform.(type) {
	case []LinearTransform:
		ctOut = make([]*Ciphertext, len(LTs))

		var maxLevel int
		for _, LT := range LTs {
			maxLevel = utils.MaxInt(maxLevel, LT.Level)
		}

		minLevel := utils.MinInt(maxLevel, ctIn.Level())

		ringQ.MulScalarBigintLvl(minLevel, ctIn.Value[0], eval.tInvModQ[minLevel], eval.buffCt.Value[0])
		ringQ.MulScalarBigintLvl(minLevel, ctIn.Value[1], eval.tInvModQ[minLevel], eval.buffCt.Value[1])

		eval.DecomposeNTT(minLevel, eval.params.PCount()-1, eval.params.PCount(), eval.buffCt.Value[1], eval.BuffDecompQP)

		for i, LT := range LTs {
			if LT.N1 == 0 {
				eval.MultiplyByDiagMatrix(eval.buffCt.Ciphertext, LT.LinearTransform, eval.BuffDecompQP, ctOut[i].Ciphertext)
			} else {
				eval.MultiplyByDiagMatrixBSGS(eval.buffCt.Ciphertext, LT.LinearTransform, eval.BuffDecompQP, ctOut[i].Ciphertext)
			}

			ringQ.MulScalarLvl(minLevel, ctOut[i].Value[0], eval.params.T(), ctOut[i].Value[0])
			ringQ.MulScalarLvl(minLevel, ctOut[i].Value[1], eval.params.T(), ctOut[i].Value[1])

			ctOut[i].Scale().Set(ctIn.Scale())
			ctOut[i].Scale().Mul(LT.Scale)
		}

	case LinearTransform:

		minLevel := utils.MinInt(LTs.Level, ctIn.Level())

		ringQ.MulScalarBigintLvl(minLevel, ctIn.Value[0], eval.tInvModQ[minLevel], eval.buffCt.Value[0])
		ringQ.MulScalarBigintLvl(minLevel, ctIn.Value[1], eval.tInvModQ[minLevel], eval.buffCt.Value[1])

		eval.DecomposeNTT(minLevel, eval.params.PCount()-1, eval.params.PCount(), eval.buffCt.Value[1], eval.BuffDecompQP)

		if LTs.N1 == 0 {
			eval.MultiplyByDiagMatrix(eval.buffCt.Ciphertext, LTs.LinearTransform, eval.BuffDecompQP, ctOut[0].Ciphertext)
		} else {
			eval.MultiplyByDiagMatrixBSGS(eval.buffCt.Ciphertext, LTs.LinearTransform, eval.BuffDecompQP, ctOut[0].Ciphertext)
		}

		ringQ.MulScalarLvl(minLevel, ctOut[0].Value[0], eval.params.T(), ctOut[0].Value[0])
		ringQ.MulScalarLvl(minLevel, ctOut[0].Value[1], eval.params.T(), ctOut[0].Value[1])

		ctOut[0].Scale().Set(ctIn.Scale())
		ctOut[0].Scale().Mul(LTs.Scale)
	}
	return
}
