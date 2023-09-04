package circuits

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// EvaluatorForLinearTransformation defines a set of common and scheme agnostic method necessary to instantiate an LinearTransformationEvaluator.
type EvaluatorForLinearTransformation interface {
	rlwe.ParameterProvider
	Rescale(ctIn, ctOut *rlwe.Ciphertext) (err error)

	// TODO: separated int
	DecomposeNTT(levelQ, levelP, nbPi int, c2 ring.Poly, c2IsNTT bool, decompQP []ringqp.Poly)
	CheckAndGetGaloisKey(galEl uint64) (evk *rlwe.GaloisKey, err error)
	GadgetProductLazy(levelQ int, cx ring.Poly, gadgetCt *rlwe.GadgetCiphertext, ct *rlwe.Element[ringqp.Poly])
	GadgetProductHoistedLazy(levelQ int, BuffQPDecompQP []ringqp.Poly, gadgetCt *rlwe.GadgetCiphertext, ct *rlwe.Element[ringqp.Poly])
	AutomorphismHoistedLazy(levelQ int, ctIn *rlwe.Ciphertext, c1DecompQP []ringqp.Poly, galEl uint64, ctQP *rlwe.Element[ringqp.Poly]) (err error)
	ModDownQPtoQNTT(levelQ, levelP int, p1Q, p1P, p2Q ring.Poly)
	AutomorphismIndex(uint64) []uint64

	GetEvaluatorBuffer() *rlwe.EvaluatorBuffers // TODO extract
}

// LinearTransformationEvaluator is an evaluator used to evaluate linear transformations on ciphertexts.
type LinearTransformationEvaluator struct {
	EvaluatorForLinearTransformation
	*rlwe.EvaluatorBuffers
}

// NewLinearTransformationEvaluator instantiates a new LinearTransformationEvaluator from an EvaluatorForLinearTransformation.
// The method is allocation free if the underlying EvaluatorForLinearTransformation returns a non-nil *rlwe.EvaluatorBuffers.
func NewLinearTransformationEvaluator(eval EvaluatorForLinearTransformation) (linTransEval *LinearTransformationEvaluator) {
	linTransEval = new(LinearTransformationEvaluator)
	linTransEval.EvaluatorForLinearTransformation = eval
	linTransEval.EvaluatorBuffers = eval.GetEvaluatorBuffer()
	if linTransEval.EvaluatorBuffers == nil {
		linTransEval.EvaluatorBuffers = rlwe.NewEvaluatorBuffers(*eval.GetRLWEParameters())
	}
	return
}

// EvaluateNew takes as input a ciphertext ctIn and a linear transformation M and evaluate and returns opOut: M(ctIn).
func (eval LinearTransformationEvaluator) EvaluateNew(ctIn *rlwe.Ciphertext, linearTransformation LinearTransformation) (opOut *rlwe.Ciphertext, err error) {
	cts, err := eval.EvaluateManyNew(ctIn, []LinearTransformation{linearTransformation})
	return cts[0], err
}

// Evaluate takes as input a ciphertext ctIn, a linear transformation M and evaluates opOut: M(ctIn).
func (eval LinearTransformationEvaluator) Evaluate(ctIn *rlwe.Ciphertext, linearTransformation LinearTransformation, opOut *rlwe.Ciphertext) (err error) {
	return eval.EvaluateMany(ctIn, []LinearTransformation{linearTransformation}, []*rlwe.Ciphertext{opOut})
}

// EvaluateManyNew takes as input a ciphertext ctIn and a list of linear transformations [M0, M1, M2, ...] and returns opOut:[M0(ctIn), M1(ctIn), M2(ctInt), ...].
func (eval LinearTransformationEvaluator) EvaluateManyNew(ctIn *rlwe.Ciphertext, linearTransformations []LinearTransformation) (opOut []*rlwe.Ciphertext, err error) {

	params := eval.GetRLWEParameters()
	opOut = make([]*rlwe.Ciphertext, len(linearTransformations))
	for i := range opOut {
		opOut[i] = rlwe.NewCiphertext(params, 1, linearTransformations[i].Level)
	}

	return opOut, eval.EvaluateMany(ctIn, linearTransformations, opOut)
}

// EvaluateMany takes as input a ciphertext ctIn, a list of linear transformations [M0, M1, M2, ...] and a list of pre-allocated receiver opOut
// and evaluates opOut: [M0(ctIn), M1(ctIn), M2(ctIn), ...]
func (eval LinearTransformationEvaluator) EvaluateMany(ctIn *rlwe.Ciphertext, linearTransformations []LinearTransformation, opOut []*rlwe.Ciphertext) (err error) {

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

// EvaluateSequentialNew takes as input a ciphertext ctIn and a list of linear transformations [M0, M1, M2, ...] and evaluates and returns opOut:...M2(M1(M0(ctIn))
func (eval LinearTransformationEvaluator) EvaluateSequentialNew(ctIn *rlwe.Ciphertext, linearTransformations []LinearTransformation) (opOut *rlwe.Ciphertext, err error) {
	params := eval.GetRLWEParameters()
	opOut = rlwe.NewCiphertext(params, 1, linearTransformations[0].Level)
	return opOut, eval.EvaluateSequential(ctIn, linearTransformations, opOut)
}

// EvaluateSequential takes as input a ciphertext ctIn and a list of linear transformations [M0, M1, M2, ...] and evaluates opOut:...M2(M1(M0(ctIn))
func (eval LinearTransformationEvaluator) EvaluateSequential(ctIn *rlwe.Ciphertext, linearTransformations []LinearTransformation, opOut *rlwe.Ciphertext) (err error) {

	if err = eval.Evaluate(ctIn, linearTransformations[0], opOut); err != nil {
		return
	}

	if err = eval.Rescale(opOut, opOut); err != nil {
		return
	}

	for i := 1; i < len(linearTransformations); i++ {
		if err = eval.Evaluate(opOut, linearTransformations[i], opOut); err != nil {
			return
		}

		if err = eval.Rescale(opOut, opOut); err != nil {
			return
		}
	}

	return
}

// MultiplyByDiagMatrix multiplies the Ciphertext "ctIn" by the plaintext matrix "matrix" and returns the result on the Ciphertext
// "opOut". Memory buffers for the decomposed ciphertext BuffDecompQP, BuffDecompQP must be provided, those are list of poly of ringQ and ringP
// respectively, each of size params.Beta().
// The naive approach is used (single hoisting and no baby-step giant-step), which is faster than MultiplyByDiagMatrixBSGS
// for matrix of only a few non-zero diagonals but uses more keys.
func (eval LinearTransformationEvaluator) MultiplyByDiagMatrix(ctIn *rlwe.Ciphertext, matrix LinearTransformation, BuffDecompQP []ringqp.Poly, opOut *rlwe.Ciphertext) (err error) {

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

	cQP := &rlwe.Element[ringqp.Poly]{}
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
func (eval LinearTransformationEvaluator) MultiplyByDiagMatrixBSGS(ctIn *rlwe.Ciphertext, matrix LinearTransformation, BuffDecompQP []ringqp.Poly, opOut *rlwe.Ciphertext) (err error) {

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
	ctInRotQP := map[int]*rlwe.Element[ringqp.Poly]{}
	for _, i := range rotN2 {
		if i != 0 {
			ctInRotQP[i] = rlwe.NewElementQP(params, 1, levelQ, levelP)
			if err = eval.AutomorphismHoistedLazy(levelQ, ctIn, BuffDecompQP, params.GaloisElement(i), ctInRotQP[i]); err != nil {
				return
			}
		}
	}

	// Accumulator inner loop
	tmp0QP := eval.BuffQP[1]
	tmp1QP := eval.BuffQP[2]

	// Accumulator outer loop
	cQP := &rlwe.Element[ringqp.Poly]{}
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
