package bgv

import (
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
	ringP := eval.params.RingP()
	ringQP := eval.params.RingQP()

	levelQ := ctIn.Level()
	levelP := len(ringP.Modulus) - 1

	ctOut.Resize(ctOut.Degree(), levelQ)
	ctOut.SetScale(ctIn.Scale())

	if n == 1 {
		if ctIn != ctOut {
			ring.CopyValuesLvl(levelQ, ctIn.Value[0], ctOut.Value[0])
			ring.CopyValuesLvl(levelQ, ctIn.Value[1], ctOut.Value[1])
		}
	} else {

		// Memory buffer for ctIn = ctIn + rot(ctIn, 2^i) in Q
		tmpc0 := eval.buffQ[0] // unused memory buffer from evaluator
		tmpc1 := eval.buffQ[1] // unused memory buffer from evaluator
		tmpc2 := eval.buffQ[2] // unused memory buffer from evaluator

		c0OutQP := eval.BuffQP[2]
		c1OutQP := eval.BuffQP[3]
		c0QP := eval.BuffQP[4]
		c1QP := eval.BuffQP[5]

		cQP := rlwe.CiphertextQP{Value: []ringqp.Poly{c0QP, c1QP}}

		tmpc0.IsNTT = true
		tmpc1.IsNTT = true
		c0QP.Q.IsNTT = true
		c1QP.Q.IsNTT = true

		tmpct := NewCiphertextAtLevelFromPoly(levelQ, [2]*ring.Poly{tmpc0, tmpc1})
		ctqp := NewCiphertextAtLevelFromPoly(levelQ, [2]*ring.Poly{c0QP.Q, c1QP.Q})

		state := false
		copy := true
		// Binary reading of the input n
		for i, j := 0, n; j > 0; i, j = i+1, j>>1 {

			// Starts by decomposing the input ciphertext
			if i == 0 {
				// If first iteration, then copies directly from the input ciphertext that hasn't been rotated
				ringQ.MulScalarBigintLvl(levelQ, ctIn.Value[1], eval.tInvModQ[levelQ], tmpc2)
				eval.DecomposeNTT(levelQ, levelP, levelP+1, tmpc2, eval.BuffDecompQP)
			} else {
				// Else copies from the rotated input ciphertext
				ringQ.MulScalarBigintLvl(levelQ, tmpc1, eval.tInvModQ[levelQ], tmpc2)
				eval.DecomposeNTT(levelQ, levelP, levelP+1, tmpc2, eval.BuffDecompQP)
			}

			// If the binary reading scans a 1
			if j&1 == 1 {

				k := n - (n & ((2 << i) - 1))
				k *= batchSize

				// If the rotation is not zero
				if k != 0 {

					// Rotate((tmpc0, tmpc1), k)
					if i == 0 {
						eval.AutomorphismHoistedNoModDown(levelQ, ctIn.Value[0], eval.BuffDecompQP, eval.params.GaloisElementForColumnRotationBy(k), cQP)
					} else {
						eval.AutomorphismHoistedNoModDown(levelQ, tmpc0, eval.BuffDecompQP, eval.params.GaloisElementForColumnRotationBy(k), cQP)
					}

					// ctOut += Rotate((tmpc0, tmpc1), k)
					if copy {
						ringQP.CopyValuesLvl(levelQ, levelP, c0QP, c0OutQP)
						ringQP.CopyValuesLvl(levelQ, levelP, c1QP, c1OutQP)
						copy = false
					} else {
						ringQP.AddLvl(levelQ, levelP, c0OutQP, c0QP, c0OutQP)
						ringQP.AddLvl(levelQ, levelP, c1OutQP, c1QP, c1OutQP)
					}
				} else {

					state = true

					// if n is not a power of two
					if n&(n-1) != 0 {
						eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c0OutQP.Q, c0OutQP.P, c0OutQP.Q) // Division by P
						eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c1OutQP.Q, c1OutQP.P, c1OutQP.Q) // Division by P

						ringQ.MulScalarLvl(levelQ, c0OutQP.Q, eval.params.T(), c0OutQP.Q)
						ringQ.MulScalarLvl(levelQ, c1OutQP.Q, eval.params.T(), c1OutQP.Q)

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
					eval.AutomorphismHoisted(levelQ, ctIn, eval.BuffDecompQP, rot, tmpct)
					ringQ.AddLvl(levelQ, tmpc0, ctIn.Value[0], tmpc0)
					ringQ.AddLvl(levelQ, tmpc1, ctIn.Value[1], tmpc1)
				} else {
					// (tmpc0, tmpc1) = Rotate((tmpc0, tmpc1), 2^i)
					eval.AutomorphismHoisted(levelQ, tmpct, eval.BuffDecompQP, rot, ctqp)
					ringQ.AddLvl(levelQ, tmpc0, c0QP.Q, tmpc0)
					ringQ.AddLvl(levelQ, tmpc1, c1QP.Q, tmpc1)
				}
			}
		}
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

type LinearTransform struct {
	rlwe.LinearTransform
	Scale uint64 // Scale is the scale at which the matrix is encoded (can be circuit dependent)
}

// NewLinearTransform allocates a new LinearTransform with zero plaintexts at the specified level.
// If BSGSRatio == 0, the LinearTransform is set to not use the BSGS approach.
// Method will panic if BSGSRatio < 0.
func NewLinearTransform(params Parameters, nonZeroDiags []int, level int, BSGSRatio float64) LinearTransform {
	vec := make(map[int]rlwe.OperandQP)
	slots := params.N() >> 1
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
		N1 = FindBestBSGSSplit(nonZeroDiags, slots, BSGSRatio)
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
			LogSlots: params.LogN() - 1,
			N1:       N1,
			Level:    level,
			Vec:      vec,
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
func (LT *LinearTransform) Encode(ecd Encoder, dMat map[int][]uint64, scale uint64) {

	enc, ok := ecd.(*encoder)
	if !ok {
		panic("encoder should be an encoderComplex128")
	}

	ringQP := enc.params.RingQP()

	levelQ := LT.Level
	levelP := enc.params.PCount() - 1

	slots := 1 << LT.LogSlots
	N1 := LT.N1

	buffT := enc.params.RingT().NewPoly()

	if N1 == 0 {
		for i := range dMat {
			idx := i
			if idx < 0 {
				idx += slots
			}

			if _, ok := LT.Vec[idx]; !ok {
				panic("error encoding on LinearTransform: input does not match the same non-zero diagonals")
			}

			enc.EncodeRingT(dMat[i], scale, buffT)
			enc.RingT2Q(levelQ, buffT, LT.Vec[idx].El().Value[0].Q)
			enc.RingT2Q(levelP, buffT, LT.Vec[idx].El().Value[0].P)

			ringQP.NTTLvl(levelQ, levelP, LT.Vec[idx].El().Value[0], LT.Vec[idx].El().Value[0])
			ringQP.MFormLvl(levelQ, levelP, LT.Vec[idx].El().Value[0], LT.Vec[idx].El().Value[0])
		}
	} else {
		index, _, _ := rlwe.BsgsIndex(dMat, slots, N1)

		for j := range index {
			for _, i := range index[j] {
				// manages inputs that have rotation between 0 and slots-1 or between -slots/2 and slots/2-1
				v, ok := dMat[j+i]
				if !ok {
					v = dMat[j+i-slots]
				}

				if _, ok := LT.Vec[j+i]; !ok {
					panic("error encoding on LinearTransform BSGS: input does not match the same non-zero diagonals")
				}

				enc.EncodeRingT(utils.RotateUint64Slots(v, -j), scale, buffT)
				enc.RingT2Q(levelQ, buffT, LT.Vec[j+i].El().Value[0].Q)
				enc.RingT2Q(levelP, buffT, LT.Vec[j+i].El().Value[0].P)

				ringQP.NTTLvl(levelQ, levelP, LT.Vec[j+i].El().Value[0], LT.Vec[j+i].El().Value[0])
				ringQP.MFormLvl(levelQ, levelP, LT.Vec[j+i].El().Value[0], LT.Vec[j+i].El().Value[0])
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
func GenLinearTransform(ecd Encoder, dMat map[int][]uint64, level int, scale uint64) LinearTransform {

	enc, ok := ecd.(*encoder)
	if !ok {
		panic("encoder should be an encoderComplex128")
	}

	params := enc.params
	vec := make(map[int]rlwe.OperandQP)
	slots := params.N() >> 1
	levelQ := level
	levelP := params.PCount() - 1
	ringQP := params.RingQP()
	buffT := params.RingT().NewPoly()
	for i := range dMat {

		idx := i
		if idx < 0 {
			idx += slots
		}
		vec[idx] = &rlwe.PlaintextQP{Value: params.RingQP().NewPolyLvl(levelQ, levelP)}

		enc.EncodeRingT(dMat[i], scale, buffT)
		enc.RingT2Q(levelQ, buffT, vec[idx].El().Value[0].Q)
		enc.RingT2Q(levelP, buffT, vec[idx].El().Value[0].P)

		ringQP.NTTLvl(levelQ, levelP, vec[idx].El().Value[0], vec[idx].El().Value[0])
		ringQP.MFormLvl(levelQ, levelP, vec[idx].El().Value[0], vec[idx].El().Value[0])
	}

	return LinearTransform{
		LinearTransform: rlwe.LinearTransform{
			LogSlots: params.LogN() - 1,
			N1:       0,
			Vec:      vec,
			Level:    level,
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
func GenLinearTransformBSGS(ecd Encoder, dMat map[int][]uint64, level int, scale uint64, BSGSRatio float64) (LT LinearTransform) {

	enc, ok := ecd.(*encoder)
	if !ok {
		panic("cannot GenLinearTransformBSGS: encoder should be an encoderComplex128")
	}

	params := enc.params

	slots := params.N() >> 1

	// N1*N2 = N
	N1 := FindBestBSGSSplit(dMat, slots, BSGSRatio)

	index, _, _ := rlwe.BsgsIndex(dMat, slots, N1)

	vec := make(map[int]rlwe.OperandQP)

	levelQ := level
	levelP := params.PCount() - 1
	ringQP := params.RingQP()

	buffT := params.RingT().NewPoly()

	for j := range index {

		for _, i := range index[j] {

			// manages inputs that have rotation between 0 and slots-1 or between -slots/2 and slots/2-1
			v, ok := dMat[j+i]
			if !ok {
				v = dMat[j+i-slots]
			}
			vec[j+i] = &rlwe.PlaintextQP{Value: params.RingQP().NewPolyLvl(levelQ, levelP)}

			enc.EncodeRingT(utils.RotateUint64Slots(v, -j), scale, buffT)
			enc.RingT2Q(levelQ, buffT, vec[j+i].El().Value[0].Q)
			enc.RingT2Q(levelP, buffT, vec[j+i].El().Value[0].P)

			ringQP.NTTLvl(levelQ, levelP, vec[j+i].El().Value[0], vec[j+i].El().Value[0])
			ringQP.MFormLvl(levelQ, levelP, vec[j+i].El().Value[0], vec[j+i].El().Value[0])
		}
	}

	return LinearTransform{
		LinearTransform: rlwe.LinearTransform{
			LogSlots: params.LogN() - 1,
			N1:       N1,
			Vec:      vec,
			Level:    level,
		},
		Scale: scale,
	}
}

// FindBestBSGSSplit finds the best N1*N2 = N for the baby-step giant-step algorithm for matrix multiplication.
func FindBestBSGSSplit(diagMatrix interface{}, maxN int, maxRatio float64) (minN int) {

	for N1 := 1; N1 < maxN; N1 <<= 1 {

		_, rotN1, rotN2 := rlwe.BsgsIndex(diagMatrix, maxN, N1)

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
			ctOut[i] = NewCiphertext(eval.params, 1, minLevel, 0)
		}

	case LinearTransform:
		ctOut = []*Ciphertext{NewCiphertext(eval.params, 1, utils.MinInt(LTs.Level, ctIn.Level()), 0)}

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

			ctOut[i].SetScale(ctIn.Scale() * LT.Scale)
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

		ctOut[0].SetScale(ctIn.Scale() * LTs.Scale)
	}
	return
}
