package bgv

import (
	"runtime"

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

type OperandQP interface {
	Degree() (deg int)
	Level() (levelQ, levelP int)
	El() rlwe.CiphertextQP
	GetDataLen64(WithMetaData bool) (dataLen int)
	WriteTo64(data []byte) (ptr int, err error)
	Decode64(data []byte) (ptr int, err error)
}

// LinearTransform is a type for linear transformations on ciphertexts.
// It stores a plaintext matrix diagonalized in diagonal form and
// can be evaluated on a ciphertext by using the evaluator.LinearTransform method.
type LinearTransform struct {
	LogSlots int
	N1       int               // N1 is the number of inner loops of the baby-step giant-step algorithm used in the evaluation (if N1 == 0, BSGS is not used).
	Level    int               // Level is the level at which the matrix is encoded (can be circuit dependent)
	Scale    uint64            // Scale is the scale at which the matrix is encoded (can be circuit dependent)
	Vec      map[int]OperandQP // Vec is the matrix, in diagonal form, where each entry of vec is an indexed non-zero diagonal.
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

// NewLinearTransform allocates a new LinearTransform with zero plaintexts at the specified level.
// If BSGSRatio == 0, the LinearTransform is set to not use the BSGS approach.
// Method will panic if BSGSRatio < 0.
func NewLinearTransform(params Parameters, nonZeroDiags []int, level int, BSGSRatio float64) LinearTransform {
	vec := make(map[int]OperandQP)
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
		index, _, _ := BsgsIndex(nonZeroDiags, slots, N1)
		for j := range index {
			for _, i := range index[j] {
				vec[j+i] = &rlwe.PlaintextQP{Value: params.RingQP().NewPolyLvl(levelQ, levelP)}
			}
		}
	} else {
		panic("BSGS ratio cannot be negative")
	}

	return LinearTransform{LogSlots: params.LogN() - 1, N1: N1, Level: level, Vec: vec}
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
		index, _, _ := BsgsIndex(dMat, slots, N1)

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
	vec := make(map[int]OperandQP)
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

	return LinearTransform{LogSlots: params.LogN() - 1, N1: 0, Vec: vec, Level: level, Scale: scale}
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

	index, _, _ := BsgsIndex(dMat, slots, N1)

	vec := make(map[int]OperandQP)

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

	return LinearTransform{LogSlots: params.LogN() - 1, N1: N1, Vec: vec, Level: level, Scale: scale}
}

// BsgsIndex returns the index map and needed rotation for the BSGS matrix-vector multiplication algorithm.
func BsgsIndex(el interface{}, slots, N1 int) (index map[int][]int, rotN1, rotN2 []int) {
	index = make(map[int][]int)
	rotN1Map := make(map[int]bool)
	rotN2Map := make(map[int]bool)
	var nonZeroDiags []int
	switch element := el.(type) {
	case map[int][]uint64:
		nonZeroDiags = make([]int, len(element))
		var i int
		for key := range element {
			nonZeroDiags[i] = key
			i++
		}
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
	case map[int]*rlwe.PlaintextQP:
		nonZeroDiags = make([]int, len(element))
		var i int
		for key := range element {
			nonZeroDiags[i] = key
			i++
		}
	case map[int]*rlwe.CiphertextQP:
		nonZeroDiags = make([]int, len(element))
		var i int
		for key := range element {
			nonZeroDiags[i] = key
			i++
		}
	case map[int]OperandQP:
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

// FindBestBSGSSplit finds the best N1*N2 = N for the baby-step giant-step algorithm for matrix multiplication.
func FindBestBSGSSplit(diagMatrix interface{}, maxN int, maxRatio float64) (minN int) {

	for N1 := 1; N1 < maxN; N1 <<= 1 {

		_, rotN1, rotN2 := BsgsIndex(diagMatrix, maxN, N1)

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
		eval.params.RingQ().MulScalarBigintLvl(minLevel, ctIn.Value[1], eval.tInvModQ[minLevel], eval.buffQ[0])
		eval.DecomposeNTT(minLevel, eval.params.PCount()-1, eval.params.PCount(), eval.buffQ[0], eval.BuffDecompQP)

		for i, LT := range LTs {
			ctOut[i] = NewCiphertext(eval.params, 1, minLevel, ctIn.Scale())

			if LT.N1 == 0 {
				eval.MultiplyByDiagMatrix(ctIn, LT, eval.BuffDecompQP, ctOut[i])
			} else {
				eval.MultiplyByDiagMatrixBSGS(ctIn, LT, eval.BuffDecompQP, ctOut[i])
			}
		}

	case LinearTransform:

		minLevel := utils.MinInt(LTs.Level, ctIn.Level())
		eval.params.RingQ().MulScalarBigintLvl(minLevel, ctIn.Value[1], eval.tInvModQ[minLevel], eval.buffQ[0])
		eval.DecomposeNTT(minLevel, eval.params.PCount()-1, eval.params.PCount(), eval.buffQ[0], eval.BuffDecompQP)

		ctOut = []*Ciphertext{NewCiphertext(eval.params, 1, minLevel, ctIn.Scale())}

		if LTs.N1 == 0 {
			eval.MultiplyByDiagMatrix(ctIn, LTs, eval.BuffDecompQP, ctOut[0])
		} else {
			eval.MultiplyByDiagMatrixBSGS(ctIn, LTs, eval.BuffDecompQP, ctOut[0])
		}
	}
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
		eval.params.RingQ().MulScalarBigintLvl(minLevel, ctIn.Value[1], eval.tInvModQ[minLevel], eval.buffQ[0])
		eval.DecomposeNTT(minLevel, eval.params.PCount()-1, eval.params.PCount(), eval.buffQ[0], eval.BuffDecompQP)

		for i, LT := range LTs {
			if LT.N1 == 0 {
				eval.MultiplyByDiagMatrix(ctIn, LT, eval.BuffDecompQP, ctOut[i])
			} else {
				eval.MultiplyByDiagMatrixBSGS(ctIn, LT, eval.BuffDecompQP, ctOut[i])
			}
		}

	case LinearTransform:
		minLevel := utils.MinInt(LTs.Level, ctIn.Level())
		eval.params.RingQ().MulScalarBigintLvl(minLevel, ctIn.Value[1], eval.tInvModQ[minLevel], eval.buffQ[0])
		eval.DecomposeNTT(minLevel, eval.params.PCount()-1, eval.params.PCount(), eval.buffQ[0], eval.BuffDecompQP)
		if LTs.N1 == 0 {
			eval.MultiplyByDiagMatrix(ctIn, LTs, eval.BuffDecompQP, ctOut[0])
		} else {
			eval.MultiplyByDiagMatrixBSGS(ctIn, LTs, eval.BuffDecompQP, ctOut[0])
		}
	}
}

// MultiplyByDiagMatrix multiplies the ciphertext "ctIn" by the plaintext matrix "matrix" and returns the result on the ciphertext
// "ctOut". Memory buffers for the decomposed ciphertext BuffDecompQP, BuffDecompQP must be provided, those are list of poly of ringQ and ringP
// respectively, each of size params.Beta().
// The naive approach is used (single hoisting and no baby-step giant-step), which is faster than MultiplyByDiagMatrixBSGS
// for matrix of only a few non-zero diagonals but uses more keys.
func (eval *evaluator) MultiplyByDiagMatrix(ctIn *Ciphertext, matrix LinearTransform, BuffDecompQP []ringqp.Poly, ctOut *Ciphertext) {

	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()
	ringQP := eval.params.RingQP()

	levelQ := utils.MinInt(ctOut.Level(), utils.MinInt(ctIn.Level(), matrix.Level))
	levelP := len(ringP.Modulus) - 1

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
	ksResQP := rlwe.CiphertextQP{Value: []ringqp.Poly{ksRes0QP, ksRes1QP}}

	ring.CopyValuesLvl(levelQ, ctIn.Value[0], eval.buffCt.Value[0])
	ring.CopyValuesLvl(levelQ, ctIn.Value[1], eval.buffCt.Value[1])
	ctInTmp0, ctInTmp1 := eval.buffCt.Value[0], eval.buffCt.Value[1]

	ringQ.MulScalarBigintLvl(levelQ, ctInTmp0, ringP.ModulusAtLevel[levelP], ct0TimesP) // P*c0
	ringQ.MulScalarBigintLvl(levelQ, ct0TimesP, eval.tInvModQ[levelQ], ct0TimesP)

	var state bool
	var cnt int
	for k := range matrix.Vec {

		k &= int((ringQ.NthRoot >> 2) - 1)

		if k == 0 {
			state = true
		} else {

			galEl := eval.params.GaloisElementForColumnRotationBy(k)

			rtk, generated := eval.Rtks.Keys[galEl]
			if !generated {
				panic("cannot MultiplyByDiagMatrix: switching key not available")
			}

			index := eval.PermuteNTTIndex[galEl]

			eval.GadgetProductHoistedNoModDown(levelQ, BuffDecompQP, rtk.GadgetCiphertext, ksResQP)
			ringQ.AddLvl(levelQ, ksRes0QP.Q, ct0TimesP, ksRes0QP.Q)
			ringQP.PermuteNTTWithIndexLvl(levelQ, levelP, ksRes0QP, index, tmp0QP)
			ringQP.PermuteNTTWithIndexLvl(levelQ, levelP, ksRes1QP, index, tmp1QP)

			if cnt == 0 {
				// keyswitch(c1_Q) = (d0_QP, d1_QP)
				ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, matrix.Vec[k].El().Value[0], tmp0QP, c0OutQP)
				ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, matrix.Vec[k].El().Value[0], tmp1QP, c1OutQP)
			} else {
				// keyswitch(c1_Q) = (d0_QP, d1_QP)
				ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, matrix.Vec[k].El().Value[0], tmp0QP, c0OutQP)
				ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, matrix.Vec[k].El().Value[0], tmp1QP, c1OutQP)
			}

			if cnt%QiOverF == QiOverF-1 {
				ringQ.ReduceLvl(levelQ, c0OutQP.Q, c0OutQP.Q)
				ringQ.ReduceLvl(levelQ, c1OutQP.Q, c1OutQP.Q)
			}

			if cnt%PiOverF == PiOverF-1 {
				ringP.ReduceLvl(levelP, c0OutQP.P, c0OutQP.P)
				ringP.ReduceLvl(levelP, c1OutQP.P, c1OutQP.P)
			}

			cnt++
		}
	}

	if cnt%QiOverF == 0 {
		ringQ.ReduceLvl(levelQ, c0OutQP.Q, c0OutQP.Q)
		ringQ.ReduceLvl(levelQ, c1OutQP.Q, c1OutQP.Q)
	}

	if cnt%PiOverF == 0 {
		ringP.ReduceLvl(levelP, c0OutQP.P, c0OutQP.P)
		ringP.ReduceLvl(levelP, c1OutQP.P, c1OutQP.P)
	}

	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c0OutQP.Q, c0OutQP.P, c0OutQP.Q) // sum(phi(c0 * P + d0_QP))/P
	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c1OutQP.Q, c1OutQP.P, c1OutQP.Q) // sum(phi(d1_QP))/P

	ringQ.MulScalarLvl(levelQ, c0OutQP.Q, eval.params.T(), c0OutQP.Q)
	ringQ.MulScalarLvl(levelQ, c1OutQP.Q, eval.params.T(), c1OutQP.Q)

	if state { // Rotation by zero
		ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.Vec[0].El().Value[0].Q, ctInTmp0, c0OutQP.Q) // ctOut += c0_Q * plaintext
		ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.Vec[0].El().Value[0].Q, ctInTmp1, c1OutQP.Q) // ctOut += c1_Q * plaintext
	}

	ctOut.SetScale(matrix.Scale * ctIn.Scale())
}

// MultiplyByDiagMatrixBSGS multiplies the ciphertext "ctIn" by the plaintext matrix "matrix" and returns the result on the ciphertext
// "ctOut". Memory buffers for the decomposed ciphertext BuffDecompQP, BuffDecompQP must be provided, those are list of poly of ringQ and ringP
// respectively, each of size params.Beta().
// The BSGS approach is used (double hoisting with baby-step giant-step), which is faster than MultiplyByDiagMatrix
// for matrix with more than a few non-zero diagonals and uses much less keys.
func (eval *evaluator) MultiplyByDiagMatrixBSGS(ctIn *Ciphertext, matrix LinearTransform, PoolDecompQP []ringqp.Poly, ctOut *Ciphertext) {

	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()
	ringQP := eval.params.RingQP()

	levelQ := utils.MinInt(ctOut.Level(), utils.MinInt(ctIn.Level(), matrix.Level))
	levelP := len(ringP.Modulus) - 1

	ctOut.Resize(ctOut.Degree(), levelQ)

	QiOverF := eval.params.QiOverflowMargin(levelQ) >> 1
	PiOverF := eval.params.PiOverflowMargin(levelP) >> 1

	isPlaintext := matrix.IsPlaintext()

	// Accumulator inner loop
	accInLoop := rlwe.CiphertextQP{Value: []ringqp.Poly{eval.BuffQP[2], eval.BuffQP[3]}}

	// Available & Used Buffers:
	// Evaluator CKKS:
	//		buffQ[0] -> accInLoop.Value[2].Q
	//		buffQ[1] -> accInLoop.Value[2].Q
	//		buffQ[2]
	//		buffCt   -> ctInMulP
	// Evaluator RLWE:
	//		BuffQ        [0]ringqp.Poly -> c2QP.Q (rlwe.Evaluator)
	//		BuffP        [0]ringqp.Poly -> c2QP.P (rlwe.Evaluator)
	//		BuffQ        [0]ringqp.Poly -> ctBuffQP.Value[0].Q
	//		BuffP        [0]ringqp.Poly -> ctBuffQP.Value[0].Q
	//		BuffQ        [1]ringqp.Poly -> ctBuffQP.Value[1].P
	//		BuffP        [1]ringqp.Poly -> ctBuffQP.Value[1].P
	//		BuffQ        [2]ringqp.Poly -> accInLoop.Value[0].Q
	//		BuffP        [2]ringqp.Poly -> accInLoop.Value[0].P
	//		BuffQ        [3]ringqp.Poly -> accInLoop.Value[1].Q
	//		BuffP        [3]ringqp.Poly -> accInLoop.Value[1].P
	//		BuffQ        [4]ringqp.Poly -> accOutLoopQP.Value[0].P
	//		BuffP        [4]ringqp.Poly -> accOutLoopQP.Value[1].P
	//		BuffInvNTT    *ring.Poly    -> cxInvNTT (rlwe.Evaluator)
	//		BuffDecompQP  []ringqp.Poly
	//

	if !isPlaintext {
		accInLoop.Value = append(accInLoop.Value, ringqp.Poly{Q: eval.buffQ[0], P: eval.buffQ[1]})
	}

	// Result in QP
	accOutLoopQP := &rlwe.CiphertextQP{Value: []ringqp.Poly{ringqp.Poly{Q: ctOut.Value[0], P: eval.BuffQP[4].Q}, ringqp.Poly{Q: ctOut.Value[1], P: eval.BuffQP[4].P}}}

	// Ciphertext buffer mod QP for gadget products without mod down
	ctBuffQP := rlwe.CiphertextQP{Value: []ringqp.Poly{eval.BuffQP[0], eval.BuffQP[1]}}

	// Computes the N2 rotations indexes of the non-zero rows of the diagonalized DFT matrix for the baby-step giang-step algorithm
	index, _, rotN2 := BsgsIndex(matrix.Vec, 1<<matrix.LogSlots, matrix.N1)

	ctInMulP := eval.buffCt

	ringQ.MulScalarBigintLvl(levelQ, ctIn.Value[0], eval.tInvModQ[levelQ], ctInMulP.Value[0])
	ringQ.MulScalarBigintLvl(levelQ, ctIn.Value[1], eval.tInvModQ[levelQ], ctInMulP.Value[1])

	// Pre-rotates ciphertext for the baby-step giant-step algorithm, does not divide by P yet
	ctInRotQP := eval.RotateHoistedNoModDownNew(levelQ, rotN2, ctInMulP.Value[0], PoolDecompQP)

	// Rotation by 0 -> pre-multiplies by P
	ringQ.MulScalarBigintLvl(levelQ, ctInMulP.Value[0], ringP.ModulusAtLevel[levelP], ctInMulP.Value[0]) // P*c0
	ringQ.MulScalarBigintLvl(levelQ, ctInMulP.Value[1], ringP.ModulusAtLevel[levelP], ctInMulP.Value[1]) // P*c1

	// OUTER LOOP
	var cnt0 int
	for j := range index {

		// INNER LOOP
		if isPlaintext {
			innerLoopBSGSDegree1(j, index[j], levelQ, levelP, *ringQP, matrix.Vec, ctInMulP, ctInRotQP, accInLoop, QiOverF, PiOverF)
		} else {
			innerLoopBSGSDegree2(j, index[j], levelQ, levelP, *ringQP, matrix.Vec, ctInMulP, ctInRotQP, accInLoop, QiOverF, PiOverF)

			// acc[2] mod QP -> acc[2]/P mod Q
			eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, accInLoop.Value[2].Q, accInLoop.Value[2].P, accInLoop.Value[2].Q)

			// Relinearize
			accInLoop.Value[2].Q.IsNTT = true
			eval.GadgetProductNoModDown(levelQ, accInLoop.Value[2].Q, eval.Rlk.Keys[0].GadgetCiphertext, ctBuffQP)

			// [acc[0] + buffRelin[0], acc[1] + buffRelin[1]] mod QP
			ringQP.AddLvl(levelQ, levelP, accInLoop.Value[0], ctBuffQP.Value[0], accInLoop.Value[0])
			ringQP.AddLvl(levelQ, levelP, accInLoop.Value[1], ctBuffQP.Value[1], accInLoop.Value[1])
		}

		// If j != 0, then rotates ((tmp0QP.Q, tmp0QP.P), (tmp1QP.Q, tmp1QP.P)) by N1*j and adds the result on ((c0QP.Q, c0QP.P), (c1QP.Q, c1QP.P))
		if j != 0 {

			// Hoisting of the ModDown of sum(sum(phi(d1) * plaintext))
			eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, accInLoop.Value[1].Q, accInLoop.Value[1].P, accInLoop.Value[1].Q) // c1 * plaintext + sum(phi(d1) * plaintext) + phi(c1) * plaintext mod Q

			galEl := eval.params.GaloisElementForColumnRotationBy(j)

			rtk, generated := eval.Rtks.Keys[galEl]
			if !generated {
				panic("cannot MultiplyByDiagMatrixBSGS: switching key not available")
			}

			rotIndex := eval.PermuteNTTIndex[galEl]

			accInLoop.Value[1].Q.IsNTT = true
			eval.GadgetProductNoModDown(levelQ, accInLoop.Value[1].Q, rtk.GadgetCiphertext, ctBuffQP) // Switchkey(P*phi(tmpRes_1)) = (d0, d1) in base QP
			ringQP.AddLvl(levelQ, levelP, ctBuffQP.Value[0], accInLoop.Value[0], ctBuffQP.Value[0])

			// Outer loop rotations
			if cnt0 == 0 {
				ringQP.PermuteNTTWithIndexLvl(levelQ, levelP, ctBuffQP.Value[0], rotIndex, accOutLoopQP.Value[0])
				ringQP.PermuteNTTWithIndexLvl(levelQ, levelP, ctBuffQP.Value[1], rotIndex, accOutLoopQP.Value[1])
			} else {
				ringQP.PermuteNTTWithIndexAndAddNoModLvl(levelQ, levelP, ctBuffQP.Value[0], rotIndex, accOutLoopQP.Value[0])
				ringQP.PermuteNTTWithIndexAndAddNoModLvl(levelQ, levelP, ctBuffQP.Value[1], rotIndex, accOutLoopQP.Value[1])
			}

			// Else directly adds on accOutLoopQP
		} else {

			if cnt0 == 0 {
				ringQP.CopyValuesLvl(levelQ, levelP, accInLoop.Value[0], accOutLoopQP.Value[0])
				ringQP.CopyValuesLvl(levelQ, levelP, accInLoop.Value[1], accOutLoopQP.Value[1])
			} else {
				ringQP.AddNoModLvl(levelQ, levelP, accOutLoopQP.Value[0], accInLoop.Value[0], accOutLoopQP.Value[0])
				ringQP.AddNoModLvl(levelQ, levelP, accOutLoopQP.Value[1], accInLoop.Value[1], accOutLoopQP.Value[1])
			}
		}

		if cnt0%QiOverF == QiOverF-1 {
			ringQ.ReduceLvl(levelQ, ctOut.Value[0], ctOut.Value[0])
			ringQ.ReduceLvl(levelQ, ctOut.Value[1], ctOut.Value[1])
		}

		if cnt0%PiOverF == PiOverF-1 {
			ringP.ReduceLvl(levelP, accOutLoopQP.Value[0].P, accOutLoopQP.Value[0].P)
			ringP.ReduceLvl(levelP, accOutLoopQP.Value[1].P, accOutLoopQP.Value[1].P)
		}

		cnt0++
	}

	if cnt0%QiOverF != 0 {
		ringQ.ReduceLvl(levelQ, ctOut.Value[0], ctOut.Value[0])
		ringQ.ReduceLvl(levelQ, ctOut.Value[1], ctOut.Value[1])
	}

	if cnt0%PiOverF != 0 {
		ringP.ReduceLvl(levelP, accOutLoopQP.Value[0].P, accOutLoopQP.Value[0].P)
		ringP.ReduceLvl(levelP, accOutLoopQP.Value[1].P, accOutLoopQP.Value[1].P)
	}

	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, ctOut.Value[0], accOutLoopQP.Value[0].P, ctOut.Value[0]) // sum(phi(c0 * P + d0_QP))/P
	eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, ctOut.Value[1], accOutLoopQP.Value[1].P, ctOut.Value[1]) // sum(phi(d1_QP))/P

	ringQ.MulScalarLvl(levelQ, ctOut.Value[0], eval.params.T(), ctOut.Value[0])
	ringQ.MulScalarLvl(levelQ, ctOut.Value[1], eval.params.T(), ctOut.Value[1])

	ctOut.SetScale(matrix.Scale * ctIn.Scale())

	ctInRotQP = nil
	runtime.GC()
}

func innerLoopBSGSDegree1(j int, index []int, levelQ, levelP int, ringQP ringqp.Ring, diagMatrix map[int]OperandQP, ctInMulP *Ciphertext, ctInRotQP map[int]rlwe.CiphertextQP, accInLoop rlwe.CiphertextQP, qiOverFlow, piOverFlow int) {

	var reduce int

	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	for _, i := range index {
		if i == 0 {
			if reduce == 0 {
				ringQ.MulCoeffsMontgomeryConstantLvl(levelQ, diagMatrix[j].El().Value[0].Q, ctInMulP.Value[0], accInLoop.Value[0].Q)
				ringQ.MulCoeffsMontgomeryConstantLvl(levelQ, diagMatrix[j].El().Value[0].Q, ctInMulP.Value[1], accInLoop.Value[1].Q)
				accInLoop.Value[0].P.Zero()
				accInLoop.Value[1].P.Zero()
			} else {
				ringQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, diagMatrix[j].El().Value[0].Q, ctInMulP.Value[0], accInLoop.Value[0].Q)
				ringQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, diagMatrix[j].El().Value[0].Q, ctInMulP.Value[1], accInLoop.Value[1].Q)
			}
		} else {
			if reduce == 0 {
				ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, diagMatrix[j+i].El().Value[0], ctInRotQP[i].Value[0], accInLoop.Value[0])
				ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, diagMatrix[j+i].El().Value[0], ctInRotQP[i].Value[1], accInLoop.Value[1])
			} else {
				ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, diagMatrix[j+i].El().Value[0], ctInRotQP[i].Value[0], accInLoop.Value[0])
				ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, diagMatrix[j+i].El().Value[0], ctInRotQP[i].Value[1], accInLoop.Value[1])
			}
		}

		if reduce%qiOverFlow == qiOverFlow-1 {
			ringQ.ReduceLvl(levelQ, accInLoop.Value[0].Q, accInLoop.Value[0].Q)
			ringQ.ReduceLvl(levelQ, accInLoop.Value[1].Q, accInLoop.Value[1].Q)
		}

		if reduce%piOverFlow == piOverFlow-1 {
			ringP.ReduceLvl(levelP, accInLoop.Value[0].P, accInLoop.Value[0].P)
			ringP.ReduceLvl(levelP, accInLoop.Value[1].P, accInLoop.Value[1].P)
		}

		reduce++
	}

	if reduce%qiOverFlow != 0 {
		ringQ.ReduceLvl(levelQ, accInLoop.Value[0].Q, accInLoop.Value[0].Q)
		ringQ.ReduceLvl(levelQ, accInLoop.Value[1].Q, accInLoop.Value[1].Q)
	}

	if reduce%piOverFlow != 0 {
		ringP.ReduceLvl(levelP, accInLoop.Value[0].P, accInLoop.Value[0].P)
		ringP.ReduceLvl(levelP, accInLoop.Value[1].P, accInLoop.Value[1].P)
	}
}

func innerLoopBSGSDegree2(j int, index []int, levelQ, levelP int, ringQP ringqp.Ring, diagMatrix map[int]OperandQP, ctInMulP *Ciphertext, ctInRotQP map[int]rlwe.CiphertextQP, accInLoop rlwe.CiphertextQP, qiOverFlow, piOverFlow int) {

	var reduce int

	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	for _, i := range index {
		if i == 0 {
			if reduce == 0 {
				tensorDegree2Q(levelQ, ringQ, diagMatrix[j].El().Value, ctInMulP.Value, accInLoop.Value)
			} else {
				tensorDegree2AndAddQ(levelQ, ringQ, diagMatrix[j].El().Value, ctInMulP.Value, accInLoop.Value)
			}
		} else {
			if reduce == 0 {
				tensorDegree2QP(levelQ, levelP, ringQP, diagMatrix[j+i].El().Value, ctInRotQP[i].Value, accInLoop.Value)
			} else {
				tensorDegree2AndAddQP(levelQ, levelP, ringQP, diagMatrix[j+i].El().Value, ctInRotQP[i].Value, accInLoop.Value)
			}
		}

		if reduce%qiOverFlow == qiOverFlow-1 {
			ringQ.ReduceLvl(levelQ, accInLoop.Value[0].Q, accInLoop.Value[0].Q)
			ringQ.ReduceLvl(levelQ, accInLoop.Value[1].Q, accInLoop.Value[1].Q)
			ringQ.ReduceLvl(levelQ, accInLoop.Value[2].Q, accInLoop.Value[2].Q)
		}

		if reduce%piOverFlow == piOverFlow-1 {
			ringP.ReduceLvl(levelP, accInLoop.Value[0].P, accInLoop.Value[0].P)
			ringP.ReduceLvl(levelP, accInLoop.Value[1].P, accInLoop.Value[1].P)
			ringP.ReduceLvl(levelP, accInLoop.Value[2].P, accInLoop.Value[2].P)
		}

		reduce++
	}

	if reduce%qiOverFlow != 0 {
		ringQ.ReduceLvl(levelQ, accInLoop.Value[0].Q, accInLoop.Value[0].Q)
		ringQ.ReduceLvl(levelQ, accInLoop.Value[1].Q, accInLoop.Value[1].Q)
		ringQ.ReduceLvl(levelQ, accInLoop.Value[2].Q, accInLoop.Value[2].Q)
	}

	if reduce%piOverFlow != 0 {
		ringP.ReduceLvl(levelP, accInLoop.Value[0].P, accInLoop.Value[0].P)
		ringP.ReduceLvl(levelP, accInLoop.Value[1].P, accInLoop.Value[1].P)
		ringP.ReduceLvl(levelP, accInLoop.Value[2].P, accInLoop.Value[2].P)
	}
}

func tensorDegree2Q(level int, r *ring.Ring, op0 []ringqp.Poly, op1 []*ring.Poly, op2 []ringqp.Poly) {
	r.MulCoeffsMontgomeryLvl(level, op0[0].Q, op1[0], op2[0].Q)
	r.MulCoeffsMontgomeryLvl(level, op0[0].Q, op1[1], op2[1].Q)
	r.MulCoeffsMontgomeryAndAddLvl(level, op0[1].Q, op1[0], op2[1].Q)
	r.MulCoeffsMontgomeryLvl(level, op0[1].Q, op1[1], op2[2].Q)
	op2[0].P.Zero()
	op2[1].P.Zero()
	op2[2].P.Zero()
}

func tensorDegree2AndAddQ(level int, r *ring.Ring, op0 []ringqp.Poly, op1 []*ring.Poly, op2 []ringqp.Poly) {
	r.MulCoeffsMontgomeryAndAddLvl(level, op0[0].Q, op1[0], op2[0].Q)
	r.MulCoeffsMontgomeryAndAddLvl(level, op0[0].Q, op1[1], op2[1].Q)
	r.MulCoeffsMontgomeryAndAddLvl(level, op0[1].Q, op1[0], op2[1].Q)
	r.MulCoeffsMontgomeryAndAddLvl(level, op0[1].Q, op1[1], op2[2].Q)
}

func tensorDegree2QP(levelQ, levelP int, r ringqp.Ring, op0, op1, op2 []ringqp.Poly) {
	r.MulCoeffsMontgomeryLvl(levelQ, levelP, op0[0], op1[0], op2[0])
	r.MulCoeffsMontgomeryLvl(levelQ, levelP, op0[0], op1[1], op2[1])
	r.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, op0[1], op1[0], op2[1])
	r.MulCoeffsMontgomeryLvl(levelQ, levelP, op0[1], op1[1], op2[2])
}

func tensorDegree2AndAddQP(levelQ, levelP int, r ringqp.Ring, op0, op1, op2 []ringqp.Poly) {
	r.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, op0[0], op1[0], op2[0])
	r.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, op0[0], op1[1], op2[1])
	r.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, op0[1], op1[0], op2[1])
	r.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, op0[1], op1[1], op2[2])
}
