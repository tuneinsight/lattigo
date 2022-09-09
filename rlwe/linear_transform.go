package rlwe

import (
	"runtime"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// LinearTransform is a type for linear transformations on ciphertexts.
// It stores a plaintext matrix diagonalized in diagonal form and
// can be evaluated on a ciphertext by using the evaluator.LinearTransform method.
type LinearTransform struct {
	LogSlots int               // Log of the number of slots of the plaintext (needed to compute the appropriate rotation keys)
	N1       int               // N1 is the number of inner loops of the baby-step giant-step algorithm used in the evaluation (if N1 == 0, BSGS is not used).
	Level    int               // Level is the level at which the matrix is encoded (can be circuit dependent)
	Vec      map[int]OperandQP // Vec is the matrix, in diagonal form, where each entry of vec is an indexed non-zero diagonal.
}

func (LT *LinearTransform) IsPlaintext() bool {
	for _, el := range LT.Vec {
		switch el.(type) {
		case *PlaintextQP:
			return true
		default:
			return false
		}
	}
	return false
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
	case map[int]*PlaintextQP:
		nonZeroDiags = make([]int, len(element))
		var i int
		for key := range element {
			nonZeroDiags[i] = key
			i++
		}
	case map[int]*CiphertextQP:
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

// MultiplyByDiagMatrix multiplies the ciphertext "ctIn" by the plaintext matrix "matrix" and returns the result on the ciphertext
// "ctOut". Memory buffers for the decomposed ciphertext BuffDecompQP, BuffDecompQP must be provided, those are list of poly of ringQ and ringP
// respectively, each of size params.Beta().
// The naive approach is used (single hoisting and no baby-step giant-step), which is faster than MultiplyByDiagMatrixBSGS
// for matrix of only a few non-zero diagonals but uses more keys.
func (eval *Evaluator) MultiplyByDiagMatrix(ctIn *Ciphertext, matrix LinearTransform, BuffDecompQP []ringqp.Poly, ctOut *Ciphertext) {

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
	ksResQP := CiphertextQP{Value: []ringqp.Poly{ksRes0QP, ksRes1QP}}

	ring.CopyValuesLvl(levelQ, ctIn.Value[0], eval.BuffCt.Value[0])
	ring.CopyValuesLvl(levelQ, ctIn.Value[1], eval.BuffCt.Value[1])
	ctInTmp0, ctInTmp1 := eval.BuffCt.Value[0], eval.BuffCt.Value[1]

	ringQ.MulScalarBigintLvl(levelQ, ctInTmp0, ringP.ModulusAtLevel[levelP], ct0TimesP) // P*c0

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

	if state { // Rotation by zero
		ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.Vec[0].El().Value[0].Q, ctInTmp0, c0OutQP.Q) // ctOut += c0_Q * plaintext
		ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.Vec[0].El().Value[0].Q, ctInTmp1, c1OutQP.Q) // ctOut += c1_Q * plaintext
	}
}

// MultiplyByDiagMatrixBSGS multiplies the ciphertext "ctIn" by the plaintext matrix "matrix" and returns the result on the ciphertext
// "ctOut". Memory buffers for the decomposed ciphertext BuffDecompQP, BuffDecompQP must be provided, those are list of poly of ringQ and ringP
// respectively, each of size params.Beta().
// The BSGS approach is used (double hoisting with baby-step giant-step), which is faster than MultiplyByDiagMatrix
// for matrix with more than a few non-zero diagonals and uses much less keys.
func (eval *Evaluator) MultiplyByDiagMatrixBSGS(ctIn *Ciphertext, matrix LinearTransform, PoolDecompQP []ringqp.Poly, ctOut *Ciphertext) {

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
	accInLoop := CiphertextQP{Value: []ringqp.Poly{eval.BuffQP[2], eval.BuffQP[3]}}

	// Available & Used Buffers:
	// Evaluator CKKS:
	//		buffQ[0] -> accInLoop.Value[2].Q
	//		buffQ[1] -> accInLoop.Value[2].Q
	//		buffQ[2]
	//		BuffCt   -> ctInMulP
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
		accInLoop.Value = append(accInLoop.Value, eval.BuffQP[5])
	}

	// Result in QP
	accOutLoopQP := &CiphertextQP{Value: []ringqp.Poly{ringqp.Poly{Q: ctOut.Value[0], P: eval.BuffQP[4].Q}, ringqp.Poly{Q: ctOut.Value[1], P: eval.BuffQP[4].P}}}

	// Ciphertext buffer mod QP for gadget products without mod down
	ctBuffQP := CiphertextQP{Value: []ringqp.Poly{eval.BuffQP[0], eval.BuffQP[1]}}

	// Computes the N2 rotations indexes of the non-zero rows of the diagonalized DFT matrix for the baby-step giang-step algorithm
	index, _, rotN2 := BsgsIndex(matrix.Vec, 1<<matrix.LogSlots, matrix.N1)

	// Pre-rotates ciphertext for the baby-step giant-step algorithm, does not divide by P yet
	ctInRotQP := make(map[int]CiphertextQP)
	for _, i := range rotN2 {
		if i != 0 {
			ctInRotQP[i] = NewCiphertextQP(eval.params, 1, levelQ, levelP)
			eval.AutomorphismHoistedNoModDown(levelQ, ctIn.Value[0], PoolDecompQP, eval.params.GaloisElementForColumnRotationBy(i), ctInRotQP[i])
		}
	}

	// Rotation by 0 -> pre-multiplies by P
	ctInMulP := eval.BuffCt
	ringQ.MulScalarBigintLvl(levelQ, ctIn.Value[0], ringP.ModulusAtLevel[levelP], ctInMulP.Value[0]) // P*c0
	ringQ.MulScalarBigintLvl(levelQ, ctIn.Value[1], ringP.ModulusAtLevel[levelP], ctInMulP.Value[1]) // P*c1

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

	ctInRotQP = nil
	runtime.GC()
}

func innerLoopBSGSDegree1(j int, index []int, levelQ, levelP int, ringQP ringqp.Ring, diagMatrix map[int]OperandQP, ctInMulP Ciphertext, ctInRotQP map[int]CiphertextQP, accInLoop CiphertextQP, qiOverFlow, piOverFlow int) {

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

func innerLoopBSGSDegree2(j int, index []int, levelQ, levelP int, ringQP ringqp.Ring, diagMatrix map[int]OperandQP, ctInMulP Ciphertext, ctInRotQP map[int]CiphertextQP, accInLoop CiphertextQP, qiOverFlow, piOverFlow int) {

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
