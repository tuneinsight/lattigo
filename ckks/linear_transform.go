package ckks

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

func (eval *evaluator) LinearTransform(vec *Ciphertext, linearTransform interface{}, rotkeys *RotationKeys) (res []*Ciphertext) {

	levelQ := vec.Level()

	eval.DecompInternal(levelQ, vec.value[1], eval.c2QiQDecomp, eval.c2QiPDecomp)

	switch element := linearTransform.(type) {
	case []*PtDiagMatrix:
		res = make([]*Ciphertext, len(element))
		for i, matrix := range element {
			eval.multiplyByDiabMatrix(vec, res[i], matrix, rotkeys, eval.c2QiQDecomp, eval.c2QiPDecomp)
		}
	case *PtDiagMatrix:
		res = []*Ciphertext{NewCiphertext(eval.params, 1, vec.Level(), vec.Scale())}
		eval.multiplyByDiabMatrix(vec, res[0], element, rotkeys, eval.c2QiQDecomp, eval.c2QiPDecomp)
	}

	return
}

func (eval *evaluator) multiplyByDiabMatrix(vec, res *Ciphertext, matrix *PtDiagMatrix, rotKeys *RotationKeys, c2QiQDecomp, c2QiPDecomp []*ring.Poly) {

	if matrix.rotOnly{
		for i := range matrix.Vec{
			eval.permuteNTTHoisted(vec, c2QiQDecomp, c2QiPDecomp, i, rotKeys, res)
		}
	}else{
		if matrix.naive {
			eval.multiplyByDiabMatrixNaive(vec, res, matrix, rotKeys, c2QiQDecomp, c2QiPDecomp)
		} else {
			eval.multiplyByDiabMatrixBSGS(vec, res, matrix, rotKeys, c2QiQDecomp, c2QiPDecomp)
		}
	}

	return
}

func (eval *evaluator) multiplyByDiabMatrixNaive(vec, res *Ciphertext, matrix *PtDiagMatrix, rotKeys *RotationKeys, c2QiQDecomp, c2QiPDecomp []*ring.Poly) {

	ringQ := eval.ringQ
	ringP := eval.ringP

	levelQ := vec.Level()
	levelP := eval.params.PiCount() - 1

	ksResP0 := eval.poolP[0]  // Key-Switch res[0] mod P
	ksResP1 := eval.poolP[1]  // Key-Switch res[1] mod P
	tmpP0 := eval.poolP[2]    // Automorphism not-inplace pool res[0] mod P
	tmpP1 := eval.poolQMul[0] // Automorphism not-inplace pool res[1] mod P
	accP0 := eval.poolP[3]    // Accumulator res[0] mod P
	accP1 := eval.poolP[4]    // Accumulator res[1] mod P

	ct0TimesP := eval.poolQ[0] // ct0 * P mod Q
	ksResQ0 := eval.poolQ[1]   // Key-Switch res[0] mod Q
	ksResQ1 := eval.poolQ[2]   // Key-Switch res[0] mod Q
	tmpQ0 := eval.poolQ[3]     // Automorphism not-inplace pool res[0] mod Q
	tmpQ1 := eval.poolQ[4]     // Automorphism not-inplace pool res[1] mod Q

	ringQ.MulScalarBigintLvl(levelQ, vec.value[0], ringP.ModulusBigint, ct0TimesP) // P*c0

	state := false
	cnt := 0
	for k := range matrix.Vec {

		k &= ((ringQ.N >> 1) - 1)

		if k == 0 {
			state = true
		} else {

			eval.keyswitchHoistedNoModDown(levelQ, c2QiQDecomp, c2QiPDecomp, rotKeys.evakeyRotColLeft[k], ksResQ0, ksResQ1, ksResP0, ksResP1)

			ringQ.AddLvl(levelQ, ksResQ0, ct0TimesP, ksResQ0) // phi(d0_Q) += phi(P*c0)

			ring.PermuteNTTWithIndexLvl(levelQ, ksResQ0, rotKeys.permuteNTTLeftIndex[k], tmpQ0) // phi(P*c0 + d0_Q)
			ring.PermuteNTTWithIndexLvl(levelQ, ksResQ1, rotKeys.permuteNTTLeftIndex[k], tmpQ1) // phi(       d1_Q)

			ring.PermuteNTTWithIndexLvl(levelP, ksResP0, rotKeys.permuteNTTLeftIndex[k], tmpP0) // phi(P*c0 + d0_P)
			ring.PermuteNTTWithIndexLvl(levelP, ksResP1, rotKeys.permuteNTTLeftIndex[k], tmpP1) // phi(       d1_P)

			plaintextQ := matrix.Vec[k][0]
			plaintextP := matrix.Vec[k][1]

			if cnt == 0 {
				// keyswitch(c1_Q) = (d0_QP, d1_QP)
				ringQ.MulCoeffsMontgomeryLvl(levelQ, plaintextQ, tmpQ0, res.value[0]) // phi(P*c0 + d0_Q) * plaintext
				ringQ.MulCoeffsMontgomeryLvl(levelQ, plaintextQ, tmpQ1, res.value[1]) // phi(d1_Q) * plaintext
				ringP.MulCoeffsMontgomery(plaintextP, tmpP0, accP0)                   // phi(d0_P) * plaintext
				ringP.MulCoeffsMontgomery(plaintextP, tmpP1, accP1)                   // phi(d1_P) * plaintext
			} else {
				// keyswitch(c1_Q) = (d0_QP, d1_QP)
				ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, plaintextQ, tmpQ0, res.value[0]) // phi(P*c0 + d0_Q) * plaintext
				ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, plaintextQ, tmpQ1, res.value[1]) // phi(d1_Q) * plaintext
				ringP.MulCoeffsMontgomeryAndAdd(plaintextP, tmpP0, accP0)                   // phi(d0_P) * plaintext
				ringP.MulCoeffsMontgomeryAndAdd(plaintextP, tmpP1, accP1)                   // phi(d1_P) * plaintext
			}

			cnt++
		}
	}

	eval.baseconverter.ModDownSplitNTTPQ(levelQ, res.value[0], accP0, res.value[0]) // sum(phi(c0 * P + d0_QP))/P
	eval.baseconverter.ModDownSplitNTTPQ(levelQ, res.value[1], accP1, res.value[1]) // sum(phi(d1_QP))/P

	if state { // Rotation by zero
		ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.Vec[0][0], vec.value[0], res.value[0]) // res += c0_Q * plaintext
		ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.Vec[0][0], vec.value[1], res.value[1]) // res += c1_Q * plaintext
	}

	res.SetScale(matrix.Scale * vec.Scale())
}

func (eval *evaluator) multiplyByDiabMatrixBSGS(vec, res *Ciphertext, matrix *PtDiagMatrix, rotKeys *RotationKeys, c2QiQDecomp, c2QiPDecomp []*ring.Poly) {

	// N1*N2 = N
	N1 := matrix.N1

	ringQ := eval.ringQ
	ringP := eval.ringP

	levelQ := vec.Level()
	levelP := eval.params.PiCount() - 1

	// Computes the rotations indexes of the non-zero rows of the diagonalized DFT matrix for the baby-step giang-step algorithm
	index := make(map[uint64][]uint64)
	rotations := []uint64{}

	for key := range matrix.Vec {

		idx1 := key / N1
		idx2 := key & (N1 - 1)

		if index[idx1] == nil {
			index[idx1] = []uint64{idx2}
		} else {
			index[idx1] = append(index[idx1], idx2)
		}

		if !utils.IsInSliceUint64(idx2, rotations) {
			rotations = append(rotations, idx2)
		}
	}

	// Pre-rotates ciphertext for the baby-step giant-step algorithm, does not divide by P yet
	vecRotQ, vecRotP := eval.rotateHoistedNoModDown(vec, rotations, eval.c2QiQDecomp, eval.c2QiPDecomp, rotKeys)

	// Accumulator inner loop
	tmpQ0 := eval.poolQMul[0] // unused memory pool from evaluator
	tmpQ1 := eval.poolQMul[1] // unused memory pool from evaluator

	// Accumulator outer loop
	tmpQ2 := eval.poolQMul[2] // unused memory pool from evaluator
	tmpQ3 := eval.poolQ[4]
	tmpP2 := eval.poolP[3]
	tmpP3 := eval.poolP[4]

	// Keyswitch accumulator
	pool2Q := eval.poolQ[1] // res(c0', c1') from evaluator keyswitch memory pool
	pool3Q := eval.poolQ[2] // res(c0', c1') from evaluator keyswitch memory pool
	pool2P := eval.poolP[1] // res(c0', c1') from evaluator keyswitch memory pool
	pool3P := eval.poolP[2] // res(c0', c1') from evaluator keyswitch memory pool

	// Do not use (used by switchKeysInPlaceNoModDown)
	// eval.PoolP[0]
	// eval.PoolQ[0]
	// eval.PoolQ[2]

	N1Rot := 0
	N2Rot := 0

	ringQ.MulScalarBigintLvl(levelQ, vec.value[0], ringP.ModulusBigint, tmpQ0) // P*c0

	for _, i := range rotations {
		if i != 0 {
			ring.PermuteNTTWithIndexLvl(levelQ, tmpQ0, rotKeys.permuteNTTLeftIndex[i], tmpQ1) // phi(P*c0)
			ringQ.AddLvl(levelQ, vecRotQ[i][0], tmpQ1, vecRotQ[i][0])                         // phi(d0_Q) += phi(P*c0)
		}
	}

	// OUTER LOOP
	cnt0 := 0
	for j := range index {

		if j != 0 {

			// INNER LOOP
			state := false
			cnt1 := 0
			for _, i := range index[j] {

				if i == 0 {
					state = true
				} else {

					N1Rot++

					plaintextQ := matrix.Vec[N1*j+uint64(i)][0]
					plaintextP := matrix.Vec[N1*j+uint64(i)][1]

					if cnt1 == 0 {
						ringQ.MulCoeffsMontgomeryLvl(levelQ, plaintextQ, vecRotQ[i][0], tmpQ0) // phi(P*c0 + d0_Q) * plaintext
						ringQ.MulCoeffsMontgomeryLvl(levelQ, plaintextQ, vecRotQ[i][1], tmpQ1) // phi(d1_Q) * plaintext
						ringP.MulCoeffsMontgomery(plaintextP, vecRotP[i][0], pool2P)           // phi(d0_P) * plaintext
						ringP.MulCoeffsMontgomery(plaintextP, vecRotP[i][1], pool3P)           // phi(d1_P) * plaintext
					} else {
						ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, plaintextQ, vecRotQ[i][0], tmpQ0) // phi(d0_Q) * plaintext
						ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, plaintextQ, vecRotQ[i][1], tmpQ1) // phi(d1_Q) * plaintext
						ringP.MulCoeffsMontgomeryAndAdd(plaintextP, vecRotP[i][0], pool2P)           // phi(d0_P) * plaintext
						ringP.MulCoeffsMontgomeryAndAdd(plaintextP, vecRotP[i][1], pool3P)           // phi(d1_P) * plaintext
					}

					cnt1++
				}
			}

			// Hoisting of the ModDown of sum(sum(phi(d0 + P*c0) * plaintext)) and sum(sum(phi(d1) * plaintext))
			eval.baseconverter.ModDownSplitNTTPQ(levelQ, tmpQ0, pool2P, tmpQ0) // sum(phi(d0) * plaintext)/P
			eval.baseconverter.ModDownSplitNTTPQ(levelQ, tmpQ1, pool3P, tmpQ1) // sum(phi(d1) * plaintext)/P

			// If i == 0
			if state {

				// If no loop before, then we copy the values on the accumulator instead of adding them
				if len(index[j]) == 1 {
					ringQ.MulCoeffsMontgomeryLvl(levelQ, matrix.Vec[N1*j][0], vec.value[0], tmpQ0) // c0 * plaintext + sum(phi(d0) * plaintext)/P + phi(c0) * plaintext mod Q
					ringQ.MulCoeffsMontgomeryLvl(levelQ, matrix.Vec[N1*j][0], vec.value[1], tmpQ1) // c1 * plaintext + sum(phi(d1) * plaintext)/P + phi(c1) * plaintext mod Q
				} else {
					ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.Vec[N1*j][0], vec.value[0], tmpQ0) // c0 * plaintext + sum(phi(d0) * plaintext)/P + phi(c0) * plaintext mod Q
					ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.Vec[N1*j][0], vec.value[1], tmpQ1) // c1 * plaintext + sum(phi(d1) * plaintext)/P + phi(c1) * plaintext mod Q
				}

				N1Rot++
			}

			eval.switchKeysInPlaceNoModDown(levelQ, tmpQ1, rotKeys.evakeyRotColLeft[N1*j], pool2Q, pool2P, pool3Q, pool3P) // Switchkey(phi(tmpRes_1)) = (d0, d1) in base QP

			// Outer loop rotations
			ring.PermuteNTTWithIndexLvl(levelQ, tmpQ0, rotKeys.permuteNTTLeftIndex[N1*j], tmpQ1) // phi(tmpRes_0)
			ringQ.AddLvl(levelQ, res.value[0], tmpQ1, res.value[0])                              // res += phi(tmpRes)

			rot := rotKeys.permuteNTTLeftIndex[N1*j]

			N2Rot++

			if cnt0 == 0 {
				ring.PermuteNTTWithIndexLvl(levelQ, pool2Q, rot, tmpQ2) // sum(phi(d0_Q))
				ring.PermuteNTTWithIndexLvl(levelQ, pool3Q, rot, tmpQ3) // sum(phi(d1_Q))
				ring.PermuteNTTWithIndexLvl(levelP, pool2P, rot, tmpP2) // sum(phi(d0_P))
				ring.PermuteNTTWithIndexLvl(levelP, pool3P, rot, tmpP3) // sum(phi(d1_P))
			} else {
				ring.PermuteNTTWithIndexAndAddNoModLvl(levelQ, pool2Q, rot, tmpQ2) // sum(phi(d0_Q))
				ring.PermuteNTTWithIndexAndAddNoModLvl(levelQ, pool3Q, rot, tmpQ3) // sum(phi(d1_Q))
				ring.PermuteNTTWithIndexAndAddNoModLvl(levelP, pool2P, rot, tmpP2) // sum(phi(d0_P))
				ring.PermuteNTTWithIndexAndAddNoModLvl(levelP, pool3P, rot, tmpP3) // sum(phi(d1_P))
			}

			if cnt0 == 7 {
				ringQ.ReduceLvl(levelQ, tmpQ2, tmpQ2)
				ringQ.ReduceLvl(levelQ, tmpQ3, tmpQ3)
				ringP.Reduce(tmpP2, tmpP2)
				ringP.Reduce(tmpP3, tmpP3)
			}

			cnt0++
		}
	}

	if cnt0 != 7 {
		ringQ.ReduceLvl(levelQ, tmpQ2, tmpQ2)
		ringQ.ReduceLvl(levelQ, tmpQ3, tmpQ3)
		ringP.Reduce(tmpP2, tmpP2)
		ringP.Reduce(tmpP3, tmpP3)
	}

	// if j == 0 (N2 rotation by zero)
	state := false
	for _, i := range index[0] {

		if i == 0 {
			state = true
		} else {

			plaintextQ := matrix.Vec[uint64(i)][0]
			plaintextP := matrix.Vec[uint64(i)][1]
			N1Rot++
			// keyswitch(c1_Q) = (d0_QP, d1_QP)
			ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, plaintextQ, vecRotQ[i][0], tmpQ2) // phi(P*c0 + d0_Q) * plaintext
			ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, plaintextQ, vecRotQ[i][1], tmpQ3) // phi(d1_Q) * plaintext
			ringP.MulCoeffsMontgomeryAndAdd(plaintextP, vecRotP[i][0], tmpP2)            // phi(d0_P) * plaintext
			ringP.MulCoeffsMontgomeryAndAdd(plaintextP, vecRotP[i][1], tmpP3)            // phi(d1_P) * plaintext
		}
	}

	eval.baseconverter.ModDownSplitNTTPQ(levelQ, tmpQ2, tmpP2, tmpQ2) // sum(phi(c0 * P + d0_QP))/P
	eval.baseconverter.ModDownSplitNTTPQ(levelQ, tmpQ3, tmpP3, tmpQ3) // sum(phi(d1_QP))/P

	ringQ.AddLvl(levelQ, res.value[0], tmpQ2, res.value[0]) // res += sum(phi(c0 * P + d0_QP))/P
	ringQ.AddLvl(levelQ, res.value[1], tmpQ3, res.value[1]) // res += sum(phi(d1_QP))/P

	if state { // Rotation by zero
		N1Rot++
		ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.Vec[0][0], vec.value[0], res.value[0]) // res += c0_Q * plaintext
		ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.Vec[0][0], vec.value[1], res.value[1]) // res += c1_Q * plaintext
	}

	res.SetScale(matrix.Scale * vec.Scale())

	vecRotQ, vecRotP = nil, nil
}
