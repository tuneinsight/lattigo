package ckks

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// Trace maps X -> sum((-1)^i * X^{i*n+1}) for 0 <= i < N
// For log(n) = logSlotStart and log(N/2) = logSlotsEnd
func (eval *evaluator) Trace(ctIn *Ciphertext, logSlotsStart, logSlotsEnd int, ctOut *Ciphertext) {
	if ctIn != ctOut {
		ctOut.Copy(ctIn)
	} else {
		ctOut = ctIn
	}
	for i := logSlotsStart; i < logSlotsEnd; i++ {
		eval.permuteNTT(ctOut, eval.params.GaloisElementForColumnRotationBy(1<<i), eval.ctxpool)
		ctPool := &Ciphertext{Ciphertext: eval.ctxpool.Ciphertext, Scale: ctOut.Scale}
		ctPool.Value = ctPool.Value[:2]
		eval.Add(ctOut, ctPool, ctOut)
	}
}

// TraceNew maps X -> sum((-1)^i * X^{i*n+1}) for 0 <= i < N and returns the result on a new ciphertext.
// For log(n) = logSlotStart and log(N/2) = logSlotsEnd
func (eval *evaluator) TraceNew(ctIn *Ciphertext, logSlotsStart, logSlotsEnd int) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1, ctIn.Level(), ctIn.Scale)
	eval.Trace(ctIn, logSlotsStart, logSlotsEnd, ctOut)
	return
}

// RotateHoistedNew takes an input Ciphertext and a list of rotations and returns a map of Ciphertext, where each element of the map is the input Ciphertext
// rotation by one element of the list. It is much faster than sequential calls to Rotate.
func (eval *evaluator) RotateHoistedNew(ctIn *Ciphertext, rotations []int) (ctOut map[int]*Ciphertext) {
	ctOut = make(map[int]*Ciphertext)
	for _, i := range rotations {
		ctOut[i] = NewCiphertext(eval.params, 1, ctIn.Level(), ctIn.Scale)
	}
	eval.RotateHoisted(ctIn, rotations, ctOut)
	return
}

// RotateHoisted takes an input Ciphertext and a list of rotations and populates a map of pre-allocated Ciphertexts,
// where each element of the map is the input Ciphertext rotation by one element of the list.
// It is much faster than sequential calls to Rotate.
func (eval *evaluator) RotateHoisted(ctIn *Ciphertext, rotations []int, ctOut map[int]*Ciphertext) {
	levelQ := ctIn.Level()
	eval.DecomposeNTT(levelQ, eval.params.PCount()-1, eval.params.PCount(), ctIn.Value[1], eval.PoolDecompQP)
	for _, i := range rotations {
		if i == 0 {
			ctOut[i].Copy(ctIn)
		} else {
			eval.PermuteNTTHoisted(levelQ, ctIn.Value[0], ctIn.Value[1], eval.PoolDecompQP, i, ctOut[i].Value[0], ctOut[i].Value[1])
		}
	}
}

// LinearTransformNew evaluates a linear transform on the ciphertext and returns the result on a new ciphertext.
// The linearTransform can either be an (ordered) list of PtDiagMatrix or a single PtDiagMatrix.
// In either case a list of ciphertext is return (the second case returnign a list of
// containing a single ciphertext. A PtDiagMatrix is a diagonalized plaintext matrix contructed with an Encoder using
// the method encoder.EncodeDiagMatrixAtLvl(*).
func (eval *evaluator) LinearTransformNew(ctIn *Ciphertext, linearTransform interface{}) (ctOut []*Ciphertext) {

	switch element := linearTransform.(type) {
	case []PtDiagMatrix:
		ctOut = make([]*Ciphertext, len(element))

		var maxLevel int
		for _, matrix := range element {
			maxLevel = utils.MaxInt(maxLevel, matrix.Level)
		}

		minLevel := utils.MinInt(maxLevel, ctIn.Level())

		eval.DecomposeNTT(minLevel, eval.params.PCount()-1, eval.params.PCount(), ctIn.Value[1], eval.PoolDecompQP)

		for i, matrix := range element {
			ctOut[i] = NewCiphertext(eval.params, 1, minLevel, ctIn.Scale)

			if matrix.Naive {
				eval.MultiplyByDiagMatrix(ctIn, matrix, eval.PoolDecompQP, ctOut[i])
			} else {
				eval.MultiplyByDiagMatrixBSGS(ctIn, matrix, eval.PoolDecompQP, ctOut[i])
			}
		}

	case PtDiagMatrix:

		minLevel := utils.MinInt(element.Level, ctIn.Level())
		eval.DecomposeNTT(minLevel, eval.params.PCount()-1, eval.params.PCount(), ctIn.Value[1], eval.PoolDecompQP)

		ctOut = []*Ciphertext{NewCiphertext(eval.params, 1, minLevel, ctIn.Scale)}

		if element.Naive {
			eval.MultiplyByDiagMatrix(ctIn, element, eval.PoolDecompQP, ctOut[0])
		} else {
			eval.MultiplyByDiagMatrixBSGS(ctIn, element, eval.PoolDecompQP, ctOut[0])
		}
	}
	return
}

// LinearTransformNew evaluates a linear transform on the pre-allocated ciphertexts.
// The linearTransform can either be an (ordered) list of PtDiagMatrix or a single PtDiagMatrix.
// In either case a list of ciphertext is return (the second case returnign a list of
// containing a single ciphertext. A PtDiagMatrix is a diagonalized plaintext matrix contructed with an Encoder using
// the method encoder.EncodeDiagMatrixAtLvl(*).
func (eval *evaluator) LinearTransform(ctIn *Ciphertext, linearTransform interface{}, ctOut []*Ciphertext) {

	switch element := linearTransform.(type) {
	case []PtDiagMatrix:
		var maxLevel int
		for _, matrix := range element {
			maxLevel = utils.MaxInt(maxLevel, matrix.Level)
		}

		minLevel := utils.MinInt(maxLevel, ctIn.Level())

		eval.DecomposeNTT(minLevel, eval.params.PCount()-1, eval.params.PCount(), ctIn.Value[1], eval.PoolDecompQP)

		for i, matrix := range element {
			if matrix.Naive {
				eval.MultiplyByDiagMatrix(ctIn, matrix, eval.PoolDecompQP, ctOut[i])
			} else {
				eval.MultiplyByDiagMatrixBSGS(ctIn, matrix, eval.PoolDecompQP, ctOut[i])
			}
		}

	case PtDiagMatrix:
		minLevel := utils.MinInt(element.Level, ctIn.Level())
		eval.DecomposeNTT(minLevel, eval.params.PCount()-1, eval.params.PCount(), ctIn.Value[1], eval.PoolDecompQP)
		if element.Naive {
			eval.MultiplyByDiagMatrix(ctIn, element, eval.PoolDecompQP, ctOut[0])
		} else {
			eval.MultiplyByDiagMatrixBSGS(ctIn, element, eval.PoolDecompQP, ctOut[0])
		}
	}
}

// InnerSumLog applies an optimized inner sum on the ciphetext (log2(n) + HW(n) rotations with double hoisting).
// The operation assumes that `ctIn` encrypts SlotCount/`batchSize` sub-vectors of size `batchSize` which it adds together (in parallel) by groups of `n`.
// It outputs in ctOut a ciphertext for which the "leftmost" sub-vector of each group is equal to the sum of the group.
// This method is faster than InnerSum when the number of rotations is large and uses log2(n) + HW(n) insteadn of 'n' keys.
func (eval *evaluator) InnerSumLog(ctIn *Ciphertext, batchSize, n int, ctOut *Ciphertext) {

	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()

	levelQ := ctIn.Level()
	levelP := len(ringP.Modulus) - 1

	//QiOverF := eval.params.QiOverflowMargin(levelQ)
	//PiOverF := eval.params.PiOverflowMargin()

	if n == 1 {
		if ctIn != ctOut {
			ring.CopyValuesLvl(levelQ, ctIn.Value[0], ctOut.Value[0])
			ring.CopyValuesLvl(levelQ, ctIn.Value[1], ctOut.Value[1])
		}
	} else {

		// Memory pool for ctIn = ctIn + rot(ctIn, 2^i) in Q
		tmpc0 := eval.poolQMul[0] // unused memory pool from evaluator
		tmpc1 := eval.poolQMul[1] // unused memory pool from evaluator
		tmpc2 := eval.poolQMul[2]

		// Accumulator outer loop for ctOut = ctOut + rot(ctIn, k) in QP
		ct0OutQ := eval.Pool[4].Q
		ct1OutQ := eval.Pool[5].Q
		ct0OutP := eval.Pool[4].P
		ct1OutP := eval.Pool[5].P

		// Memory pool for rot(ctIn, k)
		c0Q := eval.Pool[2].Q // ctOut(c0', c1') from evaluator keyswitch memory pool
		c1Q := eval.Pool[3].Q // ctOut(c0', c1') from evaluator keyswitch memory pool
		c0P := eval.Pool[2].P // ctOut(c0', c1') from evaluator keyswitch memory pool
		c1P := eval.Pool[3].P // ctOut(c0', c1') from evaluator keyswitch memory pool

		// Used by the key-switch
		// eval.Pool[0].Q
		// eval.Pool[1].Q
		// eval.Pool[0].P
		// eval.Pool[1].P

		state := false
		copy := true
		// Binary reading of the input n
		for i, j := 0, n; j > 0; i, j = i+1, j>>1 {

			// Starts by decomposing the input ciphertext
			if i == 0 {
				// If first iteration, then copies directly from the input ciphertext that hasn't been rotated
				eval.DecomposeNTT(levelQ, levelP, levelP+1, ctIn.Value[1], eval.PoolDecompQP)
			} else {
				// Else copies from the rotated input ciphertext
				tmpc1.IsNTT = true
				eval.DecomposeNTT(levelQ, levelP, levelP+1, tmpc1, eval.PoolDecompQP)
			}

			// If the binary reading scans a 1
			if j&1 == 1 {

				k := n - (n & ((2 << i) - 1))
				k *= batchSize

				// If the rotation is not zero
				if k != 0 {

					// Rotate((tmpc0, tmpc1), k)
					eval.PermuteNTTHoistedNoModDown(levelQ, eval.PoolDecompQP, k, c0Q, c1Q, c0P, c1P)

					// ctOut += Rotate((tmpc0, tmpc1), k)
					if copy {
						ring.CopyValuesLvl(levelQ, c0Q, ct0OutQ)
						ring.CopyValuesLvl(levelQ, c1Q, ct1OutQ)
						ring.CopyValues(c0P, ct0OutP)
						ring.CopyValues(c1P, ct1OutP)
						copy = false
					} else {
						ringQ.AddLvl(levelQ, ct0OutQ, c0Q, ct0OutQ)
						ringQ.AddLvl(levelQ, ct1OutQ, c1Q, ct1OutQ)
						ringP.Add(ct0OutP, c0P, ct0OutP)
						ringP.Add(ct1OutP, c1P, ct1OutP)
					}

					if i == 0 {
						ring.PermuteNTTWithIndexLvl(levelQ, ctIn.Value[0], eval.permuteNTTIndex[eval.params.GaloisElementForColumnRotationBy(k)], tmpc2)
					} else {
						ring.PermuteNTTWithIndexLvl(levelQ, tmpc0, eval.permuteNTTIndex[eval.params.GaloisElementForColumnRotationBy(k)], tmpc2)
					}

					ringQ.MulScalarBigintLvl(levelQ, tmpc2, ringP.ModulusBigint, tmpc2)
					ringQ.AddLvl(levelQ, ct0OutQ, tmpc2, ct0OutQ)

				} else {

					state = true

					// if n is not a power of two
					if n&(n-1) != 0 {
						eval.Baseconverter.ModDownQPtoQNTT(levelQ, levelP, ct0OutQ, ct0OutP, ct0OutQ) // Division by P
						eval.Baseconverter.ModDownQPtoQNTT(levelQ, levelP, ct1OutQ, ct1OutP, ct1OutQ) // Division by P

						// ctOut += (tmpc0, tmpc1)
						ringQ.AddLvl(levelQ, ct0OutQ, tmpc0, ctOut.Value[0])
						ringQ.AddLvl(levelQ, ct1OutQ, tmpc1, ctOut.Value[1])

					} else {
						ring.CopyValuesLvl(levelQ, tmpc0, ctOut.Value[0])
						ring.CopyValuesLvl(levelQ, tmpc1, ctOut.Value[1])
						ctOut.Value[0].IsNTT = true
						ctOut.Value[1].IsNTT = true
					}
				}
			}

			if !state {
				if i == 0 {
					eval.PermuteNTTHoisted(levelQ, ctIn.Value[0], ctIn.Value[1], eval.PoolDecompQP, (1<<i)*batchSize, tmpc0, tmpc1)

					ringQ.AddLvl(levelQ, tmpc0, ctIn.Value[0], tmpc0)
					ringQ.AddLvl(levelQ, tmpc1, ctIn.Value[1], tmpc1)
				} else {
					// (tmpc0, tmpc1) = Rotate((tmpc0, tmpc1), 2^i)
					eval.PermuteNTTHoisted(levelQ, tmpc0, tmpc1, eval.PoolDecompQP, (1<<i)*batchSize, c0Q, c1Q)
					ringQ.AddLvl(levelQ, tmpc0, c0Q, tmpc0)
					ringQ.AddLvl(levelQ, tmpc1, c1Q, tmpc1)
				}

			}
		}
	}
}

// InnerSum applies an naive inner sum on the ciphetext (n rotations with single hoisting).
// The operation assumes that `ctIn` encrypts SlotCount/`batchSize` sub-vectors of size `batchSize` which it adds together (in parallel) by groups of `n`.
// It outputs in ctOut a ciphertext for which the "leftmost" sub-vector of each group is equal to the sum of the group.
// This method is faster than InnerSumLog when the number of rotations is small but uses 'n' keys instead of log(n) + HW(n).
func (eval *evaluator) InnerSum(ctIn *Ciphertext, batchSize, n int, ctOut *Ciphertext) {

	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()

	levelQ := ctIn.Level()
	levelP := len(ringP.Modulus) - 1

	QiOverF := eval.params.QiOverflowMargin(levelQ) >> 1
	PiOverF := eval.params.PiOverflowMargin(levelP) >> 1

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

		// Memory pool
		tmpQ0 := eval.poolQMul[0] // unused memory pool from evaluator
		tmpQ1 := eval.poolQMul[1] // unused memory pool from evaluator

		c0P := eval.Pool[1].P // ctOut(c0', c1') from evaluator keyswitch memory pool
		c1P := eval.Pool[2].P // ctOut(c0', c1') from evaluator keyswitch memory pool

		// Basis decomposition
		eval.DecomposeNTT(levelQ, levelP, levelP+1, ctIn.Value[1], eval.PoolDecompQP)

		// Pre-rotates all [1, ..., n-1] rotations
		// Hoisted rotation without division by P
		vecRotQ, vecRotP := eval.rotateHoistedNoModDown(levelQ, rotations, eval.PoolDecompQP)

		// P*c0 -> tmpQ0
		ringQ.MulScalarBigintLvl(levelQ, ctIn.Value[0], ringP.ModulusBigint, tmpQ0)

		// Adds phi_k(P*c0) on each of the vecRotQ
		// Does not need to add on the vecRotP because mod P === 0
		for _, i := range rotations {
			if i != 0 {

				galEl := eval.params.GaloisElementForColumnRotationBy(i)

				_, generated := eval.rtks.Keys[galEl]
				if !generated {
					panic("switching key not available")
				}

				index := eval.permuteNTTIndex[galEl]

				ring.PermuteNTTWithIndexLvl(levelQ, tmpQ0, index, tmpQ1)  // phi(P*c0)
				ringQ.AddLvl(levelQ, vecRotQ[i][0], tmpQ1, vecRotQ[i][0]) // phi(d0_Q) += phi(P*c0)
			}
		}

		var reduce int
		// Sums elements [2, ..., n-1]
		for i := 1; i < n; i++ {

			j := i * batchSize

			if i == 1 {
				ring.CopyValuesLvl(levelQ, vecRotQ[j][0], tmpQ0)
				ring.CopyValuesLvl(levelQ, vecRotQ[j][1], tmpQ1)
				ring.CopyValuesLvl(levelP, vecRotP[j][0], c0P)
				ring.CopyValuesLvl(levelP, vecRotP[j][1], c1P)
			} else {
				ringQ.AddNoModLvl(levelQ, tmpQ0, vecRotQ[j][0], tmpQ0)
				ringQ.AddNoModLvl(levelQ, tmpQ1, vecRotQ[j][1], tmpQ1)
				ringP.AddNoModLvl(levelP, c0P, vecRotP[j][0], c0P)
				ringP.AddNoModLvl(levelP, c1P, vecRotP[j][1], c1P)
			}

			if reduce%QiOverF == QiOverF-1 {
				ringQ.ReduceLvl(levelQ, tmpQ0, tmpQ0)
				ringQ.ReduceLvl(levelQ, tmpQ1, tmpQ1)
			}

			if reduce%PiOverF == PiOverF-1 {
				ringP.ReduceLvl(levelP, c0P, c0P)
				ringP.ReduceLvl(levelP, c1P, c1P)
			}

			reduce++
		}

		if reduce%QiOverF != 0 {
			ringQ.ReduceLvl(levelQ, tmpQ0, tmpQ0)
			ringQ.ReduceLvl(levelQ, tmpQ1, tmpQ1)
		}

		if reduce%PiOverF != 0 {
			ringP.ReduceLvl(levelP, c0P, c0P)
			ringP.ReduceLvl(levelP, c1P, c1P)
		}

		// Division by P of sum(elements [2, ..., n-1] )
		eval.Baseconverter.ModDownQPtoQNTT(levelQ, levelP, tmpQ0, c0P, tmpQ0) // sum_{i=1, n-1}(phi(d0))/P
		eval.Baseconverter.ModDownQPtoQNTT(levelQ, levelP, tmpQ1, c1P, tmpQ1) // sum_{i=1, n-1}(phi(d1))/P

		// Adds element[1] (which did not require rotation)
		ringQ.AddLvl(levelQ, ctIn.Value[0], tmpQ0, ctOut.Value[0]) // sum_{i=1, n-1}(phi(d0))/P + ct0
		ringQ.AddLvl(levelQ, ctIn.Value[1], tmpQ1, ctOut.Value[1]) // sum_{i=1, n-1}(phi(d1))/P + ct1
	}
}

// ReplicateLog applies an optimized replication on the ciphetext (log2(n) + HW(n) rotations with double hoisting).
// It acts as the inverse of a inner sum (summing elements from left to right).
// The replication is parameterized by the size of the sub-vectors to replicate "batchSize" and
// the number of time "n" they need to be replicated.
// To ensure correctness, a gap of zero values of size batchSize * (n-1) must exist between
// two consecutive sub-vectors to replicate.
// This method is faster than Replicate when the number of rotations is large and uses log2(n) + HW(n) instead of 'n'.
func (eval *evaluator) ReplicateLog(ctIn *Ciphertext, batchSize, n int, ctOut *Ciphertext) {
	eval.InnerSumLog(ctIn, -batchSize, n, ctOut)
}

// Replicate applies naive replication on the ciphetext (n rotations with single hoisting).
// It acts as the inverse of a inner sum (summing elements from left to right).
// The replication is parameterized by the size of the sub-vectors to replicate "batchSize" and
// the number of time "n" they need to be replicated.
// To ensure correctness, a gap of zero values of size batchSize * (n-1) must exist between
// two consecutive sub-vectors to replicate.
// This method is faster than ReplicateLog when the number of rotations is small but uses 'n' keys instead of log2(n) + HW(n).
func (eval *evaluator) Replicate(ctIn *Ciphertext, batchSize, n int, ctOut *Ciphertext) {
	eval.InnerSum(ctIn, -batchSize, n, ctOut)
}

// MultiplyByDiagMatrix multiplies the ciphertext "ctIn" by the plaintext matrix "matrix" and returns the result on the ciphertext
// "ctOut". Memory pools for the decomposed ciphertext PoolDecompQ, PoolDecompP must be provided, those are list of poly of ringQ and ringP
// respectively, each of size params.Beta().
// The naive approach is used (single hoisting and no baby-step giant-step), which is faster than MultiplyByDiagMatrixBSGS
// for matrix of only a few non-zero diagonals but uses more keys.
func (eval *evaluator) MultiplyByDiagMatrix(ctIn *Ciphertext, matrix PtDiagMatrix, PoolDecompQP []rlwe.PolyQP, ctOut *Ciphertext) {

	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()

	levelQ := utils.MinInt(ctOut.Level(), utils.MinInt(ctIn.Level(), matrix.Level))
	levelP := len(ringP.Modulus) - 1

	QiOverF := eval.params.QiOverflowMargin(levelQ)
	PiOverF := eval.params.PiOverflowMargin(levelP)

	ksResP0 := eval.Pool[0].P // Key-Switch ctOut[0] mod P
	ksResP1 := eval.Pool[1].P // Key-Switch ctOut[1] mod P
	tmpP0 := eval.Pool[2].P   // Automorphism not-inplace pool res[0] mod P
	tmpP1 := eval.poolQMul[0] // Automorphism not-inplace pool res[1] mod P
	accP0 := eval.Pool[3].P   // Accumulator ctOut[0] mod P
	accP1 := eval.Pool[4].P   // Accumulator ctOut[1] mod P

	ct0TimesP := eval.Pool[0].Q // ct0 * P mod Q
	ksResQ0 := eval.Pool[1].Q   // Key-Switch ctOut[0] mod Q
	ksResQ1 := eval.Pool[2].Q   // Key-Switch ctOut[0] mod Q
	tmpQ0 := eval.Pool[3].Q     // Automorphism not-inplace pool ctOut[0] mod Q
	tmpQ1 := eval.Pool[4].Q     // Automorphism not-inplace pool ctOut[1] mod Q

	var ctInTmp0, ctInTmp1 *ring.Poly
	if ctIn != ctOut {
		ring.CopyValuesLvl(levelQ, ctIn.Value[0], eval.ctxpool.Value[0])
		ring.CopyValuesLvl(levelQ, ctIn.Value[1], eval.ctxpool.Value[1])
		ctInTmp0, ctInTmp1 = eval.ctxpool.Value[0], eval.ctxpool.Value[1]
	} else {
		ctInTmp0, ctInTmp1 = ctIn.Value[0], ctIn.Value[1]
	}

	ringQ.MulScalarBigintLvl(levelQ, ctInTmp0, ringP.ModulusBigint, ct0TimesP) // P*c0

	var state bool
	var cnt int
	for k := range matrix.Vec {

		k &= int((ringQ.N >> 1) - 1)

		if k == 0 {
			state = true
		} else {

			galEl := eval.params.GaloisElementForColumnRotationBy(k)

			rtk, generated := eval.rtks.Keys[galEl]
			if !generated {
				panic("switching key not available")
			}

			index := eval.permuteNTTIndex[galEl]

			eval.KeyswitchHoistedNoModDown(levelQ, PoolDecompQP, rtk, ksResQ0, ksResQ1, ksResP0, ksResP1)

			ringQ.AddLvl(levelQ, ksResQ0, ct0TimesP, ksResQ0) // phi(d0_Q) += phi(P*c0)

			ring.PermuteNTTWithIndexLvl(levelQ, ksResQ0, index, tmpQ0) // phi(P*c0 + d0_Q)
			ring.PermuteNTTWithIndexLvl(levelQ, ksResQ1, index, tmpQ1) // phi(       d1_Q)

			ring.PermuteNTTWithIndexLvl(levelP, ksResP0, index, tmpP0) // phi(P*c0 + d0_P)
			ring.PermuteNTTWithIndexLvl(levelP, ksResP1, index, tmpP1) // phi(       d1_P)

			plaintextQ := matrix.Vec[k].Q
			plaintextP := matrix.Vec[k].P

			if cnt == 0 {
				// keyswitch(c1_Q) = (d0_QP, d1_QP)
				ringQ.MulCoeffsMontgomeryLvl(levelQ, plaintextQ, tmpQ0, ctOut.Value[0]) // phi(P*c0 + d0_Q) * plaintext
				ringQ.MulCoeffsMontgomeryLvl(levelQ, plaintextQ, tmpQ1, ctOut.Value[1]) // phi(d1_Q) * plaintext
				ringP.MulCoeffsMontgomeryLvl(levelP, plaintextP, tmpP0, accP0)          // phi(d0_P) * plaintext
				ringP.MulCoeffsMontgomeryLvl(levelP, plaintextP, tmpP1, accP1)          // phi(d1_P) * plaintext
			} else {
				// keyswitch(c1_Q) = (d0_QP, d1_QP)
				ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, plaintextQ, tmpQ0, ctOut.Value[0]) // phi(P*c0 + d0_Q) * plaintext
				ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, plaintextQ, tmpQ1, ctOut.Value[1]) // phi(d1_Q) * plaintext
				ringP.MulCoeffsMontgomeryAndAddLvl(levelP, plaintextP, tmpP0, accP0)          // phi(d0_P) * plaintext
				ringP.MulCoeffsMontgomeryAndAddLvl(levelP, plaintextP, tmpP1, accP1)          // phi(d1_P) * plaintext
			}

			if cnt%QiOverF == QiOverF-1 {
				ringQ.ReduceLvl(levelQ, ctOut.Value[0], ctOut.Value[0])
				ringQ.ReduceLvl(levelQ, ctOut.Value[1], ctOut.Value[1])
			}

			if cnt%PiOverF == PiOverF-1 {
				ringP.ReduceLvl(levelP, accP0, accP0)
				ringP.ReduceLvl(levelP, accP1, accP1)
			}

			cnt++
		}
	}

	if cnt%QiOverF == 0 {
		ringQ.ReduceLvl(levelQ, ctOut.Value[0], ctOut.Value[0])
		ringQ.ReduceLvl(levelQ, ctOut.Value[1], ctOut.Value[1])
	}

	if cnt%PiOverF == 0 {
		ringP.ReduceLvl(levelP, accP0, accP0)
		ringP.ReduceLvl(levelP, accP1, accP1)
	}

	eval.Baseconverter.ModDownQPtoQNTT(levelQ, levelP, ctOut.Value[0], accP0, ctOut.Value[0]) // sum(phi(c0 * P + d0_QP))/P
	eval.Baseconverter.ModDownQPtoQNTT(levelQ, levelP, ctOut.Value[1], accP1, ctOut.Value[1]) // sum(phi(d1_QP))/P

	if state { // Rotation by zero
		ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.Vec[0].Q, ctInTmp0, ctOut.Value[0]) // ctOut += c0_Q * plaintext
		ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.Vec[0].Q, ctInTmp1, ctOut.Value[1]) // ctOut += c1_Q * plaintext
	}

	ctOut.Scale = matrix.Scale * ctIn.Scale
}

// MultiplyByDiagMatrixBSGS multiplies the ciphertext "ctIn" by the plaintext matrix "matrix" and returns the result on the ciphertext
// "ctOut". Memory pools for the decomposed ciphertext PoolDecompQ, PoolDecompP must be provided, those are list of poly of ringQ and ringP
// respectively, each of size params.Beta().
// The BSGS approach is used (double hoisting with baby-step giant-step), which is faster than MultiplyByDiagMatrix
// for matrix with more than a few non-zero diagonals and uses much less keys.
func (eval *evaluator) MultiplyByDiagMatrixBSGS(ctIn *Ciphertext, matrix PtDiagMatrix, PoolDecompQP []rlwe.PolyQP, ctOut *Ciphertext) {

	// N1*N2 = N
	N1 := matrix.N1

	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()

	levelQ := utils.MinInt(ctOut.Level(), utils.MinInt(ctIn.Level(), matrix.Level))
	levelP := len(ringP.Modulus) - 1

	QiOverF := eval.params.QiOverflowMargin(levelQ)
	PiOverF := eval.params.PiOverflowMargin(levelP)

	// Computes the rotations indexes of the non-zero rows of the diagonalized DFT matrix for the baby-step giang-step algorithm

	index, rotations := bsgsIndex(matrix.Vec, 1<<matrix.LogSlots, matrix.N1)

	var ctInTmp0, ctInTmp1 *ring.Poly
	if ctIn == ctOut {
		ring.CopyValuesLvl(levelQ, ctIn.Value[0], eval.ctxpool.Value[0])
		ring.CopyValuesLvl(levelQ, ctIn.Value[1], eval.ctxpool.Value[1])
		ctInTmp0, ctInTmp1 = eval.ctxpool.Value[0], eval.ctxpool.Value[1]
	} else {
		ctInTmp0, ctInTmp1 = ctIn.Value[0], ctIn.Value[1]
	}

	// Pre-rotates ciphertext for the baby-step giant-step algorithm, does not divide by P yet
	vecRotQ, vecRotP := eval.rotateHoistedNoModDown(levelQ, rotations, eval.PoolDecompQP)

	// Accumulator inner loop
	tmpQ0 := eval.Pool[1].Q
	tmpQ1 := eval.Pool[2].Q
	tmpP0 := eval.Pool[1].P
	tmpP1 := eval.Pool[2].P

	resP0 := eval.Pool[3].Q // Temporary modulus P of c0'
	resP1 := eval.Pool[3].P // Temporary modulus P of c1'

	// Keyswitch accumulator
	c0Q := eval.Pool[4].Q
	c0P := eval.Pool[4].P
	c1Q := eval.Pool[5].Q
	c1P := eval.Pool[5].P

	// Not available (used by switchKeysInPlaceNoModDown)
	// eval.Pool[0].P
	// eval.Pool[0].Q

	ringQ.MulScalarBigintLvl(levelQ, ctInTmp0, ringP.ModulusBigint, ctInTmp0) // P*c0
	ringQ.MulScalarBigintLvl(levelQ, ctInTmp1, ringP.ModulusBigint, ctInTmp1) // P*c1

	for _, i := range rotations {
		if i != 0 {

			galEl := eval.params.GaloisElementForColumnRotationBy(i)

			_, generated := eval.rtks.Keys[galEl]
			if !generated {
				panic("switching key not available")
			}

			ring.PermuteNTTWithIndexAndAddNoModLvl(levelQ, ctInTmp0, eval.permuteNTTIndex[galEl], vecRotQ[i][0]) // phi(d0_Q) += phi(P*c0)
		}
	}

	// OUTER LOOP
	var cnt0 int
	for j := range index {

		if j != 0 {

			// INNER LOOP
			var state bool
			var cnt1 int
			for _, i := range index[j] {

				if i == 0 {
					state = true
				} else {

					plaintextQ, plaintextP := matrix.Vec[N1*j+i].Q, matrix.Vec[N1*j+i].P

					if cnt1 == 0 {
						ringQ.MulCoeffsMontgomeryConstantLvl(levelQ, plaintextQ, vecRotQ[i][0], tmpQ0) // P * phi(c0 + d0_Q) * plaintext
						ringQ.MulCoeffsMontgomeryConstantLvl(levelQ, plaintextQ, vecRotQ[i][1], tmpQ1) // P * phi(d1_Q) * plaintext
						ringP.MulCoeffsMontgomeryConstantLvl(levelP, plaintextP, vecRotP[i][0], tmpP0) // P * phi(d0_P) * plaintext
						ringP.MulCoeffsMontgomeryConstantLvl(levelP, plaintextP, vecRotP[i][1], tmpP1) // P * phi(d1_P) * plaintext
					} else {
						ringQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, plaintextQ, vecRotQ[i][0], tmpQ0) // P * phi(c0 + d0_Q) * plaintext
						ringQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, plaintextQ, vecRotQ[i][1], tmpQ1) // P * phi(d1_Q) * plaintext
						ringP.MulCoeffsMontgomeryAndAddNoMod(plaintextP, vecRotP[i][0], tmpP0)                    // P * phi(d0_P) * plaintext
						ringP.MulCoeffsMontgomeryAndAddNoMod(plaintextP, vecRotP[i][1], tmpP1)                    // P * phi(d1_P) * plaintext
					}

					if cnt1%(QiOverF>>1) == (QiOverF>>1)-1 {
						ringQ.ReduceLvl(levelQ, tmpQ0, tmpQ0)
						ringQ.ReduceLvl(levelQ, tmpQ1, tmpQ1)
					}

					if cnt1%(PiOverF>>1) == (PiOverF>>1)-1 {
						ringP.ReduceLvl(levelP, tmpP0, tmpP0)
						ringP.ReduceLvl(levelP, tmpP1, tmpP1)
					}

					cnt1++
				}
			}

			if cnt1%(QiOverF>>1) != 0 {
				ringQ.ReduceLvl(levelQ, tmpQ0, tmpQ0)
				ringQ.ReduceLvl(levelQ, tmpQ1, tmpQ1)
			}

			if cnt1%(PiOverF>>1) != 0 {
				ringP.ReduceLvl(levelP, tmpP0, tmpP0)
				ringP.ReduceLvl(levelP, tmpP1, tmpP1)
			}

			// If i == 0 (rotation by 0)
			if state {

				// If no loop before, then we copy the values on the accumulator instead of adding them
				if len(index[j]) == 1 {
					ringQ.MulCoeffsMontgomeryLvl(levelQ, matrix.Vec[N1*j].Q, ctInTmp0, tmpQ0) // P * c0 * plaintext + sum(P * phi(d0) * plaintext) + P * phi(c0) * plaintext mod Q
					ringQ.MulCoeffsMontgomeryLvl(levelQ, matrix.Vec[N1*j].Q, ctInTmp1, tmpQ1) // P * c1 * plaintext + sum(P * phi(d1) * plaintext) + P * phi(c1) * plaintext mod Q
				} else {
					ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.Vec[N1*j].Q, ctInTmp0, tmpQ0) // P * c0 * plaintext + sum(P * phi(d0) * plaintext) + P * phi(c0) * plaintext mod Q
					ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.Vec[N1*j].Q, ctInTmp1, tmpQ1) // P * c1 * plaintext + sum(P * phi(d1) * plaintext) + P * phi(c1) * plaintext mod Q
				}
			}

			// Hoisting of the ModDown of sum(sum(phi(d1) * plaintext))
			eval.Baseconverter.ModDownQPtoQNTT(levelQ, levelP, tmpQ1, tmpP1, tmpQ1) // c1 * plaintext + sum(phi(d1) * plaintext) + phi(c1) * plaintext mod Q

			galEl := eval.params.GaloisElementForColumnRotationBy(N1 * j)

			rtk, generated := eval.rtks.Keys[galEl]
			if !generated {
				panic("switching key not available")
			}

			index := eval.permuteNTTIndex[galEl]

			tmpQ1.IsNTT = true
			eval.SwitchKeysInPlaceNoModDown(levelQ, tmpQ1, rtk, c0Q, c0P, c1Q, c1P) // Switchkey(P*phi(tmpRes_1)) = (d0, d1) in base QP

			// Outer loop rotations
			if cnt0 == 0 {
				ring.PermuteNTTWithIndexLvl(levelQ, c0Q, index, ctOut.Value[0]) // sum(phi(d0_Q))
				ring.PermuteNTTWithIndexLvl(levelQ, c1Q, index, ctOut.Value[1]) // sum(phi(d1_Q))
				ring.PermuteNTTWithIndexLvl(levelP, c0P, index, resP0)          // sum(phi(d0_P))
				ring.PermuteNTTWithIndexLvl(levelP, c1P, index, resP1)          // sum(phi(d1_P))

			} else {
				ring.PermuteNTTWithIndexAndAddNoModLvl(levelQ, c0Q, index, ctOut.Value[0]) // sum(phi(d0_Q))
				ring.PermuteNTTWithIndexAndAddNoModLvl(levelQ, c1Q, index, ctOut.Value[1]) // sum(phi(d1_Q))
				ring.PermuteNTTWithIndexAndAddNoModLvl(levelP, c0P, index, resP0)          // sum(phi(d0_P))
				ring.PermuteNTTWithIndexAndAddNoModLvl(levelP, c1P, index, resP1)          // sum(phi(d1_P))
			}

			// sum(sum(phi(d0 + P*c0) * plaintext))
			ring.PermuteNTTWithIndexAndAddNoModLvl(levelQ, tmpQ0, index, ctOut.Value[0]) // sum(phi(d0_Q))
			ring.PermuteNTTWithIndexAndAddNoModLvl(levelP, tmpP0, index, resP0)          // sum(phi(d1_Q))

			if cnt0%QiOverF == QiOverF-1 {
				ringQ.ReduceLvl(levelQ, ctOut.Value[0], ctOut.Value[0])
				ringQ.ReduceLvl(levelQ, ctOut.Value[1], ctOut.Value[1])
			}

			if cnt0%PiOverF == PiOverF-1 {
				ringP.ReduceLvl(levelP, resP0, resP0)
				ringP.ReduceLvl(levelP, resP1, resP1)
			}

			cnt0++
		}
	}

	if cnt0%QiOverF != 0 {
		ringQ.ReduceLvl(levelQ, ctOut.Value[0], ctOut.Value[0])
		ringQ.ReduceLvl(levelQ, ctOut.Value[1], ctOut.Value[1])
	}

	if cnt0%PiOverF != 0 {
		ringP.ReduceLvl(levelP, resP0, resP0)
		ringP.ReduceLvl(levelP, resP1, resP1)
	}

	// if j == 0 (N2 rotation by zero)
	var state bool
	var cnt1 int
	for _, i := range index[0] {

		if i == 0 {
			state = true
		} else {
			plaintextQ, plaintextP := matrix.Vec[i].Q, matrix.Vec[i].P
			// keyswitch(c1_Q) = (d0_QP, d1_QP)
			ringQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, plaintextQ, vecRotQ[i][0], ctOut.Value[0]) // phi(P*c0 + d0_Q) * plaintext
			ringQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, plaintextQ, vecRotQ[i][1], ctOut.Value[1]) // phi(d1_Q) * plaintext
			ringP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelP, plaintextP, vecRotP[i][0], resP0)          // phi(d0_P) * plaintext
			ringP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelP, plaintextP, vecRotP[i][1], resP1)          // phi(d1_P) * plaintext

			if cnt1%(QiOverF>>1) == (QiOverF>>1)-1 {
				ringQ.ReduceLvl(levelQ, ctOut.Value[0], ctOut.Value[0])
				ringQ.ReduceLvl(levelQ, ctOut.Value[1], ctOut.Value[1])
			}

			if cnt1%(PiOverF>>1) == (PiOverF>>1)-1 {
				ringP.ReduceLvl(levelP, resP0, resP0)
				ringP.ReduceLvl(levelP, resP1, resP1)
			}

			cnt1++
		}
	}

	if cnt1%(QiOverF>>1) != 0 {
		ringQ.ReduceLvl(levelQ, ctOut.Value[0], ctOut.Value[0])
		ringQ.ReduceLvl(levelQ, ctOut.Value[1], ctOut.Value[1])
	}

	if cnt1%(PiOverF>>1) != 0 {
		ringP.ReduceLvl(levelP, resP0, resP0)
		ringP.ReduceLvl(levelP, resP1, resP1)
	}

	if state { // Rotation by zero
		ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.Vec[0].Q, ctInTmp0, ctOut.Value[0]) // ctOut += c0_Q * plaintext
		ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.Vec[0].Q, ctInTmp1, ctOut.Value[1]) // ctOut += c1_Q * plaintext
	}

	eval.Baseconverter.ModDownQPtoQNTT(levelQ, levelP, ctOut.Value[0], resP0, ctOut.Value[0]) // sum(phi(c0 * P + d0_QP))/P
	eval.Baseconverter.ModDownQPtoQNTT(levelQ, levelP, ctOut.Value[1], resP1, ctOut.Value[1]) // sum(phi(d1_QP))/P

	ctOut.Scale = matrix.Scale * ctIn.Scale

	vecRotQ, vecRotP = nil, nil
}
