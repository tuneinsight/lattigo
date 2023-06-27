package rgsw

import (
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
)

// Evaluator is a type for evaluating homomorphic operations involving RGSW ciphertexts.
// It currently supports the external product between a RLWE and a RGSW ciphertext (see
// Evaluator.ExternalProduct).
type Evaluator struct {
	rlwe.Evaluator

	params rlwe.Parameters
}

// NewEvaluator creates a new Evaluator type supporting RGSW operations in addition
// to rlwe.Evaluator operations.
func NewEvaluator(params rlwe.Parameters, evk rlwe.EvaluationKeySet) *Evaluator {
	return &Evaluator{*rlwe.NewEvaluator(params, evk), params}
}

// ShallowCopy creates a shallow copy of this Evaluator in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Evaluators can be used concurrently.
func (eval Evaluator) ShallowCopy() *Evaluator {
	return &Evaluator{*eval.Evaluator.ShallowCopy(), eval.params}
}

// WithKey creates a shallow copy of the receiver Evaluator for which the new EvaluationKey is evaluationKey
// and where the temporary buffers are shared. The receiver and the returned Evaluators cannot be used concurrently.
func (eval Evaluator) WithKey(evk rlwe.EvaluationKeySet) *Evaluator {
	return &Evaluator{*eval.Evaluator.WithKey(evk), eval.params}
}

// ExternalProduct computes RLWE x RGSW -> RLWE
//
//	RLWE : (-as + m + e, a)
//	x
//	RGSW : [(-as + P*w*m1 + e, a), (-bs + e, b + P*w*m1)]
//	=
//	RLWE : (<RLWE, RGSW[0]>, <RLWE, RGSW[1]>)
func (eval Evaluator) ExternalProduct(op0 *rlwe.Ciphertext, op1 *Ciphertext, op2 *rlwe.Ciphertext) {

	levelQ, levelP := op1.LevelQ(), op1.LevelP()

	var c0QP, c1QP ringqp.Poly
	if op0 == op2 {
		c0QP, c1QP = eval.BuffQP[1], eval.BuffQP[2]
	} else {
		c0QP, c1QP = ringqp.Poly{Q: op2.Value[0], P: eval.BuffQP[1].P}, ringqp.Poly{Q: op2.Value[1], P: eval.BuffQP[2].P}
	}

	if levelP < 1 {

		// If log(Q) * (Q-1)**2 < 2^{64}-1
		if ringQ := eval.params.RingQ(); levelQ == 0 && levelP == -1 && (ringQ.SubRings[0].Modulus>>29) == 0 {
			eval.externalProduct32Bit(op0, op1, c0QP.Q, c1QP.Q)
			ringQ.AtLevel(0).IMForm(c0QP.Q, op2.Value[0])
			ringQ.AtLevel(0).IMForm(c1QP.Q, op2.Value[1])
		} else {

			eval.externalProductInPlaceSinglePAndBitDecomp(op0, op1, c0QP, c1QP)

			if levelP == 0 {
				eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c0QP.Q, c0QP.P, op2.Value[0])
				eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c1QP.Q, c1QP.P, op2.Value[1])
			} else {
				op2.Value[0].CopyValues(c0QP.Q)
				op2.Value[1].CopyValues(c1QP.Q)
			}
		}
	} else {
		eval.externalProductInPlaceMultipleP(levelQ, levelP, op0, op1, eval.BuffQP[1].Q, eval.BuffQP[1].P, eval.BuffQP[2].Q, eval.BuffQP[2].P)
		eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c0QP.Q, c0QP.P, op2.Value[0])
		eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c1QP.Q, c1QP.P, op2.Value[1])

	}
}

func (eval Evaluator) externalProduct32Bit(ct0 *rlwe.Ciphertext, rgsw *Ciphertext, c0, c1 ring.Poly) {

	// rgsw = [(-as + P*w*m1 + e, a), (-bs + e, b + P*w*m1)]
	// ct = [-cs + m0 + e, c]
	// ctOut = [<ct, rgsw[0]>, <ct, rgsw[1]>] = [ct[0] * rgsw[0][0] + ct[1] * rgsw[0][1], ct[0] * rgsw[1][0] + ct[1] * rgsw[1][1]]
	ringQ := eval.params.RingQ().AtLevel(0)
	subRing := ringQ.SubRings[0]
	pw2 := eval.params.Pow2Base()
	mask := uint64(((1 << pw2) - 1))

	cw := eval.BuffQP[0].Q.Coeffs[0]
	cwNTT := eval.BuffBitDecomp

	acc0 := c0.Coeffs[0]
	acc1 := c1.Coeffs[0]

	// (a, b) + (c0 * rgsw[0][0], c0 * rgsw[0][1])
	// (a, b) + (c1 * rgsw[1][0], c1 * rgsw[1][1])
	for i, el := range rgsw.Value {
		ringQ.INTT(ct0.Value[i], eval.BuffInvNTT)
		for j := range el.Value[0] {
			ring.MaskVec(eval.BuffInvNTT.Coeffs[0], j*pw2, mask, cw)
			if j == 0 && i == 0 {
				subRing.NTTLazy(cw, cwNTT)
				subRing.MulCoeffsLazy(el.Value[0][j][0].Q.Coeffs[0], cwNTT, acc0)
				subRing.MulCoeffsLazy(el.Value[0][j][1].Q.Coeffs[0], cwNTT, acc1)
			} else {
				subRing.NTTLazy(cw, cwNTT)
				subRing.MulCoeffsLazyThenAddLazy(el.Value[0][j][0].Q.Coeffs[0], cwNTT, acc0)
				subRing.MulCoeffsLazyThenAddLazy(el.Value[0][j][1].Q.Coeffs[0], cwNTT, acc1)
			}
		}
	}
}

func (eval Evaluator) externalProductInPlaceSinglePAndBitDecomp(ct0 *rlwe.Ciphertext, rgsw *Ciphertext, c0QP, c1QP ringqp.Poly) {

	// rgsw = [(-as + P*w*m1 + e, a), (-bs + e, b + P*w*m1)]
	// ct = [-cs + m0 + e, c]
	// ctOut = [<ct, rgsw[0]>, <ct, rgsw[1]>] = [ct[0] * rgsw[0][0] + ct[1] * rgsw[0][1], ct[0] * rgsw[1][0] + ct[1] * rgsw[1][1]]
	levelQ := rgsw.LevelQ()
	levelP := rgsw.LevelP()

	ringQP := eval.params.RingQP().AtLevel(levelQ, levelP)

	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	pw2 := eval.params.Pow2Base()
	mask := uint64(((1 << pw2) - 1))
	if mask == 0 {
		mask = 0xFFFFFFFFFFFFFFFF
	}

	decompRNS := eval.params.DecompRNS(levelQ, levelP)
	decompPw2 := eval.params.DecompPw2(levelQ, levelP)

	// (a, b) + (c0 * rgsw[k][0], c0 * rgsw[k][1])
	for k, el := range rgsw.Value {
		ringQ.INTT(ct0.Value[k], eval.BuffInvNTT)
		cw := eval.BuffQP[0].Q.Coeffs[0]
		cwNTT := eval.BuffBitDecomp
		for i := 0; i < decompRNS; i++ {
			for j := 0; j < decompPw2; j++ {
				ring.MaskVec(eval.BuffInvNTT.Coeffs[i], j*pw2, mask, cw)
				if k == 0 && i == 0 && j == 0 {

					for u, s := range ringQ.SubRings[:levelQ+1] {
						s.NTTLazy(cw, cwNTT)
						s.MulCoeffsMontgomery(el.Value[i][j][0].Q.Coeffs[u], cwNTT, c0QP.Q.Coeffs[u])
						s.MulCoeffsMontgomery(el.Value[i][j][1].Q.Coeffs[u], cwNTT, c1QP.Q.Coeffs[u])
					}

					if ringP != nil {
						for u, s := range ringP.SubRings[:levelP+1] {
							s.NTTLazy(cw, cwNTT)
							s.MulCoeffsMontgomery(el.Value[i][j][0].P.Coeffs[u], cwNTT, c0QP.P.Coeffs[u])
							s.MulCoeffsMontgomery(el.Value[i][j][1].P.Coeffs[u], cwNTT, c1QP.P.Coeffs[u])
						}
					}

				} else {

					for u, s := range ringQ.SubRings[:levelQ+1] {
						s.NTTLazy(cw, cwNTT)
						s.MulCoeffsMontgomeryThenAdd(el.Value[i][j][0].Q.Coeffs[u], cwNTT, c0QP.Q.Coeffs[u])
						s.MulCoeffsMontgomeryThenAdd(el.Value[i][j][1].Q.Coeffs[u], cwNTT, c1QP.Q.Coeffs[u])
					}

					if ringP != nil {
						for u, s := range ringP.SubRings[:levelP+1] {
							s.NTTLazy(cw, cwNTT)
							s.MulCoeffsMontgomeryThenAdd(el.Value[i][j][0].P.Coeffs[u], cwNTT, c0QP.P.Coeffs[u])
							s.MulCoeffsMontgomeryThenAdd(el.Value[i][j][1].P.Coeffs[u], cwNTT, c1QP.P.Coeffs[u])
						}
					}
				}
			}
		}
	}
}

func (eval Evaluator) externalProductInPlaceMultipleP(levelQ, levelP int, ct0 *rlwe.Ciphertext, rgsw *Ciphertext, c0OutQ, c0OutP, c1OutQ, c1OutP ring.Poly) {
	var reduce int

	ringQP := eval.params.RingQP().AtLevel(levelQ, levelP)
	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	c2QP := eval.BuffQP[0]

	c0QP := ringqp.Poly{Q: c0OutQ, P: c0OutP}
	c1QP := ringqp.Poly{Q: c1OutQ, P: c1OutP}

	decompRNS := eval.params.DecompRNS(levelQ, levelP)

	QiOverF := eval.params.QiOverflowMargin(levelQ) >> 1
	PiOverF := eval.params.PiOverflowMargin(levelP) >> 1

	var c2NTT, c2InvNTT ring.Poly

	for k, el := range rgsw.Value {

		if ct0.IsNTT {
			c2NTT = ct0.Value[k]
			c2InvNTT = eval.BuffInvNTT
			ringQ.INTT(c2NTT, c2InvNTT)
		} else {
			c2NTT = eval.BuffInvNTT
			c2InvNTT = ct0.Value[k]
			ringQ.NTT(c2InvNTT, c2NTT)
		}

		// (a, b) + (c0 * rgsw[0][0], c0 * rgsw[0][1])
		for i := 0; i < decompRNS; i++ {

			eval.DecomposeSingleNTT(levelQ, levelP, levelP+1, i, c2NTT, c2InvNTT, c2QP.Q, c2QP.P)

			if k == 0 && i == 0 {
				ringQP.MulCoeffsMontgomeryLazy(el.Value[i][0][0], c2QP, c0QP)
				ringQP.MulCoeffsMontgomeryLazy(el.Value[i][0][1], c2QP, c1QP)
			} else {
				ringQP.MulCoeffsMontgomeryLazyThenAddLazy(el.Value[i][0][0], c2QP, c0QP)
				ringQP.MulCoeffsMontgomeryLazyThenAddLazy(el.Value[i][0][1], c2QP, c1QP)
			}

			if reduce%QiOverF == QiOverF-1 {
				ringQ.Reduce(c0QP.Q, c0QP.Q)
				ringQ.Reduce(c1QP.Q, c1QP.Q)
			}

			if reduce%PiOverF == PiOverF-1 {
				ringP.Reduce(c0QP.P, c0QP.P)
				ringP.Reduce(c1QP.P, c1QP.P)
			}

			reduce++
		}
	}

	if reduce%QiOverF != 0 {
		ringQ.Reduce(c0QP.Q, c0QP.Q)
		ringQ.Reduce(c1QP.Q, c1QP.Q)
	}

	if reduce%PiOverF != 0 {
		ringP.Reduce(c0QP.P, c0QP.P)
		ringP.Reduce(c1QP.P, c1QP.P)
	}
}

// AddLazy adds op to ctOut, without modular reduction.
func AddLazy(op interface{}, ringQP ringqp.Ring, ctOut *Ciphertext) {
	switch el := op.(type) {
	case *Plaintext:

		nQ := ringQP.LevelQ() + 1
		nP := ringQP.LevelP() + 1

		if nP == 0 {
			nP = 1
		}

		s := ringQP.RingQ.SubRings[0] // Doesn't matter which one since we add without modular reduction

		for i := range ctOut.Value[0].Value {
			for j := range ctOut.Value[0].Value[i] {
				start, end := i*nP, (i+1)*nP
				if end > nQ {
					end = nQ
				}
				for k := start; k < end; k++ {
					s.AddLazy(ctOut.Value[0].Value[i][j][0].Q.Coeffs[k], el.Value[j].Coeffs[k], ctOut.Value[0].Value[i][j][0].Q.Coeffs[k])
					s.AddLazy(ctOut.Value[1].Value[i][j][1].Q.Coeffs[k], el.Value[j].Coeffs[k], ctOut.Value[1].Value[i][j][1].Q.Coeffs[k])
				}
			}
		}
	case *Ciphertext:
		for i := range el.Value[0].Value {
			for j := range el.Value[0].Value[i] {
				ringQP.AddLazy(ctOut.Value[0].Value[i][j][0], el.Value[0].Value[i][j][0], ctOut.Value[0].Value[i][j][0])
				ringQP.AddLazy(ctOut.Value[0].Value[i][j][1], el.Value[0].Value[i][j][1], ctOut.Value[0].Value[i][j][1])
				ringQP.AddLazy(ctOut.Value[1].Value[i][j][0], el.Value[1].Value[i][j][0], ctOut.Value[1].Value[i][j][0])
				ringQP.AddLazy(ctOut.Value[1].Value[i][j][1], el.Value[1].Value[i][j][1], ctOut.Value[1].Value[i][j][1])
			}
		}
	default:
		panic("cannot AddLazy: unsuported op.(type), must be either *rgsw.Plaintext or *rgsw.Ciphertext")
	}
}

// Reduce applies the modular reduction on ctIn and returns the result on ctOut.
func Reduce(ctIn *Ciphertext, ringQP ringqp.Ring, ctOut *Ciphertext) {
	for i := range ctIn.Value[0].Value {
		for j := range ctIn.Value[0].Value[i] {
			ringQP.Reduce(ctIn.Value[0].Value[i][j][0], ctOut.Value[0].Value[i][j][0])
			ringQP.Reduce(ctIn.Value[0].Value[i][j][1], ctOut.Value[0].Value[i][j][1])
			ringQP.Reduce(ctIn.Value[1].Value[i][j][0], ctOut.Value[1].Value[i][j][0])
			ringQP.Reduce(ctIn.Value[1].Value[i][j][1], ctOut.Value[1].Value[i][j][1])
		}
	}
}

// MulByXPowAlphaMinusOneLazy multiplies ctOut by (X^alpha - 1) and returns the result on ctOut.
func MulByXPowAlphaMinusOneLazy(ctIn *Ciphertext, powXMinusOne ringqp.Poly, ringQP ringqp.Ring, ctOut *Ciphertext) {
	for i := range ctIn.Value[0].Value {
		for j := range ctIn.Value[0].Value[i] {
			ringQP.MulCoeffsMontgomeryLazy(ctIn.Value[0].Value[i][j][0], powXMinusOne, ctOut.Value[0].Value[i][j][0])
			ringQP.MulCoeffsMontgomeryLazy(ctIn.Value[0].Value[i][j][1], powXMinusOne, ctOut.Value[0].Value[i][j][1])
			ringQP.MulCoeffsMontgomeryLazy(ctIn.Value[1].Value[i][j][0], powXMinusOne, ctOut.Value[1].Value[i][j][0])
			ringQP.MulCoeffsMontgomeryLazy(ctIn.Value[1].Value[i][j][1], powXMinusOne, ctOut.Value[1].Value[i][j][1])
		}
	}
}

// MulByXPowAlphaMinusOneThenAddLazy multiplies ctOut by (X^alpha - 1) and adds the result on ctOut.
func MulByXPowAlphaMinusOneThenAddLazy(ctIn *Ciphertext, powXMinusOne ringqp.Poly, ringQP ringqp.Ring, ctOut *Ciphertext) {
	for i := range ctIn.Value[0].Value {
		for j := range ctIn.Value[0].Value[i] {
			ringQP.MulCoeffsMontgomeryLazyThenAddLazy(ctIn.Value[0].Value[i][j][0], powXMinusOne, ctOut.Value[0].Value[i][j][0])
			ringQP.MulCoeffsMontgomeryLazyThenAddLazy(ctIn.Value[0].Value[i][j][1], powXMinusOne, ctOut.Value[0].Value[i][j][1])
			ringQP.MulCoeffsMontgomeryLazyThenAddLazy(ctIn.Value[1].Value[i][j][0], powXMinusOne, ctOut.Value[1].Value[i][j][0])
			ringQP.MulCoeffsMontgomeryLazyThenAddLazy(ctIn.Value[1].Value[i][j][1], powXMinusOne, ctOut.Value[1].Value[i][j][1])
		}
	}
}
