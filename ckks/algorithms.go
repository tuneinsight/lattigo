package ckks

import (
	"math/bits"
)

// PowerOf2 computes op^(2^logPow2), consuming logPow2 levels, and returns the result on opOut. Providing an evaluation
// key is necessary when logPow2 > 1.
func (evaluator *Evaluator) PowerOf2(op *Ciphertext, logPow2 uint64, evakey *EvaluationKey, opOut *Ciphertext) {

	if logPow2 == 0 {

		if op != opOut {

			opOut.Copy(op.Element())
		}

	} else {

		evaluator.MulRelin(op.Element(), op.Element(), evakey, opOut)

		evaluator.Rescale(opOut, evaluator.ckksContext.scale, opOut)

		for i := uint64(1); i < logPow2; i++ {

			evaluator.MulRelin(opOut.Element(), opOut.Element(), evakey, opOut)

			evaluator.Rescale(opOut, evaluator.ckksContext.scale, opOut)
		}
	}
}

// PowerNew computes op^degree, consuming log(degree) levels, and returns the result on a new element. Providing an evaluation
// key is necessary when degree > 2.
func (evaluator *Evaluator) PowerNew(op *Ciphertext, degree uint64, evakey *EvaluationKey) (opOut *Ciphertext) {
	opOut = NewCiphertext(evaluator.params, 1, op.Level(), op.Scale())
	evaluator.Power(op, degree, evakey, opOut)
	return
}

// Power computes op^degree, consuming log(degree) levels, and returns the result on opOut. Providing an evaluation
// key is necessary when degree > 2.
func (evaluator *Evaluator) Power(op *Ciphertext, degree uint64, evakey *EvaluationKey, opOut *Ciphertext) {

	tmpct0 := op.CopyNew()

	var logDegree, po2Degree uint64

	logDegree = uint64(bits.Len64(degree)) - 1
	po2Degree = 1 << logDegree

	evaluator.PowerOf2(tmpct0.Ciphertext(), logDegree, evakey, opOut)

	degree -= po2Degree

	for degree > 0 {

		logDegree = uint64(bits.Len64(degree)) - 1
		po2Degree = 1 << logDegree

		tmp := NewCiphertext(evaluator.params, 1, tmpct0.Level(), tmpct0.Scale())

		evaluator.PowerOf2(tmpct0.Ciphertext(), logDegree, evakey, tmp)

		evaluator.MulRelin(opOut.Element(), tmp.Element(), evakey, opOut)

		evaluator.Rescale(opOut, evaluator.ckksContext.scale, opOut)

		degree -= po2Degree
	}
}

// InverseNew computes 1/op and returns the result on a new element, iterating for n steps and consuming n levels. The algorithm requires the encrypted values to be in the range
// [-1.5 - 1.5i, 1.5 + 1.5i] or the result will be wrong. Each iteration increases the precision.
func (evaluator *Evaluator) InverseNew(op *Ciphertext, steps uint64, evakey *EvaluationKey) (opOut *Ciphertext) {

	cbar := evaluator.NegNew(op)

	evaluator.AddConst(cbar, 1, cbar)

	tmp := evaluator.AddConstNew(cbar, 1)
	opOut = tmp.CopyNew().Ciphertext()

	for i := uint64(1); i < steps; i++ {

		evaluator.MulRelin(cbar.Element(), cbar.Element(), evakey, cbar.Ciphertext())

		evaluator.Rescale(cbar, evaluator.ckksContext.scale, cbar)

		tmp = evaluator.AddConstNew(cbar, 1)

		evaluator.MulRelin(tmp.Element(), opOut.Element(), evakey, tmp.Ciphertext())

		evaluator.Rescale(tmp, evaluator.ckksContext.scale, tmp)

		opOut = tmp.CopyNew().Ciphertext()
	}

	return opOut
}
