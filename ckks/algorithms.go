package ckks

import (
	"math/bits"
)

// PowerOf2 compute ct0^(2^logPow2), consuming logPow2 levels, and returns the result on ct1. Providing an evaluation
// key is necessary when logPow2 > 1.
func (evaluator *Evaluator) PowerOf2(el0 *Ciphertext, logPow2 uint64, evakey *EvaluationKey, elOut *Ciphertext) {

	if logPow2 == 0 {

		if el0 != elOut {

			elOut.Copy(el0.Element())
		}

	} else {

		evaluator.MulRelin(el0.Element(), el0.Element(), evakey, elOut)

		evaluator.Rescale(elOut, evaluator.ckksContext.scale, elOut)

		for i := uint64(1); i < logPow2; i++ {

			evaluator.MulRelin(elOut.Element(), elOut.Element(), evakey, elOut)

			evaluator.Rescale(elOut, evaluator.ckksContext.scale, elOut)
		}
	}
}

// PowerNew compute ct0^degree, consuming log(degree) levels, and returns the result on a new element. Providing an evaluation
// key is necessary when degree > 2.
func (evaluator *Evaluator) PowerNew(op *Ciphertext, degree uint64, evakey *EvaluationKey) (opOut *Ciphertext) {
	opOut = NewCiphertext(evaluator.params, 1, op.Level(), op.Scale())
	evaluator.Power(op, degree, evakey, opOut)
	return
}

// Power compute ct0^degree, consuming log(degree) levels, and returns the result on res. Providing an evaluation
// key is necessary when degree > 2.
func (evaluator *Evaluator) Power(ct0 *Ciphertext, degree uint64, evakey *EvaluationKey, res *Ciphertext) {

	tmpct0 := ct0.CopyNew()

	var logDegree, po2Degree uint64

	logDegree = uint64(bits.Len64(degree)) - 1
	po2Degree = 1 << logDegree

	evaluator.PowerOf2(tmpct0.Ciphertext(), logDegree, evakey, res)

	degree -= po2Degree

	for degree > 0 {

		logDegree = uint64(bits.Len64(degree)) - 1
		po2Degree = 1 << logDegree

		tmp := NewCiphertext(evaluator.params, 1, tmpct0.Level(), tmpct0.Scale())

		evaluator.PowerOf2(tmpct0.Ciphertext(), logDegree, evakey, tmp)

		evaluator.MulRelin(res.Element(), tmp.Element(), evakey, res)

		evaluator.Rescale(res, evaluator.ckksContext.scale, res)

		degree -= po2Degree
	}
}

// InverseNew computes 1/ct0 and returns the result on a new element, iterating for n steps and consuming n levels. The algorithm requires the encrypted values to be in the range
// [-1.5 - 1.5i, 1.5 + 1.5i]  or the result will be  wrong. Each iteration increases the precision.
func (evaluator *Evaluator) InverseNew(ct0 *Ciphertext, steps uint64, evakey *EvaluationKey) (res *Ciphertext) {

	cbar := evaluator.NegNew(ct0)

	evaluator.AddConst(cbar, 1, cbar)

	tmp := evaluator.AddConstNew(cbar, 1)
	res = tmp.CopyNew().Ciphertext()

	for i := uint64(1); i < steps; i++ {

		evaluator.MulRelin(cbar.Element(), cbar.Element(), evakey, cbar.Ciphertext())

		evaluator.Rescale(cbar, evaluator.ckksContext.scale, cbar)

		tmp = evaluator.AddConstNew(cbar, 1)

		evaluator.MulRelin(tmp.Element(), res.Element(), evakey, tmp.Ciphertext())

		evaluator.Rescale(tmp, evaluator.ckksContext.scale, tmp)

		res = tmp.CopyNew().Ciphertext()
	}

	return res
}
