package ckks

import (
	"math/bits"
)

// PowerOf2 compute ct0^(2^logPow2), consuming logPow2 levels, and returns the result on ct1. Providing an evaluation
// key is necessary when logPow2 > 1.
func (evaluator *Evaluator) PowerOf2(ct0 CkksElement, logPow2 uint64, evakey *EvaluationKey, ct1 CkksElement) error {

	if logPow2 == 0 {

		if ct0 != ct1 {

			if err := ct0.Copy(ct1); err != nil {
				return err
			}
		}

	} else {

		if err := evaluator.MulRelin(ct0, ct0, evakey, ct1); err != nil {
			return err
		}

		evaluator.Rescale(ct1, ct1)

		for i := uint64(1); i < logPow2; i++ {
			if err := evaluator.MulRelin(ct1, ct1, evakey, ct1); err != nil {
				return err
			}
			evaluator.Rescale(ct1, ct1)
		}
	}

	return nil
}

// Power compute ct0^degree, consuming log(degree) levels, and returns the result on a new element. Providing an evaluation
// key is necessary when degree > 2.
func (evaluator *Evaluator) PowerNew(ct0 *Ciphertext, degree uint64, evakey *EvaluationKey) (res *Ciphertext) {
	res = evaluator.ckkscontext.NewCiphertext(1, ct0.Level(), ct0.Scale())
	evaluator.Power(ct0, degree, evakey, res)
	return
}

// Power compute ct0^degree, consuming log(degree) levels, and returns the result on res. Providing an evaluation
// key is necessary when degree > 2.
func (evaluator *Evaluator) Power(ct0 CkksElement, degree uint64, evakey *EvaluationKey, res CkksElement) error {

	tmpct0 := ct0.CopyNew()

	var logDegree, po2Degree uint64

	logDegree = uint64(bits.Len64(degree)) - 1
	po2Degree = 1 << logDegree

	if err := evaluator.PowerOf2(tmpct0, logDegree, evakey, res); err != nil {
		return err
	}

	degree -= po2Degree

	for degree > 0 {

		logDegree = uint64(bits.Len64(degree)) - 1
		po2Degree = 1 << logDegree

		tmp := evaluator.ckkscontext.NewCiphertext(1, tmpct0.Level(), tmpct0.Scale())

		if err := evaluator.PowerOf2(tmpct0, logDegree, evakey, tmp); err != nil {
			return err
		}

		if err := evaluator.MulRelin(res, tmp, evakey, res); err != nil {
			return err
		}

		evaluator.Rescale(res, res)

		degree -= po2Degree
	}

	return nil
}

// InverseNew computes 1/ct0 and returns the result on a new element, iterating for n steps and consuming n levels. The algorithm requires the encrypted values to be in the range
// [-1.5 - 1.5i, 1.5 + 1.5i]  or the result will be  wrong. Each iteration increases the precision.
func (evaluator *Evaluator) InverseNew(ct0 *Ciphertext, steps uint64, evakey *EvaluationKey) (res *Ciphertext, err error) {

	cbar := evaluator.NegNew(ct0)

	evaluator.AddConst(cbar, 1, cbar)

	tmp := evaluator.AddConstNew(cbar, 1)
	res = tmp.CopyNew().(*Ciphertext)

	for i := uint64(1); i < steps; i++ {

		evaluator.MulRelin(cbar, cbar, evakey, cbar)

		if err = evaluator.Rescale(cbar, cbar); err != nil {
			return nil, err
		}

		tmp = evaluator.AddConstNew(cbar, 1)

		if err := evaluator.MulRelin(tmp, res, evakey, tmp); err != nil {
			return nil, err
		}

		evaluator.Rescale(tmp, tmp)
		res = tmp.CopyNew().(*Ciphertext)
	}

	return res, nil
}
