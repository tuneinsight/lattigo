package ckks

import (
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// InverseNew computes 1/op and returns the result on a new element, iterating for n steps and consuming n levels. The algorithm requires the encrypted values to be in the range
// [-1.5 - 1.5i, 1.5 + 1.5i] or the result will be wrong. Each iteration increases the precision.
func (eval *evaluator) InverseNew(op *rlwe.Ciphertext, steps int) (opOut *rlwe.Ciphertext, err error) {

	cbar := eval.NegNew(op)

	eval.AddConst(cbar, 1, cbar)

	tmp := eval.AddConstNew(cbar, 1)
	opOut = tmp.CopyNew()

	for i := 1; i < steps; i++ {

		eval.MulRelin(cbar, cbar, cbar)

		if err = eval.Rescale(cbar, op.Scale, cbar); err != nil {
			return
		}

		tmp = eval.AddConstNew(cbar, 1)

		eval.MulRelin(tmp, opOut, tmp)

		if err = eval.Rescale(tmp, op.Scale, tmp); err != nil {
			return
		}

		opOut = tmp.CopyNew()
	}

	return
}
