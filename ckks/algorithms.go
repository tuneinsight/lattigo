package ckks

import (
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// GetGoldschmidtDivisionIterationsNumber returns the minimum number of iterations of the GoldschmidtDivision
// algorithm to get at least log2Targetprecision bits of precision, considering that the a value in the interval
// (0 < minValue < 2, 2).

// GoldschmidtDivisionNew homomorphically computes 1/x.
// input: ct: Enc(x) with values bounded in the interval (0<minvalue<2, 2) and log2Targetprecision the desired number of bits of precision after the decimal.
// output: Enc(1/x - e), where |e| <= (1-x)^2^(iters+1) -> the bit-precision doubles after each iteration.
// The method automatically estimates how many iterations are needed to achieve the desired precision, and returns an error if the input ciphertext
// does not have enough remaining level and if no bootstrapper was given.
// Note that the desired precision will never exceed log2(ct.Scale) - logN + 1.
func (eval *evaluator) GoldschmidtDivisionNew(ct *rlwe.Ciphertext, minValue, log2Targetprecision float64, btp rlwe.Bootstrapper) (ctInv *rlwe.Ciphertext, err error) {

	params := eval.params

	start := math.Log2(1 - minValue)
	var iters int
	for start+log2Targetprecision > 0.5 {
		start *= 2 // Doubles the bit-precision at each iteration
		iters++
	}

	if depth := iters * params.DefaultScaleModuliRatio(); btp == nil && depth > ct.Level() {
		return nil, fmt.Errorf("cannot GoldschmidtDivisionNew: ct.Level()=%d < depth=%d", ct.Level(), depth)
	}

	a := eval.MulNew(ct, -1)
	b := a.CopyNew()
	eval.Add(a, 2, a)
	eval.Add(b, 1, b)

	for i := 1; i < iters; i++ {

		if btp != nil && (b.Level() == btp.MinimumInputLevel() || b.Level() == params.DefaultScaleModuliRatio()-1) {
			if b, err = btp.Bootstrap(b); err != nil {
				return nil, err
			}
		}

		if btp != nil && (a.Level() == btp.MinimumInputLevel() || a.Level() == params.DefaultScaleModuliRatio()-1) {
			if a, err = btp.Bootstrap(a); err != nil {
				return nil, err
			}
		}

		eval.MulRelin(b, b, b)
		if err = eval.Rescale(b, params.DefaultScale(), b); err != nil {
			return nil, err
		}

		if btp != nil && (b.Level() == btp.MinimumInputLevel() || b.Level() == params.DefaultScaleModuliRatio()-1) {
			if b, err = btp.Bootstrap(b); err != nil {
				return nil, err
			}
		}

		tmp := eval.MulRelinNew(a, b)
		if err = eval.Rescale(tmp, params.DefaultScale(), tmp); err != nil {
			return nil, err
		}

		eval.SetScale(a, tmp.Scale)

		eval.Add(a, tmp, a)
	}

	return a, nil
}
