package ckks

import (
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// GoldschmidtDivisionNew homomorphically computes 1/x.
// input: ct: Enc(x) with values in the interval [0+minvalue, 2-minvalue] and logPrec the desired number of bits of precision.
// output: Enc(1/x - e), where |e| <= (1-x)^2^(#iterations+1) -> the bit-precision doubles after each iteration.
// The method automatically estimates how many iterations are needed to achieve the desired precision, and returns an error if the input ciphertext
// does not have enough remaining level and if no bootstrapper was given.
// This method will return an error if something goes wrong with the bootstrapping or the rescaling operations.
func (eval *Evaluator) GoldschmidtDivisionNew(ct *rlwe.Ciphertext, minValue, logPrec float64, btp rlwe.Bootstrapper) (ctInv *rlwe.Ciphertext, err error) {

	params := eval.GetParameters()

	start := math.Log2(1 - minValue)
	var iters int
	for start+logPrec > 0.5 {
		start *= 2 // Doubles the bit-precision at each iteration
		iters++
	}

	levelsPerRescaling := params.LevelsConsummedPerRescaling()

	if depth := iters * levelsPerRescaling; btp == nil && depth > ct.Level() {
		return nil, fmt.Errorf("cannot GoldschmidtDivisionNew: ct.Level()=%d < depth=%d and rlwe.Bootstrapper is nil", ct.Level(), depth)
	}

	a, err := eval.MulNew(ct, -1)
	if err != nil {
		return nil, err
	}

	b := a.CopyNew()

	if err = eval.Add(a, 2, a); err != nil {
		return nil, err
	}

	if err = eval.Add(b, 1, b); err != nil {
		return nil, err
	}

	for i := 1; i < iters; i++ {

		if btp != nil && (b.Level() == btp.MinimumInputLevel() || b.Level() == levelsPerRescaling-1) {
			if b, err = btp.Bootstrap(b); err != nil {
				return nil, err
			}
		}

		if btp != nil && (a.Level() == btp.MinimumInputLevel() || a.Level() == levelsPerRescaling-1) {
			if a, err = btp.Bootstrap(a); err != nil {
				return nil, err
			}
		}

		if err = eval.MulRelin(b, b, b); err != nil {
			return nil, err
		}

		if err = eval.Rescale(b, params.DefaultScale(), b); err != nil {
			return nil, err
		}

		if btp != nil && (b.Level() == btp.MinimumInputLevel() || b.Level() == levelsPerRescaling-1) {
			if b, err = btp.Bootstrap(b); err != nil {
				return nil, err
			}
		}

		tmp, err := eval.MulRelinNew(a, b)

		if err != nil {
			return nil, err
		}

		if err = eval.Rescale(tmp, params.DefaultScale(), tmp); err != nil {
			return nil, err
		}

		if err = eval.SetScale(a, tmp.Scale); err != nil {
			return nil, err
		}

		if err = eval.Add(a, tmp, a); err != nil {
			return nil, err
		}
	}

	return a, nil
}
