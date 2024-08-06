package ckks

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

func bigComplexToRNSScalar(r *ring.Ring, scale *big.Float, cmplx *bignum.Complex) (RNSReal, RNSImag ring.RNSScalar) {

	if scale == nil {
		scale = new(big.Float).SetFloat64(1)
	}

	real := new(big.Int)
	if cmplx[0] != nil {
		r := new(big.Float).Mul(cmplx[0], scale)

		if cmp := cmplx[0].Cmp(new(big.Float)); cmp > 0 {
			r.Add(r, new(big.Float).SetFloat64(0.5))
		} else if cmp < 0 {
			r.Sub(r, new(big.Float).SetFloat64(0.5))
		}

		r.Int(real)
	}

	imag := new(big.Int)
	if cmplx[1] != nil {
		i := new(big.Float).Mul(cmplx[1], scale)

		if cmp := cmplx[1].Cmp(new(big.Float)); cmp > 0 {
			i.Add(i, new(big.Float).SetFloat64(0.5))
		} else if cmp < 0 {
			i.Sub(i, new(big.Float).SetFloat64(0.5))
		}

		i.Int(imag)
	}

	return r.NewRNSScalarFromBigint(real), r.NewRNSScalarFromBigint(imag)
}

// Divides x by n, returns a float.
func scaleDown(coeff *big.Int, n float64) (x float64) {

	x, _ = new(big.Float).SetInt(coeff).Float64()
	x /= n

	return
}
