package ckks

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

func bigComplexToRNSScalar(r *ring.Ring, scale *big.Float, cmplx *bignum.Complex) (RNSReal, RNSImag ring.RNSScalar) {

	if scale == nil {
		scale = new(big.Float).SetFloat64(1)
	}

	real := new(big.Int)
	if cmplx[0] != nil {
		new(big.Float).Mul(cmplx[0], scale).Int(real)
	}

	imag := new(big.Int)
	if cmplx[1] != nil {
		new(big.Float).Mul(cmplx[1], scale).Int(imag)
	}

	return r.NewRNSScalarFromBigint(real), r.NewRNSScalarFromBigint(imag)
}

// Divides x by n, returns a float.
func scaleDown(coeff *big.Int, n float64) (x float64) {

	x, _ = new(big.Float).SetInt(coeff).Float64()
	x /= n

	return
}
