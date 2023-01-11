package ckks

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ring"
)

const (
	scalingPrecision = uint(128)
)

func valueToBigComplex(value interface{}, prec uint) (cmplx *ring.Complex) {

	cmplx = new(ring.Complex)

	switch value := value.(type) {
	case complex128:

		if v := real(value); v != 0 {
			cmplx[0] = new(big.Float).SetPrec(prec)
			cmplx[0].SetFloat64(v)
		}

		if v := imag(value); v != 0 {
			cmplx[1] = new(big.Float).SetPrec(prec)
			cmplx[1].SetFloat64(v)
		}

	case float64:
		return valueToBigComplex(complex(value, 0), prec)
	case int:
		return valueToBigComplex(new(big.Int).SetInt64(int64(value)), prec)
	case int64:
		return valueToBigComplex(new(big.Int).SetInt64(value), prec)
	case uint64:
		return valueToBigComplex(new(big.Int).SetUint64(value), prec)
	case *big.Float:
		cmplx[0] = new(big.Float).SetPrec(prec)
		cmplx[0].Set(value)
	case *big.Int:
		cmplx[0] = new(big.Float).SetPrec(prec)
		cmplx[0].SetInt(value)
	case *ring.Complex:
		cmplx[0] = new(big.Float).Set(value[0])
		cmplx[1] = new(big.Float).Set(value[1])
	default:
		panic(fmt.Errorf("invalid value.(type): must be int, int64, uint64, float64, complex128, *big.Int, *big.Float or *ring.Complex but is %T", value))
	}

	return
}

func bigComplexToRNSScalar(r *ring.Ring, scale *big.Float, cmplx *ring.Complex) (RNSReal, RNSImag ring.RNSScalar) {

	//fmt.Println(cmplx, scale)

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

// Divides x by n^2, returns a float
func scaleDown(coeff *big.Int, n float64) (x float64) {

	x, _ = new(big.Float).SetInt(coeff).Float64()
	x /= n

	return
}
