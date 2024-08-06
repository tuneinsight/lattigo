package ckks

import (
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// GetRootsBigComplex returns the roots e^{2*pi*i/m *j} for 0 <= j <= NthRoot
// with prec bits of precision.
func GetRootsBigComplex(NthRoot int, prec uint) (roots []*bignum.Complex) {

	roots = make([]*bignum.Complex, NthRoot+1)

	quarm := NthRoot >> 2

	Pi := bignum.Pi(prec)

	e2ipi := bignum.NewFloat(2, prec)
	e2ipi.Mul(e2ipi, Pi)
	e2ipi.Quo(e2ipi, bignum.NewFloat(float64(NthRoot), prec))

	angle := new(big.Float).SetPrec(prec)

	roots[0] = &bignum.Complex{bignum.NewFloat(1, prec), bignum.NewFloat(0, prec)}

	for i := 1; i < quarm; i++ {
		angle.Mul(e2ipi, bignum.NewFloat(float64(i), prec))
		roots[i] = &bignum.Complex{bignum.Cos(angle), nil}
	}

	for i := 1; i < quarm; i++ {
		roots[quarm-i][1] = new(big.Float).Set(roots[i].Real())
	}

	roots[quarm] = &bignum.Complex{bignum.NewFloat(0, prec), bignum.NewFloat(1, prec)}

	for i := 1; i < quarm+1; i++ {
		roots[i+1*quarm] = &bignum.Complex{new(big.Float).Neg(roots[quarm-i].Real()), new(big.Float).Set(roots[quarm-i].Imag())}
		roots[i+2*quarm] = &bignum.Complex{new(big.Float).Neg(roots[i].Real()), new(big.Float).Neg(roots[i].Imag())}
		roots[i+3*quarm] = &bignum.Complex{new(big.Float).Set(roots[quarm-i].Real()), new(big.Float).Neg(roots[quarm-i].Imag())}
	}

	roots[NthRoot] = roots[0]

	return
}

// GetRootsComplex128 returns the roots e^{2*pi*i/m *j} for 0 <= j <= NthRoot.
func GetRootsComplex128(NthRoot int) (roots []complex128) {
	roots = make([]complex128, NthRoot+1)

	quarm := NthRoot >> 2

	angle := 2 * 3.141592653589793 / float64(NthRoot)

	for i := 0; i < quarm; i++ {
		roots[i] = complex(math.Cos(angle*float64(i)), 0)
	}

	for i := 0; i < quarm; i++ {
		roots[quarm-i] += complex(0, real(roots[i]))
	}

	for i := 1; i < quarm+1; i++ {
		roots[i+1*quarm] = complex(-real(roots[quarm-i]), imag(roots[quarm-i]))
		roots[i+2*quarm] = -roots[i]
		roots[i+3*quarm] = complex(real(roots[quarm-i]), -imag(roots[quarm-i]))
	}

	roots[NthRoot] = roots[0]

	return
}

// StandardDeviation computes the scaled standard deviation of the input vector.
func StandardDeviation(vec interface{}, scale rlwe.Scale) (std float64) {

	switch vec := vec.(type) {
	case []float64:
		// We assume that the error is centered around zero
		var err, tmp, mean, n float64

		n = float64(len(vec))

		for _, c := range vec {
			mean += c
		}

		mean /= n

		for _, c := range vec {
			tmp = c - mean
			err += tmp * tmp
		}

		std = math.Sqrt(err/(n-1)) * scale.Float64()
	case []*big.Float:
		mean := new(big.Float)

		for _, c := range vec {
			mean.Add(mean, c)
		}

		mean.Quo(mean, new(big.Float).SetInt64(int64(len(vec))))

		err := new(big.Float)
		tmp := new(big.Float)
		for _, c := range vec {
			tmp.Sub(c, mean)
			tmp.Mul(tmp, tmp)
			err.Add(err, tmp)
		}

		err.Quo(err, new(big.Float).SetInt64(int64(len(vec)-1)))
		err.Sqrt(err)
		err.Mul(err, &scale.Value)

		std, _ = err.Float64()
	}

	return
}

// Complex128ToFixedPointCRT encodes a vector of complex128 on a CRT polynomial.
// The real part is put in a left N/2 coefficient and the imaginary in the right N/2 coefficients.
func Complex128ToFixedPointCRT(r *ring.Ring, values []complex128, scale float64, coeffs [][]uint64) {

	for i, v := range values {
		SingleFloat64ToFixedPointCRT(r, i, real(v), scale, coeffs)
	}

	var start int
	if r.Type() == ring.Standard {
		slots := len(values)
		for i, v := range values {
			SingleFloat64ToFixedPointCRT(r, i+slots, imag(v), scale, coeffs)
		}

		start = 2 * len(values)

	} else {
		start = len(values)
	}

	end := len(coeffs[0])
	for i := start; i < end; i++ {
		SingleFloat64ToFixedPointCRT(r, i, 0, 0, coeffs)
	}
}

// Float64ToFixedPointCRT encodes a vector of floats on a CRT polynomial.
func Float64ToFixedPointCRT(r *ring.Ring, values []float64, scale float64, coeffs [][]uint64) {

	start := len(values)
	end := len(coeffs[0])

	for i := 0; i < start; i++ {
		SingleFloat64ToFixedPointCRT(r, i, values[i], scale, coeffs)
	}

	for i := start; i < end; i++ {
		SingleFloat64ToFixedPointCRT(r, i, 0, 0, coeffs)
	}
}

// SingleFloat64ToFixedPointCRT encodes a single float64 on a CRT polynomialon in the i-th coefficient.
func SingleFloat64ToFixedPointCRT(r *ring.Ring, i int, value float64, scale float64, coeffs [][]uint64) {

	if value == 0 {
		for j := range coeffs {
			coeffs[j][i] = 0
		}

		return
	}

	var isNegative bool
	var xFlo *big.Float
	var xInt *big.Int
	tmp := new(big.Int)
	var c uint64

	isNegative = false

	if value < 0 {
		isNegative = true
		scale *= -1
	}

	value *= scale

	moduli := r.ModuliChain()[:r.Level()+1]

	if value >= 1.8446744073709552e+19 {
		xFlo = big.NewFloat(value)
		xFlo.Add(xFlo, big.NewFloat(0.5))
		xInt = new(big.Int)
		xFlo.Int(xInt)
		for j := range moduli {
			tmp.Mod(xInt, bignum.NewInt(moduli[j]))
			if isNegative {
				coeffs[j][i] = moduli[j] - tmp.Uint64()
			} else {
				coeffs[j][i] = tmp.Uint64()
			}
		}

	} else {
		brc := r.BRedConstants()

		c = uint64(value + 0.5)
		if isNegative {
			for j, qi := range moduli {
				if c > qi {
					coeffs[j][i] = qi - ring.BRedAdd(c, qi, brc[j])
				} else {
					coeffs[j][i] = qi - c
				}
			}
		} else {
			for j, qi := range moduli {
				if c > 0x1fffffffffffffff {
					coeffs[j][i] = ring.BRedAdd(c, qi, brc[j])
				} else {
					coeffs[j][i] = c
				}
			}
		}
	}
}

func ComplexArbitraryToFixedPointCRT(r *ring.Ring, values []*bignum.Complex, scale *big.Float, coeffs [][]uint64) {

	xFlo := new(big.Float)
	xInt := new(big.Int)
	tmp := new(big.Int)

	zero := new(big.Float)

	half := new(big.Float).SetFloat64(0.5)

	moduli := r.ModuliChain()[:r.Level()+1]

	var negative bool

	for i := range values {

		xFlo.Mul(scale, values[i][0])

		if values[i][0].Cmp(zero) < 0 {
			xFlo.Sub(xFlo, half)
			negative = true
		} else {
			xFlo.Add(xFlo, half)
			negative = false
		}

		xFlo.Int(xInt)

		for j := range moduli {

			Q := bignum.NewInt(moduli[j])

			tmp.Mod(xInt, Q)

			if negative {
				tmp.Add(tmp, Q)
			}

			coeffs[j][i] = tmp.Uint64()
		}
	}

	if r.Type() == ring.Standard {

		slots := len(values)

		for i := range values {

			xFlo.Mul(scale, values[i][1])

			if values[i][1].Cmp(zero) < 0 {
				xFlo.Sub(xFlo, half)
				negative = true
			} else {
				xFlo.Add(xFlo, half)
				negative = false
			}

			xFlo.Int(xInt)

			for j := range moduli {

				Q := bignum.NewInt(moduli[j])

				tmp.Mod(xInt, Q)

				if negative {
					tmp.Add(tmp, Q)
				}
				coeffs[j][i+slots] = tmp.Uint64()
			}
		}
	}
}

func BigFloatToFixedPointCRT(r *ring.Ring, values []*big.Float, scale *big.Float, coeffs [][]uint64) {

	prec := values[0].Prec()

	xFlo := bignum.NewFloat(0, prec)
	xInt := new(big.Int)
	tmp := new(big.Int)

	zero := new(big.Float)

	half := bignum.NewFloat(0.5, prec)

	moduli := r.ModuliChain()[:r.Level()+1]

	for i := range values {

		if values[i] == nil || values[i].Cmp(zero) == 0 {
			for j := range moduli {
				coeffs[j][i] = 0
			}
		} else {

			xFlo.Mul(scale, values[i])

			if values[i].Cmp(zero) < 0 {
				xFlo.Sub(xFlo, half)
			} else {
				xFlo.Add(xFlo, half)
			}

			xFlo.Int(xInt)

			for j := range moduli {

				Q := bignum.NewInt(moduli[j])

				tmp.Mod(xInt, Q)

				if values[i].Cmp(zero) < 0 {
					tmp.Add(tmp, Q)
				}

				coeffs[j][i] = tmp.Uint64()
			}
		}
	}
}
