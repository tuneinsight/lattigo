package ckks

import (
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// GetRootsbigFloat returns the roots e^{2*pi*i/m *j} for 0 <= j <= NthRoot
// with prec bits of precision.
func GetRootsbigFloat(NthRoot int, prec uint) (roots []*ring.Complex) {

	roots = make([]*ring.Complex, NthRoot+1)

	quarm := NthRoot >> 2

	var PI = new(big.Float)
	PI.SetPrec(prec)
	PI.SetString(pi)

	e2ipi := ring.NewFloat(2, prec)
	e2ipi.Mul(e2ipi, PI)
	e2ipi.Quo(e2ipi, ring.NewFloat(float64(NthRoot), prec))

	angle := new(big.Float).SetPrec(prec)

	roots[0] = &ring.Complex{ring.NewFloat(1, prec), ring.NewFloat(0, prec)}

	for i := 1; i < quarm; i++ {
		angle.Mul(e2ipi, ring.NewFloat(float64(i), prec))
		roots[i] = &ring.Complex{ring.Cos(angle), nil}
	}

	for i := 1; i < quarm; i++ {
		roots[quarm-i][1] = new(big.Float).Set(roots[i].Real())
	}

	roots[quarm] = &ring.Complex{ring.NewFloat(0, prec), ring.NewFloat(1, prec)}

	for i := 1; i < quarm+1; i++ {
		roots[i+1*quarm] = &ring.Complex{new(big.Float).Neg(roots[quarm-i].Real()), new(big.Float).Set(roots[quarm-i].Imag())}
		roots[i+2*quarm] = &ring.Complex{new(big.Float).Neg(roots[i].Real()), new(big.Float).Neg(roots[i].Imag())}
		roots[i+3*quarm] = &ring.Complex{new(big.Float).Set(roots[quarm-i].Real()), new(big.Float).Neg(roots[quarm-i].Imag())}
	}

	roots[NthRoot] = roots[0]

	return
}

// GetRootsFloat64 returns the roots e^{2*pi*i/m *j} for 0 <= j <= NthRoot.
func GetRootsFloat64(NthRoot int) (roots []complex128) {
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
func StandardDeviation(vec []float64, scale float64) (std float64) {
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

	return math.Sqrt(err/n) * scale
}

func complexToFixedPointCRT(level int, values []complex128, scale float64, ringQ *ring.Ring, coeffs [][]uint64, isRingStandard bool) {

	for i, v := range values {
		singleFloatToFixedPointCRT(level, i, real(v), scale, ringQ, coeffs)
	}

	if isRingStandard {
		slots := len(values)
		for i, v := range values {
			singleFloatToFixedPointCRT(level, i+slots, imag(v), scale, ringQ, coeffs)
		}
	}
}

func floatToFixedPointCRT(level int, values []float64, scale float64, ringQ *ring.Ring, coeffs [][]uint64) {
	for i, v := range values {
		singleFloatToFixedPointCRT(level, i, v, scale, ringQ, coeffs)
	}
}

func singleFloatToFixedPointCRT(level, i int, value float64, scale float64, ringQ *ring.Ring, coeffs [][]uint64) {

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

	moduli := ringQ.ModuliChain()

	if value > 1.8446744073709552e+19 {
		xFlo = big.NewFloat(value)
		xFlo.Add(xFlo, big.NewFloat(0.5))
		xInt = new(big.Int)
		xFlo.Int(xInt)
		for j := range moduli[:level+1] {
			tmp.Mod(xInt, ring.NewUint(moduli[j]))
			if isNegative {
				coeffs[j][i] = moduli[j] - tmp.Uint64()
			} else {
				coeffs[j][i] = tmp.Uint64()
			}
		}

	} else {
		brc := ringQ.BRedConstants()

		c = uint64(value + 0.5)
		if isNegative {
			for j, qi := range moduli[:level+1] {
				if c > qi {
					coeffs[j][i] = qi - ring.BRedAdd(c, qi, brc[j])
				} else {
					coeffs[j][i] = qi - c
				}
			}
		} else {
			for j, qi := range moduli[:level+1] {
				if c > 0x1fffffffffffffff {
					coeffs[j][i] = ring.BRedAdd(c, qi, brc[j])
				} else {
					coeffs[j][i] = c
				}
			}
		}
	}
}

func scaleUpVecExactBigFloat(values []*big.Float, scale float64, moduli []uint64, coeffs [][]uint64) {

	prec := values[0].Prec()

	xFlo := ring.NewFloat(0, prec)
	xInt := new(big.Int)
	tmp := new(big.Int)

	zero := ring.NewFloat(0, prec)

	scaleFlo := ring.NewFloat(scale, prec)
	half := ring.NewFloat(0.5, prec)

	for i := range values {

		xFlo.Mul(scaleFlo, values[i])

		if values[i].Cmp(zero) < 0 {
			xFlo.Sub(xFlo, half)
		} else {
			xFlo.Add(xFlo, half)
		}

		xFlo.Int(xInt)

		for j := range moduli {

			Q := ring.NewUint(moduli[j])

			tmp.Mod(xInt, Q)

			if values[i].Cmp(zero) < 0 {
				tmp.Add(tmp, Q)
			}

			coeffs[j][i] = tmp.Uint64()
		}
	}
}

// SliceBitReverseInPlaceComplex128 applies an in-place bit-reverse permuation on the input slice.
func SliceBitReverseInPlaceComplex128(slice []complex128, N int) {

	var bit, j int

	for i := 1; i < N; i++ {

		bit = N >> 1

		for j >= bit {
			j -= bit
			bit >>= 1
		}

		j += bit

		if i < j {
			slice[i], slice[j] = slice[j], slice[i]
		}
	}
}

// SliceBitReverseInPlaceFloat64 applies an in-place bit-reverse permuation on the input slice.
func SliceBitReverseInPlaceFloat64(slice []float64, N int) {

	var bit, j int

	for i := 1; i < N; i++ {

		bit = N >> 1

		for j >= bit {
			j -= bit
			bit >>= 1
		}

		j += bit

		if i < j {
			slice[i], slice[j] = slice[j], slice[i]
		}
	}
}

// SliceBitReverseInPlaceRingComplex applies an in-place bit-reverse permuation on the input slice.
func SliceBitReverseInPlaceRingComplex(slice []*ring.Complex, N int) {

	var bit, j int

	for i := 1; i < N; i++ {

		bit = N >> 1

		for j >= bit {
			j -= bit
			bit >>= 1
		}

		j += bit

		if i < j {
			slice[i], slice[j] = slice[j], slice[i]
		}
	}
}

// GenSwitchkeysRescalingParams generates the parameters for rescaling the switching keys
func GenSwitchkeysRescalingParams(Q, P []uint64) (params []uint64) {

	params = make([]uint64, len(Q))

	PBig := ring.NewUint(1)
	for _, pj := range P {
		PBig.Mul(PBig, ring.NewUint(pj))
	}

	tmp := ring.NewUint(0)

	for i := 0; i < len(Q); i++ {

		params[i] = tmp.Mod(PBig, ring.NewUint(Q[i])).Uint64()
		params[i] = ring.ModExp(params[i], Q[i]-2, Q[i])
		params[i] = ring.MForm(params[i], Q[i], ring.BRedConstant(Q[i]))
	}

	return
}
