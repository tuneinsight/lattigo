package bettersine

import (
	//"fmt"
	"math"
	"math/big"
)

func NewFloat(x float64) (y *big.Float) {
	y = new(big.Float)
	y.SetPrec(1000) // decimal precision
	y.SetFloat64(x)
	return
}

// Arbitrary precision computation of Cos(x)
// Iterative process with an error of ~10^{−0.60206*k} after k iterations.
// ref : Johansson, B. Tomas, An elementary algorithm to evaluate trigonometric functions to high precision, 2018
func Cos(x *big.Float) (cosx *big.Float) {
	tmp := new(big.Float)

	k := 1000 // number of iterations
	t := NewFloat(0.5)
	half := new(big.Float).Copy(t)

	for i := 1; i < k-1; i++ {
		t.Mul(t, half)
	}

	s := new(big.Float).Mul(x, t)
	s.Mul(s, x)
	s.Mul(s, t)

	four := NewFloat(4.0)

	for i := 1; i < k; i++ {
		tmp.Sub(four, s)
		s.Mul(s, tmp)
	}

	cosx = new(big.Float).Quo(s, NewFloat(2.0))
	cosx.Sub(NewFloat(1.0), cosx)
	return

}

func Sin(x *big.Float) (sinx *big.Float) {

	sinx = NewFloat(1)
	tmp := Cos(x)
	tmp.Mul(tmp, tmp)
	sinx.Add(sinx, tmp)
	sinx.Sqrt(sinx)
	return
}

func log2(x float64) float64 {
	return math.Log2(x)
}

func abs(x float64) float64 {
	return math.Abs(x)
}
