package ring

import (
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
	sinx.Sub(sinx, tmp)
	sinx.Sqrt(sinx)
	return
}

type Complex [2]*big.Float

func NewComplex(a, b *big.Float) (c *Complex) {
	c = new(Complex)
	for i := 0; i < 2; i++ {
		c[i] = new(big.Float)
		c[i].SetPrec(1000)
	}

	c[0].Set(a)
	c[1].Set(b)

	return
}

func (c *Complex) Float64() complex128 {
	a, _ := c[0].Float64()
	b, _ := c[1].Float64()

	return complex(a, b)
}

func (c *Complex) Add(a, b *Complex) {
	c[0].Add(a[0], b[0])
	c[1].Add(a[1], b[1])
}

func (c *Complex) Sub(a, b *Complex) {
	c[0].Sub(a[0], b[0])
	c[1].Sub(a[1], b[1])
}

func (c *Complex) Mul(a, b *Complex) {

	tmp0 := new(big.Float)
	tmp1 := new(big.Float)
	tmp2 := new(big.Float)
	tmp3 := new(big.Float)

	tmp0.Mul(a[0], b[0])
	tmp1.Mul(a[1], b[1])
	tmp2.Mul(a[0], b[1])
	tmp3.Mul(a[1], b[0])

	c[0].Set(tmp0)
	c[0].Sub(c[0], tmp1)

	c[1].Set(tmp2)
	c[1].Add(c[1], tmp3)
}
