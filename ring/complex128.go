package ring

import (
	"math"
	"math/big"
)

// NewFloat creates a new big.Float element with "logPrecision" bits of precision
func NewFloat(x float64, logPrecision int) (y *big.Float) {
	y = new(big.Float)
	y.SetPrec(uint(logPrecision)) // decimal precision
	y.SetFloat64(x)
	return
}

// Cos implements the arbitrary precision computation of Cos(x)
// Iterative process with an error of ~10^{−0.60206*k} after k iterations.
// ref: Johansson, B. Tomas, An elementary algorithm to evaluate trigonometric functions to high precision, 2018
func Cos(x *big.Float) (cosx *big.Float) {
	tmp := new(big.Float)

	prec := int(x.Prec())

	k := int(math.Ceil(float64(x.Prec()) / (3.3219280948873626 * 0.60206))) // number of iterations : ceil( prec(log2) / (log10* 0.60206))

	t := NewFloat(0.5, prec)
	half := new(big.Float).Copy(t)

	for i := 1; i < k-1; i++ {
		t.Mul(t, half)
	}

	s := new(big.Float).Mul(x, t)
	s.Mul(s, x)
	s.Mul(s, t)

	four := NewFloat(4.0, prec)

	for i := 1; i < k; i++ {
		tmp.Sub(four, s)
		s.Mul(s, tmp)
	}

	cosx = new(big.Float).Quo(s, NewFloat(2.0, prec))
	cosx.Sub(NewFloat(1.0, prec), cosx)
	return

}

// Complex is a type for arbitrary precision complex number
type Complex [2]*big.Float

// NewComplex creates a new arbitrary precision complex number
func NewComplex(a, b *big.Float) (c *Complex) {
	c = new(Complex)
	for i := 0; i < 2; i++ {
		c[i] = new(big.Float)
	}

	if a != nil {
		c[0].Set(a)
	}

	if b != nil {
		c[1].Set(b)
	}

	return
}

// Set sets a arbitrary precision complex number
func (c *Complex) Set(a *Complex) {
	c[0].Set(a[0])
	c[1].Set(a[1])
}

// Copy returns a new copy of the target arbitrary precision complex number
func (c *Complex) Copy() *Complex {
	return NewComplex(c[0], c[1])
}

// Real returns the real part as a big.Float
func (c *Complex) Real() *big.Float {
	return c[0]
}

// Imag returns the imaginary part as a big.Float
func (c *Complex) Imag() *big.Float {
	return c[1]
}

// Float64 returns the arbitrary precision complex number as a complex128
func (c *Complex) Float64() complex128 {
	a, _ := c[0].Float64()
	b, _ := c[1].Float64()

	return complex(a, b)
}

// Add adds two arbitrary precision complex numbers together
func (c *Complex) Add(a, b *Complex) {
	c[0].Add(a[0], b[0])
	c[1].Add(a[1], b[1])
}

// Sub subtracts two arbitrary precision complex numbers together
func (c *Complex) Sub(a, b *Complex) {
	c[0].Sub(a[0], b[0])
	c[1].Sub(a[1], b[1])
}

// ComplexMultiplier is a struct for the multiplication or division of two arbitrary precision complex numbers
type ComplexMultiplier struct {
	tmp0 *big.Float
	tmp1 *big.Float
	tmp2 *big.Float
	tmp3 *big.Float
}

// NewComplexMultiplier creates a new ComplexMultiplier
func NewComplexMultiplier() (cEval *ComplexMultiplier) {
	cEval = new(ComplexMultiplier)
	cEval.tmp0 = new(big.Float)
	cEval.tmp1 = new(big.Float)
	cEval.tmp2 = new(big.Float)
	cEval.tmp3 = new(big.Float)
	return
}

// Mul multiplies two arbitrary precision complex numbers together
func (cEval *ComplexMultiplier) Mul(a, b, c *Complex) {

	cEval.tmp0.Mul(a[0], b[0])
	cEval.tmp1.Mul(a[1], b[1])
	cEval.tmp2.Mul(a[0], b[1])
	cEval.tmp3.Mul(a[1], b[0])

	c[0].Sub(cEval.tmp0, cEval.tmp1)
	c[1].Add(cEval.tmp2, cEval.tmp3)
}

// Div divides two arbitrary precision complex numbers together
func (cEval *ComplexMultiplier) Div(a, b, c *Complex) {

	// tmp0 = (a[0] * b[0]) + (a[1] * b[1]) real part
	// tmp1 = (a[1] * b[0]) - (a[0] * b[0]) imag part
	// tmp2 = (b[0] * b[0]) + (b[1] * b[1]) denominator

	cEval.tmp0.Mul(a[0], b[0])
	cEval.tmp1.Mul(a[1], b[1])
	cEval.tmp2.Mul(a[1], b[0])
	cEval.tmp3.Mul(a[0], b[1])

	cEval.tmp0.Add(cEval.tmp0, cEval.tmp1)
	cEval.tmp1.Sub(cEval.tmp2, cEval.tmp3)

	cEval.tmp2.Mul(b[0], b[0])
	cEval.tmp3.Mul(b[1], b[1])
	cEval.tmp2.Add(cEval.tmp2, cEval.tmp3)

	c[0].Quo(cEval.tmp0, cEval.tmp2)
	c[1].Quo(cEval.tmp1, cEval.tmp2)
}
