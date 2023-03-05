package bignum

import (
	"math"
	"math/big"

	"github.com/ALTree/bigfloat"
)

const pi = "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491298336733624406566430860213949463952247371907021798609437027705392171762931767523846748184676694051320005681271452635608277857713427577896091736371787214684409012249534301465495853710507922796892589235420199561121290219608640344181598136297747713099605187072113499999983729780499510597317328160963185950244594553469083026425223082533446850352619311881710100031378387528865875332083814206171776691473035982534904287554687311595628638823537875937519577818577805321712268066130019278766111959092164201989"

// Pi returns Pi with prec bits of precision.
func Pi(prec uint) *big.Float {
	pi, _ := new(big.Float).SetPrec(prec).SetString(pi)
	return pi
}

// NewFloat creates a new big.Float element with "prec" bits of precision
func NewFloat(x interface{}, prec uint) (y *big.Float) {

	y = new(big.Float)
	y.SetPrec(prec) // decimal precision

	if x == nil {
		return
	}

	switch x := x.(type) {
	case int:
		y.SetInt64(int64(x))
	case int64:
		y.SetInt64(x)
	case uint:
		y.SetUint64(uint64(x))
	case uint64:
		y.SetUint64(x)
	case float64:
		y.SetFloat64(x)
	case *big.Int:
		y.SetInt(x)
	case *big.Float:
		y.Set(x)
	}

	return
}

// Cos is an iterative arbitrary precision computation of Cos(x)
// Iterative process with an error of ~10^{−0.60206*k} = (1/4)^k after k iterations.
// ref : Johansson, B. Tomas, An elementary algorithm to evaluate trigonometric functions to high precision, 2018
func Cos(x *big.Float) (cosx *big.Float) {
	tmp := new(big.Float)

	t := NewFloat(0.5, x.Prec())
	half := new(big.Float).Copy(t)

	for i := uint(1); i < (x.Prec()>>1)-1; i++ {
		t.Mul(t, half)
	}

	s := new(big.Float).Mul(x, t)
	s.Mul(s, x)
	s.Mul(s, t)

	four := NewFloat(4.0, x.Prec())

	for i := uint(1); i < x.Prec()>>1; i++ { // (1/4)^k = (1/2)^(2*k)
		tmp.Sub(four, s)
		s.Mul(s, tmp)
	}

	cosx = new(big.Float).Quo(s, NewFloat(2.0, x.Prec()))
	cosx.Sub(NewFloat(1.0, x.Prec()), cosx)
	return
}

func Sin(x *big.Float) (sinx *big.Float) {
	halfPi := Pi(x.Prec())
	halfPi.Quo(halfPi, new(big.Float).SetInt64(2))
	return Cos(new(big.Float).Sub(x, halfPi))
}

// Log return ln(x) with 2^precisions bits.
func Log(x *big.Float) (ln *big.Float) {
	return bigfloat.Log(x)
}

// Exp returns exp(x) with 2^precisions bits.
func Exp(x *big.Float) (exp *big.Float) {
	return bigfloat.Exp(x)
}

// Pow returns x^y
func Pow(x, y *big.Float) (pow *big.Float) {
	return bigfloat.Pow(x, y)
}

// SinH returns hyperbolic sin(x) with 2^precisions bits.
func SinH(x *big.Float) (sinh *big.Float) {
	sinh = new(big.Float).Set(x)
	sinh.Add(sinh, sinh)
	sinh.Neg(sinh)
	sinh = Exp(sinh)
	sinh.Neg(sinh)
	sinh.Add(sinh, NewFloat(1, x.Prec()))
	tmp := new(big.Float).Set(x)
	tmp.Neg(tmp)
	tmp = Exp(tmp)
	tmp.Add(tmp, tmp)
	sinh.Quo(sinh, tmp)
	return
}

// TanH returns hyperbolic tan(x) with 2^precisions bits.
func TanH(x *big.Float) (tanh *big.Float) {
	tanh = new(big.Float).Set(x)
	tanh.Add(tanh, tanh)
	tanh = Exp(tanh)
	tmp := new(big.Float).Set(tanh)
	tmp.Add(tmp, NewFloat(1, x.Prec()))
	tanh.Sub(tanh, NewFloat(1, x.Prec()))
	tanh.Quo(tanh, tmp)
	return
}

// ArithmeticGeometricMean returns the  arithmetic–geometric mean of x and y with 2^precisions bits.
func ArithmeticGeometricMean(x, y *big.Float) *big.Float {
	precision := x.Prec()
	a := new(big.Float).Set(x)
	g := new(big.Float).Set(y)
	tmp := new(big.Float)
	half := NewFloat(0.5, x.Prec())

	for i := 0; i < int(math.Log2(float64(precision))); i++ {
		tmp.Mul(a, g)
		a.Add(a, g)
		a.Mul(a, half)
		g.Sqrt(tmp)
	}

	return a
}
