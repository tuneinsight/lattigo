package rlwe

// PolynomialBasisType is a type for the polynomials basis
type PolynomialBasisType int

const (
	// Monomial : x^(a+b) = x^a * x^b
	Monomial = PolynomialBasisType(0)
	// Chebyshev : T_(a+b) = 2 * T_a * T_b - T_(|a-b|)
	Chebyshev = PolynomialBasisType(1)
)

// Polynomial is an interface to manage different types of polynomials
type Polynomial interface {
	BasisType() PolynomialBasisType
	Depth() (depth int)
	Degree() (degree int)
	OddEven() (odd, even bool)
	BSGSSplit() (giant, baby int)
	SplitBSGS(split int) (coeffq, coeffsr Polynomial)
	MarshalBinary() (data []byte, err error)
	UnmarshalBinary(data []byte) (err error)
	IsEncrypted() bool
}
