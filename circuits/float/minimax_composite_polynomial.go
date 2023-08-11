package float

import (
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// MinimaxCompositePolynomial is a struct storing P(x) = pk(x) o pk-1(x) o ... o p1(x) o p0(x).
type MinimaxCompositePolynomial []bignum.Polynomial

// NewMinimaxCompositePolynomial creates a new MinimaxCompositePolynomial from a list of coefficients.
// Coefficients are expected to be given in the Chebyshev basis.
func NewMinimaxCompositePolynomial(coeffsStr [][]string) MinimaxCompositePolynomial {
	polys := make([]bignum.Polynomial, len(coeffsStr))

	for i := range coeffsStr {

		coeffs := parseCoeffs(coeffsStr[i])

		poly := bignum.NewPolynomial(
			bignum.Chebyshev,
			coeffs,
			&bignum.Interval{
				A: *bignum.NewFloat(-1, coeffs[0].Prec()),
				B: *bignum.NewFloat(1, coeffs[0].Prec()),
			},
		)

		polys[i] = poly
	}

	return MinimaxCompositePolynomial(polys)
}

func (mcp MinimaxCompositePolynomial) MaxDepth() (depth int) {
	for i := range mcp {
		depth = utils.Max(depth, mcp[i].Depth())
	}
	return
}
