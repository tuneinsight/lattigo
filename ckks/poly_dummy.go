package ckks

import (
	"fmt"
	"math"
	"math/bits"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// BSGSCoefficients is a struct storing one or more
// polynomial in their baby-step gian-step form.
// Value: [#Poly][#Degree][#Coeff]
type BSGSCoefficients struct {
	Value [][][]complex128
}

// GetScaledBSGSCoefficients returns the pre-scaled coefficients for the baby-step giant-step polynomial evaluation.
func GetScaledBSGSCoefficients(params Parameters, level int, scale float64, pol interface{}, targetScale float64) (coeffs *BSGSCoefficients, err error) {

	var polyVec PolynomialVector

	switch pol := pol.(type) {
	case Polynomial:
		polyVec = PolynomialVector{Value: []Polynomial{pol}}
	case PolynomialVector:

		var maxDeg int
		var basis BasisType
		for _, poly := range pol.Value {
			maxDeg = utils.MaxInt(maxDeg, poly.MaxDeg)
			basis = poly.BasisType
		}

		for _, poly := range pol.Value {
			if basis != poly.BasisType {
				return nil, fmt.Errorf("polynomial basis must be the same for all polynomials in a polynomial vector")
			}

			if maxDeg != poly.MaxDeg {
				return nil, fmt.Errorf("polynomial degree must all be the same")
			}
		}

		polyVec = pol

	default:
		return nil, fmt.Errorf("EvaluatePoly: invalide pol.(type)")
	}

	return getScaledBSGSCoefficients(params, level, scale, polyVec, targetScale), nil
}

type dummyCiphertext struct {
	Level int
	Scale float64
}

func (d *dummyCiphertext) rescale(params Parameters, minScale float64) {
	var nbRescales int
	for d.Level-nbRescales >= 0 && d.Scale/float64(params.Q()[d.Level-nbRescales]) >= minScale/2 {
		d.Scale /= (float64(params.Q()[d.Level-nbRescales]))
		nbRescales++
	}

	d.Level -= nbRescales
}

type dummyPolynomialBasis struct {
	Value  map[int]*dummyCiphertext
	params Parameters
}

func newDummyPolynomialBasis(params Parameters, ct *dummyCiphertext) (p *dummyPolynomialBasis) {
	p = new(dummyPolynomialBasis)
	p.params = params
	p.Value = make(map[int]*dummyCiphertext)
	p.Value[1] = ct
	return
}

func (p *dummyPolynomialBasis) GenPower(n int, lazy bool, scale float64) {
	if p.Value[n] == nil {
		p.genPower(n, lazy, scale)
		p.Value[n].rescale(p.params, scale)
	}
}

func (p *dummyPolynomialBasis) genPower(n int, lazy bool, scale float64) (err error) {
	if p.Value[n] == nil {

		isPow2 := n&(n-1) == 0

		// Computes the index required to compute the asked ring evaluation
		var a, b int
		if isPow2 {
			a, b = n/2, n/2 //Necessary for optimal depth
		} else {
			// [Lee et al. 2020] : High-Precision and Low-Complexity Approximate Homomorphic Encryption by Error Variance Minimization
			// Maximize the number of odd terms of Chebyshev basis
			k := int(math.Ceil(math.Log2(float64(n)))) - 1
			a = (1 << k) - 1
			b = n + 1 - (1 << k)
		}

		// Recurses on the given indexes
		p.genPower(a, lazy, scale)
		p.genPower(b, lazy, scale)

		p.Value[a].rescale(p.params, scale)
		p.Value[b].rescale(p.params, scale)

		// Computes C[n] = C[a]*C[b]
		p.Value[n] = new(dummyCiphertext)
		p.Value[n].Level = utils.MinInt(p.Value[a].Level, p.Value[b].Level)
		if lazy && !isPow2 {
			p.Value[n].Scale = p.Value[a].Scale * p.Value[b].Scale
		} else {
			p.Value[n].Scale = p.Value[a].Scale * p.Value[b].Scale
			p.Value[n].rescale(p.params, scale)
		}
	}
	return
}

type dummyPolynomialEvaluator struct {
	dummyPolynomialBasis
	params     Parameters
	logDegree  int
	logSplit   int
	isOdd      bool
	isEven     bool
	bsgsCoeffs *BSGSCoefficients
}

func getScaledBSGSCoefficients(params Parameters, level int, scale float64, pol PolynomialVector, targetScale float64) (coeffs *BSGSCoefficients) {

	dummbpb := newDummyPolynomialBasis(params, &dummyCiphertext{level, scale})

	logDegree := bits.Len64(uint64(pol.Value[0].Degree()))
	logSplit := optimalSplit(logDegree)

	var odd, even bool = true, true
	for _, p := range pol.Value {
		tmp0, tmp1 := isOddOrEvenPolynomial(p.Coeffs)
		odd, even = odd && tmp0, even && tmp1
	}

	isRingStandard := params.RingType() == ring.Standard

	for i := (1 << logSplit) - 1; i > 1; i-- {
		if !(even || odd) || (i&1 == 0 && even) || (i&1 == 1 && odd) {
			dummbpb.GenPower(i, isRingStandard, targetScale)
		}
	}

	for i := logSplit; i < logDegree; i++ {
		dummbpb.GenPower(1<<i, false, targetScale)
	}

	polyEval := &dummyPolynomialEvaluator{}
	polyEval.params = params
	polyEval.dummyPolynomialBasis = *dummbpb
	polyEval.logDegree = logDegree
	polyEval.logSplit = logSplit
	polyEval.isOdd = odd
	polyEval.isEven = even
	polyEval.bsgsCoeffs = &BSGSCoefficients{Value: make([][][]complex128, len(pol.Value))}

	polyEval.recurse(dummbpb.Value[1].Level-logDegree+1, targetScale, pol)

	return polyEval.bsgsCoeffs
}

func (polyEval *dummyPolynomialEvaluator) recurse(targetLevel int, targetScale float64, pol PolynomialVector) (res *dummyCiphertext) {

	params := polyEval.params

	logSplit := polyEval.logSplit

	// Recursively computes the evaluation of the Chebyshev polynomial using a baby-set giant-step algorithm.
	if pol.Value[0].Degree() < (1 << logSplit) {

		if pol.Value[0].Lead && polyEval.logSplit > 1 && pol.Value[0].MaxDeg%(1<<(logSplit+1)) > (1<<(logSplit-1)) {

			logDegree := int(bits.Len64(uint64(pol.Value[0].Degree())))
			logSplit := logDegree >> 1

			polyEvalBis := new(dummyPolynomialEvaluator)
			polyEvalBis.params = params
			polyEvalBis.logDegree = logDegree
			polyEvalBis.logSplit = logSplit
			polyEvalBis.dummyPolynomialBasis = polyEval.dummyPolynomialBasis
			polyEvalBis.isOdd = polyEval.isOdd
			polyEvalBis.isEven = polyEval.isEven
			polyEvalBis.bsgsCoeffs = polyEval.bsgsCoeffs

			return polyEvalBis.recurse(targetLevel, targetScale, pol)
		}

		if pol.Value[0].Lead {
			targetScale *= params.QiFloat64(targetLevel)
		}

		return polyEval.evaluatePolyFromPolynomialBasis(targetScale, targetLevel, pol)
	}

	var nextPower = 1 << polyEval.logSplit
	for nextPower < (pol.Value[0].Degree()>>1)+1 {
		nextPower <<= 1
	}

	coeffsq, coeffsr := splitCoeffsPolyVector(pol, nextPower)

	XPow := polyEval.dummyPolynomialBasis.Value[nextPower]

	level := targetLevel

	var currentQi float64
	if pol.Value[0].Lead {
		currentQi = params.QiFloat64(level)
	} else {
		currentQi = params.QiFloat64(level + 1)
	}

	res = polyEval.recurse(targetLevel+1, targetScale*currentQi/XPow.Scale, coeffsq)

	res.rescale(params, params.DefaultScale())

	res.Scale *= XPow.Scale

	tmp := polyEval.recurse(res.Level, res.Scale, coeffsr)

	res.Level = utils.MinInt(tmp.Level, res.Level)

	return
}

func (polyEval *dummyPolynomialEvaluator) evaluatePolyFromPolynomialBasis(targetScale float64, level int, pol PolynomialVector) (res *dummyCiphertext) {

	X := polyEval.dummyPolynomialBasis.Value

	minimumDegreeNonZeroCoefficient := len(pol.Value[0].Coeffs) - 1

	if polyEval.isEven {
		minimumDegreeNonZeroCoefficient--
	}

	values := make([][]complex128, minimumDegreeNonZeroCoefficient+1)
	for i := range values {
		values[i] = make([]complex128, len(pol.Value))
	}

	for i, p := range pol.Value {
		if isNotNegligible(p.Coeffs[0]) {
			values[0][i] = p.Coeffs[0] * complex(targetScale, 0)
		}
	}

	for i := 1; i < pol.Value[0].Degree(); i++ {
		for j, p := range pol.Value {
			if isNotNegligible(p.Coeffs[i]) {
				values[i][j] = p.Coeffs[i] * complex(targetScale/X[i].Scale, 0)
			}
		}
	}

	polyEval.bsgsCoeffs.Value = append(polyEval.bsgsCoeffs.Value, values)

	return &dummyCiphertext{level, targetScale}
}
