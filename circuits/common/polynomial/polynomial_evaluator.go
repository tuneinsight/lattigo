package polynomial

import (
	"fmt"
	"math/bits"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// CoefficientGetter defines an interface to get the coefficients of a Polynomial.
type CoefficientGetter[T uint64 | *bignum.Complex] interface {
	// GetVectorCoefficient should return a slice []T containing the k-th coefficient
	// of each polynomial of PolynomialVector indexed by its Mapping.
	// See PolynomialVector for additional information about the Mapping.
	GetVectorCoefficient(pol PolynomialVector, k int) (values []T)
	// GetSingleCoefficient should return the k-th coefficient of Polynomial as the type T.
	GetSingleCoefficient(pol Polynomial, k int) (value T)
}

type Evaluator[T uint64 | *bignum.Complex] struct {
	schemes.Evaluator
	CoefficientGetter[T]
}

// Evaluate is a generic and scheme agnostic method to evaluate polynomials on rlwe.Ciphertexts.
func (eval Evaluator[T]) Evaluate(input interface{}, p interface{}, targetScale rlwe.Scale, levelsConsumedPerRescaling int, SimEval SimEvaluator) (opOut *rlwe.Ciphertext, err error) {

	var polyVec PolynomialVector
	switch p := p.(type) {
	case bignum.Polynomial:
		polyVec = PolynomialVector{Value: []Polynomial{{Polynomial: p, MaxDeg: p.Degree(), Lead: true, Lazy: false}}}
	case Polynomial:
		polyVec = PolynomialVector{Value: []Polynomial{p}}
	case PolynomialVector:
		polyVec = p
	default:
		return nil, fmt.Errorf("cannot Polynomial: invalid polynomial type, must be either bignum.Polynomial, polynomial.Polynomial or polynomial.PolynomialVector, but is %T", p)
	}

	var powerbasis PowerBasis
	switch input := input.(type) {
	case *rlwe.Ciphertext:
		powerbasis = NewPowerBasis(input, polyVec.Value[0].Basis)
	case PowerBasis:
		if input.Value[1] == nil {
			return nil, fmt.Errorf("cannot evaluatePolyVector: given PowerBasis.Value[1] is empty")
		}
		powerbasis = input
	default:
		return nil, fmt.Errorf("cannot evaluatePolyVector: invalid input, must be either *rlwe.Ciphertext or *PowerBasis")
	}

	if level, depth := powerbasis.Value[1].Level(), levelsConsumedPerRescaling*polyVec.Value[0].Depth(); level < depth {
		return nil, fmt.Errorf("%d levels < %d log(d) -> cannot evaluate poly", level, depth)
	}

	logDegree := bits.Len64(uint64(polyVec.Value[0].Degree()))
	logSplit := bignum.OptimalSplit(logDegree)

	var odd, even = false, false
	for _, p := range polyVec.Value {
		odd, even = odd || p.IsOdd, even || p.IsEven
	}

	// Computes all the powers of two with relinearization
	// This will recursively compute and store all powers of two up to 2^logDegree
	if err = powerbasis.GenPower(1<<(logDegree-1), false, eval); err != nil {
		return nil, err
	}

	// Computes the intermediate powers, starting from the largest, without relinearization if possible
	for i := (1 << logSplit) - 1; i > 2; i-- {
		if !(even || odd) || (i&1 == 0 && even) || (i&1 == 1 && odd) {
			if err = powerbasis.GenPower(i, polyVec.Value[0].Lazy, eval); err != nil {
				return nil, err
			}
		}
	}

	PS := polyVec.PatersonStockmeyerPolynomial(*eval.GetRLWEParameters(), powerbasis.Value[1].Level(), powerbasis.Value[1].Scale, targetScale, SimEval)

	if opOut, err = eval.EvaluatePatersonStockmeyerPolynomialVector(PS, powerbasis); err != nil {
		return nil, err
	}

	return opOut, err
}

// BabyStep is a struct storing the result of a baby-step
// of the Paterson-Stockmeyer polynomial evaluation algorithm.
type BabyStep struct {
	Degree int
	Value  *rlwe.Ciphertext
}

// EvaluatePatersonStockmeyerPolynomialVector evaluates a pre-decomposed PatersonStockmeyerPolynomialVector on a pre-computed power basis [1, X^{1}, X^{2}, ..., X^{2^{n}}, X^{2^{n+1}}, ..., X^{2^{m}}]
func (eval Evaluator[T]) EvaluatePatersonStockmeyerPolynomialVector(poly PatersonStockmeyerPolynomialVector, pb PowerBasis) (res *rlwe.Ciphertext, err error) {

	split := len(poly.Value[0].Value)

	babySteps := make([]*BabyStep, split)

	// Small steps
	for i := range babySteps {

		// eval & cg are not thread-safe
		if babySteps[split-i-1], err = eval.EvaluateBabyStep(i, poly, pb); err != nil {
			return nil, fmt.Errorf("cannot EvaluateBabyStep: %w", err)
		}
	}

	// Loops as long as there is more than one sub-polynomial
	for len(babySteps) != 1 {

		// Precomputes the ops to apply in the giant steps loop
		giantsteps := make([]int, len(babySteps))
		for i := 0; i < len(babySteps); i++ {
			if i == len(babySteps)-1 {
				giantsteps[i] = 2
			} else if babySteps[i].Degree == babySteps[i+1].Degree {
				giantsteps[i] = 1
				i++
			}
		}

		for i := 0; i < len(babySteps); i++ {

			// eval is not thread-safe
			if err = eval.EvaluateGianStep(i, giantsteps, babySteps, pb); err != nil {
				return nil, err
			}
		}

		// Discards processed sub-polynomials
		var idx int
		for i := range babySteps {
			if babySteps[i] != nil {
				babySteps[idx] = babySteps[i]
				idx++
			}
		}

		babySteps = babySteps[:idx]
	}

	if babySteps[0].Value.Degree() == 2 {
		if err = eval.Relinearize(babySteps[0].Value, babySteps[0].Value); err != nil {
			return nil, fmt.Errorf("cannot EvaluatePatersonStockmeyerPolynomial: %w", err)
		}
	}

	if err = eval.Rescale(babySteps[0].Value, babySteps[0].Value); err != nil {
		return nil, err
	}

	return babySteps[0].Value, nil
}

// EvaluateBabyStep evaluates a baby-step of the PatersonStockmeyer polynomial evaluation algorithm, i.e. the inner-product between the precomputed
// powers [1, T, T^2, ..., T^{n-1}] and the coefficients [ci0, ci1, ci2, ..., ci{n-1}].
func (eval Evaluator[T]) EvaluateBabyStep(i int, poly PatersonStockmeyerPolynomialVector, pb PowerBasis) (ct *BabyStep, err error) {

	nbPoly := len(poly.Value)

	polyVec := PolynomialVector{
		Value:   make([]Polynomial, nbPoly),
		Mapping: poly.Mapping,
	}

	// Transposes the polynomial matrix
	for j := 0; j < nbPoly; j++ {
		polyVec.Value[j] = poly.Value[j].Value[i]
	}

	level := poly.Value[0].Value[i].Level
	scale := poly.Value[0].Value[i].Scale

	ct = new(BabyStep)
	ct.Degree = poly.Value[0].Value[i].Degree()
	if ct.Value, err = eval.EvaluatePolynomialVectorFromPowerBasis(level, polyVec, pb, scale); err != nil {
		return ct, fmt.Errorf("cannot EvaluatePolynomialVectorFromPowerBasis: polynomial[%d]: %w", i, err)
	}

	return ct, nil
}

// EvaluateGianStep evaluates a giant-step of the PatersonStockmeyer polynomial evaluation algorithm, which consists
// in combining the baby-steps <[1, T, T^2, ..., T^{n-1}], [ci0, ci1, ci2, ..., ci{n-1}]> together with powers T^{2^k}.
func (eval Evaluator[T]) EvaluateGianStep(i int, giantSteps []int, babySteps []*BabyStep, pb PowerBasis) (err error) {

	// If we reach the end of the list it means we weren't able to combine
	// the last two sub-polynomials which necessarily implies that that the
	// last one has degree smaller than the previous one and that there is
	// no next polynomial to combine it with.
	// Therefore we update it's degree to the one of the previous one.
	if giantSteps[i] == 2 {
		babySteps[i].Degree = babySteps[i-1].Degree

		// If two consecutive sub-polynomials, from ascending degree order, have the
		// same degree, we combine them.
	} else if giantSteps[i] == 1 {

		even, odd := babySteps[i], babySteps[i+1]

		deg := 1 << bits.Len64(uint64(babySteps[i].Degree))

		if err = eval.EvaluateMonomial(even.Value, odd.Value, pb.Value[deg]); err != nil {
			return
		}

		odd.Degree = 2*deg - 1
		babySteps[i] = nil

		i++
	}

	return
}

// EvaluateMonomial evaluates a monomial of the form a + b * X^{pow} and writes the results in b.
func (eval Evaluator[T]) EvaluateMonomial(a, b, xpow *rlwe.Ciphertext) (err error) {

	if b.Degree() == 2 {
		if err = eval.Relinearize(b, b); err != nil {
			return fmt.Errorf("evalMonomial: %w", err)
		}
	}

	if err = eval.Rescale(b, b); err != nil {
		return fmt.Errorf("evalMonomial: %w", err)
	}

	if err = eval.Mul(b, xpow, b); err != nil {
		return fmt.Errorf("evalMonomial: %w", err)
	}

	if !a.Scale.InDelta(b.Scale, float64(rlwe.ScalePrecision-12)) {
		return fmt.Errorf("evalMonomial: scale discrepency: (rescale(b) * X^{n}).Scale = %v != a.Scale = %v", &a.Scale.Value, &b.Scale.Value)
	}

	if err = eval.Add(b, a, b); err != nil {
		return fmt.Errorf("evalMonomial: %w", err)
	}

	return
}

// EvaluatePolynomialVectorFromPowerBasis evaluates P(ct) = sum c_i * ct^{i}.
func (eval Evaluator[T]) EvaluatePolynomialVectorFromPowerBasis(targetLevel int, pol PolynomialVector, pb PowerBasis, targetScale rlwe.Scale) (res *rlwe.Ciphertext, err error) {

	// Map[int] of the powers [X^{0}, X^{1}, X^{2}, ...]
	X := pb.Value

	params := eval.GetRLWEParameters()
	mapping := pol.Mapping
	even := pol.IsEven()
	odd := pol.IsOdd()

	// Retrieve the degree of the highest degree non-zero coefficient
	// TODO: optimize for nil/zero coefficients
	minimumDegreeNonZeroCoefficient := len(pol.Value[0].Coeffs) - 1
	if even && !odd {
		minimumDegreeNonZeroCoefficient--
	}

	// Gets the maximum degree of the ciphertexts among the power basis
	// TODO: optimize for nil/zero coefficients, odd/even polynomial
	maximumCiphertextDegree := 0
	for i := pol.Value[0].Degree(); i > 0; i-- {
		if x, ok := X[i]; ok {
			maximumCiphertextDegree = utils.Max(maximumCiphertextDegree, x.Degree())
		}
	}

	// If an index slot is given (either multiply polynomials or masking)
	if mapping != nil {

		// If the degree of the poly is zero
		if minimumDegreeNonZeroCoefficient == 0 {

			// Allocates the output ciphertext
			res = rlwe.NewCiphertext(params, 1, targetLevel)
			*res.MetaData = *X[1].MetaData
			res.Scale = targetScale

			if even {

				if err = eval.Add(res, eval.GetVectorCoefficient(pol, 0), res); err != nil {
					return nil, err
				}
			}

			return
		}

		// Allocates the output ciphertext
		res = rlwe.NewCiphertext(params, maximumCiphertextDegree, targetLevel)
		*res.MetaData = *X[1].MetaData
		res.Scale = targetScale

		if even {
			if err = eval.Add(res, eval.GetVectorCoefficient(pol, 0), res); err != nil {
				return nil, err
			}
		}

		// Loops starting from the highest degree coefficient
		for key := pol.Value[0].Degree(); key > 0; key-- {
			if !(even || odd) || (key&1 == 0 && even) || (key&1 == 1 && odd) {
				if err = eval.MulThenAdd(X[key], eval.GetVectorCoefficient(pol, key), res); err != nil {
					return
				}
			}
		}

	} else {

		if minimumDegreeNonZeroCoefficient == 0 {

			res = rlwe.NewCiphertext(params, 1, targetLevel)
			*res.MetaData = *X[1].MetaData
			res.Scale = targetScale

			if even {
				if err = eval.Add(res, eval.GetSingleCoefficient(pol.Value[0], 0), res); err != nil {
					return
				}
			}

			return
		}

		res = rlwe.NewCiphertext(params, maximumCiphertextDegree, targetLevel)
		*res.MetaData = *X[1].MetaData
		res.Scale = targetScale

		if even {
			if err = eval.Add(res, eval.GetSingleCoefficient(pol.Value[0], 0), res); err != nil {
				return
			}
		}

		for key := pol.Value[0].Degree(); key > 0; key-- {
			if key != 0 && (!(even || odd) || (key&1 == 0 && even) || (key&1 == 1 && odd)) {
				// MulScalarAndAdd automatically scales c to match the scale of res.
				if err = eval.MulThenAdd(X[key], eval.GetSingleCoefficient(pol.Value[0], key), res); err != nil {
					return
				}
			}
		}
	}

	return
}
