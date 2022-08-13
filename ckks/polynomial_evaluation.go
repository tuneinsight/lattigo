package ckks

import (
	"encoding/binary"
	"fmt"
	"math"
	"math/big"
	"math/bits"
	"runtime"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// Polynomial is a struct storing the coefficients of a polynomial
// that then can be evaluated on the ciphertext
type Polynomial struct {
	BasisType
	MaxDeg int
	Coeffs []complex128
	Lead   bool
	A      float64
	B      float64
}

// BasisType is a type for the polynomials basis
type BasisType int

const (
	// Monomial : x^(a+b) = x^a * x^b
	Monomial = BasisType(0)
	// Chebyshev : T_(a+b) = 2 * T_a * T_b - T_(|a-b|)
	Chebyshev = BasisType(1)
)

// IsNegligbleThreshold : threshold under which a coefficient
// of a polynomial is ignored.
const IsNegligbleThreshold float64 = 1e-14

// Depth returns the number of levels needed to evaluate the polynomial.
func (p *Polynomial) Depth() int {
	return int(math.Ceil(math.Log2(float64(len(p.Coeffs)))))
}

// Degree returns the degree of the polynomial
func (p *Polynomial) Degree() int {
	return len(p.Coeffs) - 1
}

// NewPoly creates a new Poly from the input coefficients
func NewPoly(coeffs []complex128) (p *Polynomial) {
	c := make([]complex128, len(coeffs))
	copy(c, coeffs)
	return &Polynomial{Coeffs: c, MaxDeg: len(c) - 1, Lead: true}
}

// checkEnoughLevels checks that enough levels are available to evaluate the polynomial.
// Also checks if c is a Gaussian integer or not. If not, then one more level is needed
// to evaluate the polynomial.
func checkEnoughLevels(levels, depth int, c complex128) (err error) {

	if real(c) != float64(int64(real(c))) || imag(c) != float64(int64(imag(c))) {
		depth++
	}

	if levels < depth {
		return fmt.Errorf("%d levels < %d log(d) -> cannot evaluate", levels, depth)
	}

	return nil
}

type polynomialEvaluator struct {
	Evaluator
	Encoder
	PolynomialBasis
	slotsIndex map[int][]int
	logDegree  int
	logSplit   int
	isOdd      bool
	isEven     bool
}

// EvaluatePoly evaluates a polynomial in standard basis on the input Ciphertext in ceil(log2(deg+1)) levels.
// Returns an error if the input ciphertext does not have enough level to carry out the full polynomial evaluation.
// Returns an error if something is wrong with the scale.
// If the polynomial is given in Chebyshev basis, then a change of basis ct' = (2/(b-a)) * (ct + (-a-b)/(b-a))
// is necessary before the polynomial evaluation to ensure correctness.
// Coefficients of the polynomial with an absolute value smaller than "IsNegligbleThreshold" will automatically be set to zero
// if the polynomial is "even" or "odd" (to ensure that the even or odd property remains valid
// after the "splitCoeffs" polynomial decomposition).
// input must be either *Ciphertext or *PolynomialBasis.
// pol: a *Polynomial
// targetScale: the desired output scale. This value shouldn't differ too much from the original ciphertext scale. It can
// for example be used to correct small deviations in the ciphertext scale and reset it to the default scale.
func (eval *evaluator) EvaluatePoly(input interface{}, pol *Polynomial, targetScale float64) (opOut *Ciphertext, err error) {
	return eval.evaluatePolyVector(input, polynomialVector{Value: []*Polynomial{pol}}, targetScale)
}

type polynomialVector struct {
	Encoder    Encoder
	Value      []*Polynomial
	SlotsIndex map[int][]int
}

// EvaluatePolyVector evaluates a vector of Polynomials on the input Ciphertext in ceil(log2(deg+1)) levels.
// Returns an error if the input Ciphertext does not have enough level to carry out the full polynomial evaluation.
// Returns an error if something is wrong with the scale.
// Returns an error if polynomials are not all in the same basis.
// Returns an error if polynomials do not all have the same degree.
// If the polynomials are given in Chebyshev basis, then a change of basis ct' = (2/(b-a)) * (ct + (-a-b)/(b-a))
// is necessary before the polynomial evaluation to ensure correctness.
// Coefficients of the polynomial with an absolute value smaller than "IsNegligbleThreshold" will automatically be set to zero
// if the polynomial is "even" or "odd" (to ensure that the even or odd property remains valid
// after the "splitCoeffs" polynomial decomposition).
// input: must be either *Ciphertext or *PolynomialBasis.
// pols: a slice of up to 'n' *Polynomial ('n' being the maximum number of slots), indexed from 0 to n-1.
// encoder: an Encoder.
// slotsIndex: a map[int][]int indexing as key the polynomial to evaluate and as value the index of the slots on which to evaluate the polynomial indexed by the key.
// targetScale: the desired output scale. This value shouldn't differ too much from the original ciphertext scale. It can
// for example be used to correct small deviations in the ciphertext scale and reset it to the default scale.
//
// Example: if pols = []*Polynomial{pol0, pol1} and slotsIndex = map[int][]int:{0:[1, 2, 4, 5, 7], 1:[0, 3]},
// then pol0 will be applied to slots [1, 2, 4, 5, 7], pol1 to slots [0, 3] and the slot 6 will be zero-ed.
func (eval *evaluator) EvaluatePolyVector(input interface{}, pols []*Polynomial, encoder Encoder, slotsIndex map[int][]int, targetScale float64) (opOut *Ciphertext, err error) {
	var maxDeg int
	var basis BasisType
	for i := range pols {
		maxDeg = utils.MaxInt(maxDeg, pols[i].MaxDeg)
		basis = pols[i].BasisType
	}

	for i := range pols {
		if basis != pols[i].BasisType {
			return nil, fmt.Errorf("polynomial basis must be the same for all polynomials in a polynomial vector")
		}

		if maxDeg != pols[i].MaxDeg {
			return nil, fmt.Errorf("polynomial degree must all be the same")
		}
	}

	return eval.evaluatePolyVector(input, polynomialVector{Encoder: encoder, Value: pols, SlotsIndex: slotsIndex}, targetScale)
}

func optimalSplit(logDegree int) (logSplit int) {
	logSplit = logDegree >> 1
	a := (1 << logSplit) + (1 << (logDegree - logSplit)) + logDegree - logSplit - 3
	b := (1 << (logSplit + 1)) + (1 << (logDegree - logSplit - 1)) + logDegree - logSplit - 4
	if a > b {
		logSplit++
	}

	return
}

func (eval *evaluator) evaluatePolyVector(input interface{}, pol polynomialVector, targetScale float64) (opOut *Ciphertext, err error) {

	if pol.SlotsIndex != nil && pol.Encoder == nil {
		return nil, fmt.Errorf("cannot EvaluatePolyVector: missing Encoder input")
	}

	var monomialBasis *PolynomialBasis
	switch input := input.(type) {
	case *Ciphertext:
		monomialBasis = NewPolynomialBasis(input, pol.Value[0].BasisType)
	case *PolynomialBasis:
		if input.Value[1] == nil {
			return nil, fmt.Errorf("cannot evaluatePolyVector: given PolynomialBasis.Value[1] is empty")
		}
		monomialBasis = input
	default:
		return nil, fmt.Errorf("cannot evaluatePolyVector: invalid input, must be either *Ciphertext or *PolynomialBasis")
	}

	if err := checkEnoughLevels(monomialBasis.Value[1].Level(), pol.Value[0].Depth(), 1); err != nil {
		return nil, err
	}

	logDegree := bits.Len64(uint64(pol.Value[0].Degree()))
	logSplit := optimalSplit(logDegree)

	var odd, even bool = true, true
	for _, p := range pol.Value {
		tmp0, tmp1 := isOddOrEvenPolynomial(p.Coeffs)
		odd, even = odd && tmp0, even && tmp1
	}

	isRingStandard := eval.params.RingType() == ring.Standard

	for i := (1 << logSplit) - 1; i > 1; i-- {
		if !(even || odd) || (i&1 == 0 && even) || (i&1 == 1 && odd) {
			if err = monomialBasis.GenPower(i, isRingStandard, targetScale, eval); err != nil {
				return nil, err
			}
		}
	}

	for i := logSplit; i < logDegree; i++ {
		if err = monomialBasis.GenPower(1<<i, false, targetScale, eval); err != nil {
			return nil, err
		}
	}

	polyEval := &polynomialEvaluator{}
	polyEval.slotsIndex = pol.SlotsIndex
	polyEval.Evaluator = eval
	polyEval.Encoder = pol.Encoder
	polyEval.PolynomialBasis = *monomialBasis
	polyEval.logDegree = logDegree
	polyEval.logSplit = logSplit
	polyEval.isOdd = odd
	polyEval.isEven = even

	if opOut, err = polyEval.recurse(monomialBasis.Value[1].Level()-logDegree+1, targetScale, pol); err != nil {
		return nil, err
	}

	polyEval.Relinearize(opOut, opOut)

	if err = polyEval.Rescale(opOut, targetScale, opOut); err != nil {
		return nil, err
	}

	opOut.Scale = targetScale // solves float64 precision issues

	polyEval = nil
	runtime.GC()
	return opOut, err
}

// PolynomialBasis is a struct storing powers of a ciphertext.
type PolynomialBasis struct {
	BasisType
	Value map[int]*Ciphertext
	lazy  bool
}

// NewPolynomialBasis creates a new PolynomialBasis. It takes as input a ciphertext
// and a basistype. The struct treates the input ciphertext as a monomial X and
// can be used to generates power of this monomial X^{n} in the given BasisType.
func NewPolynomialBasis(ct *Ciphertext, basistype BasisType) (p *PolynomialBasis) {
	p = new(PolynomialBasis)
	p.Value = make(map[int]*Ciphertext)
	p.Value[1] = ct.CopyNew()
	p.BasisType = basistype
	return
}

// GenPower recursively computes X^{n}.
// If lazy = true, the final X^{n} will not be relinearized.
// Previous non-relinearized X^{n} that are required to compute the target X^{n} are automatically relinearized.
// Scale sets the threshold for rescaling (ciphertext won't be rescaled if the rescaling operation would make the scale go under this threshold).
func (p *PolynomialBasis) GenPower(n int, lazy bool, scale float64, eval Evaluator) (err error) {

	if p.Value[n] == nil {
		if err = p.genPower(n, lazy, scale, eval); err != nil {
			return
		}

		if err = eval.Rescale(p.Value[n], scale, p.Value[n]); err != nil {
			return
		}
	}

	return nil
}

func (p *PolynomialBasis) genPower(n int, lazy bool, scale float64, eval Evaluator) (err error) {
	if p.Value[n] == nil {

		isPow2 := n&(n-1) == 0

		// Computes the index required to compute the asked ring evaluation
		var a, b, c int
		if isPow2 {
			a, b = n/2, n/2 //Necessary for optimal depth
		} else {
			// [Lee et al. 2020] : High-Precision and Low-Complexity Approximate Homomorphic Encryption by Error Variance Minimization
			// Maximize the number of odd terms of Chebyshev basis
			k := int(math.Ceil(math.Log2(float64(n)))) - 1
			a = (1 << k) - 1
			b = n + 1 - (1 << k)

			if p.BasisType == Chebyshev {
				c = int(math.Abs(float64(a) - float64(b))) // Cn = 2*Ca*Cb - Cc, n = a+b and c = abs(a-b)
			}
		}

		// Recurses on the given indexes
		if err = p.genPower(a, lazy && !isPow2, scale, eval); err != nil {
			return err
		}
		if err = p.genPower(b, lazy && !isPow2, scale, eval); err != nil {
			return err
		}

		// Computes C[n] = C[a]*C[b]
		if lazy {
			if p.Value[a].Degree() == 2 {
				eval.Relinearize(p.Value[a], p.Value[a])
			}

			if p.Value[b].Degree() == 2 {
				eval.Relinearize(p.Value[b], p.Value[b])
			}

			if err = eval.Rescale(p.Value[a], scale, p.Value[a]); err != nil {
				return err
			}

			if err = eval.Rescale(p.Value[b], scale, p.Value[b]); err != nil {
				return err
			}

			p.Value[n] = eval.MulNew(p.Value[a], p.Value[b])

		} else {
			p.Value[n] = eval.MulRelinNew(p.Value[a], p.Value[b])

			if err = eval.Rescale(p.Value[n], scale, p.Value[n]); err != nil {
				return err
			}
		}

		if p.BasisType == Chebyshev {

			// Computes C[n] = 2*C[a]*C[b]
			eval.Add(p.Value[n], p.Value[n], p.Value[n])

			// Computes C[n] = 2*C[a]*C[b] - C[c]
			if c == 0 {
				eval.AddConst(p.Value[n], -1, p.Value[n])
			} else {
				// Since C[0] is not stored (but rather seen as the constant 1), only recurses on c if c!= 0
				if err = p.GenPower(c, lazy, scale, eval); err != nil {
					return err
				}
				eval.Sub(p.Value[n], p.Value[c], p.Value[n])
			}
		}
	}
	return
}

// MarshalBinary encodes the target on a slice of bytes.
func (p *PolynomialBasis) MarshalBinary() (data []byte, err error) {
	data = make([]byte, 16)
	binary.LittleEndian.PutUint64(data[0:8], uint64(len(p.Value)))
	binary.LittleEndian.PutUint64(data[8:16], uint64(p.Value[1].GetDataLen(true)))
	for key, ct := range p.Value {
		keyBytes := make([]byte, 8)
		binary.LittleEndian.PutUint64(keyBytes, uint64(key))
		data = append(data, keyBytes...)
		ctBytes, err := ct.MarshalBinary()
		if err != nil {
			return []byte{}, err
		}
		data = append(data, ctBytes...)
	}
	return
}

// UnmarshalBinary decodes a slice of bytes on the target.
func (p *PolynomialBasis) UnmarshalBinary(data []byte) (err error) {
	p.Value = make(map[int]*Ciphertext)
	nbct := int(binary.LittleEndian.Uint64(data[0:8]))
	dtLen := int(binary.LittleEndian.Uint64(data[8:16]))
	ptr := 16
	for i := 0; i < nbct; i++ {
		idx := int(binary.LittleEndian.Uint64(data[ptr : ptr+8]))
		ptr += 8
		p.Value[idx] = new(Ciphertext)
		if err = p.Value[idx].UnmarshalBinary(data[ptr : ptr+dtLen]); err != nil {
			return
		}
		ptr += dtLen
	}
	return
}

func splitCoeffs(coeffs *Polynomial, split int) (coeffsq, coeffsr *Polynomial) {

	// Splits a polynomial p such that p = q*C^degree + r.
	coeffsr = &Polynomial{}
	coeffsr.Coeffs = make([]complex128, split)
	if coeffs.MaxDeg == coeffs.Degree() {
		coeffsr.MaxDeg = split - 1
	} else {
		coeffsr.MaxDeg = coeffs.MaxDeg - (coeffs.Degree() - split + 1)
	}

	for i := 0; i < split; i++ {
		coeffsr.Coeffs[i] = coeffs.Coeffs[i]
	}

	coeffsq = &Polynomial{}
	coeffsq.Coeffs = make([]complex128, coeffs.Degree()-split+1)
	coeffsq.MaxDeg = coeffs.MaxDeg

	coeffsq.Coeffs[0] = coeffs.Coeffs[split]

	if coeffs.BasisType == Monomial {
		for i := split + 1; i < coeffs.Degree()+1; i++ {
			coeffsq.Coeffs[i-split] = coeffs.Coeffs[i]
		}
	} else if coeffs.BasisType == Chebyshev {
		for i, j := split+1, 1; i < coeffs.Degree()+1; i, j = i+1, j+1 {
			coeffsq.Coeffs[i-split] = 2 * coeffs.Coeffs[i]
			coeffsr.Coeffs[split-j] -= coeffs.Coeffs[i]
		}
	}

	if coeffs.Lead {
		coeffsq.Lead = true
	}

	coeffsq.BasisType, coeffsr.BasisType = coeffs.BasisType, coeffs.BasisType

	return
}

func splitCoeffsPolyVector(poly polynomialVector, split int) (polyq, polyr polynomialVector) {
	coeffsq := make([]*Polynomial, len(poly.Value))
	coeffsr := make([]*Polynomial, len(poly.Value))
	for i, p := range poly.Value {
		coeffsq[i], coeffsr[i] = splitCoeffs(p, split)
	}

	return polynomialVector{Value: coeffsq}, polynomialVector{Value: coeffsr}
}

func (polyEval *polynomialEvaluator) recurse(targetLevel int, targetScale float64, pol polynomialVector) (res *Ciphertext, err error) {

	params := polyEval.Evaluator.(*evaluator).params

	logSplit := polyEval.logSplit

	// Recursively computes the evaluation of the Chebyshev polynomial using a baby-set giant-step algorithm.
	if pol.Value[0].Degree() < (1 << logSplit) {

		if pol.Value[0].Lead && polyEval.logSplit > 1 && pol.Value[0].MaxDeg%(1<<(logSplit+1)) > (1<<(logSplit-1)) {

			logDegree := int(bits.Len64(uint64(pol.Value[0].Degree())))
			logSplit := logDegree >> 1

			polyEvalBis := new(polynomialEvaluator)
			polyEvalBis.Evaluator = polyEval.Evaluator
			polyEvalBis.Encoder = polyEval.Encoder
			polyEvalBis.slotsIndex = polyEval.slotsIndex
			polyEvalBis.logDegree = logDegree
			polyEvalBis.logSplit = logSplit
			polyEvalBis.PolynomialBasis = polyEval.PolynomialBasis
			polyEvalBis.isOdd = polyEval.isOdd
			polyEvalBis.isEven = polyEval.isEven

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

	XPow := polyEval.PolynomialBasis.Value[nextPower]

	level := targetLevel

	var currentQi float64
	if pol.Value[0].Lead {
		currentQi = params.QiFloat64(level)
	} else {
		currentQi = params.QiFloat64(level + 1)
	}

	if res, err = polyEval.recurse(targetLevel+1, targetScale*currentQi/XPow.Scale, coeffsq); err != nil {
		return nil, err
	}

	if res.Degree() == 2 {
		polyEval.Relinearize(res, res)
	}

	if err = polyEval.Rescale(res, params.DefaultScale(), res); err != nil {
		return nil, err
	}

	polyEval.Mul(res, XPow, res)

	var tmp *Ciphertext
	if tmp, err = polyEval.recurse(res.Level(), res.Scale, coeffsr); err != nil {
		return nil, err
	}

	polyEval.Add(res, tmp, res)

	tmp = nil

	return
}

func (polyEval *polynomialEvaluator) evaluatePolyFromPolynomialBasis(targetScale float64, level int, pol polynomialVector) (res *Ciphertext, err error) {

	X := polyEval.PolynomialBasis.Value

	params := polyEval.Evaluator.(*evaluator).params
	slotsIndex := polyEval.slotsIndex

	minimumDegreeNonZeroCoefficient := len(pol.Value[0].Coeffs) - 1

	if polyEval.isEven {
		minimumDegreeNonZeroCoefficient--
	}

	// Get the minimum non-zero degree coefficient
	maximumCiphertextDegree := 0
	for i := pol.Value[0].Degree(); i > 0; i-- {
		if x, ok := X[i]; ok {
			maximumCiphertextDegree = utils.MaxInt(maximumCiphertextDegree, x.Degree())
		}
	}

	// If an index slot is given (either multiply polynomials or masking)
	if slotsIndex != nil {

		var toEncode bool

		// Allocates temporary buffer for coefficients encoding
		values := make([]complex128, params.Slots())

		// If the degree of the poly is zero
		if minimumDegreeNonZeroCoefficient == 0 {

			// Allocates the output ciphertext
			res = NewCiphertext(params, 1, level, targetScale)

			// Looks for non-zero coefficients among the degree 0 coefficients of the polynomials
			for i, p := range pol.Value {
				if isNotNegligible(p.Coeffs[0]) {
					toEncode = true
					for _, j := range slotsIndex[i] {
						values[j] = p.Coeffs[0]
					}
				}
			}

			// If a non-zero coefficient was found, encode the values, adds on the ciphertext, and returns
			if toEncode {
				pt := NewPlaintextAtLevelFromPoly(level, res.Value[0])
				pt.Scale = res.Scale
				polyEval.EncodeSlots(values, pt, params.LogSlots())
			}

			return
		}

		// Allocates the output ciphertext
		res = NewCiphertext(params, maximumCiphertextDegree, level, targetScale)

		// Allocates a temporary plaintext to encode the values
		pt := NewPlaintextAtLevelFromPoly(level, polyEval.Evaluator.BuffCt().Value[0])

		// Looks for a non-zero coefficient among the degree zero coefficient of the polynomials
		for i, p := range pol.Value {
			if isNotNegligible(p.Coeffs[0]) {
				toEncode = true
				for _, j := range slotsIndex[i] {
					values[j] = p.Coeffs[0]
				}
			}
		}

		// If a non-zero degre coefficient was found, encode and adds the values on the output
		// ciphertext
		if toEncode {
			pt.Scale = res.Scale
			polyEval.EncodeSlots(values, pt, params.LogSlots())
			polyEval.Add(res, pt, res)
			toEncode = false
		}

		// Loops starting from the highest degree coefficient
		for key := pol.Value[0].Degree(); key > 0; key-- {

			var reset bool
			// Loops over the polynomials
			for i, p := range pol.Value {

				// Looks for a non-zero coefficient
				if isNotNegligible(p.Coeffs[key]) {
					toEncode = true

					// Resets the temporary array to zero
					// is needed if a zero coefficient
					// is at the place of a previous non-zero
					// coefficient
					if !reset {
						for j := range values {
							values[j] = 0
						}
						reset = true
					}

					// Copies the coefficient on the temporary array
					// according to the slot map index
					for _, j := range slotsIndex[i] {
						values[j] = p.Coeffs[key]
					}
				}
			}

			// If a non-zero degre coefficient was found, encode and adds the values on the output
			// ciphertext
			if toEncode {
				pt.Scale = targetScale / X[key].Scale
				polyEval.EncodeSlots(values, pt, params.LogSlots())
				polyEval.MulAndAdd(X[key], pt, res)
				toEncode = false
			}
		}

	} else {

		c := pol.Value[0].Coeffs[0]

		if minimumDegreeNonZeroCoefficient == 0 {

			res = NewCiphertext(params, 1, level, targetScale)

			if isNotNegligible(c) {
				polyEval.AddConst(res, c, res)
			}

			return
		}

		res = NewCiphertext(params, maximumCiphertextDegree, level, targetScale)

		if isNotNegligible(c) {
			polyEval.AddConst(res, c, res)
		}

		cRealFlo, cImagFlo, constScale := ring.NewFloat(0, 128), ring.NewFloat(0, 128), ring.NewFloat(0, 128)
		cRealBig, cImagBig := ring.NewUint(0), ring.NewUint(0)

		for key := pol.Value[0].Degree(); key > 0; key-- {

			c = pol.Value[0].Coeffs[key]

			if key != 0 && isNotNegligible(c) {

				cRealFlo.SetFloat64(real(c))
				cImagFlo.SetFloat64(imag(c))
				constScale.SetFloat64(targetScale / X[key].Scale)

				// Target scale * rescale-scale / power basis scale
				cRealFlo.Mul(cRealFlo, constScale)
				cImagFlo.Mul(cImagFlo, constScale)

				if cRealFlo.Sign() < 0 {
					cRealFlo.Sub(cRealFlo, new(big.Float).SetFloat64(0.5))
				} else {
					cRealFlo.Add(cRealFlo, new(big.Float).SetFloat64(0.5))
				}

				if cImagFlo.Sign() < 0 {
					cImagFlo.Sub(cImagFlo, new(big.Float).SetFloat64(0.5))
				} else {
					cImagFlo.Add(cImagFlo, new(big.Float).SetFloat64(0.5))
				}

				cRealFlo.Int(cRealBig)
				cImagFlo.Int(cImagBig)

				polyEval.MultByGaussianIntegerAndAdd(X[key], cRealBig, cImagBig, res)
			}
		}
	}

	return
}

func isNotNegligible(c complex128) bool {
	return (math.Abs(real(c)) > IsNegligbleThreshold || math.Abs(imag(c)) > IsNegligbleThreshold)
}

func isOddOrEvenPolynomial(coeffs []complex128) (odd, even bool) {
	even = true
	odd = true
	for i, c := range coeffs {
		isnotnegligible := isNotNegligible(c)
		odd = odd && !(i&1 == 0 && isnotnegligible)
		even = even && !(i&1 == 1 && isnotnegligible)
		if !odd && !even {
			break
		}
	}

	// If even or odd, then sets the expected zero coefficients to zero
	if even || odd {
		var start int
		if even {
			start = 1
		}
		for i := start; i < len(coeffs); i += 2 {
			coeffs[i] = complex(0, 0)
		}
	}

	return
}
