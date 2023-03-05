package ckks

import (
	"fmt"
	"math"
	"math/big"
	"math/bits"
	"runtime"

	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

type polynomial struct {
	*bignum.Polynomial
	Prec   uint
	MaxDeg int  // Always set to len(Coeffs)-1
	Lead   bool // Always set to true
	Lazy   bool // Flag for lazy-relinearization
}

func newPolynomial(poly *bignum.Polynomial, prec uint) (p *polynomial) {
	return &polynomial{
		Polynomial: poly,
		MaxDeg:     poly.Degree(),
		Lead:       true,
		Prec:       prec,
	}
}

type polynomialVector struct {
	Encoder    *Encoder
	Value      []*polynomial
	SlotsIndex map[int][]int
}

// checkEnoughLevels checks that enough levels are available to evaluate the polynomial.
// Also checks if c is a Gaussian integer or not. If not, then one more level is needed
// to evaluate the polynomial.
func checkEnoughLevels(levels, depth int) (err error) {

	if levels < depth {
		return fmt.Errorf("%d levels < %d log(d) -> cannot evaluate", levels, depth)
	}

	return nil
}

type polynomialEvaluator struct {
	Evaluator
	*Encoder
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
// input must be either *rlwe.Ciphertext or *PolynomialBasis.
// pol: a *Polynomial
// targetScale: the desired output scale. This value shouldn't differ too much from the original ciphertext scale. It can
// for example be used to correct small deviations in the ciphertext scale and reset it to the default scale.
func (eval *evaluator) EvaluatePoly(input interface{}, poly *bignum.Polynomial, targetScale rlwe.Scale) (opOut *rlwe.Ciphertext, err error) {
	return eval.evaluatePolyVector(input, polynomialVector{Value: []*polynomial{newPolynomial(poly, eval.params.DefaultPrecision())}}, targetScale)
}

// EvaluatePolyVector evaluates a vector of Polynomials on the input Ciphertext in ceil(log2(deg+1)) levels.
// Returns an error if the input Ciphertext does not have enough level to carry out the full polynomial evaluation.
// Returns an error if something is wrong with the scale.
// Returns an error if polynomials are not all in the same basis.
// Returns an error if polynomials do not all have the same degree.
// If the polynomials are given in Chebyshev basis, then a change of basis ct' = (2/(b-a)) * (ct + (-a-b)/(b-a))
// is necessary before the polynomial evaluation to ensure correctness.
// input: must be either *rlwe.Ciphertext or *PolynomialBasis.
// pols: a slice of up to 'n' *Polynomial ('n' being the maximum number of slots), indexed from 0 to n-1.
// encoder: an Encoder.
// slotsIndex: a map[int][]int indexing as key the polynomial to evaluate and as value the index of the slots on which to evaluate the polynomial indexed by the key.
// targetScale: the desired output scale. This value shouldn't differ too much from the original ciphertext scale. It can
// for example be used to correct small deviations in the ciphertext scale and reset it to the default scale.
//
// Example: if pols = []*Polynomial{pol0, pol1} and slotsIndex = map[int][]int:{0:[1, 2, 4, 5, 7], 1:[0, 3]},
// then pol0 will be applied to slots [1, 2, 4, 5, 7], pol1 to slots [0, 3] and the slot 6 will be zero-ed.
func (eval *evaluator) EvaluatePolyVector(input interface{}, polys []*bignum.Polynomial, encoder *Encoder, slotsIndex map[int][]int, targetScale rlwe.Scale) (opOut *rlwe.Ciphertext, err error) {
	var maxDeg int
	var basis bignum.BasisType
	for i := range polys {
		maxDeg = utils.MaxInt(maxDeg, polys[i].Degree())
		basis = polys[i].BasisType
	}

	for i := range polys {
		if basis != polys[i].BasisType {
			return nil, fmt.Errorf("polynomial basis must be the same for all polynomials in a polynomial vector")
		}

		if maxDeg != polys[i].Degree() {
			return nil, fmt.Errorf("polynomial degree must all be the same")
		}
	}

	polyvec := make([]*polynomial, len(polys))

	prec := eval.params.DefaultPrecision()
	for i := range polys {
		polyvec[i] = newPolynomial(polys[i], prec)
	}

	return eval.evaluatePolyVector(input, polynomialVector{Encoder: encoder, Value: polyvec, SlotsIndex: slotsIndex}, targetScale)
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

func (eval *evaluator) evaluatePolyVector(input interface{}, pol polynomialVector, targetScale rlwe.Scale) (opOut *rlwe.Ciphertext, err error) {

	if pol.SlotsIndex != nil && pol.Encoder == nil {
		return nil, fmt.Errorf("cannot EvaluatePolyVector: missing Encoder input")
	}

	var monomialBasis *PowerBasis
	switch input := input.(type) {
	case *rlwe.Ciphertext:
		monomialBasis = NewPowerBasis(input, pol.Value[0].Basis)
	case *PowerBasis:
		if input.Value[1] == nil {
			return nil, fmt.Errorf("cannot evaluatePolyVector: given PowerBasis.Value[1] is empty")
		}
		monomialBasis = input
	default:
		return nil, fmt.Errorf("cannot evaluatePolyVector: invalid input, must be either *rlwe.Ciphertext or *PowerBasis")
	}

	nbModuliPerRescale := eval.params.DefaultScaleModuliRatio()

	if err := checkEnoughLevels(monomialBasis.Value[1].Level(), nbModuliPerRescale*pol.Value[0].Depth()); err != nil {
		return nil, err
	}

	logDegree := bits.Len64(uint64(pol.Value[0].Degree()))
	logSplit := optimalSplit(logDegree)

	var odd, even bool = true, true
	for _, p := range pol.Value {
		odd, even = odd && p.IsOdd, even && p.IsEven
	}

	// Computes all the powers of two with relinearization
	// This will recursively compute and store all powers of two up to 2^logDegree
	if err = monomialBasis.GenPower(1<<logDegree, false, targetScale, eval); err != nil {
		return nil, err
	}

	// Computes the intermediate powers, starting from the largest, without relinearization if possible
	for i := (1 << logSplit) - 1; i > 2; i-- {
		if !(even || odd) || (i&1 == 0 && even) || (i&1 == 1 && odd) {
			if err = monomialBasis.GenPower(i, pol.Value[0].Lazy, targetScale, eval); err != nil {
				return nil, err
			}
		}
	}

	polyEval := &polynomialEvaluator{}
	polyEval.slotsIndex = pol.SlotsIndex
	polyEval.Evaluator = eval
	polyEval.Encoder = pol.Encoder
	polyEval.PowerBasis = *monomialBasis
	polyEval.logDegree = logDegree
	polyEval.logSplit = logSplit
	polyEval.isOdd = odd
	polyEval.isEven = even

	if opOut, err = polyEval.recurse(monomialBasis.Value[1].Level()-nbModuliPerRescale*(logDegree-1), targetScale, pol); err != nil {
		return nil, err
	}

	if opOut.Degree() == 2 {
		polyEval.Relinearize(opOut, opOut)
	}

	if err = polyEval.Rescale(opOut, targetScale, opOut); err != nil {
		return nil, err
	}

	opOut.Scale = targetScale

	polyEval = nil
	runtime.GC()
	return opOut, err
}

// PolynomialBasis is a struct storing powers of a ciphertext.
type PolynomialBasis struct {
	bignum.BasisType
	Value map[int]*rlwe.Ciphertext
}

// NewPolynomialBasis creates a new PolynomialBasis. It takes as input a ciphertext
// and a basistype. The struct treats the input ciphertext as a monomial X and
// can be used to generates power of this monomial X^{n} in the given BasisType.
func NewPolynomialBasis(ct *rlwe.Ciphertext, basistype bignum.BasisType) (p *PolynomialBasis) {
	p = new(PolynomialBasis)
	p.Value = make(map[int]*rlwe.Ciphertext)
	p.Value[1] = ct.CopyNew()
	p.BasisType = basistype
	return
}

// GenPower recursively computes X^{n}.
// If lazy = true, the final X^{n} will not be relinearized.
// Previous non-relinearized X^{n} that are required to compute the target X^{n} are automatically relinearized.
// Scale sets the threshold for rescaling (ciphertext won't be rescaled if the rescaling operation would make the scale go under this threshold).
func (p *PolynomialBasis) GenPower(n int, lazy bool, scale rlwe.Scale, eval Evaluator) (err error) {

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

func (p *PolynomialBasis) genPower(n int, lazy bool, scale rlwe.Scale, eval Evaluator) (err error) {

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

			if p.BasisType == bignum.Chebyshev {
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

			if err = eval.Rescale(p.Value[a], scale, p.Value[a]); err != nil {
				return err
			}

			if err = eval.Rescale(p.Value[b], scale, p.Value[b]); err != nil {
				return err
			}

			p.Value[n] = eval.MulRelinNew(p.Value[a], p.Value[b])
		}

		if p.BasisType == bignum.Chebyshev {

			// Computes C[n] = 2*C[a]*C[b]
			eval.Add(p.Value[n], p.Value[n], p.Value[n])

			// Computes C[n] = 2*C[a]*C[b] - C[c]
			if c == 0 {
				eval.Add(p.Value[n], -1, p.Value[n])
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
	binary.LittleEndian.PutUint64(data[8:16], uint64(p.Value[1].MarshalBinarySize()))
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
	p.Value = make(map[int]*rlwe.Ciphertext)
	nbct := int(binary.LittleEndian.Uint64(data[0:8]))
	dtLen := int(binary.LittleEndian.Uint64(data[8:16]))
	ptr := 16
	for i := 0; i < nbct; i++ {
		idx := int(binary.LittleEndian.Uint64(data[ptr : ptr+8]))
		ptr += 8
		p.Value[idx] = new(rlwe.Ciphertext)
		if err = p.Value[idx].UnmarshalBinary(data[ptr : ptr+dtLen]); err != nil {
			return
		}
		ptr += dtLen
	}
	return
}

// splitCoeffs splits coeffs as X^{2n} * coeffsq + coeffsr.
// This function is sensitive to the precision of the coefficients.
func splitCoeffs(coeffs *polynomial, split int) (coeffsq, coeffsr *polynomial) {

	prec := coeffs.Prec

	// Splits a polynomial p such that p = q*C^degree + r.
	coeffsr = &polynomial{Polynomial: &bignum.Polynomial{}}
	coeffsr.Coeffs = make([]*bignum.Complex, split)
	if coeffs.MaxDeg == coeffs.Degree() {
		coeffsr.MaxDeg = split - 1
	} else {
		coeffsr.MaxDeg = coeffs.MaxDeg - (coeffs.Degree() - split + 1)
	}

	for i := 0; i < split; i++ {
		if coeffs.Coeffs[i] != nil {
			coeffsr.Coeffs[i] = coeffs.Coeffs[i].Copy()
			coeffsr.Coeffs[i].SetPrec(prec)
		}

	}

	coeffsq = &polynomial{Polynomial: &bignum.Polynomial{}}
	coeffsq.Coeffs = make([]*bignum.Complex, coeffs.Degree()-split+1)
	coeffsq.MaxDeg = coeffs.MaxDeg

	if coeffs.Coeffs[split] != nil {
		coeffsq.Coeffs[0] = coeffs.Coeffs[split].Copy()
	}

	odd := coeffs.IsOdd
	even := coeffs.IsEven

	switch coeffs.BasisType {
	case bignum.Monomial:
		for i := split + 1; i < coeffs.Degree()+1; i++ {
			if coeffs.Coeffs[i] != nil && (!(even || odd) || (i&1 == 0 && even) || (i&1 == 1 && odd)) {
				coeffsq.Coeffs[i-split] = coeffs.Coeffs[i].Copy()
				coeffsr.Coeffs[i-split].SetPrec(prec)
			}
		}
	case bignum.Chebyshev:

		for i, j := split+1, 1; i < coeffs.Degree()+1; i, j = i+1, j+1 {
			if coeffs.Coeffs[i] != nil && (!(even || odd) || (i&1 == 0 && even) || (i&1 == 1 && odd)) {
				coeffsq.Coeffs[i-split] = coeffs.Coeffs[i].Copy()
				coeffsr.Coeffs[i-split].SetPrec(prec)
				coeffsq.Coeffs[i-split].Add(coeffsq.Coeffs[i-split], coeffsq.Coeffs[i-split])

				if coeffsr.Coeffs[split-j] != nil {
					coeffsr.Coeffs[split-j].Sub(coeffsr.Coeffs[split-j], coeffs.Coeffs[i])
				} else {
					coeffsr.Coeffs[split-j] = coeffs.Coeffs[i].Copy()
					coeffsr.Coeffs[split-j].SetPrec(prec)
					coeffsr.Coeffs[split-j][0].Neg(coeffsr.Coeffs[split-j][0])
					coeffsr.Coeffs[split-j][1].Neg(coeffsr.Coeffs[split-j][1])
				}
			}
		}
	}

	if coeffs.Lead {
		coeffsq.Lead = true
	}

	coeffsq.BasisType, coeffsr.BasisType = coeffs.BasisType, coeffs.BasisType
	coeffsq.IsOdd, coeffsr.IsOdd = coeffs.IsOdd, coeffs.IsOdd
	coeffsq.IsEven, coeffsr.IsEven = coeffs.IsEven, coeffs.IsEven
	coeffsq.Prec, coeffsr.Prec = prec, prec

	return
}

func splitCoeffsPolyVector(poly polynomialVector, split int) (polyq, polyr polynomialVector) {
	coeffsq := make([]*polynomial, len(poly.Value))
	coeffsr := make([]*polynomial, len(poly.Value))
	for i, p := range poly.Value {
		coeffsq[i], coeffsr[i] = splitCoeffs(p, split)
	}

	return polynomialVector{Value: coeffsq}, polynomialVector{Value: coeffsr}
}

func (polyEval *polynomialEvaluator) recurse(targetLevel int, targetScale rlwe.Scale, pol polynomialVector) (res *rlwe.Ciphertext, err error) {

	params := polyEval.Evaluator.(*evaluator).params

	logSplit := polyEval.logSplit

	nbModuliPerRescale := params.DefaultScaleModuliRatio()

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
			polyEvalBis.PowerBasis = polyEval.PowerBasis
			polyEvalBis.isOdd = polyEval.isOdd
			polyEvalBis.isEven = polyEval.isEven

			return polyEvalBis.recurse(targetLevel, targetScale, pol)
		}

		if pol.Value[0].Lead {

			targetScale = targetScale.Mul(rlwe.NewScale(params.Q()[targetLevel]))

			for i := 1; i < nbModuliPerRescale; i++ {
				targetScale = targetScale.Mul(rlwe.NewScale(params.Q()[targetLevel-i]))
			}
		}

		return polyEval.evaluatePolyFromPowerBasis(targetScale, targetLevel, pol)
	}

	var nextPower = 1 << polyEval.logSplit
	for nextPower < (pol.Value[0].Degree()>>1)+1 {
		nextPower <<= 1
	}

	coeffsq, coeffsr := splitCoeffsPolyVector(pol, nextPower)

	XPow := polyEval.PowerBasis.Value[nextPower]

	level := targetLevel

	var qi *big.Int
	if pol.Value[0].Lead {
		qi = bignum.NewInt(params.Q()[level])
		for i := 1; i < nbModuliPerRescale; i++ {
			qi.Mul(qi, bignum.NewInt(params.Q()[level-i]))
		}
	} else {
		qi = bignum.NewInt(params.Q()[level+nbModuliPerRescale])
		for i := 1; i < nbModuliPerRescale; i++ {
			qi.Mul(qi, bignum.NewInt(params.Q()[level+nbModuliPerRescale-i]))
		}
	}

	targetScale = targetScale.Mul(rlwe.NewScale(qi))
	targetScale = targetScale.Div(XPow.Scale)

	if res, err = polyEval.recurse(targetLevel+nbModuliPerRescale, targetScale, coeffsq); err != nil {
		return nil, err
	}
	if res.Degree() == 2 {
		polyEval.Relinearize(res, res)
	}

	if err = polyEval.Rescale(res, params.DefaultScale(), res); err != nil {
		return nil, err
	}

	polyEval.Mul(res, XPow, res)

	var tmp *rlwe.Ciphertext
	if tmp, err = polyEval.recurse(res.Level(), res.Scale, coeffsr); err != nil {
		return nil, err
	}

	polyEval.Add(res, tmp, res)

	tmp = nil

	return
}

func (polyEval *polynomialEvaluator) evaluatePolyFromPowerBasis(targetScale rlwe.Scale, level int, pol polynomialVector) (res *rlwe.Ciphertext, err error) {

	// Map[int] of the powers [X^{0}, X^{1}, X^{2}, ...]
	X := polyEval.PolynomialBasis.Value

	// Retrieve the number of slots
	logSlots := X[1].LogSlots
	slots := 1 << X[1].LogSlots

	params := polyEval.Evaluator.(*evaluator).params
	slotsIndex := polyEval.slotsIndex

	// Retrieve the degree of the highest degree non-zero coefficient
	// TODO: optimize for nil/zero coefficients
	minimumDegreeNonZeroCoefficient := len(pol.Value[0].Coeffs) - 1
	if polyEval.isEven && !polyEval.isOdd {
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

	// Retrieve flags for even/odd
	even := polyEval.isEven
	odd := polyEval.isOdd

	// If an index slot is given (either multiply polynomials or masking)
	if slotsIndex != nil {

		var toEncode bool

		// Allocates temporary buffer for coefficients encoding
		values := make([]*bignum.Complex, slots)

		// If the degree of the poly is zero
		if minimumDegreeNonZeroCoefficient == 0 {

			// Allocates the output ciphertext
			res = NewCiphertext(params, 1, level)
			res.Scale = targetScale
			res.LogSlots = logSlots

			// Looks for non-zero coefficients among the degree 0 coefficients of the polynomials
			for i, p := range pol.Value {
				if !isZero(p.Coeffs[0]) {
					toEncode = true
					for _, j := range slotsIndex[i] {
						values[j] = p.Coeffs[0]
					}
				}
			}

			// If a non-zero coefficient was found, encode the values, adds on the ciphertext, and returns
			if toEncode {
				pt := rlwe.NewPlaintextAtLevelFromPoly(level, res.Value[0])
				pt.LogSlots = logSlots
				pt.IsNTT = true
				pt.Scale = targetScale
				polyEval.Encode(values, pt)
			}

			return
		}

		// Allocates the output ciphertext
		res = NewCiphertext(params, maximumCiphertextDegree, level)
		res.Scale = targetScale
		res.LogSlots = logSlots

		// Allocates a temporary plaintext to encode the values
		pt := rlwe.NewPlaintextAtLevelFromPoly(level, polyEval.Evaluator.BuffCt().Value[0])
		pt.IsNTT = true
		pt.LogSlots = logSlots

		// Looks for a non-zero coefficient among the degree zero coefficient of the polynomials
		for i, p := range pol.Value {
			if !isZero(p.Coeffs[0]) {
				toEncode = true
				for _, j := range slotsIndex[i] {
					values[j] = p.Coeffs[0]
				}
			}
		}

		// If a non-zero degre coefficient was found, encode and adds the values on the output
		// ciphertext
		if toEncode {
			pt.Scale = targetScale
			polyEval.Encode(values, pt)
			polyEval.Add(res, pt, res)
			toEncode = false
		}

		// Loops starting from the highest degree coefficient
		for key := pol.Value[0].Degree(); key > 0; key-- {

			var reset bool
			// Loops over the polynomials
			for i, p := range pol.Value {

				// Looks for a non-zero coefficient
				if !isZero(p.Coeffs[key]) {
					toEncode = true

					// Resets the temporary array to zero
					// is needed if a zero coefficient
					// is at the place of a previous non-zero
					// coefficient
					if !reset {
						for j := range values {
							if values[j] != nil {
								values[j][0].SetFloat64(0)
								values[j][1].SetFloat64(0)
							}
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
				pt.Scale = targetScale.Div(X[key].Scale)
				polyEval.Encode(values, pt)
				polyEval.MulThenAdd(X[key], pt, res)
				toEncode = false
			}
		}

	} else {

		var c *bignum.Complex
		if polyEval.isEven && !isZero(pol.Value[0].Coeffs[0]) {
			c = pol.Value[0].Coeffs[0]
		}

		if minimumDegreeNonZeroCoefficient == 0 {

			res = NewCiphertext(params, 1, level)
			res.Scale = targetScale
			res.LogSlots = logSlots

			if !isZero(c) {
				polyEval.Add(res, c, res)
			}

			return
		}

		res = NewCiphertext(params, maximumCiphertextDegree, level)
		res.Scale = targetScale
		res.LogSlots = logSlots

		if c != nil {
			polyEval.Add(res, c, res)
		}

		constScale := new(big.Float).SetPrec(pol.Value[0].Prec)

		ringQ := params.RingQ().AtLevel(level)

		for key := pol.Value[0].Degree(); key > 0; key-- {

			if c = pol.Value[0].Coeffs[key]; key != 0 && !isZero(c) && (!(even || odd) || (key&1 == 0 && even) || (key&1 == 1 && odd)) {

				XScale := X[key].Scale.Value
				tgScale := targetScale.Value
				constScale.Quo(&tgScale, &XScale)

				RNSReal, RNSImag := bigComplexToRNSScalar(ringQ, constScale, bignum.ToComplex(c, pol.Value[0].Prec))

				polyEval.Evaluator.(*evaluator).evaluateWithScalar(level, X[key].Value, RNSReal, RNSImag, res.Value, ringQ.MulDoubleRNSScalarThenAdd)
			}
		}
	}

	return
}

func isZero(c *bignum.Complex) bool {
	zero := new(big.Float)
	return c == nil || (c[0].Cmp(zero) == 0 && c[1].Cmp(zero) == 0)
}
