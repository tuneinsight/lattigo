package bgv

import (
	"encoding/binary"
	"fmt"
	"math"
	"math/bits"
	"runtime"

	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// Polynomial is a struct storing the coefficients of a plaintext
// polynomial that then can be evaluated on the ciphertext.
type Polynomial struct {
	MaxDeg int
	Coeffs []uint64
	Lead   bool
}

// Depth returns the depth needed to evaluate the polynomial.
func (p *Polynomial) Depth() int {
	return int(math.Ceil(math.Log2(float64(len(p.Coeffs)))))
}

// Degree returns the degree of the polynomial.
func (p *Polynomial) Degree() int {
	return len(p.Coeffs) - 1
}

// NewPoly creates a new Poly from the input coefficients.
func NewPoly(coeffs []uint64) (p *Polynomial) {
	c := make([]uint64, len(coeffs))
	copy(c, coeffs)
	return &Polynomial{Coeffs: c, MaxDeg: len(c) - 1, Lead: true}
}

type polynomialEvaluator struct {
	Evaluator
	Encoder
	slotsIndex map[int][]int
	powerBasis map[int]*rlwe.Ciphertext
	logDegree  int
	logSplit   int
	isOdd      bool
	isEven     bool
}

// EvaluatePoly evaluates a Polynomial in standard basis on the input Ciphertext in ceil(log2(deg+1)) depth.
// input must be either *rlwe.Ciphertext or *PowerBasis.
func (eval *evaluator) EvaluatePoly(input interface{}, pol *Polynomial, targetScale rlwe.Scale) (opOut *rlwe.Ciphertext, err error) {
	return eval.evaluatePolyVector(input, polynomialVector{Value: []*Polynomial{pol}}, targetScale)
}

type polynomialVector struct {
	Encoder    Encoder
	Value      []*Polynomial
	SlotsIndex map[int][]int
}

// EvaluatePolyVector evaluates a vector of Polynomials on the input Ciphertext in ceil(log2(deg+1)) depth.
// Inputs:
// input: *rlwe.Ciphertext or *PowerBasis.
// pols: a slice of up to 'n' *Polynomial ('n' being the maximum number of slots), indexed from 0 to n-1. Returns an error if the polynomials do not all have the same degree.
// encoder: an Encoder.
// slotsIndex: a map[int][]int indexing as key the polynomial to evaluate and as value the index of the slots on which to evaluate the Polynomial indexed by the key.
//
// Example: if pols = []*Polynomial{pol0, pol1} and slotsIndex = map[int][]int:{0:[1, 2, 4, 5, 7], 1:[0, 3]},
// then pol0 will be applied to slots [1, 2, 4, 5, 7], pol1 to slots [0, 3] and the slot 6 will be zero-ed.
func (eval *evaluator) EvaluatePolyVector(input interface{}, pols []*Polynomial, encoder Encoder, slotsIndex map[int][]int, targetScale rlwe.Scale) (opOut *rlwe.Ciphertext, err error) {
	var maxDeg int
	for i := range pols {
		maxDeg = utils.MaxInt(maxDeg, pols[i].MaxDeg)
	}

	for i := range pols {
		if maxDeg != pols[i].MaxDeg {
			return nil, fmt.Errorf("cannot EvaluatePolyVector: polynomial degree must all be the same")
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

func (eval *evaluator) evaluatePolyVector(input interface{}, pol polynomialVector, targetScale rlwe.Scale) (opOut *rlwe.Ciphertext, err error) {

	if pol.SlotsIndex != nil && pol.Encoder == nil {
		return nil, fmt.Errorf("cannot evaluatePolyVector: missing Encoder input")
	}

	var powerBasis *PowerBasis
	switch input := input.(type) {
	case *rlwe.Ciphertext:

		if level, depth := input.Level(), pol.Value[0].Depth(); level < depth {
			return nil, fmt.Errorf("%d levels < %d log(d) -> cannot evaluate poly", level, depth)
		}

		powerBasis = NewPowerBasis(input)

	case *PowerBasis:
		if input.Value[1] == nil {
			return nil, fmt.Errorf("cannot evaluatePolyVector: given PowerBasis[1] is empty")
		}
		powerBasis = input
	default:
		return nil, fmt.Errorf("cannot evaluatePolyVector: invalid input, must be either *rlwe.Ciphertext or *PowerBasis")
	}

	logDegree := bits.Len64(uint64(pol.Value[0].Degree()))
	logSplit := optimalSplit(logDegree)

	var odd, even = true, true
	for _, p := range pol.Value {
		tmp0, tmp1 := isOddOrEvenPolynomial(p.Coeffs)
		odd, even = odd && tmp0, even && tmp1
	}

	for i := (1 << logSplit) - 1; i > 1; i-- {
		if !(even || odd) || (i&1 == 0 && even) || (i&1 == 1 && odd) {
			if err = powerBasis.GenPower(i, true, eval); err != nil {
				return nil, err
			}
		}
	}

	for i := logSplit; i < logDegree; i++ {
		if err = powerBasis.GenPower(1<<i, false, eval); err != nil {
			return nil, err
		}
	}

	polyEval := &polynomialEvaluator{}
	polyEval.slotsIndex = pol.SlotsIndex
	polyEval.Evaluator = eval
	polyEval.Encoder = pol.Encoder
	polyEval.powerBasis = powerBasis.Value
	polyEval.logDegree = logDegree
	polyEval.logSplit = logSplit
	polyEval.isOdd = odd
	polyEval.isEven = even

	if opOut, err = polyEval.recurse(powerBasis.Value[1].Level()-logDegree+1, targetScale, pol); err != nil {
		return
	}

	polyEval.Relinearize(opOut, opOut)

	if err = polyEval.Rescale(opOut, opOut); err != nil {
		return nil, err
	}

	polyEval = nil
	runtime.GC()
	return opOut, err
}

// PowerBasis is a struct storing powers of a ciphertext.
type PowerBasis struct {
	Value map[int]*rlwe.Ciphertext
}

// NewPowerBasis creates a new PowerBasis.
func NewPowerBasis(ct *rlwe.Ciphertext) (p *PowerBasis) {
	p = new(PowerBasis)
	p.Value = make(map[int]*rlwe.Ciphertext)
	p.Value[1] = ct.CopyNew()
	return
}

// GenPower generates the n-th power of the power basis,
// as well as all the necessary intermediate powers if
// they are not yet present.
func (p *PowerBasis) GenPower(n int, lazy bool, eval Evaluator) (err error) {

	var rescale bool
	if rescale, err = p.genPower(n, n, lazy, true, eval); err != nil {
		return
	}

	if rescale {
		if err = eval.Rescale(p.Value[n], p.Value[n]); err != nil {
			return
		}
	}

	return nil
}

func (p *PowerBasis) genPower(target, n int, lazy, rescale bool, eval Evaluator) (rescaleN bool, err error) {

	if p.Value[n] == nil {

		isPow2 := n&(n-1) == 0

		// Computes the index required to compute the required ring evaluation
		var a, b int
		if isPow2 {
			a, b = n/2, n/2 // Necessary for optimal depth
		} else {
			// Maximize the number of odd terms
			k := int(math.Ceil(math.Log2(float64(n)))) - 1
			a = (1 << k) - 1
			b = n + 1 - (1 << k)
		}

		var rescaleA, rescaleB bool

		// Recurses on the given indexes
		if rescaleA, err = p.genPower(target, a, lazy, rescale, eval); err != nil {
			return false, err
		}

		if rescaleB, err = p.genPower(target, b, lazy, rescale, eval); err != nil {
			return false, err
		}

		if p.Value[a].Degree() == 2 {
			eval.Relinearize(p.Value[a], p.Value[a])
		}

		if p.Value[b].Degree() == 2 {
			eval.Relinearize(p.Value[b], p.Value[b])
		}

		if rescaleA {
			if err = eval.Rescale(p.Value[a], p.Value[a]); err != nil {
				return false, err
			}
		}

		if rescaleB {
			if err = eval.Rescale(p.Value[b], p.Value[b]); err != nil {
				return false, err
			}
		}

		// Computes C[n] = C[a]*C[b]
		if lazy && !isPow2 {
			p.Value[n] = eval.MulNew(p.Value[a], p.Value[b])
			return true, nil
		}

		p.Value[n] = eval.MulRelinNew(p.Value[a], p.Value[b])
		if err = eval.Rescale(p.Value[n], p.Value[n]); err != nil {
			return false, err
		}

	}

	return false, nil
}

// MarshalBinary encodes the target on a slice of bytes.
func (p *PowerBasis) MarshalBinary() (data []byte, err error) {
	data = make([]byte, 8)
	binary.LittleEndian.PutUint64(data, uint64(len(p.Value)))
	for key, ct := range p.Value {

		header := make([]byte, 16)
		binary.LittleEndian.PutUint64(header[0:], uint64(key))
		binary.LittleEndian.PutUint64(header[8:], uint64(ct.BinarySize()))

		data = append(data, header...)
		ctBytes, err := ct.MarshalBinary()
		if err != nil {
			return []byte{}, err
		}
		data = append(data, ctBytes...)
	}
	return
}

// UnmarshalBinary decodes a slice of bytes on the target.
func (p *PowerBasis) UnmarshalBinary(data []byte) (err error) {
	p.Value = make(map[int]*rlwe.Ciphertext)
	nbct := int(binary.LittleEndian.Uint64(data))
	ptr := 8
	for i := 0; i < nbct; i++ {
		idx := int(binary.LittleEndian.Uint64(data[ptr : ptr+8]))
		ptr += 8
		dtLen := int(binary.LittleEndian.Uint64(data[ptr : ptr+8]))
		ptr += 8
		p.Value[idx] = &rlwe.Ciphertext{}
		if err = p.Value[idx].UnmarshalBinary(data[ptr : ptr+dtLen]); err != nil {
			fmt.Println(123)
			return
		}
		ptr += dtLen
	}
	return
}

// splitCoeffs splits a polynomial p such that p = q*C^degree + r.
func splitCoeffs(coeffs *Polynomial, split int) (coeffsq, coeffsr *Polynomial) {

	coeffsr = &Polynomial{}
	coeffsr.Coeffs = make([]uint64, split)
	if coeffs.MaxDeg == coeffs.Degree() {
		coeffsr.MaxDeg = split - 1
	} else {
		coeffsr.MaxDeg = coeffs.MaxDeg - (coeffs.Degree() - split + 1)
	}

	for i := 0; i < split; i++ {
		coeffsr.Coeffs[i] = coeffs.Coeffs[i]
	}

	coeffsq = &Polynomial{}
	coeffsq.Coeffs = make([]uint64, coeffs.Degree()-split+1)
	coeffsq.MaxDeg = coeffs.MaxDeg

	coeffsq.Coeffs[0] = coeffs.Coeffs[split]

	for i := split + 1; i < coeffs.Degree()+1; i++ {
		coeffsq.Coeffs[i-split] = coeffs.Coeffs[i]
	}

	if coeffs.Lead {
		coeffsq.Lead = true
	}

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

func (polyEval *polynomialEvaluator) recurse(targetLevel int, targetScale rlwe.Scale, pol polynomialVector) (res *rlwe.Ciphertext, err error) {

	logSplit := polyEval.logSplit

	params := polyEval.Evaluator.(*evaluator).params

	// Recursively computes the evaluation of the Chebyshev polynomial using a baby-set giant-step algorithm.
	if pol.Value[0].Degree() < (1 << logSplit) {

		if pol.Value[0].Lead && polyEval.logSplit > 1 && pol.Value[0].MaxDeg%(1<<(logSplit+1)) > (1<<(logSplit-1)) {

			logDegree := int(bits.Len64(uint64(pol.Value[0].Degree())))
			logSplit := logDegree >> 1

			polyEvalBis := new(polynomialEvaluator)
			polyEvalBis.Evaluator = polyEval.Evaluator
			polyEvalBis.Encoder = polyEval.Encoder
			polyEvalBis.logDegree = logDegree
			polyEvalBis.logSplit = logSplit
			polyEvalBis.slotsIndex = polyEval.slotsIndex
			polyEvalBis.powerBasis = polyEval.powerBasis
			polyEvalBis.isOdd = polyEval.isOdd
			polyEvalBis.isEven = polyEval.isEven

			res, err = polyEvalBis.recurse(targetLevel, targetScale, pol)

			return
		}

		if pol.Value[0].Lead {
			targetScale = targetScale.Mul(params.NewScale(params.Q()[targetLevel]))
		}

		res, err = polyEval.evaluatePolyFromPowerBasis(targetLevel, targetScale, pol)

		return
	}

	var nextPower = 1 << polyEval.logSplit
	for nextPower < (pol.Value[0].Degree()>>1)+1 {
		nextPower <<= 1
	}

	coeffsq, coeffsr := splitCoeffsPolyVector(pol, nextPower)

	XPow := polyEval.powerBasis[nextPower]

	level := targetLevel

	var currentQi uint64
	if pol.Value[0].Lead {
		currentQi = params.Q()[level]
	} else {
		currentQi = params.Q()[level+1]
	}

	// targetScale = targetScale*currentQi/XPow.Scale
	targetScale = targetScale.Mul(params.NewScale(currentQi))
	targetScale = targetScale.Div(XPow.Scale)

	if res, err = polyEval.recurse(targetLevel+1, targetScale, coeffsq); err != nil {
		return nil, err
	}

	if res.Degree() == 2 {
		polyEval.Relinearize(res, res)
	}

	if err = polyEval.Rescale(res, res); err != nil {
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

func (polyEval *polynomialEvaluator) evaluatePolyFromPowerBasis(targetLevel int, targetScale rlwe.Scale, pol polynomialVector) (res *rlwe.Ciphertext, err error) {

	X := polyEval.powerBasis

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
		values := make([]uint64, params.N())

		// If the degree of the poly is zero
		if minimumDegreeNonZeroCoefficient == 0 {

			// Allocates the output ciphertext
			res = rlwe.NewCiphertext(params.Parameters, 1, targetLevel)
			res.Scale = targetScale

			// Looks for non-zero coefficients among the degree 0 coefficients of the polynomials
			for i, p := range pol.Value {
				if p.Coeffs[0] != 0 {
					toEncode = true
					for _, j := range slotsIndex[i] {
						values[j] = p.Coeffs[0]
					}
				}
			}

			// If a non-zero coefficient was found, encode the values, adds on the ciphertext, and returns
			if toEncode {
				pt := rlwe.NewPlaintextAtLevelFromPoly(targetLevel, res.Value[0])
				pt.Scale = res.Scale
				pt.IsNTT = true
				polyEval.Encode(values, pt)
			}

			return
		}

		// Allocates the output ciphertext
		res = rlwe.NewCiphertext(params.Parameters, maximumCiphertextDegree, targetLevel)
		res.Scale = targetScale

		// Allocates a temporary plaintext to encode the values
		pt := rlwe.NewPlaintextAtLevelFromPoly(targetLevel, polyEval.Evaluator.BuffQ()[0]) // buffQ[0] is safe in this case
		pt.Scale = targetScale
		pt.IsNTT = true

		// Looks for a non-zero coefficient among the degree zero coefficient of the polynomials
		for i, p := range pol.Value {
			if p.Coeffs[0] != 0 {
				toEncode = true
				for _, j := range slotsIndex[i] {
					values[j] = p.Coeffs[0]
				}
			}
		}

		// If a non-zero degree coefficient was found, encode and adds the values on the output
		// ciphertext
		if toEncode {
			// Add would actually scale the plaintext accordingly,
			// but encoding with the correct scale is slightly faster
			pt.Scale = res.Scale
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
				if p.Coeffs[key] != 0 {
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

			// If a non-zero degree coefficient was found, encode and adds the values on the output
			// ciphertext
			if toEncode {

				// MulThenAdd would actually scale the plaintext accordingly,
				// but encoding with the correct scale is slightly faster
				pt.Scale = targetScale.Div(X[key].Scale)
				polyEval.Encode(values, pt)
				polyEval.MulThenAdd(X[key], pt, res)
				toEncode = false
			}
		}

	} else {

		c := pol.Value[0].Coeffs[0]

		if minimumDegreeNonZeroCoefficient == 0 {

			res = rlwe.NewCiphertext(params.Parameters, 1, targetLevel)
			res.Scale = targetScale

			if c != 0 {
				polyEval.AddScalar(res, c, res)
			}

			return
		}

		res = rlwe.NewCiphertext(params.Parameters, maximumCiphertextDegree, targetLevel)
		res.Scale = targetScale

		if c != 0 {
			polyEval.AddScalar(res, c, res)
		}

		for key := pol.Value[0].Degree(); key > 0; key-- {
			c = pol.Value[0].Coeffs[key]
			if key != 0 && c != 0 {
				// MulScalarAndAdd automatically scales c to match the scale of res.
				polyEval.MulScalarThenAdd(X[key], c, res)
			}
		}
	}

	return
}

func isOddOrEvenPolynomial(coeffs []uint64) (odd, even bool) {
	even = true
	odd = true
	for i, c := range coeffs {
		isnotzero := c != 0
		odd = odd && !(i&1 == 0 && isnotzero)
		even = even && !(i&1 == 1 && isnotzero)
		if !odd && !even {
			break
		}
	}

	return
}
