package bfv

import (
	"encoding/binary"
	"fmt"
	"math"
	"math/bits"
	"runtime"

	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
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
	powerBasis map[int]*Ciphertext
	logDegree  int
	logSplit   int
}

// EvaluatePoly evaluates a polynomial in standard basis on the input Ciphertext in ceil(log2(deg+1)) depth.
// input must be either *Ciphertext or *Powerbasis.
func (eval *evaluator) EvaluatePoly(input interface{}, pol *Polynomial) (opOut *Ciphertext, err error) {
	return eval.evaluatePolyVector(input, polynomialVector{Value: []*Polynomial{pol}})
}

type polynomialVector struct {
	Encoder    Encoder
	Value      []*Polynomial
	SlotsIndex map[int][]int
}

// EvaluatePolyVector evaluates a vector of Polynomials on the input Ciphertext in ceil(log2(deg+1)) depth.
// Inputs:
// input: *Ciphertext or *PowerBasis.
// pols: a slice of up to 'n' *Polynomial ('n' being the maximum number of slots), indexed from 0 to n-1. Returns an error if the polynomials do not all have the same degree.
// encoder: an Encoder.
// slotsIndex: a map[int][]int indexing as key the polynomial to evalute and as value the index of the slots on which to evaluate the polynomial indexed by the key.
//
// Example: if pols = []*Polynomial{pol0, pol1} and slotsIndex = map[int][]int:{0:[1, 2, 4, 5, 7], 1:[0, 3]},
// then pol0 will be applied to slots [1, 2, 4, 5, 7], pol1 to slots [0, 3] and the slot 6 will be zero-ed.
func (eval *evaluator) EvaluatePolyVector(input interface{}, pols []*Polynomial, encoder Encoder, slotsIndex map[int][]int) (opOut *Ciphertext, err error) {
	var maxDeg int
	for i := range pols {
		maxDeg = utils.MaxInt(maxDeg, pols[i].MaxDeg)
	}

	for i := range pols {
		if maxDeg != pols[i].MaxDeg {
			return nil, fmt.Errorf("cannot EvaluatePolyVector: polynomial degree must all be the same")
		}
	}

	return eval.evaluatePolyVector(input, polynomialVector{Encoder: encoder, Value: pols, SlotsIndex: slotsIndex})
}

func (eval *evaluator) evaluatePolyVector(input interface{}, pol polynomialVector) (opOut *Ciphertext, err error) {

	if pol.SlotsIndex != nil && pol.Encoder == nil {
		return nil, fmt.Errorf("cannot evaluatePolyVector: missing Encoder input")
	}

	var powerBasis *PowerBasis
	switch input := input.(type) {
	case *Ciphertext:
		powerBasis = NewPowerBasis(input)
	case *PowerBasis:
		if input.Value[1] == nil {
			return nil, fmt.Errorf("cannot evaluatePolyVector: given PowerBasis[1] is empty")
		}
		powerBasis = input
	default:
		return nil, fmt.Errorf("cannot evaluatePolyVector: invalid input, must be either *Ciphertext or *PowerBasis")
	}

	logDegree := bits.Len64(uint64(pol.Value[0].Degree()))
	logSplit := (logDegree >> 1) //optimalSplit(logDegree)

	var odd, even bool
	for _, p := range pol.Value {
		tmp0, tmp1 := isOddOrEvenPolynomial(p.Coeffs)
		odd, even = odd && tmp0, even && tmp1
	}

	for i := 2; i < (1 << logSplit); i++ {
		if !(even || odd) || (i&1 == 0 && even) || (i&1 == 1 && odd) {
			powerBasis.GenPower(i, eval)
		}
	}

	for i := logSplit; i < logDegree; i++ {
		powerBasis.GenPower(1<<i, eval)
	}

	polyEval := &polynomialEvaluator{}
	polyEval.slotsIndex = pol.SlotsIndex
	polyEval.Evaluator = eval
	polyEval.Encoder = pol.Encoder
	polyEval.powerBasis = powerBasis.Value
	polyEval.logDegree = logDegree
	polyEval.logSplit = logSplit

	opOut, err = polyEval.recurse(pol)

	polyEval = nil
	runtime.GC()
	return opOut, err
}

// PowerBasis is a struct storing powers of a ciphertext.
type PowerBasis struct {
	Value map[int]*Ciphertext
}

// NewPowerBasis creates a new PowerBasis.
func NewPowerBasis(ct *Ciphertext) (p *PowerBasis) {
	p = new(PowerBasis)
	p.Value = make(map[int]*Ciphertext)
	p.Value[1] = ct.CopyNew()
	return
}

// GenPower generates the n-th power of the power basis,
// as well as all the necessary intermediate powers if
// they are not yet present.
func (p *PowerBasis) GenPower(n int, eval Evaluator) {

	if p.Value[n] == nil {

		// Computes the index required to compute the required ring evaluation
		var a, b int
		if n&(n-1) == 0 {
			a, b = n/2, n/2 // Necessary for optimal depth
		} else {
			// Maximize the number of odd terms
			k := int(math.Ceil(math.Log2(float64(n)))) - 1
			a = (1 << k) - 1
			b = n + 1 - (1 << k)
		}

		// Recurses on the given indexes
		p.GenPower(a, eval)
		p.GenPower(b, eval)

		// Computes C[n] = C[a]*C[b]
		p.Value[n] = eval.MulNew(p.Value[a], p.Value[b])
		eval.Relinearize(p.Value[n], p.Value[n])
	}
}

// MarshalBinary encodes the target on a slice of bytes.
func (p *PowerBasis) MarshalBinary() (data []byte, err error) {
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
func (p *PowerBasis) UnmarshalBinary(data []byte) (err error) {
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

func (polyEval *polynomialEvaluator) recurse(pol polynomialVector) (res *Ciphertext, err error) {

	logSplit := polyEval.logSplit

	// Recursively computes the evaluation of the Chebyshev polynomial using a baby-step giant-step algorithm.
	if pol.Value[0].Degree() < (1 << logSplit) {

		if pol.Value[0].Lead && polyEval.logSplit > 1 && pol.Value[0].MaxDeg%(1<<(logSplit+1)) > (1<<(logSplit-1)) {

			logDegree := int(bits.Len64(uint64(pol.Value[0].Degree())))
			logSplit := logDegree >> 1

			polyEvalBis := new(polynomialEvaluator)
			polyEvalBis.Evaluator = polyEval.Evaluator
			polyEvalBis.slotsIndex = polyEval.slotsIndex
			polyEvalBis.Encoder = polyEval.Encoder
			polyEvalBis.logDegree = logDegree
			polyEvalBis.logSplit = logSplit
			polyEvalBis.powerBasis = polyEval.powerBasis

			return polyEvalBis.recurse(pol)
		}

		return polyEval.evaluatePolyFromPowerBasis(pol)
	}

	var nextPower = 1 << polyEval.logSplit
	for nextPower < (pol.Value[0].Degree()>>1)+1 {
		nextPower <<= 1
	}

	coeffsq, coeffsr := splitCoeffsPolyVector(pol, nextPower)

	XPow := polyEval.powerBasis[nextPower]

	if res, err = polyEval.recurse(coeffsq); err != nil {
		return nil, err
	}

	var tmp *Ciphertext
	if tmp, err = polyEval.recurse(coeffsr); err != nil {
		return nil, err
	}

	res2 := NewCiphertextLvl(polyEval.Evaluator.(*evaluator).params, 2, res.Level())
	polyEval.Mul(res, XPow, res2)
	polyEval.Relinearize(res2, res)
	polyEval.Add(res, tmp, res)

	tmp = nil

	return
}

func (polyEval *polynomialEvaluator) evaluatePolyFromPowerBasis(pol polynomialVector) (res *Ciphertext, err error) {

	X := polyEval.powerBasis
	level := X[1].Level()

	params := polyEval.Evaluator.(*evaluator).params
	slotsIndex := polyEval.slotsIndex

	minimumDegreeNonZeroCoefficient := 0

	// Get the minimum non-zero degree coefficient
	for i := pol.Value[0].Degree(); i > 0; i-- {
		for _, p := range pol.Value {
			if p.Coeffs[i] != 0 {
				minimumDegreeNonZeroCoefficient = utils.MaxInt(minimumDegreeNonZeroCoefficient, i)
				break
			}
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
			res = NewCiphertextLvl(params, 1, level)

			// Looks for non-zero coefficients among the degree-0 coefficients of the polynomials
			for i, p := range pol.Value {
				if p.Coeffs[0] != 0 {
					toEncode = true
					for _, j := range slotsIndex[i] {
						values[j] = p.Coeffs[0]
					}
				}
			}

			// If a non-zero coefficient was found, encodes the values, adds on the ciphertext, and returns
			if toEncode {
				polyEval.Encode(values, &Plaintext{&rlwe.Plaintext{Value: res.Value[0]}})
			}

			return
		}

		// Allocates the output ciphertext
		res = NewCiphertextLvl(params, 1, level)

		// Allocates a temporary plaintext to encode the values
		pt := NewPlaintextAtLevelFromPoly(level, polyEval.Evaluator.BuffPt().Value)

		// Looks for a non-zero coefficient among the degree-0 coefficient of the polynomials
		for i, p := range pol.Value {
			if p.Coeffs[0] != 0 {
				toEncode = true
				for _, j := range slotsIndex[i] {
					values[j] = p.Coeffs[0]
				}
			}
		}

		// If a non-zero degree coefficient was found, encodes and adds the values on the output
		// ciphertext
		if toEncode {
			polyEval.Encode(values, pt)
			polyEval.Add(res, pt, res)
			toEncode = false
		}

		// Loops starting from the highest-degree coefficient
		for key := pol.Value[0].Degree(); key > 0; key-- {

			var reset bool
			// Loops over the polynomials
			for i, p := range pol.Value {

				// Looks for a non-zero coefficient
				if p.Coeffs[key] != 0 {
					toEncode = true

					// Resets the temporary array to zero.
					// This is needed if a zero coefficient
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

			// If a non-zero degree coefficient was found, encodes and adds the values on the output
			// ciphertext
			if toEncode {
				polyEval.EncodeMul(values, &PlaintextMul{pt.Plaintext})
				polyEval.MulAndAdd(X[key], &PlaintextMul{pt.Plaintext}, res)
				toEncode = false
			}
		}

	} else {

		c := pol.Value[0].Coeffs[0]

		if minimumDegreeNonZeroCoefficient == 0 {

			res = NewCiphertextLvl(params, 1, level)

			if c != 0 {
				polyEval.AddScalar(res, c, res)
			}

			return
		}

		res = NewCiphertextLvl(params, 1, level)

		if c != 0 {
			polyEval.AddScalar(res, c, res)
		}

		for key := pol.Value[0].Degree(); key > 0; key-- {
			c = pol.Value[0].Coeffs[key]
			if key != 0 && c != 0 {
				polyEval.MulScalarAndAdd(X[key], c, res)
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
