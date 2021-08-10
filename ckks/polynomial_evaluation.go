package ckks

import (
	"fmt"
	"math"
	"math/bits"
)

// Poly is a struct storing the coeffients of a polynomial
// that then can be evaluated on the ciphertext
type Poly struct {
	MaxDeg int
	Coeffs []complex128
	Lead   bool
}

// Depth returns the number of levels needed to evaluate the polynomial.
func (p *Poly) Depth() int {
	return int(math.Ceil(math.Log2(float64(len(p.Coeffs)))))
}

// Degree returns the degree of the polynomial
func (p *Poly) Degree() int {
	return len(p.Coeffs) - 1
}

// NewPoly creates a new Poly from the input coefficients
func NewPoly(coeffs []complex128) (p *Poly) {

	p = new(Poly)
	p.Coeffs = make([]complex128, len(coeffs))
	copy(p.Coeffs, coeffs)
	p.MaxDeg = len(coeffs) - 1
	p.Lead = true

	return
}

// checkEnoughLevels checks that enough levels are available to evaluate the polynomial.
// Also checks if c is a gaussian integer or not. If not, then one more level is needed
// to evaluate the polynomial.
func checkEnoughLevels(levels int, pol *Poly, c complex128) (err error) {

	depth := pol.Depth()

	if real(c) != float64(int64(real(c))) || imag(c) != float64(int64(imag(c))) {
		depth++
	}

	if levels < depth {
		return fmt.Errorf("%d levels < %d log(d) -> cannot evaluate", levels, depth)
	}

	return nil
}

// EvaluatePoly evaluates a polynomial in standard basis on the input Ciphertext in ceil(log2(deg+1)) levels.
// Returns an error if the input ciphertext does not have enough level to carry out the full polynomial evaluation.
// Returns an error if something is wrong with the scale.

func (eval *evaluator) EvaluatePoly(ct0 *Ciphertext, pol *Poly, targetScale float64) (opOut *Ciphertext, err error) {

	if err := checkEnoughLevels(ct0.Level(), pol, 1); err != nil {
		return ct0, err
	}

	C := make(map[int]*Ciphertext)

	C[1] = ct0.CopyNew()

	logDegree := bits.Len64(uint64(pol.Degree()))
	logSplit := (logDegree >> 1) //optimalSplit(logDegree) //

	for i := 2; i < (1 << logSplit); i++ {
		if err = computePowerBasis(i, C, targetScale, eval); err != nil {
			return nil, err
		}
	}

	for i := logSplit; i < logDegree; i++ {
		if err = computePowerBasis(1<<i, C, targetScale, eval); err != nil {
			return nil, err
		}
	}

	opOut, err = recurse(targetScale, logSplit, logDegree, pol, C, eval)

	C = nil
	return opOut, err
}

// EvaluateCheby evaluates a polynomial in Chebyshev basis on the input Ciphertext in ceil(log2(deg+1))+1 levels.
// Returns an error if the input ciphertext does not have enough level to carry out the full polynomial evaluation.
// Returns an error if something is wrong with the scale.
// A change of basis ct' = (2/(b-a)) * (ct + (-a-b)/(b-a)) is necessary before the polynomial evaluation to ensure correctness.

func (eval *evaluator) EvaluateCheby(op *Ciphertext, cheby *ChebyshevInterpolation, targetScale float64) (opOut *Ciphertext, err error) {

	if err := checkEnoughLevels(op.Level(), &cheby.Poly, 1); err != nil {
		return op, err
	}

	C := make(map[int]*Ciphertext)

	C[1] = op.CopyNew()

	logDegree := int(bits.Len64(uint64(cheby.Degree())))
	logSplit := (logDegree >> 1) //optimalSplit(logDegree) //

	for i := 2; i < (1 << logSplit); i++ {
		if err = computePowerBasisCheby(i, C, targetScale, eval); err != nil {
			return nil, err
		}
	}

	for i := logSplit; i < logDegree; i++ {
		if err = computePowerBasisCheby(1<<i, C, targetScale, eval); err != nil {
			return nil, err
		}
	}

	opOut, err = recurseCheby(targetScale, logSplit, logDegree, &cheby.Poly, C, eval)

	C = nil

	return opOut, err
}

func computePowerBasis(n int, C map[int]*Ciphertext, scale float64, evaluator *evaluator) (err error) {

	if C[n] == nil {

		// Computes the index required to compute the asked ring evaluation
		a := int(math.Ceil(float64(n) / 2))
		b := n >> 1

		// Recurses on the given indexes
		if err = computePowerBasis(a, C, scale, evaluator); err != nil {
			return err
		}
		if err = computePowerBasis(b, C, scale, evaluator); err != nil {
			return err
		}

		// Computes C[n] = C[a]*C[b]
		C[n] = evaluator.MulRelinNew(C[a], C[b])

		if err = evaluator.Rescale(C[n], scale, C[n]); err != nil {
			return err
		}
	}

	return nil
}

func computePowerBasisCheby(n int, C map[int]*Ciphertext, scale float64, evaluator *evaluator) (err error) {

	// Given a hash table with the first three evaluations of the Chebyshev ring at x in the interval a, b:
	// C0 = 1 (actually not stored in the hash table)
	// C1 = (2*x - a - b)/(b-a)
	// C2 = 2*C1*C1 - C0
	// Evaluates the nth degree Chebyshev ring in a recursive manner, storing intermediate results in the hashtable.
	// Consumes at most ceil(sqrt(n)) levels for an evaluation at Cn.
	// Uses the following property: for a given Chebyshev ring Cn = 2*Ca*Cb - Cc, n = a+b and c = abs(a-b)

	if C[n] == nil {

		// Computes the index required to compute the asked ring evaluation
		var a, b, c int
		if n&(n-1) == 0 {
			a, b, c = n/2, n/2, 0 //Necessary for depth optimality
		} else {
			// [Lee et al. 2020] : High-Precision and Low-Complexity Approximate Homomorphic Encryption by Error Variance Minimization
			// Maximize the number of odd terms of Chebyshev basis
			k := int(math.Ceil(math.Log2(float64(n)))) - 1
			a = (1 << k) - 1
			b = n + 1 - (1 << k)
			c = int(math.Abs(float64(a) - float64(b)))
		}

		// Recurses on the given indexes
		if err = computePowerBasisCheby(a, C, scale, evaluator); err != nil {
			return err
		}
		if err = computePowerBasisCheby(b, C, scale, evaluator); err != nil {
			return err
		}

		// Since C[0] is not stored (but rather seen as the constant 1), only recurses on c if c!= 0
		if c != 0 {
			if err = computePowerBasisCheby(c, C, scale, evaluator); err != nil {
				return err
			}
		}

		// Computes C[n] = C[a]*C[b]
		//fmt.Println("Mul", C[a].Level(), C[b].Level())
		C[n] = evaluator.MulRelinNew(C[a], C[b])
		if err = evaluator.Rescale(C[n], scale, C[n]); err != nil {
			return err
		}

		// Computes C[n] = 2*C[a]*C[b]
		evaluator.Add(C[n], C[n], C[n])

		// Computes C[n] = 2*C[a]*C[b] - C[c]
		if c == 0 {
			evaluator.AddConst(C[n], -1, C[n])
		} else {
			evaluator.Sub(C[n], C[c], C[n])
		}

	}

	return nil
}

func splitCoeffs(coeffs *Poly, split int) (coeffsq, coeffsr *Poly) {

	// Splits a polynomial p such that p = q*C^degree + r.

	coeffsr = new(Poly)
	coeffsr.Coeffs = make([]complex128, split)
	if coeffs.MaxDeg == coeffs.Degree() {
		coeffsr.MaxDeg = split - 1
	} else {
		coeffsr.MaxDeg = coeffs.MaxDeg - (coeffs.Degree() - split + 1)
	}

	for i := 0; i < split; i++ {
		coeffsr.Coeffs[i] = coeffs.Coeffs[i]
	}

	coeffsq = new(Poly)
	coeffsq.Coeffs = make([]complex128, coeffs.Degree()-split+1)
	coeffsq.MaxDeg = coeffs.MaxDeg

	coeffsq.Coeffs[0] = coeffs.Coeffs[split]
	for i := split + 1; i < coeffs.Degree()+1; i++ {
		coeffsq.Coeffs[i-split] = coeffs.Coeffs[i]
	}

	if coeffs.Lead {
		coeffsq.Lead = true
	}

	return coeffsq, coeffsr
}

func splitCoeffsCheby(coeffs *Poly, split int) (coeffsq, coeffsr *Poly) {

	// Splits a Chebyshev polynomial p such that p = q*C^degree + r, where q and r are a linear combination of a Chebyshev basis.
	coeffsr = new(Poly)
	coeffsr.Coeffs = make([]complex128, split)
	if coeffs.MaxDeg == coeffs.Degree() {
		coeffsr.MaxDeg = split - 1
	} else {
		coeffsr.MaxDeg = coeffs.MaxDeg - (coeffs.Degree() - split + 1)
	}

	for i := 0; i < split; i++ {
		coeffsr.Coeffs[i] = coeffs.Coeffs[i]
	}

	coeffsq = new(Poly)
	coeffsq.Coeffs = make([]complex128, coeffs.Degree()-split+1)
	coeffsq.MaxDeg = coeffs.MaxDeg

	coeffsq.Coeffs[0] = coeffs.Coeffs[split]
	for i, j := split+1, 1; i < coeffs.Degree()+1; i, j = i+1, j+1 {
		coeffsq.Coeffs[i-split] = 2 * coeffs.Coeffs[i]
		coeffsr.Coeffs[split-j] -= coeffs.Coeffs[i]
	}

	if coeffs.Lead {
		coeffsq.Lead = true
	}

	return coeffsq, coeffsr
}

func recurse(targetScale float64, logSplit, logDegree int, coeffs *Poly, C map[int]*Ciphertext, evaluator *evaluator) (res *Ciphertext, err error) {

	// Recursively computes the evalution of the Chebyshev polynomial using a baby-set giant-step algorithm.
	if coeffs.Degree() < (1 << logSplit) {

		if coeffs.Lead && coeffs.MaxDeg > ((1<<logDegree)-(1<<(logSplit-1))) && logSplit > 1 {

			logDegree = int(bits.Len64(uint64(coeffs.Degree())))
			logSplit = logDegree >> 1

			return recurse(targetScale, logSplit, logDegree, coeffs, C, evaluator)
		}

		return evaluatePolyFromPowerBasis(targetScale, coeffs, C, evaluator)
	}

	var nextPower = 1 << logSplit
	for nextPower < (coeffs.Degree()>>1)+1 {
		nextPower <<= 1
	}

	coeffsq, coeffsr := splitCoeffs(coeffs, nextPower)

	level := C[nextPower].Level() - 1

	if coeffsq.MaxDeg >= 1<<(logDegree-1) && coeffsq.Lead {
		level++
	}

	currentQi := evaluator.params.QiFloat64(level)

	//fmt.Printf("X^%2d: %d %d %t %d\n", nextPower, coeffsq.MaxDeg, coeffsr.MaxDeg, coeffsq.MaxDeg >= 1<<(logDegree-1), level)
	//fmt.Printf("X^%2d: %f %f\n", nextPower, targetScale, targetScale* currentQi / C[nextPower].Scale())
	//fmt.Printf("X^%2d : qi %d %t %d %d\n", nextPower, level, coeffsq.Lead, coeffsq.MaxDeg, 1<<(logDegree-1))
	//fmt.Println()
	var tmp *Ciphertext
	if res, err = recurse(targetScale*currentQi/C[nextPower].Scale, logSplit, logDegree, coeffsq, C, evaluator); err != nil {
		return nil, err
	}

	if tmp, err = recurse(targetScale, logSplit, logDegree, coeffsr, C, evaluator); err != nil {
		return nil, err
	}

	if res.Level() > tmp.Level() {
		for res.Level() != tmp.Level()+1 {
			evaluator.DropLevel(res, 1)
		}
	}

	//fmt.Printf("X^%2d: (%d %f -> \n", nextPower, res.Level(), res.Scale())
	evaluator.MulRelin(res, C[nextPower], res)

	if res.Level() > tmp.Level() {
		if err = evaluator.Rescale(res, targetScale, res); err != nil {
			return nil, err
		}
		//fmt.Printf("%f = %d) + (%d %f) = ", res.Scale(), res.Level(), tmp.Level(), tmp.Scale())
		evaluator.Add(res, tmp, res)
		//fmt.Printf("(%d %f) %f\n", res.Level(), res.Scale(), res.Scale()-tmp.Scale())
	} else {
		evaluator.Add(res, tmp, res)
		if err = evaluator.Rescale(res, targetScale, res); err != nil {
			return nil, err
		}

	}

	tmp = nil

	return
}

func recurseCheby(targetScale float64, logSplit, logDegree int, coeffs *Poly, C map[int]*Ciphertext, evaluator *evaluator) (res *Ciphertext, err error) {

	// Recursively computes the evalution of the Chebyshev polynomial using a baby-set giant-step algorithm.
	if coeffs.Degree() < (1 << logSplit) {

		if coeffs.Lead && coeffs.MaxDeg > ((1<<logDegree)-(1<<(logSplit-1))) && logSplit > 1 {

			logDegree = int(bits.Len64(uint64(coeffs.Degree())))
			logSplit = logDegree >> 1

			return recurseCheby(targetScale, logSplit, logDegree, coeffs, C, evaluator)
		}

		return evaluatePolyFromPowerBasis(targetScale, coeffs, C, evaluator)
	}

	var nextPower = 1 << logSplit
	for nextPower < (coeffs.Degree()>>1)+1 {
		nextPower <<= 1
	}

	coeffsq, coeffsr := splitCoeffsCheby(coeffs, nextPower)

	level := C[nextPower].Level() - 1

	if coeffsq.MaxDeg >= 1<<(logDegree-1) && coeffsq.Lead {
		level++
	}

	currentQi := evaluator.params.QiFloat64(level)

	//fmt.Printf("X^%2d: %d %d %t %d\n", nextPower, coeffsq.MaxDeg, coeffsr.MaxDeg, coeffsq.MaxDeg >= 1<<(logDegree-1), level)
	//fmt.Printf("X^%2d: %f %f\n", nextPower, targetScale, targetScale* currentQi / C[nextPower].Scale())
	//fmt.Printf("X^%2d : qi %d %t %d %d\n", nextPower, level, coeffsq.Lead, coeffsq.MaxDeg, 1<<(logDegree-1))
	//fmt.Println()

	if res, err = recurseCheby(targetScale*currentQi/C[nextPower].Scale, logSplit, logDegree, coeffsq, C, evaluator); err != nil {
		return nil, err
	}

	var tmp *Ciphertext
	if tmp, err = recurseCheby(targetScale, logSplit, logDegree, coeffsr, C, evaluator); err != nil {
		return nil, err
	}

	if res.Level() > tmp.Level() {
		for res.Level() != tmp.Level()+1 {
			evaluator.DropLevel(res, 1)
		}
	}

	//fmt.Printf("X^%2d: (%d %f -> \n", nextPower, res.Level(), res.Scale())
	evaluator.MulRelin(res, C[nextPower], res)

	if res.Level() > tmp.Level() {
		if err = evaluator.Rescale(res, targetScale, res); err != nil {
			return nil, err
		}
		//fmt.Printf("%f = %d) + (%d %f) = ", res.Scale(), res.Level(), tmp.Level(), tmp.Scale())
		evaluator.Add(res, tmp, res)
		//fmt.Printf("(%d %f) %f\n", res.Level(), res.Scale(), res.Scale()-tmp.Scale())
	} else {
		evaluator.Add(res, tmp, res)
		if err = evaluator.Rescale(res, targetScale, res); err != nil {
			return nil, err
		}
	}

	tmp = nil

	return

}

func evaluatePolyFromPowerBasis(targetScale float64, coeffs *Poly, C map[int]*Ciphertext, evaluator *evaluator) (res *Ciphertext, err error) {

	if coeffs.Degree() == 0 {

		res = NewCiphertext(evaluator.params, 1, C[1].Level(), targetScale)

		if math.Abs(real(coeffs.Coeffs[0])) > 1e-14 || math.Abs(imag(coeffs.Coeffs[0])) > 1e-14 {
			evaluator.AddConst(res, coeffs.Coeffs[0], res)
		}

		return
	}

	currentQi := evaluator.params.QiFloat64(C[coeffs.Degree()].Level())

	ctScale := targetScale * currentQi

	//fmt.Printf("%d %f\n", coeffs.MaxDeg, targetScale)
	//fmt.Println("current Qi", evaluator.params.Qi[C[coeffs.Degree()].Level()])
	//fmt.Println(coeffs.Degree(), C[coeffs.Degree()].Level())

	res = NewCiphertext(evaluator.params, 1, C[coeffs.Degree()].Level(), ctScale)

	if math.Abs(real(coeffs.Coeffs[0])) > 1e-14 || math.Abs(imag(coeffs.Coeffs[0])) > 1e-14 {
		evaluator.AddConst(res, coeffs.Coeffs[0], res)
	}

	for key := coeffs.Degree(); key > 0; key-- {

		if key != 0 && (math.Abs(real(coeffs.Coeffs[key])) > 1e-14 || math.Abs(imag(coeffs.Coeffs[key])) > 1e-14) {

			// Target scale * rescale-scale / power basis scale
			constScale := targetScale * currentQi / C[key].Scale

			cReal := int64(real(coeffs.Coeffs[key]) * constScale)
			cImag := int64(imag(coeffs.Coeffs[key]) * constScale)

			evaluator.MultByGaussianIntegerAndAdd(C[key], cReal, cImag, res)
		}
	}

	if err = evaluator.Rescale(res, targetScale, res); err != nil {
		return nil, err
	}

	return
}
