package ckks

import (
	"math"
	"math/bits"
)

type poly struct {
	maxDeg uint64
	coeffs []complex128
}

func (p *poly) degree() uint64 {
	return uint64(len(p.coeffs) - 1)
}

// EvaluateChebyFast evaluates the input Chebyshev polynomial with the input ciphertext.
// Faster than EvaluateChebyEco but consumes ceil(log(deg)) + 2 levels.
func (eval *evaluator) EvaluateChebyFast(op *Ciphertext, cheby *ChebyshevInterpolation, evakey *EvaluationKey) (opOut *Ciphertext) {

	C := make(map[uint64]*Ciphertext)

	C[1] = op.CopyNew().Ciphertext()

	eval.MultByConst(C[1], 2/(cheby.b-cheby.a), C[1])
	eval.AddConst(C[1], (-cheby.a-cheby.b)/(cheby.b-cheby.a), C[1])
	eval.Rescale(C[1], eval.ckksContext.scale, C[1])

	M := uint64(bits.Len64(cheby.degree()))
	L := uint64(M >> 1)

	for i := uint64(2); i <= (1 << L); i++ {
		computePowerBasisCheby(i, C, eval, evakey)
	}

	for i := L + 1; i < M; i++ {
		computePowerBasisCheby(1<<i, C, eval, evakey)
	}

	return recurseCheby(cheby.degree(), L, M, cheby.Poly(), C, eval, evakey)
}

// EvaluateChebyEco evaluates the input Chebyshev polynomial on the input ciphertext.
// It is slower than EvaluateChebyFast but consumes one less level (ceil(log(deg)) + 1 levels).
func (eval *evaluator) EvaluateChebyEco(op *Ciphertext, cheby *ChebyshevInterpolation, evakey *EvaluationKey) (opOut *Ciphertext) {

	C := make(map[uint64]*Ciphertext)

	C[1] = op.CopyNew().Ciphertext()

	eval.MultByConst(C[1], 2/(cheby.b-cheby.a), C[1])
	eval.AddConst(C[1], (-cheby.a-cheby.b)/(cheby.b-cheby.a), C[1])
	eval.Rescale(C[1], eval.ckksContext.scale, C[1])

	M := uint64(bits.Len64(cheby.degree() - 1))
	L := uint64(1)

	for i := uint64(2); i <= (1 << L); i++ {
		computePowerBasisCheby(i, C, eval, evakey)
	}

	for i := L + 1; i < M; i++ {
		computePowerBasisCheby(1<<i, C, eval, evakey)
	}

	return recurseCheby(cheby.degree(), L, M, cheby.Poly(), C, eval, evakey)
}

// EvaluateChebyEcoSpecial evaluates the input Chebyshev polynomial with the input ciphertext.
// Slower than EvaluateChebyFast but consumes ceil(log(deg)) + 1 levels.
func (evaluator *evaluator) EvaluateChebyEcoSpecial(ct *Ciphertext, n complex128, cheby *ChebyshevInterpolation, evakey *EvaluationKey) (res *Ciphertext) {

	C := make(map[uint64]*Ciphertext)

	C[1] = ct.CopyNew().Ciphertext()

	evaluator.MultByConst(C[1], 2/((cheby.b-cheby.a)*n), C[1])
	evaluator.AddConst(C[1], (-cheby.a-cheby.b)/(cheby.b-cheby.a), C[1])
	evaluator.Rescale(C[1], evaluator.params.Scale, C[1])

	M := uint64(bits.Len64(cheby.degree() - 1))
	L := uint64(1)

	for i := uint64(2); i <= (1 << L); i++ {
		computePowerBasisCheby(i, C, evaluator, evakey)
	}

	for i := L + 1; i < M; i++ {
		computePowerBasisCheby(1<<i, C, evaluator, evakey)
	}

	return recurseCheby(cheby.degree(), L, M, cheby.Poly(), C, evaluator, evakey)
}

// EvaluateChebyFastSpecial evaluates the input Chebyshev polynomial with the input ciphertext.
// Slower than EvaluateChebyFast but consumes ceil(log(deg)) + 1 levels.
func (evaluator *evaluator) EvaluateChebyFastSpecial(ct *Ciphertext, n complex128, cheby *ChebyshevInterpolation, evakey *EvaluationKey) (res *Ciphertext) {

	C := make(map[uint64]*Ciphertext)

	C[1] = ct.CopyNew().Ciphertext()

	evaluator.MultByConst(C[1], 2/((cheby.b-cheby.a)*n), C[1])
	evaluator.AddConst(C[1], (-cheby.a-cheby.b)/(cheby.b-cheby.a), C[1])
	evaluator.Rescale(C[1], evaluator.params.Scale, C[1])

	M := uint64(bits.Len64(cheby.degree() - 1))
	L := uint64(1)

	for i := uint64(2); i <= (1 << L); i++ {
		computePowerBasisCheby(i, C, evaluator, evakey)
	}

	for i := L + 1; i < M; i++ {
		computePowerBasisCheby(1<<i, C, evaluator, evakey)
	}

	return recurseCheby(cheby.degree(), L, M, cheby.Poly(), C, evaluator, evakey)
}

func computePowerBasisCheby(n uint64, C map[uint64]*Ciphertext, evaluator *evaluator, evakey *EvaluationKey) {

	// Given a hash table with the first three evaluations of the Chebyshev ring at x in the interval a, b:
	// C0 = 1 (actually not stored in the hash table)
	// C1 = (2*x - a - b)/(b-a)
	// C2 = 2*C1*C1 - C0
	// Evaluates the nth degree Chebyshev ring in a recursive manner, storing intermediate results in the hashtable.
	// Consumes at most ceil(sqrt(n)) levels for an evaluation at Cn.
	// Uses the following property: for a given Chebyshev ring Cn = 2*Ca*Cb - Cc, n = a+b and c = abs(a-b)

	if C[n] == nil {

		// Computes the index required to compute the asked ring evaluation
		a := uint64(math.Ceil(float64(n) / 2))
		b := n >> 1
		c := uint64(math.Abs(float64(a) - float64(b)))

		// Recurses on the given indexes
		computePowerBasisCheby(a, C, evaluator, evakey)
		computePowerBasisCheby(b, C, evaluator, evakey)

		// Since C[0] is not stored (but rather seen as the constant 1), only recurses on c if c!= 0
		if c != 0 {
			computePowerBasisCheby(c, C, evaluator, evakey)
		}

		// Computes C[n] = C[a]*C[b]
		C[n] = evaluator.MulRelinNew(C[a], C[b], evakey)

		evaluator.Rescale(C[n], evaluator.ckksContext.scale, C[n])

		// Computes C[n] = 2*C[a]*C[b]
		evaluator.Add(C[n], C[n], C[n])

		// Computes C[n] = 2*C[a]*C[b] - C[c]
		if c == 0 {
			evaluator.AddConst(C[n], -1, C[n])
		} else {
			evaluator.Sub(C[n], C[c], C[n])
		}
	}
}

func splitCoeffsCheby(coeffs *poly, degree, maxDegree uint64) (coeffsq, coeffsr *poly) {

	// Splits a Chebyshev polynomial p such that p = q*C^degree + r, where q and r are a linear combination of a Chebyshev basis.

	coeffsr = new(poly)
	coeffsr.maxDeg = coeffs.maxDeg - degree
	coeffsr.coeffs = make([]complex128, degree)
	for i := uint64(0); i < degree; i++ {
		coeffsr.coeffs[i] = coeffs.coeffs[i]
	}

	coeffsq = new(poly)
	coeffsq.coeffs = make([]complex128, maxDegree-degree+1)
	coeffsq.maxDeg = coeffs.maxDeg
	coeffsq.coeffs[0] = coeffs.coeffs[degree]
	for i := uint64(degree + 1); i < maxDegree+1; i++ {
		coeffsq.coeffs[i-degree] = 2 * coeffs.coeffs[i]
		coeffsr.coeffs[2*degree-i] -= coeffs.coeffs[i]
	}

	return coeffsq, coeffsr
}

func recurseCheby(maxDegree, L, M uint64, coeffs *poly, C map[uint64]*Ciphertext, evaluator *evaluator, evakey *EvaluationKey) (res *Ciphertext) {

	// Recursively computes the evalution of the Chebyshev polynomial using a baby-set giant-step algorithm.
	if maxDegree <= (1 << L) {
		return evaluatePolyFromChebyBasis(coeffs, C, evaluator, evakey)
	}

	for 1<<(M-1) > maxDegree {
		M--
	}

	coeffsq, coeffsr := splitCoeffsCheby(coeffs, 1<<(M-1), maxDegree)

	res = recurseCheby(maxDegree-(1<<(M-1)), L, M-1, coeffsq, C, evaluator, evakey)

	var tmp *Ciphertext
	tmp = recurseCheby((1<<(M-1))-1, L, M-1, coeffsr, C, evaluator, evakey)

	evaluator.MulRelin(res, C[1<<(M-1)], evakey, res)

	evaluator.Add(res, tmp, res)
	evaluator.Rescale(res, evaluator.ckksContext.scale, res)

	return res

}

func evaluatePolyFromChebyBasis(coeffs *poly, C map[uint64]*Ciphertext, evaluator *evaluator, evakey *EvaluationKey) (res *Ciphertext) {

	res = NewCiphertext(evaluator.params, 1, C[1].Level(), 1.718160618321759e+10)

	for key := coeffs.degree(); key > 0; key-- {
		if key != 0 && (math.Abs(real(coeffs.coeffs[key])) > 1e-15 || math.Abs(imag(coeffs.coeffs[key])) > 1e-15) {
			evaluator.MultByConstAndAdd(C[key], coeffs.coeffs[key], res)
		}
	}

	if math.Abs(real(coeffs.coeffs[0])) > 1e-15 || math.Abs(imag(coeffs.coeffs[0])) > 1e-15 {
		evaluator.AddConst(res, coeffs.coeffs[0], res)
	}

	evaluator.Rescale(res, evaluator.ckksContext.scale, res)

	return
}
