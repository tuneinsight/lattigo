package ckks

import (
	"math"
	"math/bits"
)

// EvaluateCheby evaluates the input Chebyshev polynomial with the input ciphertext.
func (evaluator *Evaluator) EvaluateCheby(ct *Ciphertext, cheby *ChebyshevInterpolation, evakey *EvaluationKey) (res *Ciphertext) {

	C := make(map[uint64]*Ciphertext)

	C[1] = ct.CopyNew().Ciphertext()

	evaluator.MultByConst(C[1], 2/(cheby.b-cheby.a), C[1])
	evaluator.AddConst(C[1], (-cheby.a-cheby.b)/(cheby.b-cheby.a), C[1])
	evaluator.Rescale(C[1], evaluator.ckkscontext.scale, C[1])

	C[2] = evaluator.MulRelinNew(C[1], C[1], evakey)
	evaluator.Rescale(C[2], evaluator.ckkscontext.scale, C[2])

	evaluator.Add(C[2], C[2], C[2])
	evaluator.AddConst(C[2], -1, C[2])

	M := uint64(bits.Len64(cheby.degree - 1))
	L := uint64(M >> 1)

	for i := uint64(3); i <= (1 << L); i++ {
		computePowerBasisCheby(i, C, evaluator, evakey)
	}

	for i := L + 1; i < M; i++ {
		computePowerBasisCheby(1<<i, C, evaluator, evakey)
	}

	return recurseCheby(cheby.degree, L, M, cheby.coeffs, C, evaluator, evakey)
}

func computePowerBasisCheby(n uint64, C map[uint64]*Ciphertext, evaluator *Evaluator, evakey *EvaluationKey) {

	// Given a hash table with the first three evaluations of the Chebyshev ring at x in the interval a, b:
	// C0 = 1 (actually not stored in the hash table)
	// C1 = (2*x - a - b)/(b-a)
	// C2 = 2*C1*C1 - C0
	// Evaluates the nth degree Chebyshev ring in a recursive manner, storing intermediate results in the hashtable.
	// Consumes at most ceil(sqrt(n)) levels for an evaluation at Cn.
	// Uses the following property : for a given Chebyshev ring Cn = 2*Ca*Cb - Cc, n = a+b and c = abs(a-b)

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

		evaluator.Rescale(C[n], evaluator.ckkscontext.scale, C[n])

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

func recurseCheby(maxDegree, L, M uint64, coeffs map[uint64]complex128, C map[uint64]*Ciphertext, evaluator *Evaluator, evakey *EvaluationKey) (res *Ciphertext) {

	// Recursively computes the evalution of the Chebyshev polynomial using a baby-set giant-step algorithm.

	if maxDegree <= (1 << L) {
		return evaluateRecurse(coeffs, C, evaluator, evakey)
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

	evaluator.Rescale(res, evaluator.ckkscontext.scale, res)

	return res

}

func splitCoeffsCheby(coeffs map[uint64]complex128, degree, maxDegree uint64) (coeffsq, coeffsr map[uint64]complex128) {

	// Splits a Chebyshev polynomial p such that p = q*C^degree + r, where q and r are a linear combination of a Chebyshev basis.

	coeffsr = make(map[uint64]complex128)
	coeffsq = make(map[uint64]complex128)

	for i := uint64(0); i < degree; i++ {
		coeffsr[i] = coeffs[i]
	}

	coeffsq[0] = coeffs[degree]

	for i := uint64(degree + 1); i < maxDegree+1; i++ {
		coeffsq[i-degree] = 2 * coeffs[i]
		coeffsr[2*degree-i] -= coeffs[i]
	}

	return coeffsq, coeffsr
}
