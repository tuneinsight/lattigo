package ckks

import (
	"math"
	"math/bits"
	"math/cmplx"
)

func chebyshevNodesU(n int, a, b complex128) (u []complex128) {
	u = make([]complex128, n)
	var x, y complex128
	for k := 1; k < n+1; k++ {
		x = 0.5 * (a + b)
		y = 0.5 * (b - a)
		u[n-k] = x + y*complex(math.Cos((float64(k)-0.5)*(3.141592653589793/float64(n))), 0)
	}
	return
}

func chebyshevNodesX(u []complex128, a, b complex128) (x []complex128) {
	x = make([]complex128, len(u))
	for i := 0; i < len(u); i++ {
		x[i] = 0.5*(b-a)*u[i] + 0.5*(a+b)
	}
	return
}

func average(y []complex128) (avg complex128) {
	avg = 0
	for i := 0; i < len(y); i++ {
		avg += y[i]
	}
	avg /= complex(float64(len(y)), 0)

	return
}

func evaluate_cheby(coeffs []complex128, x complex128, a, b complex128) (y complex128) {
	var u, Tprev, Tnext, T complex128
	u = (2*x - a - b) / (b - a)
	Tprev = 1
	T = u
	y = coeffs[0]
	for i := 1; i < len(coeffs); i++ {
		y = y + T*coeffs[i]
		Tnext = 2*u*T - Tprev
		Tprev = T
		T = Tnext
	}
	return
}

func evaluate(degree int, x, a, b complex128) (T complex128) {
	if degree == 0 {
		return 1
	}
	var u, Tprev, Tnext complex128
	u = (2*x - a - b) / (b - a)
	Tprev = 1
	T = u
	for i := 1; i < degree; i++ {
		Tnext = 2*u*T - Tprev
		Tprev = T
		T = Tnext
	}
	return
}

func chebyCoeffs(u, y []complex128, a, b complex128) (coeffs []complex128) {
	n := len(y)
	coeffs = make([]complex128, n)
	var tmp []complex128
	for i := 0; i < n; i++ {
		tmp = make([]complex128, n)
		for j := 0; j < n; j++ {
			tmp[j] = y[j] * evaluate(i, u[j], -1, 1)
		}
		if i != 0 {
			coeffs[i] = 2 * average(tmp)
		} else {
			coeffs[i] = average(tmp)
		}
	}

	return
}

type ChebyshevInterpolation struct {
	coeffs map[uint64]complex128
	degree uint64
	a      complex128
	b      complex128
}

// Approximate computes a Chebyshev approximation of the input function, for the tange [-a, b] of degree degree.
// To be used in conjonction with the function EvaluateCheby.
func Approximate(function func(complex128) complex128, a, b complex128, degree int) (cheby *ChebyshevInterpolation) {

	cheby = new(ChebyshevInterpolation)
	cheby.coeffs = make(map[uint64]complex128)
	cheby.a = a
	cheby.b = b
	cheby.degree = uint64(degree)

	u := chebyshevNodesU(degree+1, -1, 1)
	x := chebyshevNodesX(u, a, b)

	y := make([]complex128, len(x))
	for i := range x {
		y[i] = function(x[i])
	}

	coeffs := chebyCoeffs(u, y, a, b)

	for i := range coeffs {
		cheby.coeffs[uint64(i)] = coeffs[i]
	}

	return
}

// Given a hash table with the first three evaluations of the Chebyshev ring at x in the interval a, b:
// C0 = 1 (actually not stored in the hash table)
// C1 = (2*x - a - b)/(b-a)
// C2 = 2*C1*C1 - C0
// Evaluates the nth degree Chebyshev ring in a recursive manner, storing intermediate results in the hashtable.
// Consumes at most ceil(sqrt(n)) levels for an evaluation at Cn.
// Uses the following property : for a given Chebyshev ring Cn = 2*Ca*Cb - Cc, n = a+b and c = abs(a-b)
func computePowerBasis(n uint64, C map[uint64]*Ciphertext, evaluator *Evaluator, evakey *EvaluationKey) (err error) {

	if C[n] == nil {

		// Computes the index required to compute the asked ring evaluation
		a := uint64(math.Ceil(float64(n) / 2))
		b := n >> 1
		c := uint64(math.Abs(float64(a) - float64(b)))

		// Recurses on the given indexes
		computePowerBasis(a, C, evaluator, evakey)
		computePowerBasis(b, C, evaluator, evakey)

		// Since C[0] is not stored (but rather seen as the constant 1), only recurses on c if c!= 0
		if c != 0 {
			computePowerBasis(c, C, evaluator, evakey)
		}

		// Computes C[n] = C[a]*C[b]
		if C[n], err = evaluator.MulRelinNew(C[a], C[b], evakey); err != nil {
			return err
		}

		evaluator.Rescale(C[n], evaluator.ckkscontext.scale, C[n])

		// Computes C[n] = 2*C[a]*C[b]
		evaluator.Add(C[n], C[n], C[n])

		// Computes C[n] = 2*C[a]*C[b] - C[c]
		if c == 0 {
			evaluator.AddConst(C[n], -1, C[n])
		} else {

			if err = evaluator.Sub(C[n], C[c], C[n]); err != nil {
				return err
			}
		}

		return nil
	}

	return nil
}

func evaluateRecurse(coeffs map[uint64]complex128, C map[uint64]*Ciphertext, evaluator *Evaluator, evakey *EvaluationKey) (res *Ciphertext, err error) {

	res = evaluator.ckkscontext.NewCiphertext(1, C[1].Level(), C[1].Scale())

	if math.Abs(real(coeffs[0])) > 1e-15 || math.Abs(imag(coeffs[0])) > 1e-15 {
		if err = evaluator.AddConst(res, coeffs[0], res); err != nil {
			return nil, err
		}
	}

	for key := range coeffs {
		if key != 0 && (math.Abs(real(coeffs[key])) > 1e-15 || math.Abs(imag(coeffs[key])) > 1e-15) {

			if err = evaluator.MultByConstAndAdd(C[key], coeffs[key], res); err != nil {
				return nil, err
			}
		}
	}

	if err = evaluator.Rescale(res, evaluator.ckkscontext.scale, res); err != nil {
		return nil, err
	}

	return
}

func recurse(maxDegree, L, M uint64, coeffs map[uint64]complex128, C map[uint64]*Ciphertext, evaluator *Evaluator, evakey *EvaluationKey) (res *Ciphertext, err error) {

	if maxDegree <= 1<<L {

		return evaluateRecurse(coeffs, C, evaluator, evakey)

	} else {

		for 1<<(M-1) > maxDegree{
			M--
		}

		coeffsq, coeffsr := split(coeffs, 1<<(M-1), maxDegree)

		if res, err = recurse(maxDegree-(1<<(M-1)), L, M-1, coeffsq, C, evaluator, evakey); err != nil {
			return nil, err
		}

		var tmp *Ciphertext
		if tmp, err = recurse((1<<(M-1))-1, L, M-1, coeffsr, C, evaluator, evakey); err != nil {
			return nil, err
		}

		if err = evaluator.MulRelin(res, C[1<<(M-1)], evakey, res); err != nil {
			return nil, err
		}

		if err = evaluator.Add(res, tmp, res); err != nil {
			return nil, err
		}

		if err = evaluator.Rescale(res, evaluator.ckkscontext.scale, res); err != nil {
			return nil, err
		}

		return res, nil
	}
}

func split(coeffs map[uint64]complex128, degree, maxDegree uint64) (coeffsq, coeffsr map[uint64]complex128) {

	// Splits a Chebyshev polynomial p such that p = q*C^degree + r, where q and r are a linear combination of a Chebyshev basis.

	coeffsr = make(map[uint64]complex128)
	coeffsq = make(map[uint64]complex128)

	for i := uint64(0); i < degree; i++ {
		if coeffs[i] != 0 {
			coeffsr[i] = coeffs[i]
		}
	}

	for i := uint64(degree); i < maxDegree+1; i++ {
		if coeffs[i] != 0 {
			if i != degree {
				coeffsq[i-degree] = 2 * coeffs[i]
				coeffsr[2*degree-i] -= coeffs[i]
			} else {
				coeffsr[i-degree] = coeffs[i]
			}
		}
	}

	return coeffsq, coeffsr
}

func (evaluator *Evaluator) EvaluateCheby(ct *Ciphertext, cheby *ChebyshevInterpolation, evakey *EvaluationKey) (res *Ciphertext, err error) {

	C := make(map[uint64]*Ciphertext)

	C[1] = ct.CopyNew().Ciphertext()

	evaluator.MultConst(C[1], 2/(cheby.b-cheby.a), C[1])
	evaluator.AddConst(C[1], (-cheby.a-cheby.b)/(cheby.b-cheby.a), C[1])
	evaluator.Rescale(C[1], evaluator.ckkscontext.scale, C[1])

	C[2], _ = evaluator.MulRelinNew(C[1], C[1], evakey)
	evaluator.Rescale(C[2], evaluator.ckkscontext.scale, C[2])

	evaluator.Add(C[2], C[2], C[2])
	evaluator.AddConst(C[2], -1, C[2])

	M := uint64(bits.Len64(cheby.degree - 1))
	L := uint64(M >> 1)

	for i := uint64(3); i < (1<<L)+1; i++ {
		if err = computePowerBasis(i, C, evaluator, evakey); err != nil {
			return nil, err
		}
	}

	for i := L + 1; i < M; i++ {
		if err = computePowerBasis(1<<i, C, evaluator, evakey); err != nil {
			return nil, err
		}
	}

	return recurse(cheby.degree, L, M, cheby.coeffs, C, evaluator, evakey)
}

func exp2pi(x complex128) complex128 {
	return cmplx.Exp(2 * 3.141592653589793 * complex(0, 1) * x)
}

func sin2pi2pi(x complex128) complex128 {
	return cmplx.Sin(6.283185307179586*(1.0/1.0)*x) * (1.0 / 6.283185307179586)
}
