package ckks

import (
	//"fmt"
	"math"
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

type chebyshevinterpolation struct {
	coeffs []complex128
	a      complex128
	b      complex128
}

// Approximate computes a Chebyshev approximation of the input function, for the tange [-a, b] of degree degree.
// To be used in conjonction with the function EvaluateCheby.
func Approximate(function func(complex128) complex128, a, b complex128, degree int) (cheby *chebyshevinterpolation) {

	cheby = new(chebyshevinterpolation)
	cheby.coeffs = make([]complex128, degree+1)
	cheby.a = a
	cheby.b = b

	u := chebyshevNodesU(degree, -1, 1)
	x := chebyshevNodesX(u, a, b)

	y := make([]complex128, len(x))
	for i := range x {
		y[i] = function(x[i])
	}

	cheby.coeffs = chebyCoeffs(u, y, a, b)

	return
}

// Given a hash table with the first three evaluations of the Chebyshev ring at x in the interval a, b:
// C0 = 1 (actually not stored in the hash table)
// C1 = (2*x - a - b)/(b-a)
// C2 = 2*C1*C1 - C0
// Evaluates the nth degree Chebyshev ring in a recursive manner, storing intermediate results in the hashtable.
// Consumes at most ceil(sqrt(n)) levels for an evaluation at Cn.
// Uses the following property : for a given Chebyshev ring Cn = 2*Ca*Cb - Cc, n = a+b and c = abs(a-b)
func evaluateCheby(n uint64, C map[uint64]*Ciphertext, evaluator *Evaluator, evakey *EvaluationKey) (err error) {

	if C[n] == nil {

		// Computes the index required to compute the asked ring evaluation
		a := uint64(math.Ceil(float64(n) / 2))
		b := n >> 1
		c := uint64(math.Abs(float64(a) - float64(b)))

		// Recurses on the given indexes
		if err = evaluateCheby(a, C, evaluator, evakey); err != nil {
			return err
		}
		if err = evaluateCheby(b, C, evaluator, evakey); err != nil {
			return err
		}

		// Since C[0] is not stored (but rather seen as the constant 1), only recurses on c if c!= 0
		if c != 0 {
			evaluateCheby(c, C, evaluator, evakey)
		}

		// Computes C[n] = C[a]*C[b]
		if C[n], err = evaluator.MulRelinNew(C[a].Element(), C[b].Element(), evakey); err != nil {
			return err
		}

		if err = evaluator.Rescale(C[n], C[n]); err != nil {
			return err
		}

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

// EvaluateCheby evaluates a chebyshev approximation in log(n) + 1 (+1 if 2/(b-a) is not a gaussian integer) levels.
func (evaluator *Evaluator) EvaluateCheby(ct *Ciphertext, cheby *chebyshevinterpolation, evakey *EvaluationKey) (res *Ciphertext, err error) {

	a := cheby.a
	b := cheby.b
	ChebyCoeffs := cheby.coeffs

	C := make(map[uint64]*Ciphertext)

	// C0 = 1, so we treat it as a constant
	// Computes C1 and C2 which are required for the rest of the recursive computation of the Chebyshev ring

	C[1] = ct.CopyNew().Ciphertext()
	evaluator.MultConst(C[1], 2/(b-a), C[1])
	evaluator.AddConst(C[1], (-a-b)/(b-a), C[1])

	if C[1].Scale() > ct.Scale() {
		evaluator.Rescale(C[1], C[1])
	}

	C[2], _ = evaluator.MulRelinNew(C[1].Element(), C[1].Element(), evakey)

	if err = evaluator.Rescale(C[2], C[2]); err != nil {
		return nil, err
	}

	evaluator.Add(C[2], C[2], C[2])
	evaluator.AddConst(C[2], -1, C[2])

	res = C[1].CopyNew().Ciphertext() // res = C[1]

	if err = evaluator.MultConst(res, ChebyCoeffs[1], res); err != nil { // res = A[1]*C[1]
		return nil, err
	}

	if err = evaluator.AddConst(res, ChebyCoeffs[0], res); err != nil { // res = A[0] + A[1]*C[1]
		return nil, err
	}

	for i := uint64(2); i < uint64(len(ChebyCoeffs)); i++ {

		// Evaluates the C[i] Chebyshev ring
		if err = evaluateCheby(i, C, evaluator, evakey); err != nil {
			return nil, err
		}

		if err = evaluator.MultByConstAndAdd(C[i], ChebyCoeffs[i], res); err != nil { // res = A[0] + A[1]*C[1] + ... + A[i]*C[i]
			return nil, err
		}

	}

	// We only rescale at the end to save computation
	if err = evaluator.Rescale(res, res); err != nil {
		return nil, err
	}

	return res, nil

}

func exp2pi(x complex128) complex128 {
	return cmplx.Exp(2 * 3.141592653589793 * complex(0, 1) * x)
}

func sin2pi2pi(x complex128) complex128 {
	return cmplx.Sin(2*3.141592653589793*x) / (2 * 3.141592653589793)
}
