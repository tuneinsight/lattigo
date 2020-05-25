package ckks

import (
	//"fmt"
	"math"
	"math/bits"
)

type poly struct {
	maxDeg uint64
	coeffs []complex128
	lead   bool
}

func (p *poly) degree() uint64 {
	return uint64(len(p.coeffs) - 1)
}

func optimalSplit(logDegree uint64) (logSplit uint64) {
	logSplit = logDegree >> 1
	a := (1 << logSplit) + (1 << (logDegree - logSplit)) + logDegree - logSplit - 3
	b := (1 << (logSplit + 1)) + (1 << (logDegree - logSplit - 1)) + logDegree - logSplit - 4
	if a > b {
		logSplit++
	}

	return
}

// EvaluateChebyFast evaluates the input Chebyshev polynomial with the input ciphertext.
// Faster than EvaluateChebyEco but consumes ceil(log(deg)) + 2 levels.
func (eval *evaluator) EvaluateCheby(op *Ciphertext, cheby *ChebyshevInterpolation, evakey *EvaluationKey) (opOut *Ciphertext) {

	C := make(map[uint64]*Ciphertext)

	C[1] = op.CopyNew().Ciphertext()

	eval.MultByConst(C[1], 2/(cheby.b-cheby.a), C[1])
	eval.AddConst(C[1], (-cheby.a-cheby.b)/(cheby.b-cheby.a), C[1])
	eval.Rescale(C[1], eval.ckksContext.scale, C[1])

	return eval.evalCheby(cheby, C, evakey)
}

// EvaluateChebyFastSpecial evaluates the input Chebyshev polynomial with the input ciphertext.
// Slower than EvaluateChebyFast but consumes ceil(log(deg)) + 1 levels.
func (eval *evaluator) EvaluateChebySpecial(ct *Ciphertext, n complex128, cheby *ChebyshevInterpolation, evakey *EvaluationKey) (res *Ciphertext) {

	C := make(map[uint64]*Ciphertext)

	C[1] = ct.CopyNew().Ciphertext()

	eval.MultByConst(C[1], 2/((cheby.b-cheby.a)*n), C[1])
	eval.AddConst(C[1], (-cheby.a-cheby.b)/(cheby.b-cheby.a), C[1])
	eval.Rescale(C[1], eval.params.Scale, C[1])

	return eval.evalCheby(cheby, C, evakey)
}

func (eval *evaluator) evalCheby(cheby *ChebyshevInterpolation, C map[uint64]*Ciphertext, evakey *EvaluationKey) (res *Ciphertext) {

	logDegree := uint64(bits.Len64(cheby.degree()))
	logSplit := (logDegree >> 1) //optimalSplit(logDegree) //

	for i := uint64(2); i < (1 << logSplit); i++ {
		computePowerBasisCheby(i, C, eval, evakey)
	}

	for i := logSplit; i < logDegree; i++ {
		computePowerBasisCheby(1<<i, C, eval, evakey)
	}

	/*
		for i := range C {
			fmt.Println(i, C[i].Level(), C[i].Scale())
		}
	*/

	return recurseCheby(eval.ckksContext.scale, logSplit, logDegree, cheby.Poly(), C, eval, evakey)

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
		//fmt.Println("Mul", C[a].Level(), C[b].Level())
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

func computeSmallPoly(split uint64, coeffs *poly) (polyList []*poly) {

	if coeffs.degree() < (1 << split) {
		return []*poly{coeffs}
	}

	var nextPower = uint64(1 << split)
	for nextPower < (coeffs.degree()>>1)+1 {
		nextPower <<= 1
	}

	coeffsq, coeffsr := splitCoeffsCheby(coeffs, nextPower)

	a := computeSmallPoly(split, coeffsq)
	b := computeSmallPoly(split, coeffsr)

	return append(a, b...)
}

func splitCoeffsCheby(coeffs *poly, split uint64) (coeffsq, coeffsr *poly) {

	// Splits a Chebyshev polynomial p such that p = q*C^degree + r, where q and r are a linear combination of a Chebyshev basis.
	coeffsr = new(poly)
	coeffsr.coeffs = make([]complex128, split)
	if coeffs.maxDeg == coeffs.degree() {
		coeffsr.maxDeg = split - 1
	} else {
		coeffsr.maxDeg = coeffs.maxDeg - (coeffs.degree() - split + 1)
	}

	for i := uint64(0); i < split; i++ {
		coeffsr.coeffs[i] = coeffs.coeffs[i]
	}

	coeffsq = new(poly)
	coeffsq.coeffs = make([]complex128, coeffs.degree()-split+1)
	coeffsq.maxDeg = coeffs.maxDeg

	coeffsq.coeffs[0] = coeffs.coeffs[split]
	for i, j := split+1, uint64(1); i < coeffs.degree()+1; i, j = i+1, j+1 {
		coeffsq.coeffs[i-split] = 2 * coeffs.coeffs[i]
		coeffsr.coeffs[split-j] -= coeffs.coeffs[i]
	}

	if coeffs.lead {
		coeffsq.lead = true
	}

	return coeffsq, coeffsr
}

func recurseCheby(target_scale float64, logSplit, logDegree uint64, coeffs *poly, C map[uint64]*Ciphertext, evaluator *evaluator, evakey *EvaluationKey) (res *Ciphertext) {

	// Recursively computes the evalution of the Chebyshev polynomial using a baby-set giant-step algorithm.
	if coeffs.degree() < (1 << logSplit) {

		if coeffs.lead && logSplit > 1 && coeffs.degree() > 1<<(logSplit-1) {

			logDegree = uint64(bits.Len64(coeffs.degree()))
			logSplit = optimalSplit(logDegree)

			return recurseCheby(target_scale, logSplit, logDegree, coeffs, C, evaluator, evakey)
		}

		return evaluatePolyFromChebyBasis(target_scale, coeffs, C, evaluator)
	}

	var nextPower = uint64(1 << logSplit)
	for nextPower < (coeffs.degree()>>1)+1 {
		nextPower <<= 1
	}

	coeffsq, coeffsr := splitCoeffsCheby(coeffs, nextPower)

	level := C[nextPower].Level() - 1

	if coeffsq.maxDeg >= 1<<(logDegree-1) && coeffsq.lead {
		level++
	}

	current_qi := float64(evaluator.params.Qi[level])

	//fmt.Printf("X^%2d: %d %d %t %d\n", nextPower, coeffsq.maxDeg, coeffsr.maxDeg, coeffsq.maxDeg >= 1<<(logDegree-1), level)
	//fmt.Printf("X^%2d: %f %f\n", nextPower, target_scale, target_scale* current_qi / C[nextPower].Scale())
	//fmt.Printf("X^%2d : qi %d %t %d %d\n", nextPower, level, coeffsq.lead, coeffsq.maxDeg, 1<<(logDegree-1))
	//fmt.Println()

	res = recurseCheby(target_scale*current_qi/C[nextPower].Scale(), logSplit, logDegree, coeffsq, C, evaluator, evakey)

	var tmp *Ciphertext
	tmp = recurseCheby(target_scale, logSplit, logDegree, coeffsr, C, evaluator, evakey)

	if res.Level() > tmp.Level() {
		for res.Level() != tmp.Level()+1 {
			evaluator.DropLevel(res, 1)
		}
	}

	//fmt.Printf("X^%2d: (%d %f -> ", nextPower, res.Level(), res.Scale())
	evaluator.MulRelin(res, C[nextPower], evakey, res)

	if res.Level() > tmp.Level() {
		evaluator.Rescale(res, evaluator.ckksContext.scale, res)
		//fmt.Printf("%f = %d) + (%d %f) = ", res.Scale(), res.Level(), tmp.Level(), tmp.Scale())
		evaluator.Add(res, tmp, res)
		//fmt.Printf("(%d %f) %f\n", res.Level(), res.Scale(), res.Scale()-tmp.Scale())
	} else {
		evaluator.Add(res, tmp, res)
		evaluator.Rescale(res, evaluator.ckksContext.scale, res)
	}

	return

}

func evaluatePolyFromChebyBasis(target_scale float64, coeffs *poly, C map[uint64]*Ciphertext, evaluator *evaluator) (res *Ciphertext) {

	//fmt.Println(coeffs)

	if coeffs.degree() == 0 {

		res = NewCiphertext(evaluator.params, 1, C[1].Level(), target_scale)

		if math.Abs(real(coeffs.coeffs[0])) > 1e-14 || math.Abs(imag(coeffs.coeffs[0])) > 1e-14 {
			evaluator.AddConst(res, coeffs.coeffs[0], res)
		}

		return
	}

	current_qi := float64(evaluator.params.Qi[C[coeffs.degree()].Level()])

	ct_scale := target_scale * current_qi

	//fmt.Printf("%d %f\n", coeffs.maxDeg, target_scale)
	//fmt.Println("current Qi", evaluator.params.Qi[C[coeffs.degree()].Level()])
	//fmt.Println(coeffs.degree(), C[coeffs.degree()].Level())

	res = NewCiphertext(evaluator.params, 1, C[coeffs.degree()].Level(), ct_scale)

	if math.Abs(real(coeffs.coeffs[0])) > 1e-14 || math.Abs(imag(coeffs.coeffs[0])) > 1e-14 {
		evaluator.AddConst(res, coeffs.coeffs[0], res)
	}

	for key := coeffs.degree(); key > 0; key-- {

		if key != 0 && (math.Abs(real(coeffs.coeffs[key])) > 1e-14 || math.Abs(imag(coeffs.coeffs[key])) > 1e-14) {

			// Target scale * rescale-scale / power basis scale
			const_scale := target_scale * current_qi / C[key].Scale()

			cReal := int64(real(coeffs.coeffs[key]) * const_scale)
			cImag := int64(imag(coeffs.coeffs[key]) * const_scale)

			evaluator.multByGaussianIntegerAndAdd(C[key], cReal, cImag, res)
		}
	}

	evaluator.Rescale(res, evaluator.ckksContext.scale, res)

	return
}

/*
func evaluateInLine(target_scale float64, L uint64, coeffs *poly, C map[uint64]*Ciphertext, evaluator *evaluator, evakey *EvaluationKey) *Ciphertext {

	cf := computeSmallPoly(L, coeffs)

	lvl := C[1].Level() - 2
	fmt.Println(lvl, C[1].Level())

	for i := lvl; i > lvl-5; i-- {
		fmt.Println(i, evaluator.params.Qi[i])
	}
	fmt.Println()

	for i := range C {
		fmt.Printf("C[%d] : %d %f\n", i, C[i].Level(), C[i].Scale())
	}
	fmt.Println()

	target_scale03 := target_scale
	res03 := evaluatePolyFromChebyBasis(target_scale03, cf[9], C, evaluator)

	target_scale07 := target_scale03 * float64(evaluator.params.Qi[lvl-1]) / C[4].Scale()
	res07 := evaluatePolyFromChebyBasis(target_scale07, cf[8], C, evaluator)

	target_scale11 := target_scale * float64(evaluator.params.Qi[lvl-2]) / C[8].Scale()
	res11 := evaluatePolyFromChebyBasis(target_scale11, cf[7], C, evaluator)

	target_scale15 := target_scale11 * float64(evaluator.params.Qi[lvl-1]) / C[4].Scale()
	res15 := evaluatePolyFromChebyBasis(target_scale15, cf[6], C, evaluator)

	target_scale19 := target_scale * float64(evaluator.params.Qi[lvl-3]) / C[16].Scale()
	res19 := evaluatePolyFromChebyBasis(target_scale19, cf[5], C, evaluator)

	target_scale23 := target_scale19 * float64(evaluator.params.Qi[lvl-1]) / C[4].Scale()
	res23 := evaluatePolyFromChebyBasis(target_scale23, cf[4], C, evaluator)

	target_scale27 := target_scale19 * float64(evaluator.params.Qi[lvl-2]) / C[8].Scale()
	res27 := evaluatePolyFromChebyBasis(target_scale27, cf[3], C, evaluator)

	target_scale31 := target_scale27 * float64(evaluator.params.Qi[lvl-1]) / C[4].Scale()
	res31 := evaluatePolyFromChebyBasis(target_scale31, cf[2], C, evaluator)

	target_scale35 := target_scale * float64(evaluator.params.Qi[lvl-3]) / C[32].Scale()
	res35 := evaluatePolyFromChebyBasis(target_scale35, cf[1], C, evaluator)

	target_scale38 := target_scale35 * float64(evaluator.params.Qi[lvl]) / C[4].Scale()
	res38 := evaluatePolyFromChebyBasis(target_scale38, cf[0], C, evaluator)

	fmt.Printf("res03 : %d %f\n", res03.Level(), res03.Scale())
	fmt.Printf("res07 : %d %f\n", res07.Level(), res07.Scale())
	fmt.Printf("res11 : %d %f\n", res11.Level(), res11.Scale())
	fmt.Printf("res15 : %d %f\n", res15.Level(), res15.Scale())
	fmt.Printf("res19 : %d %f\n", res19.Level(), res19.Scale())
	fmt.Printf("res23 : %d %f\n", res23.Level(), res23.Scale())
	fmt.Printf("res27 : %d %f\n", res27.Level(), res27.Scale())
	fmt.Printf("res31 : %d %f\n", res31.Level(), res31.Scale())
	fmt.Printf("res35 : %d %f\n", res35.Level(), res35.Scale())
	fmt.Printf("res38 : %d %f\n", res38.Level(), res38.Scale())
	fmt.Println()

	// 3 + 7*C[4]
	evaluator.MulRelin(res07, C[4], evakey, res07)
	evaluator.Rescale(res07, evaluator.ckksContext.scale, res07)
	fmt.Printf("03 + 07 : %d %d %f %f\n", res03.Level(), res07.Level(), res03.Scale(), res07.Scale())
	evaluator.Add(res07, res03, res07)

	// 11 + 15*C[4]
	evaluator.MulRelin(res15, C[4], evakey, res15)
	evaluator.Rescale(res15, evaluator.ckksContext.scale, res15)
	fmt.Printf("11 + 15 : %d %d %f %f\n", res11.Level(), res15.Level(), res11.Scale(), res15.Scale())
	evaluator.Add(res15, res11, res15)

	// 3 + 7*C[4] + (11 + 15*C[4])*C[8]
	evaluator.MulRelin(res15, C[8], evakey, res15)
	evaluator.Rescale(res15, evaluator.ckksContext.scale, res15)
	fmt.Printf("07 + 15 : %d %d %f %f\n", res07.Level(), res15.Level(), res07.Scale(), res15.Scale())
	evaluator.Add(res15, res07, res15)
	fmt.Println()

	// 19 + 23*C[4]
	evaluator.MulRelin(res23, C[4], evakey, res23)
	evaluator.Rescale(res23, evaluator.ckksContext.scale, res23)
	fmt.Printf("19 + 23 : %d %d %f %f\n", res19.Level(), res23.Level(), res19.Scale(), res23.Scale())
	evaluator.Add(res23, res19, res23)

	// 27 + 31*C[4]
	evaluator.MulRelin(res31, C[4], evakey, res31)
	evaluator.Rescale(res31, evaluator.ckksContext.scale, res31)
	fmt.Printf("27 + 31 : %d %d %f %f\n", res27.Level(), res31.Level(), res27.Scale(), res31.Scale())
	evaluator.Add(res31, res27, res31)

	// 19 + 23*C[4] + (27 + 31*C[4])*C[8]
	evaluator.MulRelin(res31, C[8], evakey, res31)
	evaluator.Rescale(res31, evaluator.ckksContext.scale, res31)
	fmt.Printf("23 + 31 : %d %d %f %f\n", res23.Level(), res31.Level(), res23.Scale(), res31.Scale())
	evaluator.Add(res31, res23, res31)
	fmt.Println()

	// 3 + 7*C[4] + (11 + 15*C[4])*C[8] + (19 + 23*C[4] + (27 + 31*C[4])*C[8])*C[16]
	evaluator.MulRelin(res31, C[16], evakey, res31)
	evaluator.Rescale(res31, evaluator.ckksContext.scale, res31)
	fmt.Printf("15 + 31 : %d %d %f %f\n", res15.Level(), res31.Level(), res15.Scale(), res31.Scale())
	evaluator.Add(res31, res15, res31)

	// 35 + 38*C[4]
	evaluator.MulRelin(res38, C[4], evakey, res38)
	evaluator.Rescale(res38, evaluator.ckksContext.scale, res38)
	fmt.Printf("35 + 38 : %d %d %f %f\n", res35.Level(), res38.Level(), res35.Scale(), res38.Scale())
	evaluator.Add(res38, res35, res38)

	// 3 + 7*C[4] + (11 + 15*C[4])*C[8] + (19 + 23*C[4] + (27 + 31*C[4])*C[8])*C[16] + (35 + 38*C[4])*C[32]
	evaluator.MulRelin(res38, C[32], evakey, res38)
	evaluator.Rescale(res38, evaluator.ckksContext.scale, res38)
	fmt.Printf("31 + 38 : %d %d %f %f\n", res31.Level(), res38.Level(), res31.Scale(), res38.Scale())
	evaluator.Add(res38, res31, res38)

	return res38
}
*/
