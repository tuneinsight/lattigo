package advanced

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"math"
	"math/cmplx"
)

// SineType is the type of function used during the bootstrapping
// for the homomorphic modular reduction
type SineType uint64

func sin2pi2pi(x complex128) complex128 {
	return cmplx.Sin(6.283185307179586 * x) // 6.283185307179586
}

func cos2pi(x complex128) complex128 {
	return cmplx.Cos(6.283185307179586 * x)
}

// Sin and Cos are the two proposed functions for SineType
const (
	Sin  = SineType(0) // Standard Chebyshev approximation of (1/2pi) * sin(2pix)
	Cos1 = SineType(1) // Special approximation (Han and Ki) of pow((1/2pi), 1/2^r) * cos(2pi(x-0.25)/2^r)
	Cos2 = SineType(2) // Standard Chebyshev approximation of pow((1/2pi), 1/2^r) * cos(2pi(x-0.25)/2^r)
)

// EvalModParameters a struct for the paramters of the EvalMod step
// of the bootstrapping
type EvalModParameters struct {
	Q             uint64   // Q to reduce by during EvalMod
	LevelStart    int      // Starting level of EvalMod
	ScalingFactor float64  // Scaling factor used during EvalMod
	SineType      SineType // Chose betwenn [Sin(2*pi*x)] or [cos(2*pi*x/r) with double angle formula]
	MessageRatio  float64  // Ratio between Q0 and m, i.e. Q[0]/|m|
	K             int      // K parameter (interpolation in the range -K to K)
	SineDeg       int      // Degree of the interpolation
	DoubleAngle   int      // Number of rescale and double angle formula (only applies for cos)
	ArcSineDeg    int      // Degree of the Taylor arcsine composed with f(2*pi*x) (if zero then not used)
}

// QDiff return Q/ClosestedPow2
// This is the error introduced by the approximate division by Q
func (evm *EvalModParameters) QDiff() float64 {
	return float64(evm.Q) / math.Exp2(math.Round(math.Log2(float64(evm.Q))))
}

// EvalModPoly is a struct storing the EvalModParameters with
// the polynomials.
type EvalModPoly struct {
	EvalModParameters
	ScFac       float64
	Sqrt2Pi     float64
	SinePoly    ckks.ChebyshevInterpolation
	ArcSinePoly *ckks.Poly
}

// GenPoly generates an EvalModPoly fromt the EvalModParameters.
func (evm *EvalModParameters) GenPoly() EvalModPoly {

	var arcSinePoly *ckks.Poly
	var sinePoly *ckks.ChebyshevInterpolation
	var sqrt2pi float64

	scFac := math.Exp2(float64(evm.DoubleAngle))

	qDiff := float64(evm.Q) / math.Exp2(math.Round(math.Log2(float64(evm.Q))))

	if evm.ArcSineDeg > 0 {

		sqrt2pi = 1.0

		coeffs := make([]complex128, evm.ArcSineDeg+1)

		coeffs[1] = 0.15915494309189535 * complex(qDiff, 0)

		for i := 3; i < evm.ArcSineDeg+1; i += 2 {

			coeffs[i] = coeffs[i-2] * complex(float64(i*i-4*i+4)/float64(i*i-i), 0)

		}

		arcSinePoly = ckks.NewPoly(coeffs)

	} else {
		sqrt2pi = math.Pow(0.15915494309189535*qDiff, 1.0/scFac)
	}

	if evm.SineType == Sin {

		if evm.DoubleAngle != 0 {
			panic("cannot user double angle with SineType == Sin")
		}

		sinePoly = ckks.Approximate(sin2pi2pi, -complex(float64(evm.K), 0), complex(float64(evm.K), 0), evm.SineDeg)

	} else if evm.SineType == Cos1 {

		sinePoly = new(ckks.ChebyshevInterpolation)
		sinePoly.Coeffs = ApproximateCos(evm.K, evm.SineDeg, evm.MessageRatio, int(evm.DoubleAngle))
		sinePoly.MaxDeg = sinePoly.Degree()
		sinePoly.A = complex(float64(-evm.K)/scFac, 0)
		sinePoly.B = complex(float64(evm.K)/scFac, 0)
		sinePoly.Lead = true

	} else if evm.SineType == Cos2 {
		sinePoly = ckks.Approximate(cos2pi, -complex(float64(evm.K)/scFac, 0), complex(float64(evm.K)/scFac, 0), evm.SineDeg)
	} else {
		panic("invalid SineType")
	}

	for i := range sinePoly.Coeffs {
		sinePoly.Coeffs[i] *= complex(sqrt2pi, 0)
	}

	return EvalModPoly{EvalModParameters: *evm, ScFac: scFac, Sqrt2Pi: sqrt2pi, ArcSinePoly: arcSinePoly, SinePoly: *sinePoly}
}

// Depth returns the depth of the SineEval. If true, then also
// counts the double angle formula.
func (evm *EvalModParameters) Depth() int {
	depth := int(math.Ceil(math.Log2(float64(evm.SineDeg + 1))))
	depth += evm.DoubleAngle
	depth += int(math.Ceil(math.Log2(float64(evm.ArcSineDeg + 1))))
	return depth
}

// EvalMod does :
//
//	1) Delta * (Q/Delta * I(X) + m(X)) (Delta = scaling factor, I(X) integer poly, m(X) message)
//	2) Delta * (I(X) + Delta/Q * m(X)) (divide by Q/Delta)
//	3) Delta * (Delta/Q * m(X)) (x mod 1)
//	4) Delta * (m(X)) (multiply back by Q/Delta)
//
// Since Q is not a power of two, but Delta is, then does an approximate division by the closest
// power of two to Q instead. Hence, it assumes that the input plaintext is already scaled by
// the correcting factor Q/2^{round(log(Q))}.
//
// !! Assumes that the input is normalized by 1/K for K the range of the approximation.
//
// Scaling back error correction by 2^{round(log(Q))}/Q afterward is included in the polynomial
func (eval *evaluator) EvalMod(ct *ckks.Ciphertext, evalModPoly EvalModPoly) *ckks.Ciphertext {

	// Stores default scales
	prevScaleCt := ct.Scale

	// Normalize the modular reduction to mod by 1 (division by Q)
	ct.Scale = evalModPoly.ScalingFactor

	var err error

	// Compute the scales that the ciphertext should have before the double angle
	// formula such that after it it has the scale it had before the polynomial
	// evaluation
	targetScale := ct.Scale
	for i := 0; i < evalModPoly.DoubleAngle; i++ {
		targetScale = math.Sqrt(targetScale * eval.params.QiFloat64(evalModPoly.LevelStart-evalModPoly.SinePoly.Depth()-evalModPoly.DoubleAngle+i+1))
	}

	// Division by 1/2^r and change of variable for the Chebysehev evaluation
	if evalModPoly.SineType == Cos1 || evalModPoly.SineType == Cos2 {
		eval.AddConst(ct, -0.5/(complex(evalModPoly.ScFac, 0)*(evalModPoly.SinePoly.B-evalModPoly.SinePoly.A)), ct)
	}

	// Chebyshev evaluation
	if ct, err = eval.EvaluateCheby(ct, &evalModPoly.SinePoly, targetScale); err != nil {
		panic(err)
	}

	// Double angle
	sqrt2pi := evalModPoly.Sqrt2Pi
	for i := 0; i < evalModPoly.DoubleAngle; i++ {
		sqrt2pi *= sqrt2pi
		eval.MulRelin(ct, ct, ct)
		eval.Add(ct, ct, ct)
		eval.AddConst(ct, -sqrt2pi, ct)
		if err := eval.Rescale(ct, targetScale, ct); err != nil {
			panic(err)
		}
	}

	// ArcSine
	if evalModPoly.ArcSinePoly != nil {
		if ct, err = eval.EvaluatePoly(ct, evalModPoly.ArcSinePoly, ct.Scale); err != nil {
			panic(err)
		}
	}

	// Multiplies back by q
	ct.Scale = prevScaleCt
	return ct
}
