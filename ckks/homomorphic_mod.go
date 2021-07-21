package ckks

import (
	"github.com/ldsec/lattigo/v2/ckks/bettersine"
	"math"
)

// SineType is the type of function used during the bootstrapping
// for the homomorphic modular reduction
type SineType uint64

// Sin and Cos are the two proposed functions for SineType
const (
	Sin  = SineType(0) // Standard Chebyshev approximation of (1/2pi) * sin(2pix)
	Cos1 = SineType(1) // Special approximation (Han and Ki) of pow((1/2pi), 1/2^r) * cos(2pi(x-0.25)/2^r)
	Cos2 = SineType(2) // Standard Chebyshev approximation of pow((1/2pi), 1/2^r) * cos(2pi(x-0.25)/2^r)
)

// EvalModParameters a struct for the paramters of the EvalMod step
// of the bootstrapping
type EvalModParameters struct {
	Q             uint64   // Q0 to reduce by during EvalMod
	LevelStart    int      // Starting level of EvalMod
	ScalingFactor float64  // Scaling factor used during EvalMod
	SineType      SineType // Chose betwenn [Sin(2*pi*x)] or [cos(2*pi*x/r) with double angle formula]
	MessageRatio  float64  // Ratio between Q0 and m, i.e. Q[0]/|m|
	K             int      // K parameter (interpolation in the range -K to K)
	SineDeg       int      // Degree of the interpolation
	DoubleAngle   int      // Number of rescale and double angle formula (only applies for cos)
	ArcSineDeg    int      // Degree of the Taylor arcsine composed with f(2*pi*x) (if zero then not used)
}

// EvalModPoly is a struct storing the EvalModParameters with
// the polynomials.
type EvalModPoly struct {
	EvalModParameters
	ScFac       float64
	Sqrt2Pi     float64
	SinePoly    ChebyshevInterpolation
	ArcSinePoly *Poly
}

// GenPoly generates an EvalModPoly fromt the EvalModParameters.
func (evm *EvalModParameters) GenPoly() EvalModPoly {

	var arcSinePoly *Poly
	var sinePoly *ChebyshevInterpolation
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

		arcSinePoly = NewPoly(coeffs)

	} else {
		sqrt2pi = math.Pow(0.15915494309189535*qDiff, 1.0/scFac)
	}

	if evm.SineType == Sin {

		sinePoly = Approximate(sin2pi2pi, -complex(float64(evm.K)/scFac, 0), complex(float64(evm.K)/scFac, 0), evm.SineDeg)

	} else if evm.SineType == Cos1 {

		sinePoly = new(ChebyshevInterpolation)
		sinePoly.coeffs = bettersine.Approximate(evm.K, evm.SineDeg, evm.MessageRatio, int(evm.DoubleAngle))
		sinePoly.maxDeg = sinePoly.Degree()
		sinePoly.a = complex(float64(-evm.K)/scFac, 0)
		sinePoly.b = complex(float64(evm.K)/scFac, 0)
		sinePoly.lead = true

	} else if evm.SineType == Cos2 {
		sinePoly = Approximate(cos2pi, -complex(float64(evm.K)/scFac, 0), complex(float64(evm.K)/scFac, 0), evm.SineDeg)
	} else {
		panic("Bootstrapper -> invalid sineType")
	}

	for i := range sinePoly.coeffs {
		sinePoly.coeffs[i] *= complex(sqrt2pi, 0)
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

// EvalMod : homomorphic modular reduction by 1 in the range -K to K
//
// Assumes that ct is scaled by 1/(2^r * K) * q/2^{round(log(q))}
//
// (1/(2^r * K)) : Chebyshev change of basis and double angle formula
// (q/2^{round(log(q))}) : correcting factor for q not being a power of two
//
// Scaling back by 2^{round(log(q))}/q afterward is included in the polynomial
func (eval *evaluator) EvalMod(ct *Ciphertext, evalModPoly EvalModPoly) *Ciphertext {

	// Stores default scales
	prevScaleCt := ct.Scale
	prevScaleEval := eval.scale

	// Normalize the modular reduction to mod by 1
	ct.Scale = evalModPoly.ScalingFactor
	eval.scale = evalModPoly.ScalingFactor

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
		eval.AddConst(ct, -0.5/(complex(evalModPoly.ScFac, 0)*(evalModPoly.SinePoly.b-evalModPoly.SinePoly.a)), ct)
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

	ct.Scale = prevScaleCt
	eval.scale = prevScaleEval

	return ct
}
