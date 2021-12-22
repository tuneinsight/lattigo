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
	Cos1 = SineType(1) // Special approximation (Han and Ki) of pow((1/2pi), 1/2^r) * cos(2pi(x-0.25)/2^r), this method requires a minimum degree of 2*(K-1).
	Cos2 = SineType(2) // Standard Chebyshev approximation of pow((1/2pi), 1/2^r) * cos(2pi(x-0.25)/2^r)
)

// EvalModLiteral a struct for the paramters of the EvalMod step
// of the bootstrapping
type EvalModLiteral struct {
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
func (evm *EvalModLiteral) QDiff() float64 {
	return float64(evm.Q) / math.Exp2(math.Round(math.Log2(float64(evm.Q))))
}

// EvalModPoly is a struct storing the EvalModLiteral with
// the polynomials.
type EvalModPoly struct {
	levelStart    int
	scalingFactor float64
	sineType      SineType
	messageRatio  float64
	doubleAngle   int
	qDiff         float64
	scFac         float64
	sqrt2Pi       float64
	sinePoly      *ckks.Polynomial
	arcSinePoly   *ckks.Polynomial
}

// LevelStart returns the starting level of the EvalMod.
func (evp *EvalModPoly) LevelStart() int {
	return evp.levelStart
}

// ScalingFactor returns scaling factor used during the EvalMod.
func (evp *EvalModPoly) ScalingFactor() float64 {
	return evp.scalingFactor
}

// ScFac returns 1/2^r where r is the number of double angle evaluation.
func (evp *EvalModPoly) ScFac() float64 {
	return evp.scFac
}

// MessageRatio returns the pre-set ratio Q[0]/|m|.
func (evp *EvalModPoly) MessageRatio() float64 {
	return evp.messageRatio
}

// A returns the left bound of the sine approximation (scaled by 1/2^r).
func (evp *EvalModPoly) A() float64 {
	return evp.sinePoly.A
}

// B returns the right bound of the sine approximation (scaled by 1/2^r).
func (evp *EvalModPoly) B() float64 {
	return evp.sinePoly.B
}

// K return the sine approximation range.
func (evp *EvalModPoly) K() float64 {
	return evp.sinePoly.B * evp.scFac
}

// QDiff return Q/ClosestedPow2
// This is the error introduced by the approximate division by Q.
func (evp *EvalModPoly) QDiff() float64 {
	return evp.qDiff
}

// NewEvalModPolyFromLiteral generates an EvalModPoly fromt the EvalModLiteral.
func NewEvalModPolyFromLiteral(evm EvalModLiteral) EvalModPoly {

	var arcSinePoly *ckks.Polynomial
	var sinePoly *ckks.Polynomial
	var sqrt2pi float64

	scFac := math.Exp2(float64(evm.DoubleAngle))

	qDiff := evm.QDiff()

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

		sinePoly = ckks.Approximate(sin2pi2pi, -float64(evm.K), float64(evm.K), evm.SineDeg)

	} else if evm.SineType == Cos1 {

		sinePoly = new(ckks.Polynomial)
		sinePoly.Coeffs = ApproximateCos(evm.K, evm.SineDeg, evm.MessageRatio, int(evm.DoubleAngle))
		sinePoly.MaxDeg = sinePoly.Degree()
		sinePoly.A = float64(-evm.K) / scFac
		sinePoly.B = float64(evm.K) / scFac
		sinePoly.Lead = true
		sinePoly.Basis = ckks.ChebyshevBasis

	} else if evm.SineType == Cos2 {
		sinePoly = ckks.Approximate(cos2pi, -float64(evm.K)/scFac, float64(evm.K)/scFac, evm.SineDeg)
	} else {
		panic("invalid SineType")
	}

	for i := range sinePoly.Coeffs {
		sinePoly.Coeffs[i] *= complex(sqrt2pi, 0)
	}

	return EvalModPoly{
		levelStart:    evm.LevelStart,
		scalingFactor: evm.ScalingFactor,
		sineType:      evm.SineType,
		messageRatio:  evm.MessageRatio,
		doubleAngle:   evm.DoubleAngle,
		qDiff:         qDiff,
		scFac:         scFac,
		sqrt2Pi:       sqrt2pi,
		arcSinePoly:   arcSinePoly,
		sinePoly:      sinePoly}
}

// Depth returns the depth of the SineEval. If true, then also
// counts the double angle formula.
func (evm *EvalModLiteral) Depth() int {
	depth := int(math.Ceil(math.Log2(float64(evm.SineDeg + 1))))
	depth += evm.DoubleAngle
	depth += int(math.Ceil(math.Log2(float64(evm.ArcSineDeg + 1))))
	return depth
}
