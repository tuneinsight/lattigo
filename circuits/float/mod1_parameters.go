package float

import (
	"encoding/json"
	"fmt"
	"math"
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/circuits/float/cosine"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// Mod1Type is the type of function/approximation used to evaluate x mod 1.
type Mod1Type uint64

// Sin and Cos are the two proposed functions for Mod1Type.
// These trigonometric functions offer a good approximation of the function x mod 1 when the values are close to the origin.
const (
	CosDiscrete   = Mod1Type(0) // Special approximation (Han and Ki) of pow((1/2pi), 1/2^r) * cos(2pi(x-0.25)/2^r); this method requires a minimum degree of 2*(K-1).
	SinContinuous = Mod1Type(1) // Standard Chebyshev approximation of (1/2pi) * sin(2pix) on the full interval
	CosContinuous = Mod1Type(2) // Standard Chebyshev approximation of pow((1/2pi), 1/2^r) * cos(2pi(x-0.25)/2^r) on the full interval
)

// Mod1ParametersLiteral a struct for the parameters of the mod 1 procedure.
// The x mod 1 procedure goal is to homomorphically evaluate a modular reduction by Q[0] (the first prime of the moduli chain) on the encrypted plaintext.
// This struct is consumed by `NewMod1ParametersLiteralFromLiteral` to generate the `Mod1ParametersLiteral` struct, which notably stores
// the coefficient of the polynomial approximating the function x mod Q[0].
type Mod1ParametersLiteral struct {
	LevelStart      int      // Starting level of x mod 1
	LogScale        int      // Log2 of the scaling factor used during x mod 1
	Mod1Type        Mod1Type // Chose between [Sin(2*pi*x)] or [cos(2*pi*x/r) with double angle formula]
	LogMessageRatio int      // Log2 of the ratio between Q0 and m, i.e. Q[0]/|m|
	K               int      // K parameter (interpolation in the range -K to K)
	SineDegree      int      // Degree of the interpolation
	DoubleAngle     int      // Number of rescale and double angle formula (only applies for cos and is ignored if sin is used)
	ArcSineDegree   int      // Degree of the Taylor arcsine composed with f(2*pi*x) (if zero then not used)
}

// MarshalBinary returns a JSON representation of the the target Mod1ParametersLiteral struct on a slice of bytes.
// See `Marshal` from the `encoding/json` package.
func (evm Mod1ParametersLiteral) MarshalBinary() (data []byte, err error) {
	return json.Marshal(evm)
}

// UnmarshalBinary reads a JSON representation on the target Mod1ParametersLiteral struct.
// See `Unmarshal` from the `encoding/json` package.
func (evm *Mod1ParametersLiteral) UnmarshalBinary(data []byte) (err error) {
	return json.Unmarshal(data, evm)
}

// Depth returns the depth required to evaluate x mod 1.
func (evm Mod1ParametersLiteral) Depth() (depth int) {

	if evm.Mod1Type == CosDiscrete { // this method requires a minimum degree of 2*K-1.
		depth += int(bits.Len64(uint64(utils.Max(evm.SineDegree, 2*evm.K-1))))
	} else {
		depth += int(bits.Len64(uint64(evm.SineDegree)))
	}

	if evm.Mod1Type != SinContinuous {
		depth += evm.DoubleAngle
	}

	depth += int(bits.Len64(uint64(evm.ArcSineDegree)))
	return depth
}

// Mod1Parameters is a struct storing the parameters and polynomials approximating the function x mod Q[0] (the first prime of the moduli chain).
type Mod1Parameters struct {
	levelStart      int
	LogDefaultScale int
	Mod1Type        Mod1Type
	LogMessageRatio int
	doubleAngle     int
	qDiff           float64
	scFac           float64
	sqrt2Pi         float64
	sinePoly        bignum.Polynomial
	arcSinePoly     *bignum.Polynomial
	k               float64
}

// LevelStart returns the starting level of the x mod 1.
func (evp Mod1Parameters) LevelStart() int {
	return evp.levelStart
}

// ScalingFactor returns scaling factor used during the x mod 1.
func (evp Mod1Parameters) ScalingFactor() rlwe.Scale {
	return rlwe.NewScale(math.Exp2(float64(evp.LogDefaultScale)))
}

// ScFac returns 1/2^r where r is the number of double angle evaluation.
func (evp Mod1Parameters) ScFac() float64 {
	return evp.scFac
}

// MessageRatio returns the pre-set ratio Q[0]/|m|.
func (evp Mod1Parameters) MessageRatio() float64 {
	return float64(uint(1 << evp.LogMessageRatio))
}

// K return the sine approximation range.
func (evp Mod1Parameters) K() float64 {
	return evp.k * evp.scFac
}

// QDiff return Q[0]/ClosetPow2
// This is the error introduced by the approximate division by Q[0].
func (evp Mod1Parameters) QDiff() float64 {
	return evp.qDiff
}

// NewMod1ParametersFromLiteral generates an Mod1Parameters struct from the Mod1ParametersLiteral struct.
// The Mod1Parameters struct is to instantiates a Mod1Evaluator, which homomorphically evaluates x mod 1.
func NewMod1ParametersFromLiteral(params ckks.Parameters, evm Mod1ParametersLiteral) (Mod1Parameters, error) {

	var arcSinePoly *bignum.Polynomial
	var sinePoly bignum.Polynomial
	var sqrt2pi float64

	doubleAngle := evm.DoubleAngle
	if evm.Mod1Type == SinContinuous {
		doubleAngle = 0
	}

	scFac := math.Exp2(float64(doubleAngle))

	K := float64(evm.K) / scFac

	Q := params.Q()[0]
	qDiff := float64(Q) / math.Exp2(math.Round(math.Log2(float64(Q))))

	if evm.ArcSineDegree > 0 {

		sqrt2pi = 1.0

		coeffs := make([]complex128, evm.ArcSineDegree+1)

		coeffs[1] = 0.15915494309189535 * complex(qDiff, 0)

		for i := 3; i < evm.ArcSineDegree+1; i += 2 {
			coeffs[i] = coeffs[i-2] * complex(float64(i*i-4*i+4)/float64(i*i-i), 0)
		}

		p := bignum.NewPolynomial(bignum.Monomial, coeffs, nil)

		arcSinePoly = &p
		arcSinePoly.IsEven = false

		for i := range arcSinePoly.Coeffs {
			if i&1 == 0 {
				arcSinePoly.Coeffs[i] = nil
			}
		}

	} else {
		sqrt2pi = math.Pow(0.15915494309189535*qDiff, 1.0/scFac)
	}

	switch evm.Mod1Type {
	case SinContinuous:

		sinePoly = bignum.ChebyshevApproximation(sin2pi, bignum.Interval{
			Nodes: evm.SineDegree,
			A:     *new(big.Float).SetPrec(cosine.EncodingPrecision).SetFloat64(-K),
			B:     *new(big.Float).SetPrec(cosine.EncodingPrecision).SetFloat64(K),
		})
		sinePoly.IsEven = false

		for i := range sinePoly.Coeffs {
			if i&1 == 0 {
				sinePoly.Coeffs[i] = nil
			}
		}

	case CosDiscrete:
		sinePoly = bignum.NewPolynomial(bignum.Chebyshev, cosine.ApproximateCos(evm.K, evm.SineDegree, float64(uint(1<<evm.LogMessageRatio)), int(evm.DoubleAngle)), [2]float64{-K, K})
		sinePoly.IsOdd = false

		for i := range sinePoly.Coeffs {
			if i&1 == 1 {
				sinePoly.Coeffs[i] = nil
			}
		}

	case CosContinuous:
		sinePoly = bignum.ChebyshevApproximation(cos2pi, bignum.Interval{
			Nodes: evm.SineDegree,
			A:     *new(big.Float).SetPrec(cosine.EncodingPrecision).SetFloat64(-K),
			B:     *new(big.Float).SetPrec(cosine.EncodingPrecision).SetFloat64(K),
		})
		sinePoly.IsOdd = false

		for i := range sinePoly.Coeffs {
			if i&1 == 1 {
				sinePoly.Coeffs[i] = nil
			}
		}

	default:
		return Mod1Parameters{}, fmt.Errorf("invalid Mod1Type")
	}

	sqrt2piBig := new(big.Float).SetFloat64(sqrt2pi)
	for i := range sinePoly.Coeffs {
		if sinePoly.Coeffs[i] != nil {
			sinePoly.Coeffs[i][0].Mul(sinePoly.Coeffs[i][0], sqrt2piBig)
			sinePoly.Coeffs[i][1].Mul(sinePoly.Coeffs[i][1], sqrt2piBig)
		}
	}

	return Mod1Parameters{
		levelStart:      evm.LevelStart,
		LogDefaultScale: evm.LogScale,
		Mod1Type:        evm.Mod1Type,
		LogMessageRatio: evm.LogMessageRatio,
		doubleAngle:     doubleAngle,
		qDiff:           qDiff,
		scFac:           scFac,
		sqrt2Pi:         sqrt2pi,
		arcSinePoly:     arcSinePoly,
		sinePoly:        sinePoly,
		k:               K,
	}, nil
}

func sin2pi(x *big.Float) (y *big.Float) {
	y = new(big.Float).Set(x)
	y.Mul(y, new(big.Float).SetFloat64(2))
	y.Mul(y, bignum.Pi(x.Prec()))
	return bignum.Sin(y)
}

func cos2pi(x *big.Float) (y *big.Float) {
	y = new(big.Float).Set(x)
	y.Mul(y, new(big.Float).SetFloat64(2))
	y.Mul(y, bignum.Pi(x.Prec()))
	y = bignum.Cos(y)
	return y
}
