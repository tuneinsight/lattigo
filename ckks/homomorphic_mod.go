package ckks

import (
	"encoding/json"
	"math"
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/ckks/cosine"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// SineType is the type of function used during the bootstrapping
// for the homomorphic modular reduction
type SineType uint64

func sin2pi(x *bignum.Complex) (y *bignum.Complex) {
	y = bignum.NewComplex().Set(x)
	y[0].Mul(y[0], new(big.Float).SetFloat64(2))
	y[0].Mul(y[0], bignum.Pi(x.Prec()))
	y[0] = bignum.Sin(y[0])
	return
}

func cos2pi(x *bignum.Complex) (y *bignum.Complex) {
	y = bignum.NewComplex().Set(x)
	y[0].Mul(y[0], new(big.Float).SetFloat64(2))
	y[0].Mul(y[0], bignum.Pi(x.Prec()))
	y[0] = bignum.Cos(y[0])
	return y
}

// Sin and Cos are the two proposed functions for SineType.
// These trigonometric functions offer a good approximation of the function x mod 1 when the values are close to the origin.
const (
	CosDiscrete   = SineType(0) // Special approximation (Han and Ki) of pow((1/2pi), 1/2^r) * cos(2pi(x-0.25)/2^r); this method requires a minimum degree of 2*(K-1).
	SinContinuous = SineType(1) // Standard Chebyshev approximation of (1/2pi) * sin(2pix) on the full interval
	CosContinuous = SineType(2) // Standard Chebyshev approximation of pow((1/2pi), 1/2^r) * cos(2pi(x-0.25)/2^r) on the full interval
)

// EvalModLiteral a struct for the parameters of the EvalMod procedure.
// The EvalMod procedure goal is to homomorphically evaluate a modular reduction by Q[0] (the first prime of the moduli chain) on the encrypted plaintext.
// This struct is consumed by `NewEvalModPolyFromLiteral` to generate the `EvalModPoly` struct, which notably stores
// the coefficient of the polynomial approximating the function x mod Q[0].
type EvalModLiteral struct {
	LevelStart        int      // Starting level of EvalMod
	LogPlaintextScale int      // Log2 of the scaling factor used during EvalMod
	SineType          SineType // Chose between [Sin(2*pi*x)] or [cos(2*pi*x/r) with double angle formula]
	LogMessageRatio   int      // Log2 of the ratio between Q0 and m, i.e. Q[0]/|m|
	K                 int      // K parameter (interpolation in the range -K to K)
	SineDegree        int      // Degree of the interpolation
	DoubleAngle       int      // Number of rescale and double angle formula (only applies for cos and is ignored if sin is used)
	ArcSineDegree     int      // Degree of the Taylor arcsine composed with f(2*pi*x) (if zero then not used)
}

// MarshalBinary returns a JSON representation of the the target EvalModLiteral struct on a slice of bytes.
// See `Marshal` from the `encoding/json` package.
func (evm EvalModLiteral) MarshalBinary() (data []byte, err error) {
	return json.Marshal(evm)
}

// UnmarshalBinary reads a JSON representation on the target EvalModLiteral struct.
// See `Unmarshal` from the `encoding/json` package.
func (evm *EvalModLiteral) UnmarshalBinary(data []byte) (err error) {
	return json.Unmarshal(data, evm)
}

// EvalModPoly is a struct storing the parameters and polynomials approximating the function x mod Q[0] (the first prime of the moduli chain).
type EvalModPoly struct {
	levelStart        int
	LogPlaintextScale int
	sineType          SineType
	LogMessageRatio   int
	doubleAngle       int
	qDiff             float64
	scFac             float64
	sqrt2Pi           float64
	sinePoly          bignum.Polynomial
	arcSinePoly       *bignum.Polynomial
	k                 float64
}

// LevelStart returns the starting level of the EvalMod.
func (evp EvalModPoly) LevelStart() int {
	return evp.levelStart
}

// ScalingFactor returns scaling factor used during the EvalMod.
func (evp EvalModPoly) ScalingFactor() rlwe.Scale {
	return rlwe.NewScale(math.Exp2(float64(evp.LogPlaintextScale)))
}

// ScFac returns 1/2^r where r is the number of double angle evaluation.
func (evp EvalModPoly) ScFac() float64 {
	return evp.scFac
}

// MessageRatio returns the pre-set ratio Q[0]/|m|.
func (evp EvalModPoly) MessageRatio() float64 {
	return float64(uint(1 << evp.LogMessageRatio))
}

// K return the sine approximation range.
func (evp EvalModPoly) K() float64 {
	return evp.k * evp.scFac
}

// QDiff return Q[0]/ClosetPow2
// This is the error introduced by the approximate division by Q[0].
func (evp EvalModPoly) QDiff() float64 {
	return evp.qDiff
}

// NewEvalModPolyFromLiteral generates an EvalModPoly struct from the EvalModLiteral struct.
// The EvalModPoly struct is used by the `EvalModNew` method from the `Evaluator`, which
// homomorphically evaluates x mod Q[0] (the first prime of the moduli chain) on the ciphertext.
func NewEvalModPolyFromLiteral(params Parameters, evm EvalModLiteral) EvalModPoly {

	var arcSinePoly *bignum.Polynomial
	var sinePoly bignum.Polynomial
	var sqrt2pi float64

	doubleAngle := evm.DoubleAngle
	if evm.SineType == SinContinuous {
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

	switch evm.SineType {
	case SinContinuous:

		sinePoly = bignum.ChebyshevApproximation(sin2pi, bignum.Interval{
			Nodes: evm.SineDegree,
			A:     *new(big.Float).SetPrec(cosine.PlaintextPrecision).SetFloat64(-K),
			B:     *new(big.Float).SetPrec(cosine.PlaintextPrecision).SetFloat64(K),
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
			A:     *new(big.Float).SetPrec(cosine.PlaintextPrecision).SetFloat64(-K),
			B:     *new(big.Float).SetPrec(cosine.PlaintextPrecision).SetFloat64(K),
		})
		sinePoly.IsOdd = false

		for i := range sinePoly.Coeffs {
			if i&1 == 1 {
				sinePoly.Coeffs[i] = nil
			}
		}

	default:
		panic("invalid SineType")
	}

	sqrt2piBig := new(big.Float).SetFloat64(sqrt2pi)
	for i := range sinePoly.Coeffs {
		if sinePoly.Coeffs[i] != nil {
			sinePoly.Coeffs[i][0].Mul(sinePoly.Coeffs[i][0], sqrt2piBig)
			sinePoly.Coeffs[i][1].Mul(sinePoly.Coeffs[i][1], sqrt2piBig)
		}
	}

	return EvalModPoly{
		levelStart:        evm.LevelStart,
		LogPlaintextScale: evm.LogPlaintextScale,
		sineType:          evm.SineType,
		LogMessageRatio:   evm.LogMessageRatio,
		doubleAngle:       doubleAngle,
		qDiff:             qDiff,
		scFac:             scFac,
		sqrt2Pi:           sqrt2pi,
		arcSinePoly:       arcSinePoly,
		sinePoly:          sinePoly,
		k:                 K,
	}
}

// Depth returns the depth of the SineEval.
func (evm EvalModLiteral) Depth() (depth int) {

	if evm.SineType == CosDiscrete { // this method requires a minimum degree of 2*K-1.
		depth += int(bits.Len64(uint64(utils.Max(evm.SineDegree, 2*evm.K-1))))
	} else {
		depth += int(bits.Len64(uint64(evm.SineDegree)))
	}

	if evm.SineType != SinContinuous {
		depth += evm.DoubleAngle
	}

	depth += int(bits.Len64(uint64(evm.ArcSineDegree)))
	return depth
}

// EvalModNew applies a homomorphic mod Q on a vector scaled by Delta, scaled down to mod 1 :
//
//  1. Delta * (Q/Delta * I(X) + m(X)) (Delta = scaling factor, I(X) integer poly, m(X) message)
//  2. Delta * (I(X) + Delta/Q * m(X)) (divide by Q/Delta)
//  3. Delta * (Delta/Q * m(X)) (x mod 1)
//  4. Delta * (m(X)) (multiply back by Q/Delta)
//
// Since Q is not a power of two, but Delta is, then does an approximate division by the closest
// power of two to Q instead. Hence, it assumes that the input plaintext is already scaled by
// the correcting factor Q/2^{round(log(Q))}.
//
// !! Assumes that the input is normalized by 1/K for K the range of the approximation.
//
// Scaling back error correction by 2^{round(log(Q))}/Q afterward is included in the polynomial
func (eval Evaluator) EvalModNew(ct *rlwe.Ciphertext, evalModPoly EvalModPoly) *rlwe.Ciphertext {

	if ct.Level() < evalModPoly.LevelStart() {
		panic("ct.Level() < evalModPoly.LevelStart")
	}

	if ct.Level() > evalModPoly.LevelStart() {
		eval.DropLevel(ct, ct.Level()-evalModPoly.LevelStart())
	}

	// Stores default scales
	prevScaleCt := ct.PlaintextScale

	// Normalize the modular reduction to mod by 1 (division by Q)
	ct.PlaintextScale = evalModPoly.ScalingFactor()

	var err error

	// Compute the scales that the ciphertext should have before the double angle
	// formula such that after it it has the scale it had before the polynomial
	// evaluation

	Qi := eval.Parameters().Q()

	targetScale := ct.PlaintextScale
	for i := 0; i < evalModPoly.doubleAngle; i++ {
		targetScale = targetScale.Mul(rlwe.NewScale(Qi[evalModPoly.levelStart-evalModPoly.sinePoly.Depth()-evalModPoly.doubleAngle+i+1]))
		targetScale.Value.Sqrt(&targetScale.Value)
	}

	// Division by 1/2^r and change of variable for the Chebyshev evaluation
	if evalModPoly.sineType == CosDiscrete || evalModPoly.sineType == CosContinuous {
		offset := new(big.Float).Sub(&evalModPoly.sinePoly.B, &evalModPoly.sinePoly.A)
		offset.Mul(offset, new(big.Float).SetFloat64(evalModPoly.scFac))
		offset.Quo(new(big.Float).SetFloat64(-0.5), offset)
		eval.Add(ct, offset, ct)
	}

	// Chebyshev evaluation
	if ct, err = eval.Polynomial(ct, evalModPoly.sinePoly, rlwe.NewScale(targetScale)); err != nil {
		panic(err)
	}

	// Double angle
	sqrt2pi := evalModPoly.sqrt2Pi
	for i := 0; i < evalModPoly.doubleAngle; i++ {
		sqrt2pi *= sqrt2pi
		eval.MulRelin(ct, ct, ct)
		eval.Add(ct, ct, ct)
		eval.Add(ct, -sqrt2pi, ct)
		if err := eval.Rescale(ct, rlwe.NewScale(targetScale), ct); err != nil {
			panic(err)
		}
	}

	// ArcSine
	if evalModPoly.arcSinePoly != nil {
		if ct, err = eval.Polynomial(ct, *evalModPoly.arcSinePoly, ct.PlaintextScale); err != nil {
			panic(err)
		}
	}

	// Multiplies back by q
	ct.PlaintextScale = prevScaleCt
	return ct
}
