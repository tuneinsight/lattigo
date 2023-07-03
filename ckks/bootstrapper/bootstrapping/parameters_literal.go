package bootstrapping

import (
	"encoding/json"
	"fmt"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// ParametersLiteral is a struct to parameterize the bootstrapping parameters.
// The `ParametersLiteral` struct an unchecked struct that is given to the method `NewParametersFromLiteral` to validate them
// and create the bootstrapping `Parameter` struct, which is used to instantiate a `Bootstrapper`.
// This struct contains only optional fields.
// The default bootstrapping (with no optional field) has
// - Depth 4 for CoeffsToSlots
// - Depth 8 for EvalMod
// - Depth 3 for SlotsToCoeffs
// for a total depth of 15 and a bit consumption of 821
// A precision, for complex values with both real and imaginary parts uniformly distributed in -1, 1 of
// - 27.25 bits for H=192
// - 23.8 bits for H=32768,
// And a failure probability of 2^{-138.7} for 2^{15} slots.
//
// =====================================
// Optional fields (with default values)
// =====================================
//
// LogSlots: the maximum number of slots of the ciphertext. Default value: LogN-1.
//
// CoeffsToSlotsFactorizationDepthAndLogPlaintextScales: the scaling factor and distribution of the moduli for the SlotsToCoeffs (homomorphic encoding) step.
//
//	Default value is [][]int{min(4, max(LogSlots, 1)) * 56}.
//	This is a double slice where the first dimension is the index of the prime to be used, and the second dimension the scaling factors to be used: [level][scaling].
//	For example: [][]int{{45}, {46}, {47}} means that the CoeffsToSlots step will use three levels, each with one prime. Primes are consumed in reverse order,
//	so in this example the first matrix will use the prime of 47 bits, the second the prime of 46 bits, and so on.
//	Non standard parameterization can include multiple scaling factors for a same prime, for example [][]int{{30}, {30, 30}} will use two levels for three matrices.
//	The first two matrices will consume a prime of 30 + 30 bits, and have a scaling factor which prime^(1/2), and the third matrix will consume the second prime of 30 bits.
//
// SlotsToCoeffsFactorizationDepthAndLogPlaintextScales: the scaling factor and distribution of the moduli for the CoeffsToSlots (homomorphic decoding) step.
//
//	Parameterization is identical to C2SLogPlaintextScale. and the default value is [][]int{min(3, max(LogSlots, 1)) * 39}.
//
// EvalModLogPlaintextScale: the scaling factor used during the EvalMod step (all primes will have this bit-size).
//
//	Default value is 60.
//
// EphemeralSecretWeight: the Hamming weight of the ephemeral secret, by default set to 32, which ensure over 128-bit security for an evaluation key of modulus 121 bits.
//
//	The user can set this value to 0 to use the regular bootstrapping circuit without the ephemeral secret encapsulation.
//	Be aware that doing so will impact the security, precision, and failure probability of the bootstrapping circuit.
//	See https://eprint.iacr.org/2022/024 for more information.
//
// Iterations: the number of iterations of the bootstrapping (see https://eprint.iacr.org/2022/1167).
//
// By treating the bootstrapping as a blackbox bootstrapping with precision n, we can construct a bootstrapping of precision ~k*n by iteration.
//
// Here is an example for a two iterations bootstrapping of an input message mod [logq0=55, logq1=45] with scaling factor 2^{90}:
//
// INPUT:
// 1) The input is a ciphertext encrypting [2^{90} * M]_{q0, q1}
// BOOTSTRAP N°1
// 2) Rescale  [M^{90}]_{q0, q1} to [M^{90}/q1]_{q0} (ensure that M^{90}/q1 ~ q0/messageratio)
// 3) Bootsrap [M^{90}/q1]_{q0} to [M^{90}/q1 + e^{90 - logprec}/q1]_{q0, q1, q2, ...}
// 4) Scale up [M^{90}/q1 + e^{90 - logprec}/q1]_{q0, q1, q2, ...} to [M^{d} * M + e^{d - logprec}]_{q0, q1, q2, ...}
// BOOTSTRAP N°2
// 5) Subtract [M^{d}]_{q0, q1} to [M^{d} + e^{d - logprec}]_{q0, q1, q2, ...} to get [e^{d - logprec}]_{q0, q1}
// 6) Scale up [e^{90 - logprec}]_{q0, q1} to [e^{d}]_{q0, q1}
// 7) Rescale  [e^{90}]_{q0, q1} to [{90}/q1]_{q0}
// 8) Bootsrap [e^{90}/q1]_{q0} to [e^{90}/q1 + e'{90 - logprec}/q1]_{q0, q1, q2, ...}
// 9) Scale up [e^{90}/q1 + e'{90 - logprec}/q0]_{q0, q1, q2, ...} by round(q1/2^{logprec}) to get [e^{90-logprec} + e'{90 - 2logprec}/q1]_{q0, q1, q2, ...}
// 10) Subtract [e^{d - logprec} + e'{d - 2logprec}]_{q0, q1, q2, ...} to [M^{d} + e^{d - logprec}]_{q0, q1, q2, ...} to get [M^{d} * M + e'^{d - 2logprec}]_{q0, q1, q2, ...}
// 11) Go back to step 5 for more iterations until 2^{k * logprec} >= 2^{90}
//
// This example can be generalized to input message of any scaling factor and desired output precision by substituting q1 by some product of primes.
//
// Notes:
//   - Although the rescaling of 2) and 7) are approximate, we can ignore them and treat them as being part of the bootstrapping error
//   - As long as round(q1/2^{k*logprec}) >= 2^{logprec}, for k the iteration number, we are guaranteed that the error due to the approximate scale up of step 8) is smaller than 2^{logprec}
//   - The gain in precision for each iteration is proportional to round(q1/2^{k*logprec})
//   - If round(q1/2^{k * logprec}) < 2^{logprec}, where k is the iteration number, then the gain in precision will be less than the expected logprec.
//     This can happen during the last iteration when q1/2^{k * logprec} < 1, and gets rounded to 1 or 0.
//     To solve this issue, we can reduce logprec for the last iterations, but this increases the number of iterations, or reserve a prime of size at least 2^{logprec} to get
//     a proper scaling by q1/2^{k * logprec} (i.e. not a rounded scaling).
//   - If the input ciphertext is at level 0, we must reserve a prime because everything happens within Q[0] and we have no other prime to use for rescaling.
//
// Each subsequent BOOTSTRAP after the first one must be provided with [2]{expected bit precision of the previous iteration, bit-size of the prime to consume for the scale up}
//
// Note that it isn't required to actually provide the expected bit precision of the first iteration as it is usually greater than half the scale
//
// LogMessageRatio: the log of expected ratio Q[0]/|m|, by default set to 8 (ratio of 256.0).
//
//		This ratio directly impacts the precision of the bootstrapping.
//		The homomorphic modular reduction x mod 1 is approximated with by sin(2*pi*x)/(2*pi), which is a good approximation
//		when x is close to the origin. Thus a large message ratio (i.e. 2^8) implies that x is small with respect to Q, and thus close to the origin.
//		When using a small ratio (i.e. 2^4), for example if ct.PlaintextScale is close to Q[0] is small or if |m| is large, the ArcSine degree can be set to
//	 a non zero value (i.e. 5 or 7). This will greatly improve the precision of the bootstrapping, at the expense of slightly increasing its depth.
//
// SineType: the type of approximation for the modular reduction polynomial. By default set to ckks.CosDiscrete.
//
// K: the range of the approximation interval, by default set to 16.
//
// SineDeg: the degree of the polynomial approximation of the modular reduction polynomial. By default set to 30.
//
// DoubleAngle: the number of double angle evaluation. By default set to 3.
//
// ArcSineDeg: the degree of the ArcSine Taylor polynomial, by default set to 0.
type ParametersLiteral struct {
	LogSlots                                             *int                  // Default: LogN-1
	CoeffsToSlotsFactorizationDepthAndLogPlaintextScales [][]int               // Default: [][]int{min(4, max(LogSlots, 1)) * 56}
	SlotsToCoeffsFactorizationDepthAndLogPlaintextScales [][]int               // Default: [][]int{min(3, max(LogSlots, 1)) * 39}
	EvalModLogPlaintextScale                             *int                  // Default: 60
	EphemeralSecretWeight                                *int                  // Default: 32
	IterationsParameters                                 *IterationsParameters // Default: nil (default starting level of 0 and 1 iteration)
	SineType                                             ckks.SineType         // Default: ckks.CosDiscrete
	LogMessageRatio                                      *int                  // Default: 8
	K                                                    *int                  // Default: 16
	SineDegree                                           *int                  // Default: 30
	DoubleAngle                                          *int                  // Default: 3
	ArcSineDegree                                        *int                  // Default: 0
}

const (
	// DefaultCoeffsToSlotsFactorizationDepth is the default factorization depth CoeffsToSlots step.
	DefaultCoeffsToSlotsFactorizationDepth = 4
	// DefaultSlotsToCoeffsFactorizationDepth is the default factorization depth SlotsToCoeffs step.
	DefaultSlotsToCoeffsFactorizationDepth = 3
	// DefaultCoeffsToSlotsLogPlaintextScale is the default scaling factors for the CoeffsToSlots step.
	DefaultCoeffsToSlotsLogPlaintextScale = 56
	// DefaultSlotsToCoeffsLogPlaintextScale is the default scaling factors for the SlotsToCoeffs step.
	DefaultSlotsToCoeffsLogPlaintextScale = 39
	// DefaultEvalModLogPlaintextScale is the default scaling factor for the EvalMod step.
	DefaultEvalModLogPlaintextScale = 60
	// DefaultEphemeralSecretWeight is the default Hamming weight of the ephemeral secret.
	DefaultEphemeralSecretWeight = 32
	// DefaultIterations is the default number of bootstrapping iterations.
	DefaultIterations = 1
	// DefaultSineType is the default function and approximation technique for the homomorphic modular reduction polynomial.
	DefaultSineType = ckks.CosDiscrete
	// DefaultLogMessageRatio is the default ratio between Q[0] and |m|.
	DefaultLogMessageRatio = 8
	// DefaultK is the default interval [-K+1, K-1] for the polynomial approximation of the homomorphic modular reduction.
	DefaultK = 16
	// DefaultSineDeg is the default degree for the polynomial approximation of the homomorphic modular reduction.
	DefaultSineDegree = 30
	// DefaultDoubleAngle is the default number of double iterations for the homomorphic modular reduction.
	DefaultDoubleAngle = 3
	// DefaultArcSineDeg is the default degree of the arcsine polynomial for the homomorphic modular reduction.
	DefaultArcSineDegree = 0
)

type IterationsParameters struct {
	BootstrappingPrecision []float64
	ReservedPrimeBitSize   int
}

// MarshalBinary returns a JSON representation of the the target ParametersLiteral struct on a slice of bytes.
// See `Marshal` from the `encoding/json` package.
func (p *ParametersLiteral) MarshalBinary() (data []byte, err error) {
	return json.Marshal(p)
}

// UnmarshalBinary reads a JSON representation on the target ParametersLiteral struct.
// See `Unmarshal` from the `encoding/json` package.
func (p *ParametersLiteral) UnmarshalBinary(data []byte) (err error) {
	return json.Unmarshal(data, p)
}

// GetLogSlots returns the LogSlots field of the target ParametersLiteral.
// The default value LogN-1 is returned is the field is nil.
func (p *ParametersLiteral) GetLogSlots(LogN int) (LogSlots int, err error) {
	if v := p.LogSlots; v == nil {
		LogSlots = LogN - 1

	} else {
		LogSlots = *v

		if LogSlots < 1 || LogSlots > LogN-1 {
			return LogSlots, fmt.Errorf("field LogSlots cannot be smaller than 1 or greater than LogN-1")
		}
	}

	return
}

// GetCoeffsToSlotsFactorizationDepthAndLogPlaintextScales returns a copy of the CoeffsToSlotsFactorizationDepthAndLogPlaintextScales field of the target ParametersLiteral.
// The default value constructed from DefaultC2SFactorization and DefaultC2SLogPlaintextScale is returned if the field is nil.
func (p *ParametersLiteral) GetCoeffsToSlotsFactorizationDepthAndLogPlaintextScales(LogSlots int) (CoeffsToSlotsFactorizationDepthAndLogPlaintextScales [][]int, err error) {
	if p.CoeffsToSlotsFactorizationDepthAndLogPlaintextScales == nil {
		CoeffsToSlotsFactorizationDepthAndLogPlaintextScales = make([][]int, utils.Min(DefaultCoeffsToSlotsFactorizationDepth, utils.Max(LogSlots, 1)))
		for i := range CoeffsToSlotsFactorizationDepthAndLogPlaintextScales {
			CoeffsToSlotsFactorizationDepthAndLogPlaintextScales[i] = []int{DefaultCoeffsToSlotsLogPlaintextScale}
		}
	} else {
		var depth int
		for _, level := range p.CoeffsToSlotsFactorizationDepthAndLogPlaintextScales {
			for range level {
				depth++
				if depth > LogSlots {
					return nil, fmt.Errorf("field CoeffsToSlotsFactorizationDepthAndLogPlaintextScales cannot contain parameters for a depth > LogSlots")
				}
			}
		}
		CoeffsToSlotsFactorizationDepthAndLogPlaintextScales = p.CoeffsToSlotsFactorizationDepthAndLogPlaintextScales
	}
	return
}

// GetSlotsToCoeffsFactorizationDepthAndLogPlaintextScales returns a copy of the SlotsToCoeffsFactorizationDepthAndLogPlaintextScales field of the target ParametersLiteral.
// The default value constructed from DefaultS2CFactorization and DefaultS2CLogPlaintextScale is returned if the field is nil.
func (p *ParametersLiteral) GetSlotsToCoeffsFactorizationDepthAndLogPlaintextScales(LogSlots int) (SlotsToCoeffsFactorizationDepthAndLogPlaintextScales [][]int, err error) {
	if p.SlotsToCoeffsFactorizationDepthAndLogPlaintextScales == nil {
		SlotsToCoeffsFactorizationDepthAndLogPlaintextScales = make([][]int, utils.Min(DefaultSlotsToCoeffsFactorizationDepth, utils.Max(LogSlots, 1)))
		for i := range SlotsToCoeffsFactorizationDepthAndLogPlaintextScales {
			SlotsToCoeffsFactorizationDepthAndLogPlaintextScales[i] = []int{DefaultSlotsToCoeffsLogPlaintextScale}
		}
	} else {
		var depth int
		for _, level := range p.SlotsToCoeffsFactorizationDepthAndLogPlaintextScales {
			for range level {
				depth++
				if depth > LogSlots {
					return nil, fmt.Errorf("field SlotsToCoeffsFactorizationDepthAndLogPlaintextScales cannot contain parameters for a depth > LogSlots")
				}
			}
		}
		SlotsToCoeffsFactorizationDepthAndLogPlaintextScales = p.SlotsToCoeffsFactorizationDepthAndLogPlaintextScales
	}
	return
}

// GetEvalModLogPlaintextScale returns the EvalModLogPlaintextScale field of the target ParametersLiteral.
// The default value DefaultEvalModLogPlaintextScale is returned is the field is nil.
func (p *ParametersLiteral) GetEvalModLogPlaintextScale() (EvalModLogPlaintextScale int, err error) {
	if v := p.EvalModLogPlaintextScale; v == nil {
		EvalModLogPlaintextScale = DefaultEvalModLogPlaintextScale

	} else {
		EvalModLogPlaintextScale = *v

		if EvalModLogPlaintextScale < 0 || EvalModLogPlaintextScale > 60 {
			return EvalModLogPlaintextScale, fmt.Errorf("field EvalModLogPlaintextScale cannot be smaller than 0 or greater than 60")
		}
	}

	return
}

// GetIterationsParameters returns the IterationsParmaeters field of the target ParametersLiteral.
// The default value is nil.
func (p *ParametersLiteral) GetIterationsParameters() (Iterations *IterationsParameters, err error) {

	if v := p.IterationsParameters; v == nil {
		return nil, nil
	} else {

		if len(v.BootstrappingPrecision) < 1 {
			return nil, fmt.Errorf("field BootstrappingPrecision of IterationsParameters must be greater than 0")
		}

		for _, prec := range v.BootstrappingPrecision {
			if prec == 0 {
				return nil, fmt.Errorf("field BootstrappingPrecision of IterationsParameters cannot be 0")
			}
		}

		if v.ReservedPrimeBitSize > 61 {
			return nil, fmt.Errorf("field ReservedPrimeBitSize of IterationsParameters cannot be larger than 61")
		}

		return v, nil
	}
}

// GetSineType returns the SineType field of the target ParametersLiteral.
// The default value DefaultSineType is returned is the field is nil.
func (p *ParametersLiteral) GetSineType() (SineType ckks.SineType) {
	return p.SineType
}

// GetArcSineDegree returns the ArcSineDegree field of the target ParametersLiteral.
// The default value DefaultArcSineDegree is returned is the field is nil.
func (p *ParametersLiteral) GetArcSineDegree() (ArcSineDegree int, err error) {
	if v := p.ArcSineDegree; v == nil {
		ArcSineDegree = 0
	} else {
		ArcSineDegree = *v

		if ArcSineDegree < 0 {
			return ArcSineDegree, fmt.Errorf("field ArcSineDegree cannot be negative")
		}
	}

	return
}

// GetLogMessageRatio returns the LogMessageRatio field of the target ParametersLiteral.
// The default value DefaultLogMessageRatio is returned is the field is nil.
func (p *ParametersLiteral) GetLogMessageRatio() (LogMessageRatio int, err error) {
	if v := p.LogMessageRatio; v == nil {
		LogMessageRatio = DefaultLogMessageRatio
	} else {
		LogMessageRatio = *v

		if LogMessageRatio < 0 {
			return LogMessageRatio, fmt.Errorf("field LogMessageRatio cannot be negative")
		}
	}

	return
}

// GetK returns the K field of the target ParametersLiteral.
// The default value DefaultK is returned is the field is nil.
func (p *ParametersLiteral) GetK() (K int, err error) {
	if v := p.K; v == nil {
		K = DefaultK
	} else {
		K = *v

		if K < 0 {
			return K, fmt.Errorf("field K cannot be negative")
		}
	}

	return
}

// GetDoubleAngle returns the DoubleAngle field of the target ParametersLiteral.
// The default value DefaultDoubleAngle is returned is the field is nil.
func (p *ParametersLiteral) GetDoubleAngle() (DoubleAngle int, err error) {

	if v := p.DoubleAngle; v == nil {

		switch p.GetSineType() {
		case ckks.SinContinuous:
			DoubleAngle = 0
		default:
			DoubleAngle = DefaultDoubleAngle
		}

	} else {
		DoubleAngle = *v

		if DoubleAngle < 0 {
			return DoubleAngle, fmt.Errorf("field DoubleAngle cannot be negative")
		}
	}
	return
}

// GetSineDegree returns the SineDegree field of the target ParametersLiteral.
// The default value DefaultSineDegree is returned is the field is nil.
func (p *ParametersLiteral) GetSineDegree() (SineDegree int, err error) {
	if v := p.SineDegree; v == nil {
		SineDegree = DefaultSineDegree
	} else {
		SineDegree = *v

		if SineDegree < 0 {
			return SineDegree, fmt.Errorf("field SineDegree cannot be negative")
		}
	}
	return
}

// GetEphemeralSecretWeight returns the EphemeralSecretWeight field of the target ParametersLiteral.
// The default value DefaultEphemeralSecretWeight is returned is the field is nil.
func (p *ParametersLiteral) GetEphemeralSecretWeight() (EphemeralSecretWeight int, err error) {
	if v := p.EphemeralSecretWeight; v == nil {
		EphemeralSecretWeight = DefaultEphemeralSecretWeight
	} else {
		EphemeralSecretWeight = *v

		if EphemeralSecretWeight < 0 {
			return EphemeralSecretWeight, fmt.Errorf("field EphemeralSecretWeight cannot be negative")
		}
	}
	return
}

// BitConsumption returns the expected consumption in bits of
// bootstrapping circuit of the target ParametersLiteral.
// The value is rounded up and thus will overestimate the value by up to 1 bit.
func (p *ParametersLiteral) BitConsumption(LogSlots int) (logQ int, err error) {

	var C2SLogPlaintextScale [][]int
	if C2SLogPlaintextScale, err = p.GetCoeffsToSlotsFactorizationDepthAndLogPlaintextScales(LogSlots); err != nil {
		return
	}

	for i := range C2SLogPlaintextScale {
		for _, logQi := range C2SLogPlaintextScale[i] {
			logQ += logQi
		}
	}

	var S2CLogPlaintextScale [][]int
	if S2CLogPlaintextScale, err = p.GetSlotsToCoeffsFactorizationDepthAndLogPlaintextScales(LogSlots); err != nil {
		return
	}

	for i := range S2CLogPlaintextScale {
		for _, logQi := range S2CLogPlaintextScale[i] {
			logQ += logQi
		}
	}

	var SineDegree int
	if SineDegree, err = p.GetSineDegree(); err != nil {
		return
	}

	var EvalModLogPlaintextScale int
	if EvalModLogPlaintextScale, err = p.GetEvalModLogPlaintextScale(); err != nil {
		return
	}

	var DoubleAngle int
	if DoubleAngle, err = p.GetDoubleAngle(); err != nil {
		return
	}

	var ArcSineDegree int
	if ArcSineDegree, err = p.GetArcSineDegree(); err != nil {
		return
	}

	var Iterations *IterationsParameters
	if Iterations, err = p.GetIterationsParameters(); err != nil {
		return
	}

	var ReservedPrimeBitSize int
	if Iterations != nil {
		ReservedPrimeBitSize = Iterations.ReservedPrimeBitSize
	}

	logQ += 1 + EvalModLogPlaintextScale*(bits.Len64(uint64(SineDegree))+DoubleAngle+bits.Len64(uint64(ArcSineDegree))) + ReservedPrimeBitSize

	return
}
