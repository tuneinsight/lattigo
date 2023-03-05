package bootstrapping

import (
	"encoding/json"
	"fmt"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/ckks/advanced"
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
// CoeffsToSlotsFactorizationDepthAndLogScales: the scaling factor and distribution of the moduli for the SlotsToCoeffs (homomorphic encoding) step.
//
//	Default value is [][]int{min(4, max(LogSlots, 1)) * 56}.
//	This is a double slice where the first dimension is the index of the prime to be used, and the second dimension the scaling factors to be used: [level][scaling].
//	For example: [][]int{{45}, {46}, {47}} means that the CoeffsToSlots step will use three levels, each with one prime. Primes are consumed in reverse order,
//	so in this example the first matrix will use the prime of 47 bits, the second the prime of 46 bits, and so on.
//	Non standard parameterization can include multiple scaling factors for a same prime, for example [][]int{{30}, {30, 30}} will use two levels for three matrices.
//	The first two matrices will consume a prime of 30 + 30 bits, and have a scaling factor which prime^(1/2), and the third matrix will consume the second prime of 30 bits.
//
// SlotsToCoeffsFactorizationDepthAndLogScales: the scaling factor and distribution of the moduli for the CoeffsToSlots (homomorphic decoding) step.
//
//	Parameterization is identical to C2SLogScale. and the default value is [][]int{min(3, max(LogSlots, 1)) * 39}.
//
// EvalModLogScale: the scaling factor used during the EvalMod step (all primes will have this bit-size).
//
//	Default value is 60.
//
// EphemeralSecretWeight: the Hamming weight of the ephemeral secret, by default set to 32, which ensure over 128-bit security for an evaluation key of modulus 121 bits.
//
//	The user can set this value to 0 to use the regular bootstrapping circuit without the ephemeral secret encapsulation.
//	Be aware that doing so will impact the security, precision, and failure probability of the bootstrapping circuit.
//	See https://eprint.iacr.org/2022/024 for more information.
//
// LogMessageRatio: the log of expected ratio Q[0]/|m|, by default set to 8 (ratio of 256.0).
//
//		This ratio directly impacts the precision of the bootstrapping.
//		The homomorphic modular reduction x mod 1 is approximated with by sin(2*pi*x)/(2*pi), which is a good approximation
//		when x is close to the origin. Thus a large message ratio (i.e. 2^8) implies that x is small with respect to Q, and thus close to the origin.
//		When using a small ratio (i.e. 2^4), for example if ct.Scale is close to Q[0] is small or if |m| is large, the ArcSine degree can be set to
//	 a non zero value (i.e. 5 or 7). This will greatly improve the precision of the bootstrapping, at the expense of slightly increasing its depth.
//
// SineType: the type of approximation for the modular reduction polynomial. By default set to advanced.CosDiscrete.
//
// K: the range of the approximation interval, by default set to 16.
//
// SineDeg: the degree of the polynomial approximation of the modular reduction polynomial. By default set to 30.
//
// DoubleAngle: the number of double angle evaluation. By default set to 3.
//
// ArcSineDeg: the degree of the ArcSine Taylor polynomial, by default set to 0.
type ParametersLiteral struct {
	LogSlots                                    *int              // Default: LogN-1
	CoeffsToSlotsFactorizationDepthAndLogScales [][]int           // Default: [][]int{min(4, max(LogSlots, 1)) * 56}
	SlotsToCoeffsFactorizationDepthAndLogScales [][]int           // Default: [][]int{min(3, max(LogSlots, 1)) * 39}
	EvalModLogScale                             *int              // Default: 60
	EphemeralSecretWeight                       *int              // Default: 32
	Iterations                                  *int              // Default: 1
	SineType                                    advanced.SineType // Default: advanced.CosDiscrete
	LogMessageRatio                             *int              // Default: 8
	K                                           *int              // Default: 16
	SineDegree                                  *int              // Default: 30
	DoubleAngle                                 *int              // Default: 3
	ArcSineDegree                               *int              // Default: 0
}

const (
	// DefaultCoeffsToSlotsFactorizationDepth is the default factorization depth CoeffsToSlots step.
	DefaultCoeffsToSlotsFactorizationDepth = 4
	// DefaultSlotsToCoeffsFactorizationDepth is the default factorization depth SlotsToCoeffs step.
	DefaultSlotsToCoeffsFactorizationDepth = 3
	// DefaultCoeffsToSlotsLogScale is the default scaling factors for the CoeffsToSlots step.
	DefaultCoeffsToSlotsLogScale = 56
	// DefaultSlotsToCoeffsLogScale is the default scaling factors for the SlotsToCoeffs step.
	DefaultSlotsToCoeffsLogScale = 39
	// DefaultEvalModLogScale is the default scaling factor for the EvalMod step.
	DefaultEvalModLogScale = 60
	// DefaultEphemeralSecretWeight is the default Hamming weight of the ephemeral secret.
	DefaultEphemeralSecretWeight = 32
	// DefaultIterations is the default number of bootstrapping iterations.
	DefaultIterations = 1
	// DefaultIterationsLogScale is the default scaling factor for the additional prime consumed per additional bootstrapping iteration above 1.
	DefaultIterationsLogScale = 25
	// DefaultSineType is the default function and approximation technique for the homomorphic modular reduction polynomial.
	DefaultSineType = advanced.CosDiscrete
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

// GetCoeffsToSlotsFactorizationDepthAndLogScales returns a copy of the CoeffsToSlotsFactorizationDepthAndLogScales field of the target ParametersLiteral.
// The default value constructed from DefaultC2SFactorization and DefaultC2SLogScale is returned if the field is nil.
func (p *ParametersLiteral) GetCoeffsToSlotsFactorizationDepthAndLogScales(LogSlots int) (CoeffsToSlotsFactorizationDepthAndLogScales [][]int) {
	if p.CoeffsToSlotsFactorizationDepthAndLogScales == nil {
		CoeffsToSlotsFactorizationDepthAndLogScales = make([][]int, utils.MinInt(DefaultCoeffsToSlotsFactorizationDepth, utils.MaxInt(LogSlots, 1)))
		for i := range CoeffsToSlotsFactorizationDepthAndLogScales {
			CoeffsToSlotsFactorizationDepthAndLogScales[i] = []int{DefaultCoeffsToSlotsLogScale}
		}
	} else {
		CoeffsToSlotsFactorizationDepthAndLogScales = p.CoeffsToSlotsFactorizationDepthAndLogScales
	}
	return
}

// GetSlotsToCoeffsFactorizationDepthAndLogScales returns a copy of the SlotsToCoeffsFactorizationDepthAndLogScales field of the target ParametersLiteral.
// The default value constructed from DefaultS2CFactorization and DefaultS2CLogScale is returned if the field is nil.
func (p *ParametersLiteral) GetSlotsToCoeffsFactorizationDepthAndLogScales(LogSlots int) (SlotsToCoeffsFactorizationDepthAndLogScales [][]int) {
	if p.SlotsToCoeffsFactorizationDepthAndLogScales == nil {
		SlotsToCoeffsFactorizationDepthAndLogScales = make([][]int, utils.MinInt(DefaultSlotsToCoeffsFactorizationDepth, utils.MaxInt(LogSlots, 1)))
		for i := range SlotsToCoeffsFactorizationDepthAndLogScales {
			SlotsToCoeffsFactorizationDepthAndLogScales[i] = []int{DefaultSlotsToCoeffsLogScale}
		}
	} else {
		SlotsToCoeffsFactorizationDepthAndLogScales = p.SlotsToCoeffsFactorizationDepthAndLogScales
	}
	return
}

// GetEvalModLogScale returns the EvalModLogScale field of the target ParametersLiteral.
// The default value DefaultEvalModLogScale is returned is the field is nil.
func (p *ParametersLiteral) GetEvalModLogScale() (EvalModLogScale int, err error) {
	if v := p.EvalModLogScale; v == nil {
		EvalModLogScale = DefaultEvalModLogScale

	} else {
		EvalModLogScale = *v

		if EvalModLogScale < 0 || EvalModLogScale > 60 {
			return EvalModLogScale, fmt.Errorf("field EvalModLogScale cannot be smaller than 0 or greater than 60")
		}
	}

	return
}

// GetIterations returns the Iterations field of the target ParametersLiteral.
// The default value DefaultIterations is returned is the field is nil.
func (p *ParametersLiteral) GetIterations() (Iterations int, err error) {
	if v := p.Iterations; v == nil {
		Iterations = DefaultIterations
	} else {
		Iterations = *v

		if Iterations < 1 || Iterations > 2 {
			return Iterations, fmt.Errorf("field Iterations cannot be smaller than 1 or greater than 2")
		}
	}

	return
}

// GetSineType returns the SineType field of the target ParametersLiteral.
// The default value DefaultSineType is returned is the field is nil.
func (p *ParametersLiteral) GetSineType() (SineType advanced.SineType) {
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
		case advanced.SinContinuous:
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
func (p *ParametersLiteral) BitComsumption(LogSlots int) (logQ int, err error) {

	C2SLogScale := p.GetCoeffsToSlotsFactorizationDepthAndLogScales(LogSlots)
	for i := range C2SLogScale {
		for _, logQi := range C2SLogScale[i] {
			logQ += logQi
		}
	}

	S2CLogScale := p.GetSlotsToCoeffsFactorizationDepthAndLogScales(LogSlots)
	for i := range S2CLogScale {
		for _, logQi := range S2CLogScale[i] {
			logQ += logQi
		}
	}

	var SineDegree int
	if SineDegree, err = p.GetSineDegree(); err != nil {
		return
	}

	var EvalModLogScale int
	if EvalModLogScale, err = p.GetEvalModLogScale(); err != nil {
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

	var Iterations int
	if Iterations, err = p.GetIterations(); err != nil {
		return
	}

	logQ += 1 + EvalModLogScale*(bits.Len64(uint64(SineDegree))+DoubleAngle+bits.Len64(uint64(ArcSineDegree))) + (Iterations-1)*DefaultIterationsLogScale

	return
}
