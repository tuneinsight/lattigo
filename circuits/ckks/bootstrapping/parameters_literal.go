package bootstrapping

import (
	"encoding/json"
	"fmt"
	"math"
	"math/bits"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils"
)

// ParametersLiteral is a struct to parameterize the bootstrapping parameters.
// The ParametersLiteral struct an unchecked struct that is given to the method [NewParametersFromLiteral] to validate them
// and create the bootstrapping [Parameters] struct, which is used to instantiate a Bootstrapper.
// This struct contains only optional fields.
// The default bootstrapping (with no optional field) has
// - Depth 4 for CoeffsToSlots
// - Depth 8 for EvalMod
// - Depth 3 for SlotsToCoeffs
// for a total depth of 15 and a bit consumption of 821
// A precision, for complex values with both real and imaginary parts uniformly distributed in -1, 1 of
// - 27.9 bits (27.4 L2) for H=192
// - 24.4 bits (23.9 L2) for H=32768,
// And a failure probability of 2^{-138.7} for 2^{15} slots.
//
// =====================================
// Optional fields (with default values)
// =====================================
//
// LogN: the log2 of the ring degree of the bootstrapping parameters. The default value is 16.
//
// LogP: the log2 of the auxiliary primes during the key-switching operation of the bootstrapping parameters.
// The default value is [61]*max(1, floor(sqrt(#Qi))).
//
// Xs: the distribution of the secret-key used to generate the bootstrapping evaluation keys.
// The default value is ring.Ternary{H: 192}.
//
// Xe: the distribution of the error sampled to generate the bootstrapping evaluation keys.
// The default value is rlwe.DefaultXe.
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
// IterationsParameters : by treating the bootstrapping as a black box with precision logprec, we can construct a bootstrapping of precision ~k*logprec by iteration (see https://eprint.iacr.org/2022/1167).
//   - BootstrappingPrecision: []float64, the list of iterations (after the initial bootstrapping) given by the expected precision of each previous iteration.
//   - ReservedPrimeBitSize: the size of the reserved prime for the scaling after the initial bootstrapping.
//
// For example: &bootstrapping.IterationsParameters{BootstrappingPrecision: []float64{16}, ReservedPrimeBitSize: 16} will define a two iteration bootstrapping (the first iteration being the initial bootstrapping)
// with a additional prime close to 2^{16} reserved for the scaling of the error during the second iteration.
//
// Here is an example for a two iterations bootstrapping of an input message mod [logq0=55, logq1=45] with scaling factor 2^{90}:
//
// INPUT:
//  1. The input is a ciphertext encrypting [2^{90} * M]_{q0, q1}
//
// ITERATION N°0
//  2. Rescale  [M^{90}]_{q0, q1} to [M^{90}/q1]_{q0} (ensure that M^{90}/q1 ~ q0/message-ratio by additional scaling if necessary)
//  3. Bootstrap [M^{90}/q1]_{q0} to [M^{90}/q1 + e^{90 - logprec}/q1]_{q0, q1, q2, ...}
//  4. Scale up [M^{90}/q1 + e^{90 - logprec}/q1]_{q0, q1, q2, ...} to [M^{d} + e^{d - logprec}]_{q0, q1, q2, ...}
//
// ITERATION N°1
//  5. Subtract [M^{d}]_{q0, q1} to [M^{d} + e^{d - logprec}]_{q0, q1, q2, ...} to get [e^{d - logprec}]_{q0, q1}
//  6. Scale up [e^{90 - logprec}]_{q0, q1} by 2^{logprec} to get [e^{d}]_{q0, q1}
//  7. Rescale  [e^{90}]_{q0, q1} to [{90}/q1]_{q0}
//  8. Bootstrap [e^{90}/q1]_{q0} to [e^{90}/q1 + e'^{90 - logprec}/q1]_{q0, q1, q2, ...}
//  9. Scale up [e^{90}/q1 + e'^{90 - logprec}/q0]_{q0, q1, q2, ...} by round(q1/2^{logprec}) to get [e^{90-logprec} + e'^{90 - 2logprec}]_{q0, q1, q2, ...}
//  10. Subtract [e^{d - logprec} + e'^{d - 2logprec}]_{q0, q1, q2, ...} to [M^{d} + e^{d - logprec}]_{q0, q1, q2, ...} to get [M^{d} + e'^{d - 2logprec}]_{q0, q1, q2, ...}
//  11. Go back to step 5 for more iterations until 2^{k * logprec} >= 2^{90}
//
// This example can be generalized to input messages of any scaling factor and desired output precision by increasing the input scaling factor and substituting q1 by a larger product of primes.
//
// Notes:
//   - The bootstrapping precision cannot exceed the original input ciphertext precision.
//   - Although the rescalings of 2) and 7) are approximate, we can ignore them and treat them as being part of the bootstrapping error
//   - As long as round(q1/2^{k*logprec}) >= 2^{logprec}, for k the iteration number, we are guaranteed that the error due to the approximate scale up of step 8) is smaller than 2^{logprec}
//   - The gain in precision for each iteration is proportional to min(round(q1/2^{k*logprec}), 2^{logprec})
//   - If round(q1/2^{k * logprec}) < 2^{logprec}, where k is the iteration number, then the gain in precision will be less than the expected logprec.
//     This can happen during the last iteration when q1/2^{k * logprec} < 1, and gets rounded to 1 or 0.
//     To solve this issue, we can reduce logprec for the last iterations, but this increases the number of iterations, or reserve a prime of size at least 2^{logprec} to get
//     a proper scaling by q1/2^{k * logprec} (i.e. not a integer rounded scaling).
//   - If the input ciphertext is at level 0, we must reserve a prime because everything happens within Q[0] and we have no other prime to use for rescaling.
//
// LogMessageRatio: the log of expected ratio Q[0]/|m|, by default set to 8 (ratio of 256.0).
//
//		This ratio directly impacts the precision of the bootstrapping.
//		The homomorphic modular reduction x mod 1 is approximated with by sin(2*pi*x)/(2*pi), which is a good approximation
//		when x is close to the origin. Thus a large message ratio (i.e. 2^8) implies that x is small with respect to Q, and thus close to the origin.
//		When using a small ratio (i.e. 2^4), for example if ct.PlaintextScale is close to Q[0] is small or if |m| is large, the Mod1InvDegree can be set to
//	 a non zero value (i.e. 5 or 7). This will greatly improve the precision of the bootstrapping, at the expense of slightly increasing its depth.
//
// Mod1Type: the type of approximation for the modular reduction polynomial. By default set to mod1.CosDiscrete.
//
// K: the range of the approximation interval, by default set to 16.
//
// Mod1Degree: the degree of f: x mod 1. By default set to 30.
//
// DoubleAngle: the number of double angle evaluation. By default set to 3.
//
// Mod1InvDegree: the degree of the f^-1: (x mod 1)^-1, by default set to 0.
type ParametersLiteral struct {
	LogN                                        *int                        // Default: 16
	LogP                                        []int                       // Default: 61 * max(1, floor(sqrt(#Qi)))
	Xs                                          ring.DistributionParameters // Default: ring.Ternary{H: 192}
	Xe                                          ring.DistributionParameters // Default: rlwe.DefaultXe
	LogSlots                                    *int                        // Default: LogN-1
	CoeffsToSlotsFactorizationDepthAndLogScales [][]int                     // Default: [][]int{min(4, max(LogSlots, 1)) * 56}
	SlotsToCoeffsFactorizationDepthAndLogScales [][]int                     // Default: [][]int{min(3, max(LogSlots, 1)) * 39}
	EvalModLogScale                             *int                        // Default: 60
	EphemeralSecretWeight                       *int                        // Default: 32
	IterationsParameters                        *IterationsParameters       // Default: nil (default starting level of 0 and 1 iteration)
	Mod1Type                                    mod1.Type                   // Default: mod1.CosDiscrete
	LogMessageRatio                             *int                        // Default: 8
	K                                           *int                        // Default: 16
	Mod1Degree                                  *int                        // Default: 30
	DoubleAngle                                 *int                        // Default: 3
	Mod1InvDegree                               *int                        // Default: 0
}

type CircuitOrder int

const (
	ModUpThenEncode = CircuitOrder(0) // ScaleDown -> ModUp -> CoeffsToSlots -> EvalMod -> SlotsToCoeffs.
	DecodeThenModUp = CircuitOrder(1) // SlotsToCoeffs -> ScaleDown -> ModUp -> CoeffsToSlots -> EvalMod.
	Custom          = CircuitOrder(2) // Custom order (e.g. partial bootstrapping), disables checks.
)

const (
	// DefaultLogN is the default ring degree for the bootstrapping.
	DefaultLogN = 16
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
	// DefaultMod1Type is the default function and approximation technique for the homomorphic modular reduction polynomial.
	DefaultMod1Type = mod1.CosDiscrete
	// DefaultLogMessageRatio is the default ratio between Q[0] and |m|.
	DefaultLogMessageRatio = 8
	// DefaultK is the default interval [-K+1, K-1] for the polynomial approximation of the homomorphic modular reduction.
	DefaultK = 16
	// DefaultMod1Degree is the default degree for the polynomial approximation of the homomorphic modular reduction.
	DefaultMod1Degree = 30
	// DefaultDoubleAngle is the default number of double iterations for the homomorphic modular reduction.
	DefaultDoubleAngle = 3
	// DefaultMod1InvDegree is the default degree of the f^-1: (x mod 1)^-1 polynomial for the homomorphic modular reduction.
	DefaultMod1InvDegree = 0
)

var (
	// DefaultXs is the default secret distribution of the bootstrapping parameters.
	DefaultXs = ring.Ternary{H: 192}
	// DefaultXe is the default error distribution of the bootstrapping parameters.
	DefaultXe = rlwe.DefaultXe
)

type IterationsParameters struct {
	BootstrappingPrecision []float64
	ReservedPrimeBitSize   int
}

// MarshalBinary returns a JSON representation of the target ParametersLiteral struct on a slice of bytes.
// See Marshal from the [encoding/json] package.
func (p ParametersLiteral) MarshalBinary() (data []byte, err error) {
	return json.Marshal(p)
}

// UnmarshalBinary reads a JSON representation on the target ParametersLiteral struct.
// See Unmarshal from the [encoding/json] package.
func (p *ParametersLiteral) UnmarshalBinary(data []byte) (err error) {
	return json.Unmarshal(data, p)
}

// GetLogN returns the LogN field of the target [ParametersLiteral].
// The default value DefaultLogN is returned if the field is nil.
func (p ParametersLiteral) GetLogN() (LogN int) {
	if v := p.LogN; v == nil {
		LogN = DefaultLogN
	} else {
		LogN = *v
	}

	return
}

// GetDefaultXs returns the Xs field of the target [ParametersLiteral].
// The default value DefaultXs is returned if the field is nil.
func (p ParametersLiteral) GetDefaultXs() (Xs ring.DistributionParameters) {
	if v := p.Xs; v == nil {
		Xs = DefaultXs
	} else {
		Xs = v
	}

	return
}

// GetDefaultXe returns the Xe field of the target [ParametersLiteral].
// The default value DefaultXe is returned if the field is nil.
func (p ParametersLiteral) GetDefaultXe() (Xe ring.DistributionParameters) {
	if v := p.Xe; v == nil {
		Xe = DefaultXe
	} else {
		Xe = v
	}

	return
}

// GetLogP returns the list of bit-size of the primes Pi (extended primes for the key-switching)
// according to the number of #Qi (ciphertext primes).
// The default value is 61 * max(1, floor(sqrt(#Qi))).
func (p ParametersLiteral) GetLogP(NumberOfQi int) (LogP []int) {
	if v := p.LogP; v == nil {
		LogP = make([]int, utils.Max(1, int(math.Sqrt(float64(NumberOfQi)))))
		for i := range LogP {
			LogP[i] = 61
		}
	} else {
		LogP = v
	}

	return
}

// GetLogSlots returns the LogSlots field of the target [ParametersLiteral].
// The default value LogN-1 is returned if the field is nil.
func (p ParametersLiteral) GetLogSlots() (LogSlots int, err error) {

	LogN := p.GetLogN()

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

// GetCoeffsToSlotsFactorizationDepthAndLogScales returns a copy of the CoeffsToSlotsFactorizationDepthAndLogScales field of the target [ParametersLiteral].
// The default value constructed from [DefaultSlotsToCoeffsFactorizationDepth] and [DefaultSlotsToCoeffsLogScale] is returned if the field is nil.
func (p ParametersLiteral) GetCoeffsToSlotsFactorizationDepthAndLogScales(LogSlots int) (CoeffsToSlotsFactorizationDepthAndLogScales [][]int, err error) {
	if p.CoeffsToSlotsFactorizationDepthAndLogScales == nil {
		CoeffsToSlotsFactorizationDepthAndLogScales = make([][]int, utils.Min(DefaultCoeffsToSlotsFactorizationDepth, utils.Max(LogSlots, 1)))
		for i := range CoeffsToSlotsFactorizationDepthAndLogScales {
			CoeffsToSlotsFactorizationDepthAndLogScales[i] = []int{DefaultCoeffsToSlotsLogScale}
		}
	} else {
		var depth int
		for _, level := range p.CoeffsToSlotsFactorizationDepthAndLogScales {
			for range level {
				depth++
				if depth > LogSlots {
					return nil, fmt.Errorf("field CoeffsToSlotsFactorizationDepthAndLogScales cannot contain parameters for a depth > LogSlots")
				}
			}
		}
		CoeffsToSlotsFactorizationDepthAndLogScales = p.CoeffsToSlotsFactorizationDepthAndLogScales
	}
	return
}

// GetSlotsToCoeffsFactorizationDepthAndLogScales returns a copy of the SlotsToCoeffsFactorizationDepthAndLogScales field of the target [ParametersLiteral].
// The default value constructed from [DefaultSlotsToCoeffsFactorizationDepth] and [DefaultSlotsToCoeffsLogScale] is returned if the field is nil.
func (p ParametersLiteral) GetSlotsToCoeffsFactorizationDepthAndLogScales(LogSlots int) (SlotsToCoeffsFactorizationDepthAndLogScales [][]int, err error) {
	if p.SlotsToCoeffsFactorizationDepthAndLogScales == nil {
		SlotsToCoeffsFactorizationDepthAndLogScales = make([][]int, utils.Min(DefaultSlotsToCoeffsFactorizationDepth, utils.Max(LogSlots, 1)))
		for i := range SlotsToCoeffsFactorizationDepthAndLogScales {
			SlotsToCoeffsFactorizationDepthAndLogScales[i] = []int{DefaultSlotsToCoeffsLogScale}
		}
	} else {
		var depth int
		for _, level := range p.SlotsToCoeffsFactorizationDepthAndLogScales {
			for range level {
				depth++
				if depth > LogSlots {
					return nil, fmt.Errorf("field SlotsToCoeffsFactorizationDepthAndLogScales cannot contain parameters for a depth > LogSlots")
				}
			}
		}
		SlotsToCoeffsFactorizationDepthAndLogScales = p.SlotsToCoeffsFactorizationDepthAndLogScales
	}
	return
}

// GetEvalMod1LogScale returns the EvalModLogScale field of the target [ParametersLiteral].
// The default value [DefaultEvalModLogScale] is returned if the field is nil.
func (p ParametersLiteral) GetEvalMod1LogScale() (EvalModLogScale int, err error) {
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

// GetIterationsParameters returns the [IterationsParameters] field of the target [ParametersLiteral].
// The default value is nil.
func (p ParametersLiteral) GetIterationsParameters() (Iterations *IterationsParameters, err error) {

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

// GetLogMessageRatio returns the []LogMessageRatio field of the target [ParametersLiteral].
// The default value [DefaultLogMessageRatio] is returned if the field is nil.
func (p ParametersLiteral) GetLogMessageRatio() (LogMessageRatio int, err error) {
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

// GetK returns the K field of the target [ParametersLiteral].
// The default value [DefaultK] is returned if the field is nil.
func (p ParametersLiteral) GetK() (K int, err error) {
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

// GetMod1Type returns the Mod1Type field of the target ParametersLiteral.
// The default value DefaultMod1Type is returned if the field is nil.
func (p ParametersLiteral) GetMod1Type() (Mod1Type mod1.Type) {
	return p.Mod1Type
}

// GetDoubleAngle returns the DoubleAngle field of the target [ParametersLiteral].
// The default value [DefaultDoubleAngle] is returned if the field is nil.
func (p ParametersLiteral) GetDoubleAngle() (DoubleAngle int, err error) {

	if v := p.DoubleAngle; v == nil {

		switch p.GetMod1Type() {
		case mod1.SinContinuous:
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

// GetMod1Degree returns the Mod1Degree field of the target [ParametersLiteral].
// The default value [DefaultMod1Degree] is returned if the field is nil.
func (p ParametersLiteral) GetMod1Degree() (Mod1Degree int, err error) {
	if v := p.Mod1Degree; v == nil {
		Mod1Degree = DefaultMod1Degree
	} else {
		Mod1Degree = *v

		if Mod1Degree < 0 {
			return Mod1Degree, fmt.Errorf("field Mod1Degree cannot be negative")
		}
	}
	return
}

// GetMod1InvDegree returns the Mod1InvDegree field of the target [ParametersLiteral].
// The default value [DefaultMod1InvDegree] is returned if the field is nil.
func (p ParametersLiteral) GetMod1InvDegree() (Mod1InvDegree int, err error) {
	if v := p.Mod1InvDegree; v == nil {
		Mod1InvDegree = DefaultMod1InvDegree
	} else {
		Mod1InvDegree = *v

		if Mod1InvDegree < 0 {
			return Mod1InvDegree, fmt.Errorf("field Mod1InvDegree cannot be negative")
		}
	}

	return
}

// GetEphemeralSecretWeight returns the EphemeralSecretWeight field of the target [ParametersLiteral].
// The default value [DefaultEphemeralSecretWeight] is returned if the field is nil.
func (p ParametersLiteral) GetEphemeralSecretWeight() (EphemeralSecretWeight int, err error) {
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
// bootstrapping circuit of the target [ParametersLiteral].
// The value is rounded up and thus will overestimate the value by up to 1 bit.
func (p ParametersLiteral) BitConsumption(LogSlots int) (logQ int, err error) {

	var C2SLogPlaintextScale [][]int
	if C2SLogPlaintextScale, err = p.GetCoeffsToSlotsFactorizationDepthAndLogScales(LogSlots); err != nil {
		return
	}

	for i := range C2SLogPlaintextScale {
		for _, logQi := range C2SLogPlaintextScale[i] {
			logQ += logQi
		}
	}

	var S2CLogPlaintextScale [][]int
	if S2CLogPlaintextScale, err = p.GetSlotsToCoeffsFactorizationDepthAndLogScales(LogSlots); err != nil {
		return
	}

	for i := range S2CLogPlaintextScale {
		for _, logQi := range S2CLogPlaintextScale[i] {
			logQ += logQi
		}
	}

	var Mod1Degree int
	if Mod1Degree, err = p.GetMod1Degree(); err != nil {
		return
	}

	var EvalModLogPlaintextScale int
	if EvalModLogPlaintextScale, err = p.GetEvalMod1LogScale(); err != nil {
		return
	}

	var DoubleAngle int
	if DoubleAngle, err = p.GetDoubleAngle(); err != nil {
		return
	}

	var Mod1InvDegree int
	if Mod1InvDegree, err = p.GetMod1InvDegree(); err != nil {
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

	logQ += 1 + EvalModLogPlaintextScale*(bits.Len64(uint64(Mod1Degree))+DoubleAngle+bits.Len64(uint64(Mod1InvDegree))) + ReservedPrimeBitSize

	return
}
