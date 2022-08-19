package bootstrapping

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v3/ckks"
	"github.com/tuneinsight/lattigo/v3/ckks/advanced"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// ParametersLiteral is a struct to parameterize the bootstrapping parameters.
// It contains mandatory and optional fields (with default values)
//
// ================
// Mandatory fields
// ================
//
// C2SLogScale: the scaling factor and distribution of the moduli for the Coeffs To Slots (homomorphic encoding) step.
//              This is a double slice where the first dimension is the index of the prime to be used, and the second
//              dimension the scaling factors to be used: [level][scaling].
//              For example: [][]int{{45}, {46}, {47}} means that the Coeffs To Slots step will use three levels, each
//              with one prime. Primes are consumed in reverse order, so in this example the first matrix will use the
//              prime of 47 bits, the second the prime of 46 bits, and so on.
//              Non standard parameterization can include multiple scaling factors for a same prime, for example
//              [][]int{{30}, {30, 30}} will use two levels for three matrices. The first two matrices will consume a prime
//              of 30 + 30 bits, and have a scaling factor which prime^(1/2), and the third matrix will consume the second
//              prime of 30 bits.
//
// S2CLogScale: the scaling factor and distribution of the moduli for the Coeffs To Slots (homomorphic encoding) step.
//              Parameterization is identical to C2SLogScale.
//
// EvalModLogScale: the scaling factor used during the EvalMod step (all primes will have this bit-size).
//
// =====================================
// Optional fields (with default values)
// =====================================
//
// The default value of the following optional fields result in a depth 8 EvalMod and ensure a failure probability of 2^{-138.7} for 2^{15} slots.
//
// EphemeralSecretWeight: the Hamming weight of the ephemeral secret, by default set to 32, which ensure over 128-bit security
//                        for an evaluation key of modulus 121 bits. The user can set this value to -1 to use the regular
//                        bootstrapping circuit without the ephemeral secret encapsulation. Be aware that doing so will impact
//                        the security, precision, and failure probability of the bootstrapping circuit.
//                        See https://eprint.iacr.org/2022/024 for more informations.
//
// LogMsgRatio: the log of expected ratio Q[0]/|m|, by default set to 8 (ratio of 256.0). This field zero value is obtained by setting it to -1.
//
// SineType: the type of approximation for the modular reduction polynomial. By default set to advanced.CosDiscret.
//
// K: the range of the approximation interval, by default set to 16. This field zero value is obtained by setting it to -1.
//
// SineDeg: the degree of the polynomial approximation of the modular reduction polynomial. By default set to 30.
//
// DoubleAngle: the number of double angle evaluation. By default set to 3. This field zero value is obtained by setting it to -1.
//
// ArcSineDeg: the degree of the ArcSine Taylor polynomial, by default set to 0.
type ParametersLiteral struct {

	// Mandatory fields:
	C2SLogScale     [][]int
	S2CLogScale     [][]int
	EvalModLogScale int

	// Optional fields:
	EphemeralSecretWeight int               // Default: 32
	SineType              advanced.SineType // Default: advanced.CosDiscret
	LogMsgRatio           int               // Default: 8
	K                     int               // Default: 16
	SineDeg               int               // Default: 30
	DoubleAngle           int               // Default: 3
	ArcSineDeg            int               // Default: 0
}

func NewParametersFromLiteral(paramsCKKS ckks.ParametersLiteral, paramsBootstrap ParametersLiteral) (ckks.ParametersLiteral, Parameters, error) {

	// Slots To Coeffs params
	S2CScalingFactor := make([][]float64, len(paramsBootstrap.S2CLogScale))
	for i := range S2CScalingFactor {
		S2CScalingFactor[i] = make([]float64, len(paramsBootstrap.S2CLogScale[i]))
	}

	S2CParams := advanced.EncodingMatrixLiteral{
		LinearTransformType: advanced.SlotsToCoeffs,
		RepackImag2Real:     true,
		LevelStart:          len(paramsCKKS.LogQ) - 1 + len(paramsBootstrap.S2CLogScale),
		BSGSRatio:           2.0,
		ScalingFactor:       S2CScalingFactor,
	}

	// Eval Mod Params
	if paramsBootstrap.EvalModLogScale > 60 {
		return ckks.ParametersLiteral{}, Parameters{}, fmt.Errorf("NewParametersFromLiteral: EvalModLogScale cannot be greater than 60")
	}

	EvalModParams := advanced.EvalModLiteral{
		ScalingFactor: float64(uint64(1 << paramsBootstrap.EvalModLogScale)),
		SineType:      paramsBootstrap.SineType,
		ArcSineDeg:    paramsBootstrap.ArcSineDeg,
	}

	if paramsBootstrap.LogMsgRatio == 0 {
		EvalModParams.MessageRatio = 256.0
	} else if paramsBootstrap.LogMsgRatio == -1 {
		EvalModParams.MessageRatio = 1.0
	} else {
		EvalModParams.MessageRatio = float64(uint64(1 << paramsBootstrap.LogMsgRatio))
	}

	if paramsBootstrap.K == 0 {
		EvalModParams.K = 16
	} else if paramsBootstrap.K == -1 {
		EvalModParams.K = 0
	} else {
		EvalModParams.K = paramsBootstrap.K
	}

	if paramsBootstrap.DoubleAngle == 0 {
		EvalModParams.DoubleAngle = 3
	} else if paramsBootstrap.DoubleAngle == -1 {
		EvalModParams.DoubleAngle = 0
	} else {
		EvalModParams.DoubleAngle = paramsBootstrap.DoubleAngle
	}

	if paramsBootstrap.SineDeg == 0 {
		EvalModParams.SineDeg = 30
	} else {
		EvalModParams.SineDeg = paramsBootstrap.SineDeg
	}

	var EphemeralSecretWeight int
	if paramsBootstrap.EphemeralSecretWeight == 0 {
		EphemeralSecretWeight = 32
	} else if paramsBootstrap.EphemeralSecretWeight == -1 {
		EphemeralSecretWeight = 0
	} else {
		EphemeralSecretWeight = paramsBootstrap.EphemeralSecretWeight
	}

	// Coeffs To Slots params
	EvalModParams.LevelStart = S2CParams.LevelStart + EvalModParams.Depth()

	C2SScalingFactor := make([][]float64, len(paramsBootstrap.C2SLogScale))
	for i := range C2SScalingFactor {
		C2SScalingFactor[i] = make([]float64, len(paramsBootstrap.C2SLogScale[i]))
	}

	C2SParams := advanced.EncodingMatrixLiteral{
		LinearTransformType: advanced.CoeffsToSlots,
		RepackImag2Real:     true,
		LevelStart:          EvalModParams.LevelStart + len(paramsBootstrap.C2SLogScale),
		BSGSRatio:           2.0,
		ScalingFactor:       C2SScalingFactor,
	}

	LogQ := make([]int, len(paramsCKKS.LogQ))
	copy(LogQ, paramsCKKS.LogQ)

	for i := range paramsBootstrap.S2CLogScale {
		var qi int
		for j := range paramsBootstrap.S2CLogScale[i] {
			qi += paramsBootstrap.S2CLogScale[i][j]
		}

		if qi+paramsCKKS.LogDefaultScale < 61 {
			qi += paramsCKKS.LogDefaultScale
		}

		LogQ = append(LogQ, qi)
	}

	for i := 0; i < EvalModParams.Depth(); i++ {
		LogQ = append(LogQ, paramsBootstrap.EvalModLogScale)
	}

	for i := range paramsBootstrap.C2SLogScale {
		var qi int
		for j := range paramsBootstrap.C2SLogScale[i] {
			qi += paramsBootstrap.C2SLogScale[i][j]
		}
		LogQ = append(LogQ, qi)
	}

	LogP := make([]int, len(paramsCKKS.LogP))
	copy(LogP, paramsCKKS.LogP)

	return ckks.ParametersLiteral{
			LogN:            paramsCKKS.LogN,
			LogQ:            LogQ,
			LogP:            LogP,
			LogSlots:        paramsCKKS.LogSlots,
			LogDefaultScale: paramsCKKS.LogDefaultScale,
			Sigma:           paramsCKKS.Sigma,
			H:               paramsCKKS.H,
		},
		Parameters{
			EphemeralSecretWeight:   EphemeralSecretWeight,
			SlotsToCoeffsParameters: S2CParams,
			EvalModParameters:       EvalModParams,
			CoeffsToSlotsParameters: C2SParams,
		}, nil
}

// Parameters is a struct for the default bootstrapping parameters
type Parameters struct {
	SlotsToCoeffsParameters advanced.EncodingMatrixLiteral
	EvalModParameters       advanced.EvalModLiteral
	CoeffsToSlotsParameters advanced.EncodingMatrixLiteral
	EphemeralSecretWeight   int // Hamming weight of the ephemeral secret. If 0, no ephemeral secret is used during the bootstrapping.
}

// Depth returns the depth of the bootstrapping circuit.
func (p *Parameters) Depth() (depth int) {
	return p.SlotsToCoeffsParameters.Depth(true) + p.EvalModParameters.Depth() + p.CoeffsToSlotsParameters.Depth(true)
}

// MarshalBinary encode the target Parameters on a slice of bytes.
func (p *Parameters) MarshalBinary() (data []byte, err error) {
	data = []byte{}

	var tmp []byte
	if tmp, err = p.SlotsToCoeffsParameters.MarshalBinary(); err != nil {
		return nil, err
	}

	data = append(data, uint8(len(tmp)))
	data = append(data, tmp...)

	if tmp, err = p.EvalModParameters.MarshalBinary(); err != nil {
		return nil, err
	}

	data = append(data, uint8(len(tmp)))
	data = append(data, tmp...)

	if tmp, err = p.CoeffsToSlotsParameters.MarshalBinary(); err != nil {
		return nil, err
	}

	data = append(data, uint8(len(tmp)))
	data = append(data, tmp...)

	tmp = make([]byte, 4)
	tmp[0] = uint8(p.EphemeralSecretWeight >> 24)
	tmp[1] = uint8(p.EphemeralSecretWeight >> 16)
	tmp[2] = uint8(p.EphemeralSecretWeight >> 8)
	tmp[3] = uint8(p.EphemeralSecretWeight >> 0)
	data = append(data, tmp...)
	return
}

// UnmarshalBinary decodes a slice of bytes on the target Parameters.
func (p *Parameters) UnmarshalBinary(data []byte) (err error) {

	pt := 0
	dLen := int(data[pt])

	if err := p.SlotsToCoeffsParameters.UnmarshalBinary(data[pt+1 : pt+dLen+1]); err != nil {
		return err
	}

	pt += dLen
	pt++
	dLen = int(data[pt])

	if err := p.EvalModParameters.UnmarshalBinary(data[pt+1 : pt+dLen+1]); err != nil {
		return err
	}

	pt += dLen
	pt++
	dLen = int(data[pt])

	if err := p.CoeffsToSlotsParameters.UnmarshalBinary(data[pt+1 : pt+dLen+1]); err != nil {
		return err
	}

	pt += dLen
	pt++

	p.EphemeralSecretWeight = int(data[pt])<<24 | int(data[pt+1])<<16 | int(data[pt+2])<<8 | int(data[pt+3])

	return
}

// RotationsForBootstrapping returns the list of rotations performed during the Bootstrapping operation.
func (p *Parameters) RotationsForBootstrapping(params ckks.Parameters) (rotations []int) {

	logN := params.LogN()
	logSlots := params.LogSlots()

	// List of the rotation key values to needed for the bootstrapp
	rotations = []int{}

	//SubSum rotation needed X -> Y^slots rotations
	for i := logSlots; i < logN-1; i++ {
		if !utils.IsInSliceInt(1<<i, rotations) {
			rotations = append(rotations, 1<<i)
		}
	}

	p.CoeffsToSlotsParameters.LogN = logN
	p.SlotsToCoeffsParameters.LogN = logN

	p.CoeffsToSlotsParameters.LogSlots = logSlots
	p.SlotsToCoeffsParameters.LogSlots = logSlots

	rotations = append(rotations, p.CoeffsToSlotsParameters.Rotations()...)
	rotations = append(rotations, p.SlotsToCoeffsParameters.Rotations()...)

	return
}
