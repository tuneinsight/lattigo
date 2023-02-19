package bootstrapping

import (
	"encoding/binary"
	"fmt"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ckks/advanced"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// ParametersLiteral is a struct to parameterize the bootstrapping parameters.
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
// C2SLogScale: the scaling factor and distribution of the moduli for the CoeffsToSlots (homomorphic encoding) step.
//
//	 Default value is [][]int{{56}, {56}, {56}, {56}}
//
//		This is a double slice where the first dimension is the index of the prime to be used, and the second
//		dimension the scaling factors to be used: [level][scaling].
//		For example: [][]int{{45}, {46}, {47}} means that the CoeffsToSlots step will use three levels, each
//		with one prime. Primes are consumed in reverse order, so in this example the first matrix will use the
//		prime of 47 bits, the second the prime of 46 bits, and so on.
//		Non standard parameterization can include multiple scaling factors for a same prime, for example
//		[][]int{{30}, {30, 30}} will use two levels for three matrices. The first two matrices will consume a prime
//		of 30 + 30 bits, and have a scaling factor which prime^(1/2), and the third matrix will consume the second
//		prime of 30 bits.
//
// S2CLogScale: the scaling factor and distribution of the moduli for the CoeffsToSlots (homomorphic encoding) step.
//
//		Parameterization is identical to C2SLogScale.
//
//	 Default value is [][]int{{39}, {39}, {39}}
//
// EvalModLogScale: the scaling factor used during the EvalMod step (all primes will have this bit-size).
//
//	Default value is 60.
//
// EphemeralSecretWeight: the Hamming weight of the ephemeral secret, by default set to 32, which ensure over 128-bit security
//
//	for an evaluation key of modulus 121 bits. The user can set this value to -1 to use the regular
//	bootstrapping circuit without the ephemeral secret encapsulation. Be aware that doing so will impact
//	the security, precision, and failure probability of the bootstrapping circuit.
//	See https://eprint.iacr.org/2022/024 for more information.
//
// LogMessageRatio: the log of expected ratio Q[0]/|m|, by default set to 8 (ratio of 256.0). This field zero value is obtained by setting it to -1.
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
	C2SLogScale           [][]int           // Default: [][]int{{56}, {56}, {56}, {56}}
	S2CLogScale           [][]int           // Default: [][]int{{39}, {39}, {39}}
	EvalModLogScale       int               // Default. 60
	EphemeralSecretWeight int               // Default: 32 | set to -1 not use sparse-secret encapsulation
	Iterations            int               // Default:  1
	SineType              advanced.SineType // Default: advanced.CosDiscret
	LogMessageRatio       int               // Default: 8  | set to -1 to set to 0
	K                     int               // Default: 16 | set to -1 to set to 0
	SineDeg               int               // Default: 30
	DoubleAngle           int               // Default: 3  | set to -1 to not use double angle
	ArcSineDeg            int               // Default: 0
}

func (p *ParametersLiteral) MarshalBinarySize() (dataLen int) {
	dataLen++ // Has C2SLogScale
	if p.C2SLogScale != nil {
		dataLen += len(p.C2SLogScale)
		for i := range p.C2SLogScale {
			dataLen++                        // #Levels
			dataLen += len(p.C2SLogScale[i]) // #log
		}
	}

	dataLen++ // Has S2CLogScale
	if p.S2CLogScale != nil {
		dataLen += len(p.S2CLogScale)
		for i := range p.S2CLogScale {
			dataLen++                        // #Levels
			dataLen += len(p.S2CLogScale[i]) // #log
		}
	}

	dataLen++    //EvalModLogScale
	dataLen += 4 // EphemeralSecretWeight
	dataLen++    // Iterations
	dataLen++    // SyneType
	dataLen++    // LogMessageRatio
	dataLen += 2 // K
	dataLen++    //SineDeg
	dataLen++    //DoubleAngle
	dataLen++    //ArcSineDeg

	return
}

// MarshalBinary encodes the target ParametersLiteral on a slice of bytes.
func (p *ParametersLiteral) MarshalBinary() (data []byte, err error) {

	data = make([]byte, p.MarshalBinarySize())

	var ptr int

	if p.C2SLogScale != nil {
		data[ptr] = 1
		ptr++

		data[ptr] = uint8(len(p.C2SLogScale))
		ptr++
		for i := range p.C2SLogScale {
			data[ptr] = uint8(len(p.C2SLogScale[i]))
			ptr++
			for j := range p.C2SLogScale[i] {
				data[ptr] = uint8(p.C2SLogScale[i][j])
				ptr++
			}
		}

	} else {
		ptr++
	}

	if p.S2CLogScale != nil {
		data[ptr] = 1
		ptr++

		data[ptr] = uint8(len(p.S2CLogScale))
		ptr++
		for i := range p.S2CLogScale {
			data[ptr] = uint8(len(p.S2CLogScale[i]))
			ptr++
			for j := range p.S2CLogScale[i] {
				data[ptr] = uint8(p.S2CLogScale[i][j])
				ptr++
			}
		}

	} else {
		ptr++
	}

	data[ptr] = uint8(p.EvalModLogScale)
	ptr++
	binary.LittleEndian.PutUint32(data[ptr:], uint32(p.EphemeralSecretWeight))
	ptr += 4

	data[ptr] = uint8(p.Iterations)
	ptr++

	data[ptr] = uint8(p.SineType)
	ptr++

	data[ptr] = uint8(p.LogMessageRatio)
	ptr++

	binary.LittleEndian.PutUint16(data[ptr:], uint16(p.K))
	ptr += 2

	data[ptr] = uint8(p.SineDeg)
	ptr++

	data[ptr] = uint8(p.DoubleAngle)
	ptr++

	data[ptr] = uint8(p.ArcSineDeg)

	return
}

// UnmarshalBinary decodes a slice of bytes on the target ParametersLiteral.
func (p *ParametersLiteral) UnmarshalBinary(data []byte) (err error) {
	var ptr int

	if data[ptr] == 1 {

		ptr++

		p.C2SLogScale = make([][]int, data[ptr])
		ptr++

		for i := range p.C2SLogScale {

			p.C2SLogScale[i] = make([]int, data[ptr])
			ptr++

			for j := range p.C2SLogScale[i] {
				p.C2SLogScale[i][j] = int(data[ptr])
				ptr++
			}
		}

	} else {
		ptr++
	}

	if data[ptr] == 1 {

		ptr++

		p.S2CLogScale = make([][]int, data[ptr])
		ptr++

		for i := range p.S2CLogScale {

			p.S2CLogScale[i] = make([]int, data[ptr])
			ptr++

			for j := range p.S2CLogScale[i] {
				p.S2CLogScale[i][j] = int(data[ptr])
				ptr++
			}
		}

	} else {
		ptr++
	}

	p.EvalModLogScale = int(data[ptr])
	ptr++

	p.EphemeralSecretWeight = int(binary.LittleEndian.Uint32(data[ptr:]))
	ptr += 4

	if p.EphemeralSecretWeight == 0xFFFFFFFF {
		p.EphemeralSecretWeight = -1
	}

	p.Iterations = int(data[ptr])
	ptr++

	p.SineType = advanced.SineType(data[ptr])
	ptr++

	p.LogMessageRatio = int(data[ptr])
	ptr++

	if p.LogMessageRatio == 0xFF {
		p.LogMessageRatio = -1
	}

	p.K = int(binary.LittleEndian.Uint16(data[ptr:]))
	ptr += 2

	if p.K == 0xFFFF {
		p.K = -1
	}

	p.SineDeg = int(data[ptr])
	ptr++

	p.DoubleAngle = int(data[ptr])
	ptr++

	if p.DoubleAngle == 0xFF {
		p.DoubleAngle = -1
	}

	p.ArcSineDeg = int(data[ptr])

	return
}

// BitComsumption returns the expected consumption in bits of
// bootstrapping circuit of the target ParametersLiteral.
// The value is rounded up and thus will overestimate the value by up to 1 bit.
func (p *ParametersLiteral) BitComsumption() (logQ int) {

	if p.C2SLogScale == nil {
		logQ += 4 * 56
	} else {
		for i := range p.C2SLogScale {
			for _, logQi := range p.C2SLogScale[i] {
				logQ += logQi
			}
		}
	}

	if p.S2CLogScale == nil {
		logQ += 3 * 39
	} else {
		for i := range p.S2CLogScale {
			for _, logQi := range p.S2CLogScale[i] {
				logQ += logQi
			}
		}
	}

	EvalModLogScale := p.EvalModLogScale
	if EvalModLogScale < 1 {
		EvalModLogScale = 60
	}

	SineDeg := p.SineDeg

	if SineDeg < 1 {
		SineDeg = 30
	}

	DoubleAngle := p.DoubleAngle

	if DoubleAngle == -1 {
		DoubleAngle = 0
	} else if DoubleAngle == 0 {
		DoubleAngle = 3
	}

	ArcSineDeg := p.ArcSineDeg

	Iterations := p.Iterations
	if Iterations < 1 {
		Iterations = 1
	}

	logQ += 1 + EvalModLogScale*(bits.Len64(uint64(SineDeg))+DoubleAngle+bits.Len64(uint64(ArcSineDeg))) + (Iterations-1)*25

	return
}

// NewParametersFromLiteral takes as input a ckks.ParametersLiteral and a bootstrapping.ParametersLiteral structs and returns the
// appropriate ckks.ParametersLiteral for the bootstrapping circuit as well as the instantiated bootstrapping.Parameters.
// The returned ckks.ParametersLiteral contains allocated primes.
func NewParametersFromLiteral(paramsCKKS ckks.ParametersLiteral, paramsBootstrap ParametersLiteral) (ckks.ParametersLiteral, Parameters, error) {

	if paramsCKKS.RingType != ring.Standard {
		return ckks.ParametersLiteral{}, Parameters{}, fmt.Errorf("NewParametersFromLiteral: invalid ring.RingType: must be ring.Standard")
	}

	if paramsBootstrap.C2SLogScale == nil {
		paramsBootstrap.C2SLogScale = [][]int{{56}, {56}, {56}, {56}}
	}

	if paramsBootstrap.S2CLogScale == nil {
		paramsBootstrap.S2CLogScale = [][]int{{39}, {39}, {39}}
	}

	if paramsBootstrap.EvalModLogScale == 0 {
		paramsBootstrap.EvalModLogScale = 60
	}

	// Slots To Coeffs params
	S2CLevels := make([]int, len(paramsBootstrap.S2CLogScale))
	for i := range S2CLevels {
		S2CLevels[i] = len(paramsBootstrap.S2CLogScale[i])
	}

	// Number of iterations of the bootstrapping
	// Each iteration consumes an additional prime
	if paramsBootstrap.Iterations < 1 {
		paramsBootstrap.Iterations = 1
	}

	S2CParams := advanced.EncodingMatrixLiteral{
		LinearTransformType: advanced.SlotsToCoeffs,
		RepackImag2Real:     true,
		LevelStart:          len(paramsCKKS.LogQ) - 1 + len(paramsBootstrap.S2CLogScale) + paramsBootstrap.Iterations - 1,
		LogBSGSRatio:        1,
		Levels:              S2CLevels,
	}

	// Eval Mod Params
	if paramsBootstrap.EvalModLogScale > 60 {
		return ckks.ParametersLiteral{}, Parameters{}, fmt.Errorf("NewParametersFromLiteral: EvalModLogScale cannot be greater than 60")
	}

	EvalModParams := advanced.EvalModLiteral{
		LogScale:   paramsBootstrap.EvalModLogScale,
		SineType:   paramsBootstrap.SineType,
		ArcSineDeg: paramsBootstrap.ArcSineDeg,
	}

	if paramsBootstrap.LogMessageRatio == 0 {
		EvalModParams.LogMessageRatio = 8
	} else if paramsBootstrap.LogMessageRatio == -1 {
		EvalModParams.LogMessageRatio = 1
	} else {
		EvalModParams.LogMessageRatio = paramsBootstrap.LogMessageRatio
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

	C2SLevels := make([]int, len(paramsBootstrap.C2SLogScale))
	for i := range C2SLevels {
		C2SLevels[i] = len(paramsBootstrap.C2SLogScale[i])
	}

	C2SParams := advanced.EncodingMatrixLiteral{
		LinearTransformType: advanced.CoeffsToSlots,
		RepackImag2Real:     true,
		LevelStart:          EvalModParams.LevelStart + len(paramsBootstrap.C2SLogScale),
		LogBSGSRatio:        1,
		Levels:              C2SLevels,
	}

	LogQ := make([]int, len(paramsCKKS.LogQ))
	copy(LogQ, paramsCKKS.LogQ)

	for i := 0; i < paramsBootstrap.Iterations-1; i++ {
		LogQ = append(LogQ, 25)
	}

	for i := range paramsBootstrap.S2CLogScale {
		var qi int
		for j := range paramsBootstrap.S2CLogScale[i] {
			qi += paramsBootstrap.S2CLogScale[i][j]
		}

		if qi+paramsCKKS.LogScale < 61 {
			qi += paramsCKKS.LogScale
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

	Q, P, err := rlwe.GenModuli(paramsCKKS.LogN, LogQ, LogP)

	if err != nil {
		return ckks.ParametersLiteral{}, Parameters{}, err
	}

	return ckks.ParametersLiteral{
			LogN:     paramsCKKS.LogN,
			Q:        Q,
			P:        P,
			LogSlots: paramsCKKS.LogSlots,
			LogScale: paramsCKKS.LogScale,
			Sigma:    paramsCKKS.Sigma,
			H:        paramsCKKS.H,
		},
		Parameters{
			EphemeralSecretWeight:   EphemeralSecretWeight,
			SlotsToCoeffsParameters: S2CParams,
			EvalModParameters:       EvalModParams,
			CoeffsToSlotsParameters: C2SParams,
			Iterations:              paramsBootstrap.Iterations,
		}, nil
}

// Parameters is a struct for the default bootstrapping parameters
type Parameters struct {
	SlotsToCoeffsParameters advanced.EncodingMatrixLiteral
	EvalModParameters       advanced.EvalModLiteral
	CoeffsToSlotsParameters advanced.EncodingMatrixLiteral
	Iterations              int
	EphemeralSecretWeight   int // Hamming weight of the ephemeral secret. If 0, no ephemeral secret is used during the bootstrapping.
}

// DepthCoeffsToSlots returns the depth of the Coeffs to Slots of the CKKS bootstrapping.
func (p *Parameters) DepthCoeffsToSlots() (depth int) {
	return p.SlotsToCoeffsParameters.Depth(true)
}

// DepthEvalMod returns the depth of the EvalMod step of the CKKS bootstrapping.
func (p *Parameters) DepthEvalMod() (depth int) {
	return p.EvalModParameters.Depth()
}

// DepthSlotsToCoeffs returns the depth of the Slots to Coeffs step of the CKKS bootstrapping.
func (p *Parameters) DepthSlotsToCoeffs() (depth int) {
	return p.CoeffsToSlotsParameters.Depth(true)
}

// Depth returns the depth of the full bootstrapping circuit.
func (p *Parameters) Depth() (depth int) {
	return p.DepthCoeffsToSlots() + p.DepthEvalMod() + p.DepthSlotsToCoeffs()
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

	data = append(data, uint8(p.Iterations))

	tmp = make([]byte, 4)
	binary.LittleEndian.PutUint32(tmp, uint32(p.EphemeralSecretWeight))
	data = append(data, tmp...)
	return
}

// UnmarshalBinary decodes a slice of bytes on the target Parameters.
func (p *Parameters) UnmarshalBinary(data []byte) (err error) {

	var ptr int
	dLen := int(data[ptr])
	ptr++

	if err := p.SlotsToCoeffsParameters.UnmarshalBinary(data[ptr : ptr+dLen]); err != nil {
		return err
	}

	ptr += dLen

	dLen = int(data[ptr])
	ptr++

	if err := p.EvalModParameters.UnmarshalBinary(data[ptr : ptr+dLen]); err != nil {
		return err
	}

	ptr += dLen

	dLen = int(data[ptr])
	ptr++

	if err := p.CoeffsToSlotsParameters.UnmarshalBinary(data[ptr : ptr+dLen]); err != nil {
		return err
	}

	ptr += dLen

	p.Iterations = int(data[ptr])
	ptr++

	p.EphemeralSecretWeight = int(binary.LittleEndian.Uint32(data[ptr:]))

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
