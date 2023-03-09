package bootstrapping

import (
	"encoding/json"
	"fmt"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ckks/advanced"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// Parameters is a struct for the default bootstrapping parameters
type Parameters struct {
	SlotsToCoeffsParameters advanced.HomomorphicDFTMatrixLiteral
	EvalModParameters       advanced.EvalModLiteral
	CoeffsToSlotsParameters advanced.HomomorphicDFTMatrixLiteral
	Iterations              int
	EphemeralSecretWeight   int // Hamming weight of the ephemeral secret. If 0, no ephemeral secret is used during the bootstrapping.
}

// NewParametersFromLiteral takes as input a ckks.ParametersLiteral and a bootstrapping.ParametersLiteral structs and returns the
// appropriate ckks.ParametersLiteral for the bootstrapping circuit as well as the instantiated bootstrapping.Parameters.
// The returned ckks.ParametersLiteral contains allocated primes.
func NewParametersFromLiteral(paramsCKKS ckks.ParametersLiteral, paramsBootstrap ParametersLiteral) (ckks.ParametersLiteral, Parameters, error) {

	var err error

	if paramsCKKS.RingType != ring.Standard {
		return ckks.ParametersLiteral{}, Parameters{}, fmt.Errorf("NewParametersFromLiteral: invalid ring.RingType: must be ring.Standard")
	}

	CoeffsToSlotsFactorizationDepthAndLogScales := paramsBootstrap.GetCoeffsToSlotsFactorizationDepthAndLogScales()
	SlotsToCoeffsFactorizationDepthAndLogScales := paramsBootstrap.GetSlotsToCoeffsFactorizationDepthAndLogScales()

	// Slots To Coeffs params
	SlotsToCoeffsLevels := make([]int, len(SlotsToCoeffsFactorizationDepthAndLogScales))
	for i := range SlotsToCoeffsLevels {
		SlotsToCoeffsLevels[i] = len(SlotsToCoeffsFactorizationDepthAndLogScales[i])
	}

	var Iterations int
	if Iterations, err = paramsBootstrap.GetIterations(); err != nil {
		return ckks.ParametersLiteral{}, Parameters{}, err
	}

	S2CParams := advanced.HomomorphicDFTMatrixLiteral{
		Type:            advanced.Decode,
		RepackImag2Real: true,
		LevelStart:      len(paramsCKKS.LogQ) - 1 + len(SlotsToCoeffsFactorizationDepthAndLogScales) + Iterations - 1,
		LogBSGSRatio:    1,
		Levels:          SlotsToCoeffsLevels,
	}

	var EvalModLogScale int
	if EvalModLogScale, err = paramsBootstrap.GetEvalModLogScale(); err != nil {
		return ckks.ParametersLiteral{}, Parameters{}, err
	}

	SineType := paramsBootstrap.GetSineType()

	var ArcSineDegree int
	if ArcSineDegree, err = paramsBootstrap.GetArcSineDegree(); err != nil {
		return ckks.ParametersLiteral{}, Parameters{}, err
	}

	var LogMessageRatio int
	if LogMessageRatio, err = paramsBootstrap.GetLogMessageRatio(); err != nil {
		return ckks.ParametersLiteral{}, Parameters{}, err
	}

	var K int
	if K, err = paramsBootstrap.GetK(); err != nil {
		return ckks.ParametersLiteral{}, Parameters{}, err
	}

	var DoubleAngle int
	if DoubleAngle, err = paramsBootstrap.GetDoubleAngle(); err != nil {
		return ckks.ParametersLiteral{}, Parameters{}, err
	}

	var SineDegree int
	if SineDegree, err = paramsBootstrap.GetSineDegree(); err != nil {
		return ckks.ParametersLiteral{}, Parameters{}, err
	}

	EvalModParams := advanced.EvalModLiteral{
		LogScale:        EvalModLogScale,
		SineType:        SineType,
		SineDegree:      SineDegree,
		DoubleAngle:     DoubleAngle,
		K:               K,
		LogMessageRatio: LogMessageRatio,
		ArcSineDegree:   ArcSineDegree,
	}

	var EphemeralSecretWeight int
	if EphemeralSecretWeight, err = paramsBootstrap.GetEphemeralSecretWeight(); err != nil {
		return ckks.ParametersLiteral{}, Parameters{}, err
	}

	// Coeffs To Slots params
	EvalModParams.LevelStart = S2CParams.LevelStart + EvalModParams.Depth()

	CoeffsToSlotsLevels := make([]int, len(CoeffsToSlotsFactorizationDepthAndLogScales))
	for i := range CoeffsToSlotsLevels {
		CoeffsToSlotsLevels[i] = len(CoeffsToSlotsFactorizationDepthAndLogScales[i])
	}

	C2SParams := advanced.HomomorphicDFTMatrixLiteral{
		Type:            advanced.Encode,
		RepackImag2Real: true,
		LevelStart:      EvalModParams.LevelStart + len(CoeffsToSlotsFactorizationDepthAndLogScales),
		LogBSGSRatio:    1,
		Levels:          CoeffsToSlotsLevels,
	}

	LogQ := make([]int, len(paramsCKKS.LogQ))
	copy(LogQ, paramsCKKS.LogQ)

	for i := 0; i < Iterations-1; i++ {
		LogQ = append(LogQ, DefaultIterationsLogScale)
	}

	for i := range SlotsToCoeffsFactorizationDepthAndLogScales {
		var qi int
		for j := range SlotsToCoeffsFactorizationDepthAndLogScales[i] {
			qi += SlotsToCoeffsFactorizationDepthAndLogScales[i][j]
		}

		if qi+paramsCKKS.LogScale < 61 {
			qi += paramsCKKS.LogScale
		}

		LogQ = append(LogQ, qi)
	}

	for i := 0; i < EvalModParams.Depth(); i++ {
		LogQ = append(LogQ, EvalModLogScale)
	}

	for i := range CoeffsToSlotsFactorizationDepthAndLogScales {
		var qi int
		for j := range CoeffsToSlotsFactorizationDepthAndLogScales[i] {
			qi += CoeffsToSlotsFactorizationDepthAndLogScales[i][j]
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
			Xe:       paramsCKKS.Xe,
			Xs:       paramsCKKS.Xs,
		},
		Parameters{
			EphemeralSecretWeight:   EphemeralSecretWeight,
			SlotsToCoeffsParameters: S2CParams,
			EvalModParameters:       EvalModParams,
			CoeffsToSlotsParameters: C2SParams,
			Iterations:              Iterations,
		}, nil
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

// MarshalBinary returns a JSON representation of the the target Parameters struct on a slice of bytes.
// See `Marshal` from the `encoding/json` package.
func (p *Parameters) MarshalBinary() (data []byte, err error) {
	return json.Marshal(p)
}

// UnmarshalBinary reads a JSON representation on the target Parameters struct.
// See `Unmarshal` from the `encoding/json` package.
func (p *Parameters) UnmarshalBinary(data []byte) (err error) {
	return json.Unmarshal(data, p)
}

// GaloisElements returns the list of Galois elements required to evaluate the bootstrapping.
func (p *Parameters) GaloisElements(params ckks.Parameters) (galEls []uint64) {

	logN := params.LogN()
	logSlots := params.LogSlots()

	// List of the rotation key values to needed for the bootstrapp
	keys := make(map[uint64]bool)

	//SubSum rotation needed X -> Y^slots rotations
	for i := logSlots; i < logN-1; i++ {
		keys[params.GaloisElementForColumnRotationBy(1<<i)] = true
	}

	p.CoeffsToSlotsParameters.LogN = logN
	p.SlotsToCoeffsParameters.LogN = logN

	p.CoeffsToSlotsParameters.LogSlots = logSlots
	p.SlotsToCoeffsParameters.LogSlots = logSlots

	for _, galEl := range p.CoeffsToSlotsParameters.GaloisElements(params) {
		keys[galEl] = true
	}

	for _, galEl := range p.SlotsToCoeffsParameters.GaloisElements(params) {
		keys[galEl] = true
	}

	galEls = make([]uint64, len(keys))

	var i int
	for key := range keys {
		galEls[i] = key
		i++
	}

	return
}
