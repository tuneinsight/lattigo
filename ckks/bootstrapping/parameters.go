package bootstrapping

import (
	"encoding/json"
	"fmt"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// Parameters is a struct for the default bootstrapping parameters
type Parameters struct {
	SlotsToCoeffsParameters ckks.HomomorphicDFTMatrixLiteral
	EvalModParameters       ckks.EvalModLiteral
	CoeffsToSlotsParameters ckks.HomomorphicDFTMatrixLiteral
	Iterations              int
	EphemeralSecretWeight   int // Hamming weight of the ephemeral secret. If 0, no ephemeral secret is used during the bootstrapping.
}

// NewParametersFromLiteral takes as input a ckks.ParametersLiteral and a bootstrapping.ParametersLiteral structs and returns the
// appropriate ckks.ParametersLiteral for the bootstrapping circuit as well as the instantiated bootstrapping.Parameters.
// The returned ckks.ParametersLiteral contains allocated primes.
func NewParametersFromLiteral(ckksLit ckks.ParametersLiteral, btpLit ParametersLiteral) (ckks.ParametersLiteral, Parameters, error) {

	var err error

	if ckksLit.RingType != ring.Standard {
		return ckks.ParametersLiteral{}, Parameters{}, fmt.Errorf("NewParametersFromLiteral: invalid ring.RingType: must be ring.Standard")
	}

	var LogSlots int
	if LogSlots, err = btpLit.GetLogSlots(ckksLit.LogN); err != nil {
		return ckks.ParametersLiteral{}, Parameters{}, err
	}

	var CoeffsToSlotsFactorizationDepthAndLogDefaultScales [][]int
	if CoeffsToSlotsFactorizationDepthAndLogDefaultScales, err = btpLit.GetCoeffsToSlotsFactorizationDepthAndLogDefaultScales(LogSlots); err != nil {
		return ckks.ParametersLiteral{}, Parameters{}, err
	}

	var SlotsToCoeffsFactorizationDepthAndLogDefaultScales [][]int
	if SlotsToCoeffsFactorizationDepthAndLogDefaultScales, err = btpLit.GetSlotsToCoeffsFactorizationDepthAndLogDefaultScales(LogSlots); err != nil {
		return ckks.ParametersLiteral{}, Parameters{}, err
	}

	// Slots To Coeffs params
	SlotsToCoeffsLevels := make([]int, len(SlotsToCoeffsFactorizationDepthAndLogDefaultScales))
	for i := range SlotsToCoeffsLevels {
		SlotsToCoeffsLevels[i] = len(SlotsToCoeffsFactorizationDepthAndLogDefaultScales[i])
	}

	var Iterations int
	if Iterations, err = btpLit.GetIterations(); err != nil {
		return ckks.ParametersLiteral{}, Parameters{}, err
	}

	S2CParams := ckks.HomomorphicDFTMatrixLiteral{
		Type:            ckks.Decode,
		LogSlots:        LogSlots,
		RepackImag2Real: true,
		LevelStart:      len(ckksLit.LogQ) - 1 + len(SlotsToCoeffsFactorizationDepthAndLogDefaultScales) + Iterations - 1,
		LogBSGSRatio:    1,
		Levels:          SlotsToCoeffsLevels,
	}

	var EvalModLogDefaultScale int
	if EvalModLogDefaultScale, err = btpLit.GetEvalModLogDefaultScale(); err != nil {
		return ckks.ParametersLiteral{}, Parameters{}, err
	}

	SineType := btpLit.GetSineType()

	var ArcSineDegree int
	if ArcSineDegree, err = btpLit.GetArcSineDegree(); err != nil {
		return ckks.ParametersLiteral{}, Parameters{}, err
	}

	var LogMessageRatio int
	if LogMessageRatio, err = btpLit.GetLogMessageRatio(); err != nil {
		return ckks.ParametersLiteral{}, Parameters{}, err
	}

	var K int
	if K, err = btpLit.GetK(); err != nil {
		return ckks.ParametersLiteral{}, Parameters{}, err
	}

	var DoubleAngle int
	if DoubleAngle, err = btpLit.GetDoubleAngle(); err != nil {
		return ckks.ParametersLiteral{}, Parameters{}, err
	}

	var SineDegree int
	if SineDegree, err = btpLit.GetSineDegree(); err != nil {
		return ckks.ParametersLiteral{}, Parameters{}, err
	}

	EvalModParams := ckks.EvalModLiteral{
		LogDefaultScale: EvalModLogDefaultScale,
		SineType:        SineType,
		SineDegree:      SineDegree,
		DoubleAngle:     DoubleAngle,
		K:               K,
		LogMessageRatio: LogMessageRatio,
		ArcSineDegree:   ArcSineDegree,
	}

	var EphemeralSecretWeight int
	if EphemeralSecretWeight, err = btpLit.GetEphemeralSecretWeight(); err != nil {
		return ckks.ParametersLiteral{}, Parameters{}, err
	}

	// Coeffs To Slots params
	EvalModParams.LevelStart = S2CParams.LevelStart + EvalModParams.Depth()

	CoeffsToSlotsLevels := make([]int, len(CoeffsToSlotsFactorizationDepthAndLogDefaultScales))
	for i := range CoeffsToSlotsLevels {
		CoeffsToSlotsLevels[i] = len(CoeffsToSlotsFactorizationDepthAndLogDefaultScales[i])
	}

	C2SParams := ckks.HomomorphicDFTMatrixLiteral{
		Type:            ckks.Encode,
		LogSlots:        LogSlots,
		RepackImag2Real: true,
		LevelStart:      EvalModParams.LevelStart + len(CoeffsToSlotsFactorizationDepthAndLogDefaultScales),
		LogBSGSRatio:    1,
		Levels:          CoeffsToSlotsLevels,
	}

	LogQ := make([]int, len(ckksLit.LogQ))
	copy(LogQ, ckksLit.LogQ)

	for i := 0; i < Iterations-1; i++ {
		LogQ = append(LogQ, DefaultIterationsLogDefaultScale)
	}

	for i := range SlotsToCoeffsFactorizationDepthAndLogDefaultScales {
		var qi int
		for j := range SlotsToCoeffsFactorizationDepthAndLogDefaultScales[i] {
			qi += SlotsToCoeffsFactorizationDepthAndLogDefaultScales[i][j]
		}

		if qi+ckksLit.LogDefaultScale < 61 {
			qi += ckksLit.LogDefaultScale
		}

		LogQ = append(LogQ, qi)
	}

	for i := 0; i < EvalModParams.Depth(); i++ {
		LogQ = append(LogQ, EvalModLogDefaultScale)
	}

	for i := range CoeffsToSlotsFactorizationDepthAndLogDefaultScales {
		var qi int
		for j := range CoeffsToSlotsFactorizationDepthAndLogDefaultScales[i] {
			qi += CoeffsToSlotsFactorizationDepthAndLogDefaultScales[i][j]
		}
		LogQ = append(LogQ, qi)
	}

	LogP := make([]int, len(ckksLit.LogP))
	copy(LogP, ckksLit.LogP)

	Q, P, err := rlwe.GenModuli(ckksLit.LogN+1, LogQ, LogP)

	if err != nil {
		return ckks.ParametersLiteral{}, Parameters{}, err
	}

	return ckks.ParametersLiteral{
			LogN:            ckksLit.LogN,
			Q:               Q,
			P:               P,
			LogDefaultScale: ckksLit.LogDefaultScale,
			Xe:              ckksLit.Xe,
			Xs:              ckksLit.Xs,
		},
		Parameters{
			EphemeralSecretWeight:   EphemeralSecretWeight,
			SlotsToCoeffsParameters: S2CParams,
			EvalModParameters:       EvalModParams,
			CoeffsToSlotsParameters: C2SParams,
			Iterations:              Iterations,
		}, nil
}

// LogMaxDimensions returns the log plaintext dimensions of the target Parameters.
func (p *Parameters) LogMaxDimensions() ring.Dimensions {
	return ring.Dimensions{Rows: 0, Cols: p.SlotsToCoeffsParameters.LogSlots}
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

// MarshalBinary returns a JSON representation of the bootstrapping Parameters struct.
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

	// List of the rotation key values to needed for the bootstrapp
	keys := make(map[uint64]bool)

	//SubSum rotation needed X -> Y^slots rotations
	for i := p.LogMaxDimensions().Cols; i < logN-1; i++ {
		keys[params.GaloisElement(1<<i)] = true
	}

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
