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

	var hasLogQ bool
	var residualLevel int

	if len(ckksLit.LogQ)+len(ckksLit.LogQ) == 0 && len(ckksLit.Q)+len(ckksLit.P) != 0 {
		residualLevel = len(ckksLit.Q) - 1

	} else if len(ckksLit.LogQ)+len(ckksLit.LogQ) != 0 && len(ckksLit.Q)+len(ckksLit.P) == 0 {
		hasLogQ = true
		residualLevel = len(ckksLit.LogQ) - 1

	} else {
		return ckks.ParametersLiteral{}, Parameters{}, fmt.Errorf("cannot NewParametersFromLiteral: must specify (LogQ, LogP) or (Q, P) but not a mix of both")
	}

	var LogSlots int
	if LogSlots, err = btpLit.GetLogSlots(ckksLit.LogN); err != nil {
		return ckks.ParametersLiteral{}, Parameters{}, err
	}

	var CoeffsToSlotsFactorizationDepthAndLogPlaintextScales [][]int
	if CoeffsToSlotsFactorizationDepthAndLogPlaintextScales, err = btpLit.GetCoeffsToSlotsFactorizationDepthAndLogPlaintextScales(LogSlots); err != nil {
		return ckks.ParametersLiteral{}, Parameters{}, err
	}

	var SlotsToCoeffsFactorizationDepthAndLogPlaintextScales [][]int
	if SlotsToCoeffsFactorizationDepthAndLogPlaintextScales, err = btpLit.GetSlotsToCoeffsFactorizationDepthAndLogPlaintextScales(LogSlots); err != nil {
		return ckks.ParametersLiteral{}, Parameters{}, err
	}

	// Slots To Coeffs params
	SlotsToCoeffsLevels := make([]int, len(SlotsToCoeffsFactorizationDepthAndLogPlaintextScales))
	for i := range SlotsToCoeffsLevels {
		SlotsToCoeffsLevels[i] = len(SlotsToCoeffsFactorizationDepthAndLogPlaintextScales[i])
	}

	var Iterations int
	if Iterations, err = btpLit.GetIterations(); err != nil {
		return ckks.ParametersLiteral{}, Parameters{}, err
	}

	S2CParams := ckks.HomomorphicDFTMatrixLiteral{
		Type:            ckks.Decode,
		LogSlots:        LogSlots,
		RepackImag2Real: true,
		LevelStart:      residualLevel + len(SlotsToCoeffsFactorizationDepthAndLogPlaintextScales) + Iterations - 1,
		LogBSGSRatio:    1,
		Levels:          SlotsToCoeffsLevels,
	}

	var EvalModLogPlaintextScale int
	if EvalModLogPlaintextScale, err = btpLit.GetEvalModLogPlaintextScale(); err != nil {
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
		LogPlaintextScale: EvalModLogPlaintextScale,
		SineType:          SineType,
		SineDegree:        SineDegree,
		DoubleAngle:       DoubleAngle,
		K:                 K,
		LogMessageRatio:   LogMessageRatio,
		ArcSineDegree:     ArcSineDegree,
	}

	var EphemeralSecretWeight int
	if EphemeralSecretWeight, err = btpLit.GetEphemeralSecretWeight(); err != nil {
		return ckks.ParametersLiteral{}, Parameters{}, err
	}

	// Coeffs To Slots params
	EvalModParams.LevelStart = S2CParams.LevelStart + EvalModParams.Depth()

	CoeffsToSlotsLevels := make([]int, len(CoeffsToSlotsFactorizationDepthAndLogPlaintextScales))
	for i := range CoeffsToSlotsLevels {
		CoeffsToSlotsLevels[i] = len(CoeffsToSlotsFactorizationDepthAndLogPlaintextScales[i])
	}

	C2SParams := ckks.HomomorphicDFTMatrixLiteral{
		Type:            ckks.Encode,
		LogSlots:        LogSlots,
		RepackImag2Real: true,
		LevelStart:      EvalModParams.LevelStart + len(CoeffsToSlotsFactorizationDepthAndLogPlaintextScales),
		LogBSGSRatio:    1,
		Levels:          CoeffsToSlotsLevels,
	}

	LogQBootstrappingCircuit := []int{}

	for i := 0; i < Iterations-1; i++ {
		LogQBootstrappingCircuit = append(LogQBootstrappingCircuit, DefaultIterationsLogPlaintextScale)
	}

	for i := range SlotsToCoeffsFactorizationDepthAndLogPlaintextScales {
		var qi int
		for j := range SlotsToCoeffsFactorizationDepthAndLogPlaintextScales[i] {
			qi += SlotsToCoeffsFactorizationDepthAndLogPlaintextScales[i][j]
		}

		if qi+ckksLit.LogPlaintextScale < 61 {
			qi += ckksLit.LogPlaintextScale
		}

		LogQBootstrappingCircuit = append(LogQBootstrappingCircuit, qi)
	}

	for i := 0; i < EvalModParams.Depth(); i++ {
		LogQBootstrappingCircuit = append(LogQBootstrappingCircuit, EvalModLogPlaintextScale)
	}

	for i := range CoeffsToSlotsFactorizationDepthAndLogPlaintextScales {
		var qi int
		for j := range CoeffsToSlotsFactorizationDepthAndLogPlaintextScales[i] {
			qi += CoeffsToSlotsFactorizationDepthAndLogPlaintextScales[i][j]
		}
		LogQBootstrappingCircuit = append(LogQBootstrappingCircuit, qi)
	}

	var Q, P []uint64
	// Specific moduli are given in the residual parameters
	if hasLogQ {

		if Q, P, err = rlwe.GenModuli(ckksLit.LogN+1, append(ckksLit.LogQ, LogQBootstrappingCircuit...), ckksLit.LogP); err != nil {
			return ckks.ParametersLiteral{}, Parameters{}, fmt.Errorf("cannot NewParametersFromLiteral: %w", err)
		}

		// Only the bit-size of the moduli are given in the residual parameters
	} else {

		// Extracts all the different primes
		primesHave := map[uint64]bool{}

		for _, qi := range ckksLit.Q {
			primesHave[qi] = true
		}

		for _, pj := range ckksLit.P {
			primesHave[pj] = true
		}

		// Maps the number of primes per bit size
		primesBitLenNew := map[int]int{}
		for _, logqi := range LogQBootstrappingCircuit {
			primesBitLenNew[logqi]++
		}

		// Map to store [bit-size][]primes
		primesNew := map[int][]uint64{}

		// For each bit-size
		for logqi, k := range primesBitLenNew {

			// Creates a new prime generator
			g := ring.NewNTTFriendlyPrimesGenerator(uint64(logqi), 1<<(ckksLit.LogN+1))

			// Populates the list with primes that aren't yet in primesHave
			primes := make([]uint64, k)
			var i int
			for i < k {

				for {
					qi, err := g.NextAlternatingPrime()

					if err != nil {
						return ckks.ParametersLiteral{}, Parameters{}, fmt.Errorf("cannot NewParametersFromLiteral: NextAlternatingPrime for 2^{%d} +/- k*2N + 1: %w", logqi, err)

					}

					if _, ok := primesHave[qi]; !ok {
						primes[i] = qi
						i++
						break
					}
				}
			}

			primesNew[logqi] = primes
		}

		Q = make([]uint64, len(ckksLit.Q))
		copy(Q, ckksLit.Q)

		// Appends to the residual modli
		for _, qi := range LogQBootstrappingCircuit {
			Q = append(Q, primesNew[qi][0])
			primesNew[qi] = primesNew[qi][1:]
		}

		P = make([]uint64, len(ckksLit.P))
		copy(P, ckksLit.P)
	}

	return ckks.ParametersLiteral{
			LogN:              ckksLit.LogN,
			Q:                 Q,
			P:                 P,
			LogPlaintextScale: ckksLit.LogPlaintextScale,
			Xe:                ckksLit.Xe,
			Xs:                ckksLit.Xs,
		},
		Parameters{
			EphemeralSecretWeight:   EphemeralSecretWeight,
			SlotsToCoeffsParameters: S2CParams,
			EvalModParameters:       EvalModParams,
			CoeffsToSlotsParameters: C2SParams,
			Iterations:              Iterations,
		}, nil
}

// PlaintextLogDimensions returns the log plaintext dimensions of the target Parameters.
func (p *Parameters) PlaintextLogDimensions() [2]int {
	return [2]int{0, p.SlotsToCoeffsParameters.LogSlots}
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
	for i := p.PlaintextLogDimensions()[1]; i < logN-1; i++ {
		keys[params.GaloisElement(1<<i)] = true
	}

	for _, galEl := range p.CoeffsToSlotsParameters.GaloisElements(params) {
		keys[galEl] = true
	}

	for _, galEl := range p.SlotsToCoeffsParameters.GaloisElements(params) {
		keys[galEl] = true
	}

	galEls = make([]uint64, len(keys)+1)

	var i int
	for key := range keys {
		galEls[i] = key
		i++
	}

	galEls[len(galEls)-1] = params.GaloisElementInverse()

	return
}
