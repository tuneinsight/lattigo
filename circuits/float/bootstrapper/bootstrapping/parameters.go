package bootstrapping

import (
	"encoding/json"
	"fmt"

	"github.com/google/go-cmp/cmp"

	"github.com/tuneinsight/lattigo/v4/circuits/float"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// Parameters is a struct storing the parameters
// of the bootstrapping circuit.
type Parameters struct {
	ckks.Parameters
	SlotsToCoeffsParameters float.DFTMatrixLiteral
	Mod1ParametersLiteral   float.Mod1ParametersLiteral
	CoeffsToSlotsParameters float.DFTMatrixLiteral
	IterationsParameters    *IterationsParameters
	EphemeralSecretWeight   int // Hamming weight of the ephemeral secret. If 0, no ephemeral secret is used during the bootstrapping.
}

// NewParametersFromLiteral instantiates a bootstrapping.Parameters from the residual ckks.Parameters and
// a bootstrapping.ParametersLiteral struct.
//
// The residualParameters corresponds to the ckks.Parameters that are left after the bootstrapping circuit is evaluated.
// These are entirely independant of the bootstrapping parameters with one exception: the ciphertext primes Qi must be
// congruent to 1 mod 2N of the bootstrapping parameters (note that the auxiliary primes Pi do not need to be).
// This is required because the primes Qi of the residual parameters and the bootstrapping parameters are the same between
// the two sets of parameters.
//
// The user can ensure that this condition is met by setting the appropriate LogNThRoot in the ckks.ParametersLiteral before
// instantiating them.
//
// The method NewParametersFromLiteral will automatically allocate the ckks.Parameters of the bootstrapping circuit based on
// the provided residualParameters and the information given in the bootstrapping.ParametersLiteral.
func NewParametersFromLiteral(residualParameters ckks.Parameters, btpLit ParametersLiteral) (Parameters, error) {

	var err error

	// Retrieve the LogN of the bootstrapping circuit
	LogN := btpLit.GetLogN()

	// Retrive the NthRoot
	var NthRoot uint64
	switch residualParameters.RingType() {
	case ring.ConjugateInvariant:

		// If ConjugateInvariant, then the bootstrapping LogN must be at least 1 greater
		// than the residualParameters LogN
		if LogN <= residualParameters.LogN() {
			return Parameters{}, fmt.Errorf("cannot NewParametersFromLiteral: LogN of bootstrapping parameters must be greater than LogN of residual parameters if ringtype is ConjugateInvariant")
		}

		// Takes the greatest NthRoot between the residualParameters NthRoot and the bootstrapphg NthRoot
		NthRoot = utils.Max(uint64(residualParameters.N()<<2), uint64(2<<LogN))

	default:

		// The LogN of the bootstrapping parameters cannot be smaller than the LogN of the residualParameters.
		if LogN < residualParameters.LogN() {
			return Parameters{}, fmt.Errorf("cannot NewParametersFromLiteral: LogN of bootstrapping parameters must be greater or equal to LogN of residual parameters")
		}

		// Takes the greatest NthRoot between the residualParameters NthRoot and the bootstrapphg NthRoot
		NthRoot = utils.Max(uint64(residualParameters.N()<<1), uint64(2<<LogN))
	}

	// Checks that all primes Qi of the residualParameters are congruent to 1 mod NthRoot of the bootstrapping parameters.
	for i, q := range residualParameters.Q() {
		if q&(NthRoot-1) != 1 {
			return Parameters{}, fmt.Errorf("cannot NewParametersFromLiteral: Q[%d]=%d != 1 mod NthRoot=%d", i, q, NthRoot)
		}
	}

	// Retrieves the variable LogSlots, which is used to instantiates the encoding/decoding matrices.
	var LogSlots int
	if LogSlots, err = btpLit.GetLogSlots(); err != nil {
		return Parameters{}, err
	}

	// Retrieves the factorization depth and scaling factor of the encoding matrix
	var CoeffsToSlotsFactorizationDepthAndLogScales [][]int
	if CoeffsToSlotsFactorizationDepthAndLogScales, err = btpLit.GetCoeffsToSlotsFactorizationDepthAndLogScales(LogSlots); err != nil {
		return Parameters{}, err
	}

	// Retrieves the factorization depth and scaling factor of the decoding matrix
	var SlotsToCoeffsFactorizationDepthAndLogScales [][]int
	if SlotsToCoeffsFactorizationDepthAndLogScales, err = btpLit.GetSlotsToCoeffsFactorizationDepthAndLogScales(LogSlots); err != nil {
		return Parameters{}, err
	}

	// Slots To Coeffs params
	SlotsToCoeffsLevels := make([]int, len(SlotsToCoeffsFactorizationDepthAndLogScales))
	for i := range SlotsToCoeffsLevels {
		SlotsToCoeffsLevels[i] = len(SlotsToCoeffsFactorizationDepthAndLogScales[i])
	}

	// Number of bootstrapping iterations
	var iterParams *IterationsParameters
	if iterParams, err = btpLit.GetIterationsParameters(); err != nil {
		return Parameters{}, err
	}

	// Boolean if there is a reserved prime for the bootstrapping iterations
	var hasReservedIterationPrime int
	if iterParams != nil && iterParams.ReservedPrimeBitSize > 0 {
		hasReservedIterationPrime = 1
	}

	// SlotsToCoeffs parameters (homomorphic decoding)
	S2CParams := float.DFTMatrixLiteral{
		Type:            float.HomomorphicDecode,
		LogSlots:        LogSlots,
		RepackImag2Real: true,
		LevelStart:      residualParameters.MaxLevel() + len(SlotsToCoeffsFactorizationDepthAndLogScales) + hasReservedIterationPrime,
		LogBSGSRatio:    1,
		Levels:          SlotsToCoeffsLevels,
	}

	// Scaling factor of the homomorphic modular reduction x mod 1
	var EvalMod1LogScale int
	if EvalMod1LogScale, err = btpLit.GetEvalMod1LogScale(); err != nil {
		return Parameters{}, err
	}

	// Type of polynomial approximation of x mod 1
	SineType := btpLit.GetSineType()

	// Degree of the taylor series of arc sine
	var ArcSineDegree int
	if ArcSineDegree, err = btpLit.GetArcSineDegree(); err != nil {
		return Parameters{}, err
	}

	// Log2 ratio between Q[0] and |m| (i.e. gap between the message and Q[0])
	var LogMessageRatio int
	if LogMessageRatio, err = btpLit.GetLogMessageRatio(); err != nil {
		return Parameters{}, err
	}

	// Interval [-K+1, K-1] of the polynomial approximation of x mod 1
	var K int
	if K, err = btpLit.GetK(); err != nil {
		return Parameters{}, err
	}

	// Number of double angle evaluation if x mod 1 is approximated with cos
	var DoubleAngle int
	if DoubleAngle, err = btpLit.GetDoubleAngle(); err != nil {
		return Parameters{}, err
	}

	// Degree of the polynomial approximation of x mod 1
	var SineDegree int
	if SineDegree, err = btpLit.GetSineDegree(); err != nil {
		return Parameters{}, err
	}

	// Parameters of the homomorphic modular reduction x mod 1
	Mod1ParametersLiteral := float.Mod1ParametersLiteral{
		LogScale:        EvalMod1LogScale,
		SineType:        SineType,
		SineDegree:      SineDegree,
		DoubleAngle:     DoubleAngle,
		K:               K,
		LogMessageRatio: LogMessageRatio,
		ArcSineDegree:   ArcSineDegree,
	}

	// Hamming weight of the ephemeral secret key to which the ciphertext is
	// switched to during the ModUp step.
	var EphemeralSecretWeight int
	if EphemeralSecretWeight, err = btpLit.GetEphemeralSecretWeight(); err != nil {
		return Parameters{}, err
	}

	// Coeffs To Slots params
	Mod1ParametersLiteral.LevelStart = S2CParams.LevelStart + Mod1ParametersLiteral.Depth()

	CoeffsToSlotsLevels := make([]int, len(CoeffsToSlotsFactorizationDepthAndLogScales))
	for i := range CoeffsToSlotsLevels {
		CoeffsToSlotsLevels[i] = len(CoeffsToSlotsFactorizationDepthAndLogScales[i])
	}

	// Parameters of the CoeffsToSlots (homomorphic encoding)
	C2SParams := float.DFTMatrixLiteral{
		Type:            float.HomomorphicEncode,
		LogSlots:        LogSlots,
		RepackImag2Real: true,
		LevelStart:      Mod1ParametersLiteral.LevelStart + len(CoeffsToSlotsFactorizationDepthAndLogScales),
		LogBSGSRatio:    1,
		Levels:          CoeffsToSlotsLevels,
	}

	// List of the prime-size of all primes required by the bootstrapping circuit.
	LogQBootstrappingCircuit := []int{}

	// appends the reserved prime first for multiple iteraiton, if any
	if hasReservedIterationPrime == 1 {
		LogQBootstrappingCircuit = append(LogQBootstrappingCircuit, iterParams.ReservedPrimeBitSize)
	}

	// Appends all other primes in reverse order of the circuit
	for i := range SlotsToCoeffsFactorizationDepthAndLogScales {
		var qi int
		for j := range SlotsToCoeffsFactorizationDepthAndLogScales[i] {
			qi += SlotsToCoeffsFactorizationDepthAndLogScales[i][j]
		}

		if qi+residualParameters.LogDefaultScale() < 61 {
			qi += residualParameters.LogDefaultScale()
		}

		LogQBootstrappingCircuit = append(LogQBootstrappingCircuit, qi)
	}

	for i := 0; i < Mod1ParametersLiteral.Depth(); i++ {
		LogQBootstrappingCircuit = append(LogQBootstrappingCircuit, EvalMod1LogScale)
	}

	for i := range CoeffsToSlotsFactorizationDepthAndLogScales {
		var qi int
		for j := range CoeffsToSlotsFactorizationDepthAndLogScales[i] {
			qi += CoeffsToSlotsFactorizationDepthAndLogScales[i][j]
		}
		LogQBootstrappingCircuit = append(LogQBootstrappingCircuit, qi)
	}

	var Q, P []uint64

	// Extracts all the different primes Qi that are
	// in the residualParameters
	primesHave := map[uint64]bool{}
	for _, qi := range residualParameters.Q() {
		primesHave[qi] = true
	}

	// Maps the number of primes per bit size
	primesBitLenNew := map[int]int{}
	for _, logqi := range LogQBootstrappingCircuit {
		primesBitLenNew[logqi]++
	}

	// Retrieve the number of primes #Pi of the bootstrapping circuit
	// and adds them to the list of bit-size
	LogP := btpLit.GetLogP(C2SParams.LevelStart + 1)
	for _, logpi := range LogP {
		primesBitLenNew[logpi]++
	}

	// Map to store [bit-size][]primes
	primesNew := map[int][]uint64{}

	// For each bit-size sample a pair-wise coprime prime
	for logqi, k := range primesBitLenNew {

		// Creates a new prime generator
		g := ring.NewNTTFriendlyPrimesGenerator(uint64(logqi), NthRoot)

		// Populates the list with primes that aren't yet in primesHave
		primes := make([]uint64, k)
		var i int
		for i < k {

			for {
				qi, err := g.NextAlternatingPrime()

				if err != nil {
					return Parameters{}, fmt.Errorf("cannot NewParametersFromLiteral: NextAlternatingPrime for 2^{%d} +/- k*2N + 1: %w", logqi, err)

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

	// Constructs the set of primes Qi
	Q = make([]uint64, len(residualParameters.Q()))
	copy(Q, residualParameters.Q())

	// Appends to the residual modli
	for _, qi := range LogQBootstrappingCircuit {
		Q = append(Q, primesNew[qi][0])
		primesNew[qi] = primesNew[qi][1:]
	}

	// Constructs the set of primes Pi
	P = make([]uint64, len(LogP))
	for i, logpi := range LogP {
		P[i] = primesNew[logpi][0]
		primesNew[logpi] = primesNew[logpi][1:]
	}

	// Instantiates the ckks.Parameters of the bootstrapping circuit.
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            LogN,
		Q:               Q,
		P:               P,
		LogDefaultScale: residualParameters.LogDefaultScale(),
		Xs:              btpLit.GetDefaultXs(),
		Xe:              btpLit.GetDefaultXe(),
	})

	if err != nil {
		return Parameters{}, err
	}

	return Parameters{
		Parameters:              params,
		EphemeralSecretWeight:   EphemeralSecretWeight,
		SlotsToCoeffsParameters: S2CParams,
		Mod1ParametersLiteral:   Mod1ParametersLiteral,
		CoeffsToSlotsParameters: C2SParams,
		IterationsParameters:    iterParams,
	}, nil
}

func (p Parameters) Equal(other *Parameters) (res bool) {
	res = p.Parameters.Equal(&other.Parameters)

	res = res && p.EphemeralSecretWeight == other.EphemeralSecretWeight
	res = res && cmp.Equal(p.SlotsToCoeffsParameters, other.SlotsToCoeffsParameters)
	res = res && cmp.Equal(p.Mod1ParametersLiteral, other.Mod1ParametersLiteral)
	res = res && cmp.Equal(p.CoeffsToSlotsParameters, other.CoeffsToSlotsParameters)
	res = res && cmp.Equal(p.IterationsParameters, other.IterationsParameters)
	return
}

// LogMaxDimensions returns the log plaintext dimensions of the target Parameters.
func (p Parameters) LogMaxDimensions() ring.Dimensions {
	return ring.Dimensions{Rows: 0, Cols: p.SlotsToCoeffsParameters.LogSlots}
}

// LogMaxSlots returns the log of the maximum number of slots.
func (p Parameters) LogMaxSlots() int {
	return p.SlotsToCoeffsParameters.LogSlots
}

// DepthCoeffsToSlots returns the depth of the Coeffs to Slots of the CKKS bootstrapping.
func (p Parameters) DepthCoeffsToSlots() (depth int) {
	return p.SlotsToCoeffsParameters.Depth(true)
}

// DepthEvalMod returns the depth of the EvalMod step of the CKKS bootstrapping.
func (p Parameters) DepthEvalMod() (depth int) {
	return p.Mod1ParametersLiteral.Depth()
}

// DepthSlotsToCoeffs returns the depth of the Slots to Coeffs step of the CKKS bootstrapping.
func (p Parameters) DepthSlotsToCoeffs() (depth int) {
	return p.CoeffsToSlotsParameters.Depth(true)
}

// Depth returns the depth of the full bootstrapping circuit.
func (p Parameters) Depth() (depth int) {
	return p.DepthCoeffsToSlots() + p.DepthEvalMod() + p.DepthSlotsToCoeffs()
}

// MarshalBinary returns a JSON representation of the bootstrapping Parameters struct.
// See `Marshal` from the `encoding/json` package.
func (p Parameters) MarshalBinary() (data []byte, err error) {
	return json.Marshal(p)
}

// UnmarshalBinary reads a JSON representation on the target Parameters struct.
// See `Unmarshal` from the `encoding/json` package.
func (p *Parameters) UnmarshalBinary(data []byte) (err error) {
	return json.Unmarshal(data, p)
}

func (p Parameters) MarshalJSON() (data []byte, err error) {
	return json.Marshal(struct {
		Parameters              ckks.Parameters
		SlotsToCoeffsParameters float.DFTMatrixLiteral
		Mod1ParametersLiteral   float.Mod1ParametersLiteral
		CoeffsToSlotsParameters float.DFTMatrixLiteral
		IterationsParameters    *IterationsParameters
		EphemeralSecretWeight   int
	}{
		Parameters:              p.Parameters,
		SlotsToCoeffsParameters: p.SlotsToCoeffsParameters,
		Mod1ParametersLiteral:   p.Mod1ParametersLiteral,
		CoeffsToSlotsParameters: p.CoeffsToSlotsParameters,
		IterationsParameters:    p.IterationsParameters,
		EphemeralSecretWeight:   p.EphemeralSecretWeight,
	})
}

func (p *Parameters) UnmarshalJSON(data []byte) (err error) {
	var params struct {
		Parameters              ckks.Parameters
		SlotsToCoeffsParameters float.DFTMatrixLiteral
		Mod1ParametersLiteral   float.Mod1ParametersLiteral
		CoeffsToSlotsParameters float.DFTMatrixLiteral
		IterationsParameters    *IterationsParameters
		EphemeralSecretWeight   int
	}

	if err = json.Unmarshal(data, &params); err != nil {
		return
	}

	p.Parameters = params.Parameters
	p.SlotsToCoeffsParameters = params.SlotsToCoeffsParameters
	p.Mod1ParametersLiteral = params.Mod1ParametersLiteral
	p.CoeffsToSlotsParameters = params.CoeffsToSlotsParameters
	p.IterationsParameters = params.IterationsParameters
	p.EphemeralSecretWeight = params.EphemeralSecretWeight

	return
}

// GaloisElements returns the list of Galois elements required to evaluate the bootstrapping.
func (p Parameters) GaloisElements(params ckks.Parameters) (galEls []uint64) {

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

	keys[params.GaloisElementForComplexConjugation()] = true

	return utils.GetSortedKeys(keys)
}
