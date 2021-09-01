package bootstrapping

import (
	"fmt"
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ckks/advanced"
	"github.com/ldsec/lattigo/v2/rlwe"
	"math"
)

// Bootstrapper is a struct to stores a memory pool the plaintext matrices
// the polynomial approximation and the keys for the bootstrapping.
type Bootstrapper struct {
	advanced.Evaluator
	Parameters
	*Key
	params ckks.Parameters

	dslots    int // Number of plaintext slots after the re-encoding
	logdslots int

	encoder ckks.Encoder // Encoder

	evalModPoly advanced.EvalModPoly
	stcMatrices advanced.EncodingMatrix
	ctsMatrices advanced.EncodingMatrix

	rotKeyIndex []int // a list of the required rotation keys
}

// Key is a type for a CKKS bootstrapping key, wich regroups the necessary public relinearization
// and rotation keys (i.e., an EvaluationKey).
type Key rlwe.EvaluationKey

// NewBootstrapper creates a new Bootstrapper.
func NewBootstrapper(params ckks.Parameters, btpParams Parameters, btpKey Key) (btp *Bootstrapper, err error) {

	if btpParams.EvalModParameters.SineType == advanced.Sin && btpParams.EvalModParameters.DoubleAngle != 0 {
		return nil, fmt.Errorf("cannot use double angle formul for SineType = Sin -> must use SineType = Cos")
	}

	if btpParams.CoeffsToSlotsParameters.LevelStart-btpParams.CoeffsToSlotsParameters.Depth(true) != btpParams.EvalModParameters.LevelStart {
		return nil, fmt.Errorf("starting level and depth of CoeffsToSlotsParameters inconsistent starting level of SineEvalParameters")
	}

	if btpParams.EvalModParameters.LevelStart-btpParams.EvalModParameters.Depth() != btpParams.SlotsToCoeffsParameters.LevelStart {
		return nil, fmt.Errorf("starting level and depth of SineEvalParameters inconsistent starting level of CoeffsToSlotsParameters")
	}

	btp = new(Bootstrapper)

	btp.params = params
	btp.Parameters = btpParams

	btp.dslots = params.Slots()
	btp.logdslots = params.LogSlots()
	if params.LogSlots() < params.MaxLogSlots() {
		btp.dslots <<= 1
		btp.logdslots++
	}

	btp.encoder = ckks.NewEncoder(params)

	btp.evalModPoly = advanced.NewEvalModPolyFromLiteral(btpParams.EvalModParameters)
	btp.genDFTMatrices()

	btp.Key = &Key{btpKey.Rlk, btpKey.Rtks}
	if err = btp.CheckKeys(); err != nil {
		return nil, fmt.Errorf("invalid bootstrapping key: %w", err)
	}

	btp.Evaluator = advanced.NewEvaluator(params, rlwe.EvaluationKey{Rlk: btpKey.Rlk, Rtks: btpKey.Rtks})

	return btp, nil
}

// CheckKeys checks if all the necessary keys are present
func (btp *Bootstrapper) CheckKeys() (err error) {

	if btp.Rlk == nil {
		return fmt.Errorf("relinearization key is nil")
	}

	if btp.Rtks == nil {
		return fmt.Errorf("rotation key is nil")
	}

	rotMissing := []int{}
	for _, i := range btp.rotKeyIndex {
		galEl := btp.params.GaloisElementForColumnRotationBy(int(i))
		if _, generated := btp.Rtks.Keys[galEl]; !generated {
			rotMissing = append(rotMissing, i)
		}
	}

	if len(rotMissing) != 0 {
		return fmt.Errorf("rotation key(s) missing: %d", rotMissing)
	}

	return nil
}

func (btp *Bootstrapper) genDFTMatrices() {

	K := btp.evalModPoly.K()
	n := float64(btp.params.N())
	scFac := btp.evalModPoly.ScFac()
	ctsDepth := float64(btp.CoeffsToSlotsParameters.Depth(false))
	stcDepth := float64(btp.SlotsToCoeffsParameters.Depth(false))

	// Correcting factor for approximate division by Q
	// The second correcting factor for approximate multiplication by Q is included in the coefficients of the EvalMod polynomials
	qDiff := btp.evalModPoly.QDiff()

	// CoeffsToSlots vectors
	// Change of variable for the evaluation of the Chebyshev polynomial + cancelling factor for the DFT and SubSum + evantual scaling factor for the double angle formula
	coeffsToSlotsDiffScale := complex(math.Pow(1.0/(K*n*scFac*qDiff), 1.0/ctsDepth), 0)
	btp.ctsMatrices = advanced.NewHomomorphicEncodingMatrixFromLiteral(btp.CoeffsToSlotsParameters, btp.encoder, btp.params.LogN(), btp.params.LogSlots(), coeffsToSlotsDiffScale)

	// SlotsToCoeffs vectors
	// Rescaling factor to set the final ciphertext to the desired scale
	slotsToCoeffsDiffScale := complex(math.Pow(btp.params.Scale()/(btp.evalModPoly.ScalingFactor()/btp.evalModPoly.MessageRatio()), 1.0/stcDepth), 0)
	btp.stcMatrices = advanced.NewHomomorphicEncodingMatrixFromLiteral(btp.SlotsToCoeffsParameters, btp.encoder, btp.params.LogN(), btp.params.LogSlots(), slotsToCoeffsDiffScale)

	// List of the rotation key values to needed for the bootstrapp
	btp.rotKeyIndex = []int{}
	btp.rotKeyIndex = append(btp.rotKeyIndex, btp.params.RotationsForTrace(btp.params.LogSlots(), btp.params.MaxLogSlots())...)
	btp.rotKeyIndex = append(btp.rotKeyIndex, btp.CoeffsToSlotsParameters.Rotations(btp.params.LogN(), btp.params.LogSlots())...)
	btp.rotKeyIndex = append(btp.rotKeyIndex, btp.SlotsToCoeffsParameters.Rotations(btp.params.LogN(), btp.params.LogSlots())...)
}
