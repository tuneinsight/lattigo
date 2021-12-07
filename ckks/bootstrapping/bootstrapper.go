package bootstrapping

import (
	"fmt"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ckks/advanced"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// Bootstrapper is a struct to stores a memory pool the plaintext matrices
// the polynomial approximation and the keys for the bootstrapping.
type Bootstrapper struct {
	advanced.Evaluator
	*bootstrapperBase
}

type bootstrapperBase struct {
	Parameters
	params ckks.Parameters

	dslots    int // Number of plaintext slots after the re-encoding
	logdslots int

	evalModPoly advanced.EvalModPoly
	stcMatrices advanced.EncodingMatrix
	ctsMatrices advanced.EncodingMatrix
}

// NewBootstrapper creates a new Bootstrapper.
func NewBootstrapper(params ckks.Parameters, btpParams Parameters, btpKey rlwe.EvaluationKey) (btp *Bootstrapper, err error) {

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
	btp.bootstrapperBase = newBootstrapperBase(params, btpParams, btpKey)

	if err = btp.bootstrapperBase.CheckKeys(btpKey); err != nil {
		return nil, fmt.Errorf("invalid bootstrapping key: %w", err)
	}

	btp.Evaluator = advanced.NewEvaluator(params, btpKey)

	return
}

// ShallowCopy creates a shallow copy of this Bootstrapper in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Bootstrapper can be used concurrently.
func (btp *Bootstrapper) ShallowCopy() *Bootstrapper {
	return &Bootstrapper{
		Evaluator:        btp.Evaluator.ShallowCopy(),
		bootstrapperBase: btp.bootstrapperBase,
	}
}

// CheckKeys checks if all the necessary keys are present in the instantiated Bootstrapper
func (bb *bootstrapperBase) CheckKeys(btpKey rlwe.EvaluationKey) (err error) {

	if btpKey.Rlk == nil {
		return fmt.Errorf("relinearization key is nil")
	}

	if btpKey.Rtks == nil {
		return fmt.Errorf("rotation key is nil")
	}

	rotKeyIndex := []int{}
	rotKeyIndex = append(rotKeyIndex, bb.params.RotationsForTrace(bb.params.LogSlots(), bb.params.MaxLogSlots())...)
	rotKeyIndex = append(rotKeyIndex, bb.CoeffsToSlotsParameters.Rotations(bb.params.LogN(), bb.params.LogSlots())...)
	rotKeyIndex = append(rotKeyIndex, bb.SlotsToCoeffsParameters.Rotations(bb.params.LogN(), bb.params.LogSlots())...)

	rotMissing := []int{}
	for _, i := range rotKeyIndex {
		galEl := bb.params.GaloisElementForColumnRotationBy(int(i))
		if _, generated := btpKey.Rtks.Keys[galEl]; !generated {
			rotMissing = append(rotMissing, i)
		}
	}

	if len(rotMissing) != 0 {
		return fmt.Errorf("rotation key(s) missing: %d", rotMissing)
	}

	return nil
}

func newBootstrapperBase(params ckks.Parameters, btpParams Parameters, btpKey rlwe.EvaluationKey) (bb *bootstrapperBase) {
	bb = new(bootstrapperBase)
	bb.params = params
	bb.Parameters = btpParams

	bb.dslots = params.Slots()
	bb.logdslots = params.LogSlots()
	if params.LogSlots() < params.MaxLogSlots() {
		bb.dslots <<= 1
		bb.logdslots++
	}

	bb.evalModPoly = advanced.NewEvalModPolyFromLiteral(btpParams.EvalModParameters)

	scFac := bb.evalModPoly.ScFac()
	K := bb.evalModPoly.K() / scFac
	n := float64(2 * params.Slots())

	// Correcting factor for approximate division by Q
	// The second correcting factor for approximate multiplication by Q is included in the coefficients of the EvalMod polynomials
	qDiff := bb.evalModPoly.QDiff()

	encoder := ckks.NewEncoder(bb.params)

	// CoeffsToSlots vectors
	// Change of variable for the evaluation of the Chebyshev polynomial + cancelling factor for the DFT and SubSum + evantual scaling factor for the double angle formula
	bb.CoeffsToSlotsParameters.Scaling = 1.0 / (K * n * scFac * qDiff)
	bb.ctsMatrices = advanced.NewHomomorphicEncodingMatrixFromLiteral(bb.CoeffsToSlotsParameters, encoder)

	// SlotsToCoeffs vectors
	// Rescaling factor to set the final ciphertext to the desired scale
	bb.SlotsToCoeffsParameters.Scaling = bb.params.DefaultScale() / (bb.evalModPoly.ScalingFactor() / bb.evalModPoly.MessageRatio())
	bb.stcMatrices = advanced.NewHomomorphicEncodingMatrixFromLiteral(bb.SlotsToCoeffsParameters, encoder)

	encoder = nil

	return
}
