package bootstrapping

import (
	"fmt"
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
	"math"
)

// Bootstrapper is a struct to stores a memory pool the plaintext matrices
// the polynomial approximation and the keys for the bootstrapping.
type Bootstrapper struct {
	ckks.Evaluator
	Parameters
	*Key
	params ckks.Parameters

	dslots    int // Number of plaintext slots after the re-encoding
	logdslots int

	encoder ckks.Encoder // Encoder

	evalModPoly ckks.EvalModPoly
	stcMatrices ckks.EncodingMatrices
	ctsMatrices ckks.EncodingMatrices

	rotKeyIndex []int // a list of the required rotation keys
}

// Key is a type for a CKKS bootstrapping key, wich regroups the necessary public relinearization
// and rotation keys (i.e., an EvaluationKey).
type Key rlwe.EvaluationKey

// NewBootstrapper creates a new Bootstrapper.
func NewBootstrapper(params ckks.Parameters, btpParams Parameters, btpKey Key) (btp *Bootstrapper, err error) {

	if btpParams.EvalModParameters.SineType == ckks.Sin && btpParams.EvalModParameters.DoubleAngle != 0 {
		return nil, fmt.Errorf("cannot use double angle formul for SineType = Sin -> must use SineType = Cos")
	}

	if btpParams.CoeffsToSlotsParameters.LevelStart-btpParams.CoeffsToSlotsParameters.Depth(true) != btpParams.EvalModParameters.LevelStart {
		return nil, fmt.Errorf("starting level and depth of CoeffsToSlotsParameters inconsistent starting level of SineEvalParameters")
	}

	if btpParams.EvalModParameters.LevelStart-btpParams.EvalModParameters.Depth() != btpParams.SlotsToCoeffsParameters.LevelStart {
		return nil, fmt.Errorf("starting level and depth of SineEvalParameters inconsistent starting level of CoeffsToSlotsParameters")
	}

	btp = newBootstrapper(params, btpParams)

	btp.Key = &Key{btpKey.Rlk, btpKey.Rtks}
	if err = btp.CheckKeys(); err != nil {
		return nil, fmt.Errorf("invalid bootstrapping key: %w", err)
	}
	btp.Evaluator = btp.Evaluator.WithKey(rlwe.EvaluationKey{Rlk: btpKey.Rlk, Rtks: btpKey.Rtks})

	return btp, nil
}

// newBootstrapper is a constructor of "dummy" bootstrapper to enable the generation of bootstrapping-related constants
// without providing a bootstrapping key. To be replaced by a proper factorization of the bootstrapping pre-computations.
func newBootstrapper(params ckks.Parameters, btpParams Parameters) (btp *Bootstrapper) {
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
	btp.Evaluator = ckks.NewEvaluator(params, rlwe.EvaluationKey{})

	btp.evalModPoly = btpParams.EvalModParameters.GenPoly()
	btp.genDFTMatrices()

	return btp
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

	a := real(btp.evalModPoly.SinePoly.A)
	b := real(btp.evalModPoly.SinePoly.B)
	n := float64(btp.params.N())

	// Correcting factor for approximate division by Q
	// The second correcting factor for approximate multiplication by Q is included in the coefficients of the EvalMod polynomials
	qDiff := btp.params.QiFloat64(0) / math.Exp2(math.Round(math.Log2(btp.params.QiFloat64(0))))

	// CoeffsToSlots vectors
	// Change of variable for the evaluation of the Chebyshev polynomial + cancelling factor for the DFT and SubSum + evantual scaling factor for the double angle formula
	coeffsToSlotsDiffScale := complex(math.Pow(2.0/((b-a)*n*btp.evalModPoly.ScFac*qDiff), 1.0/float64(btp.CoeffsToSlotsParameters.Depth(false))), 0)
	btp.ctsMatrices = btp.encoder.GenHomomorphicEncodingMatrices(btp.CoeffsToSlotsParameters, coeffsToSlotsDiffScale)

	// SlotsToCoeffs vectors
	// Rescaling factor to set the final ciphertext to the desired scale
	slotsToCoeffsDiffScale := complex(math.Pow(btp.params.Scale()/(btp.evalModPoly.ScalingFactor/btp.evalModPoly.MessageRatio), 1.0/float64(btp.SlotsToCoeffsParameters.Depth(false))), 0)
	btp.stcMatrices = btp.encoder.GenHomomorphicEncodingMatrices(btp.SlotsToCoeffsParameters, slotsToCoeffsDiffScale)

	// List of the rotation key values to needed for the bootstrapp
	btp.rotKeyIndex = []int{}

	//SubSum rotation needed X -> Y^slots rotations
	for i := btp.params.LogSlots(); i < btp.params.MaxLogSlots(); i++ {
		if !utils.IsInSliceInt(1<<i, btp.rotKeyIndex) {
			btp.rotKeyIndex = append(btp.rotKeyIndex, 1<<i)
		}
	}

	// Coeffs to Slots rotations
	for _, pVec := range btp.ctsMatrices.Matrices {
		btp.rotKeyIndex = addEncodingMatrixRotationsToList(pVec, btp.rotKeyIndex, btp.params.Slots(), false)
	}

	// Slots to Coeffs rotations
	for i, pVec := range btp.stcMatrices.Matrices {
		btp.rotKeyIndex = addEncodingMatrixRotationsToList(pVec, btp.rotKeyIndex, btp.params.Slots(), (i == 0) && (btp.params.LogSlots() < btp.params.MaxLogSlots()))
	}
}

// AddMatrixRotToList adds the rotations neede to evaluate pVec to the list rotations
func addEncodingMatrixRotationsToList(pVec ckks.PtDiagMatrix, rotations []int, slots int, repack bool) []int {

	if pVec.Naive {
		for j := range pVec.Vec {
			if !utils.IsInSliceInt(j, rotations) {
				rotations = append(rotations, j)
			}
		}
	} else {
		var index int
		for j := range pVec.Vec {

			N1 := pVec.N1

			index = ((j / N1) * N1)

			if repack {
				// Sparse repacking, occurring during the first DFT matrix of the CoeffsToSlots.
				index &= 2*slots - 1
			} else {
				// Other cases
				index &= slots - 1
			}

			if index != 0 && !utils.IsInSliceInt(index, rotations) {
				rotations = append(rotations, index)
			}

			index = j & (N1 - 1)

			if index != 0 && !utils.IsInSliceInt(index, rotations) {
				rotations = append(rotations, index)
			}
		}
	}

	return rotations
}
