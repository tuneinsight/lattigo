package ckks

import (
	"fmt"
	"math"
	"math/cmplx"

	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// Bootstrapper is a struct to stores a memory pool the plaintext matrices
// the polynomial approximation and the keys for the bootstrapping.
type Bootstrapper struct {
	*evaluator
	BootstrappingParameters
	*BootstrappingKey
	params Parameters

	dslots    int // Number of plaintext slots after the re-encoding
	logdslots int

	encoder Encoder // Encoder

	evalModPoly EvalModPoly
	stcMatrices EncodingMatrices
	ctsMatrices EncodingMatrices

	rotKeyIndex []int // a list of the required rotation keys
}

func sin2pi2pi(x complex128) complex128 {
	return cmplx.Sin(6.283185307179586*x) / 6.283185307179586
}

func cos2pi(x complex128) complex128 {
	return cmplx.Cos(6.283185307179586 * x)
}

// NewBootstrapper creates a new Bootstrapper.
func NewBootstrapper(params Parameters, btpParams BootstrappingParameters, btpKey BootstrappingKey) (btp *Bootstrapper, err error) {

	if btpParams.SineType == SineType(Sin) && btpParams.DoubleAngle != 0 {
		return nil, fmt.Errorf("cannot use double angle formul for SineType = Sin -> must use SineType = Cos")
	}

	if btpParams.CoeffsToSlotsParameters.LevelStart-btpParams.CoeffsToSlotsParameters.Depth(true) != btpParams.EvalModParameters.LevelStart {
		return nil, fmt.Errorf("starting level and depth of CoeffsToSlotsParameters inconsistent starting level of SineEvalParameters")
	}

	if btpParams.EvalModParameters.LevelStart-btpParams.EvalModParameters.Depth() != btpParams.SlotsToCoeffsParameters.LevelStart {
		return nil, fmt.Errorf("starting level and depth of SineEvalParameters inconsistent starting level of CoeffsToSlotsParameters")
	}

	btp = newBootstrapper(params, btpParams)

	btp.BootstrappingKey = &BootstrappingKey{btpKey.Rlk, btpKey.Rtks}
	if err = btp.CheckKeys(); err != nil {
		return nil, fmt.Errorf("invalid bootstrapping key: %w", err)
	}
	btp.evaluator = btp.evaluator.WithKey(rlwe.EvaluationKey{Rlk: btpKey.Rlk, Rtks: btpKey.Rtks}).(*evaluator)

	return btp, nil
}

// newBootstrapper is a constructor of "dummy" bootstrapper to enable the generation of bootstrapping-related constants
// without providing a bootstrapping key. To be replaced by a proper factorization of the bootstrapping pre-computations.
func newBootstrapper(params Parameters, btpParams BootstrappingParameters) (btp *Bootstrapper) {
	btp = new(Bootstrapper)

	btp.params = params
	btp.BootstrappingParameters = btpParams

	btp.dslots = params.Slots()
	btp.logdslots = params.LogSlots()
	if params.LogSlots() < params.MaxLogSlots() {
		btp.dslots <<= 1
		btp.logdslots++
	}

	btp.encoder = NewEncoder(params)
	btp.evaluator = NewEvaluator(params, rlwe.EvaluationKey{}).(*evaluator) // creates an evaluator without keys for genDFTMatrices

	btp.evalModPoly = btpParams.EvalModParameters.GenPoly()
	btp.genDFTMatrices()

	btp.ctxpool = NewCiphertext(params, 1, params.MaxLevel(), 0)

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

	a := real(btp.evalModPoly.SinePoly.a)
	b := real(btp.evalModPoly.SinePoly.b)
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
	slotsToCoeffsDiffScale := complex(math.Pow(btp.params.Scale()/(btp.evalModPoly.ScalingFactor/btp.MessageRatio), 1.0/float64(btp.SlotsToCoeffsParameters.Depth(false))), 0)
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
func addEncodingMatrixRotationsToList(pVec PtDiagMatrix, rotations []int, slots int, repack bool) []int {

	if pVec.naive {
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
