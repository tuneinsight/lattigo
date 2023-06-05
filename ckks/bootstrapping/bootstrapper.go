package bootstrapping

import (
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// Bootstrapper is a struct to store a memory buffer with the plaintext matrices,
// the polynomial approximation, and the keys for the bootstrapping.
type Bootstrapper struct {
	*ckks.Evaluator
	*bootstrapperBase
}

type bootstrapperBase struct {
	Parameters
	*EvaluationKeySet
	params ckks.Parameters

	dslots    int // Number of plaintext slots after the re-encoding
	logdslots int

	evalModPoly ckks.EvalModPoly
	stcMatrices ckks.HomomorphicDFTMatrix
	ctsMatrices ckks.HomomorphicDFTMatrix

	q0OverMessageRatio float64
}

// EvaluationKeySet is a type for a CKKS bootstrapping key, which
// regroups the necessary public relinearization and rotation keys.
type EvaluationKeySet struct {
	*rlwe.EvaluationKeySet
	EvkDtS *rlwe.EvaluationKey
	EvkStD *rlwe.EvaluationKey
}

// NewBootstrapper creates a new Bootstrapper.
func NewBootstrapper(params ckks.Parameters, btpParams Parameters, btpKeys *EvaluationKeySet) (btp *Bootstrapper, err error) {

	if btpParams.EvalModParameters.SineType == ckks.SinContinuous && btpParams.EvalModParameters.DoubleAngle != 0 {
		return nil, fmt.Errorf("cannot use double angle formul for SineType = Sin -> must use SineType = Cos")
	}

	if btpParams.EvalModParameters.SineType == ckks.CosDiscrete && btpParams.EvalModParameters.SineDegree < 2*(btpParams.EvalModParameters.K-1) {
		return nil, fmt.Errorf("SineType 'ckks.CosDiscrete' uses a minimum degree of 2*(K-1) but EvalMod degree is smaller")
	}

	if btpParams.CoeffsToSlotsParameters.LevelStart-btpParams.CoeffsToSlotsParameters.Depth(true) != btpParams.EvalModParameters.LevelStart {
		return nil, fmt.Errorf("starting level and depth of CoeffsToSlotsParameters inconsistent starting level of SineEvalParameters")
	}

	if btpParams.EvalModParameters.LevelStart-btpParams.EvalModParameters.Depth() != btpParams.SlotsToCoeffsParameters.LevelStart {
		return nil, fmt.Errorf("starting level and depth of SineEvalParameters inconsistent starting level of CoeffsToSlotsParameters")
	}

	btp = new(Bootstrapper)
	btp.bootstrapperBase = newBootstrapperBase(params, btpParams, btpKeys)

	if err = btp.bootstrapperBase.CheckKeys(btpKeys); err != nil {
		return nil, fmt.Errorf("invalid bootstrapping key: %w", err)
	}

	btp.EvaluationKeySet = btpKeys

	btp.Evaluator = ckks.NewEvaluator(params, btpKeys)

	return
}

// GenEvaluationKeySetNew generates a new bootstrapping EvaluationKeySet, which contain:
//
//	EvaluationKeySet: struct compliant to the interface rlwe.EvaluationKeySetInterface.
//	EvkDtS: *rlwe.EvaluationKey
//	EvkStD: *rlwe.EvaluationKey
func GenEvaluationKeySetNew(btpParams Parameters, ckksParams ckks.Parameters, sk *rlwe.SecretKey) *EvaluationKeySet {

	kgen := ckks.NewKeyGenerator(ckksParams)

	evk := rlwe.NewEvaluationKeySet()

	evk.RelinearizationKey = kgen.GenRelinearizationKeyNew(sk)

	for _, galEl := range btpParams.GaloisElements(ckksParams) {
		evk.GaloisKeys[galEl] = kgen.GenGaloisKeyNew(galEl, sk)
	}

	evk.GaloisKeys[ckksParams.GaloisElementInverse()] = kgen.GenGaloisKeyNew(ckksParams.GaloisElementInverse(), sk)

	EvkDtS, EvkStD := btpParams.GenEncapsulationEvaluationKeysNew(ckksParams, sk)

	return &EvaluationKeySet{
		EvaluationKeySet: evk,
		EvkDtS:           EvkDtS,
		EvkStD:           EvkStD,
	}
}

// GenEncapsulationEvaluationKeysNew generates the low level encapsulation EvaluationKeys for the bootstrapping.
func (p *Parameters) GenEncapsulationEvaluationKeysNew(params ckks.Parameters, skDense *rlwe.SecretKey) (EvkDtS, EvkStD *rlwe.EvaluationKey) {

	if p.EphemeralSecretWeight == 0 {
		return
	}

	paramsSparse, _ := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		LogN: params.LogN(),
		Q:    params.Q()[:1],
		P:    params.P()[:1],
	})

	kgenSparse := rlwe.NewKeyGenerator(paramsSparse)
	kgenDense := rlwe.NewKeyGenerator(params.Parameters)
	skSparse := kgenSparse.GenSecretKeyWithHammingWeightNew(p.EphemeralSecretWeight)

	return kgenDense.GenEvaluationKeyNew(skDense, skSparse), kgenDense.GenEvaluationKeyNew(skSparse, skDense)
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
func (bb *bootstrapperBase) CheckKeys(btpKeys *EvaluationKeySet) (err error) {

	if _, err = btpKeys.GetRelinearizationKey(); err != nil {
		return
	}

	for _, galEl := range bb.GaloisElements(bb.params) {
		if _, err = btpKeys.GetGaloisKey(galEl); err != nil {
			return
		}
	}

	if btpKeys.EvkDtS == nil && bb.Parameters.EphemeralSecretWeight != 0 {
		return fmt.Errorf("rlwe.EvaluationKey key dense to sparse is nil")
	}

	if btpKeys.EvkStD == nil && bb.Parameters.EphemeralSecretWeight != 0 {
		return fmt.Errorf("rlwe.EvaluationKey key sparse to dense is nil")
	}

	return
}

func newBootstrapperBase(params ckks.Parameters, btpParams Parameters, btpKey *EvaluationKeySet) (bb *bootstrapperBase) {
	bb = new(bootstrapperBase)
	bb.params = params
	bb.Parameters = btpParams

	bb.logdslots = btpParams.PlaintextLogDimensions()[1]
	bb.dslots = 1 << bb.logdslots
	if maxLogSlots := params.PlaintextLogDimensions()[1]; bb.dslots < maxLogSlots {
		bb.dslots <<= 1
		bb.logdslots++
	}

	bb.evalModPoly = ckks.NewEvalModPolyFromLiteral(params, btpParams.EvalModParameters)

	scFac := bb.evalModPoly.ScFac()
	K := bb.evalModPoly.K() / scFac

	// Correcting factor for approximate division by Q
	// The second correcting factor for approximate multiplication by Q is included in the coefficients of the EvalMod polynomials
	qDiff := bb.evalModPoly.QDiff()

	Q0 := params.Q()[0]

	// Q0/|m|
	bb.q0OverMessageRatio = math.Exp2(math.Round(math.Log2(float64(Q0) / bb.evalModPoly.MessageRatio())))

	// If the scale used during the EvalMod step is smaller than Q0, then we cannot increase the scale during
	// the EvalMod step to get a free division by MessageRatio, and we need to do this division (totally or partly)
	// during the CoeffstoSlots step
	qDiv := bb.evalModPoly.ScalingFactor().Float64() / math.Exp2(math.Round(math.Log2(float64(Q0))))

	// Sets qDiv to 1 if there is enough room for the division to happen using scale manipulation.
	if qDiv > 1 {
		qDiv = 1
	}

	encoder := ckks.NewEncoder(bb.params)

	// CoeffsToSlots vectors
	// Change of variable for the evaluation of the Chebyshev polynomial + cancelling factor for the DFT and SubSum + eventual scaling factor for the double angle formula

	if bb.CoeffsToSlotsParameters.Scaling == nil {
		bb.CoeffsToSlotsParameters.Scaling = new(big.Float).SetFloat64(qDiv / (K * scFac * qDiff))
	} else {
		bb.CoeffsToSlotsParameters.Scaling.Mul(bb.CoeffsToSlotsParameters.Scaling, new(big.Float).SetFloat64(qDiv/(K*scFac*qDiff)))
	}

	bb.ctsMatrices = ckks.NewHomomorphicDFTMatrixFromLiteral(bb.CoeffsToSlotsParameters, encoder)

	// SlotsToCoeffs vectors
	// Rescaling factor to set the final ciphertext to the desired scale

	if bb.SlotsToCoeffsParameters.Scaling == nil {
		bb.SlotsToCoeffsParameters.Scaling = new(big.Float).SetFloat64(bb.params.PlaintextScale().Float64() / (bb.evalModPoly.ScalingFactor().Float64() / bb.evalModPoly.MessageRatio()) * qDiff)
	} else {
		bb.SlotsToCoeffsParameters.Scaling.Mul(bb.SlotsToCoeffsParameters.Scaling, new(big.Float).SetFloat64(bb.params.PlaintextScale().Float64()/(bb.evalModPoly.ScalingFactor().Float64()/bb.evalModPoly.MessageRatio())*qDiff))
	}

	bb.stcMatrices = ckks.NewHomomorphicDFTMatrixFromLiteral(bb.SlotsToCoeffsParameters, encoder)

	encoder = nil

	return
}
