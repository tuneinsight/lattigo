package bootstrapping

import (
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/circuits/float"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// Bootstrapper is a struct to store a memory buffer with the plaintext matrices,
// the polynomial approximation, and the keys for the bootstrapping.
type Bootstrapper struct {
	*ckks.Evaluator
	*float.DFTEvaluator
	*float.Mod1Evaluator
	*bootstrapperBase
}

type bootstrapperBase struct {
	Parameters
	*EvaluationKeySet
	params ckks.Parameters

	dslots    int // Number of plaintext slots after the re-encoding
	logdslots int

	mod1Parameters float.Mod1Parameters
	stcMatrices    float.DFTMatrix
	ctsMatrices    float.DFTMatrix

	q0OverMessageRatio float64
}

// EvaluationKeySet is a type for a CKKS bootstrapping key, which
// regroups the necessary public relinearization and rotation keys.
type EvaluationKeySet struct {
	*rlwe.MemEvaluationKeySet
	EvkDtS *rlwe.EvaluationKey
	EvkStD *rlwe.EvaluationKey
}

// NewBootstrapper creates a new Bootstrapper.
func NewBootstrapper(btpParams Parameters, btpKeys *EvaluationKeySet) (btp *Bootstrapper, err error) {

	if btpParams.Mod1ParametersLiteral.SineType == float.SinContinuous && btpParams.Mod1ParametersLiteral.DoubleAngle != 0 {
		return nil, fmt.Errorf("cannot use double angle formul for SineType = Sin -> must use SineType = Cos")
	}

	if btpParams.Mod1ParametersLiteral.SineType == float.CosDiscrete && btpParams.Mod1ParametersLiteral.SineDegree < 2*(btpParams.Mod1ParametersLiteral.K-1) {
		return nil, fmt.Errorf("SineType 'ckks.CosDiscrete' uses a minimum degree of 2*(K-1) but EvalMod degree is smaller")
	}

	if btpParams.CoeffsToSlotsParameters.LevelStart-btpParams.CoeffsToSlotsParameters.Depth(true) != btpParams.Mod1ParametersLiteral.LevelStart {
		return nil, fmt.Errorf("starting level and depth of CoeffsToSlotsParameters inconsistent starting level of SineEvalParameters")
	}

	if btpParams.Mod1ParametersLiteral.LevelStart-btpParams.Mod1ParametersLiteral.Depth() != btpParams.SlotsToCoeffsParameters.LevelStart {
		return nil, fmt.Errorf("starting level and depth of SineEvalParameters inconsistent starting level of CoeffsToSlotsParameters")
	}

	params := btpParams.Parameters

	btp = new(Bootstrapper)
	if btp.bootstrapperBase, err = newBootstrapperBase(params, btpParams, btpKeys); err != nil {
		return
	}

	if err = btp.bootstrapperBase.CheckKeys(btpKeys); err != nil {
		return nil, fmt.Errorf("invalid bootstrapping key: %w", err)
	}

	btp.EvaluationKeySet = btpKeys

	btp.Evaluator = ckks.NewEvaluator(params, btpKeys)

	btp.DFTEvaluator = float.NewDFTEvaluator(params, btp.Evaluator)

	btp.Mod1Evaluator = float.NewMod1Evaluator(btp.Evaluator, float.NewPolynomialEvaluator(params, btp.Evaluator), btp.bootstrapperBase.mod1Parameters)

	return
}

// GenEvaluationKeySetNew generates a new bootstrapping EvaluationKeySet, which contain:
//
//	EvaluationKeySet: struct compliant to the interface rlwe.EvaluationKeySetInterface.
//	EvkDtS: *rlwe.EvaluationKey
//	EvkStD: *rlwe.EvaluationKey
func (p Parameters) GenEvaluationKeySetNew(sk *rlwe.SecretKey) *EvaluationKeySet {

	ringQ := p.Parameters.RingQ()
	ringP := p.Parameters.RingP()

	if sk.Value.Q.N() != ringQ.N() {
		panic(fmt.Sprintf("invalid secret key: secret key ring degree = %d does not match bootstrapping parameters ring degree = %d", sk.Value.Q.N(), ringQ.N()))
	}

	params := p.Parameters

	skExtended := rlwe.NewSecretKey(params)
	buff := ringQ.NewPoly()

	// Extends basis Q0 -> QL
	rlwe.ExtendBasisSmallNormAndCenterNTTMontgomery(ringQ, ringQ, sk.Value.Q, buff, skExtended.Value.Q)

	// Extends basis Q0 -> P
	rlwe.ExtendBasisSmallNormAndCenterNTTMontgomery(ringQ, ringP, sk.Value.Q, buff, skExtended.Value.P)

	kgen := ckks.NewKeyGenerator(params)

	EvkDtS, EvkStD := p.GenEncapsulationEvaluationKeysNew(skExtended)

	rlk := kgen.GenRelinearizationKeyNew(skExtended)
	gks := kgen.GenGaloisKeysNew(append(p.GaloisElements(params), params.GaloisElementForComplexConjugation()), skExtended)

	evk := rlwe.NewMemEvaluationKeySet(rlk, gks...)
	return &EvaluationKeySet{
		MemEvaluationKeySet: evk,
		EvkDtS:              EvkDtS,
		EvkStD:              EvkStD,
	}
}

// GenEncapsulationEvaluationKeysNew generates the low level encapsulation EvaluationKeys for the bootstrapping.
func (p Parameters) GenEncapsulationEvaluationKeysNew(skDense *rlwe.SecretKey) (EvkDtS, EvkStD *rlwe.EvaluationKey) {

	params := p.Parameters

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

	EvkDtS = kgenDense.GenEvaluationKeyNew(skDense, skSparse)
	EvkStD = kgenDense.GenEvaluationKeyNew(skSparse, skDense)
	return
}

// ShallowCopy creates a shallow copy of this Bootstrapper in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Bootstrapper can be used concurrently.
func (btp Bootstrapper) ShallowCopy() *Bootstrapper {
	return &Bootstrapper{
		Evaluator:        btp.Evaluator.ShallowCopy(),
		bootstrapperBase: btp.bootstrapperBase,
		//DFTEvaluator: btp.DFTEvaluator.ShallowCopy(),
		//Mod1Evaluator: btp.Mod1Evaluator.ShallowCopy(),

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

func newBootstrapperBase(params ckks.Parameters, btpParams Parameters, btpKey *EvaluationKeySet) (bb *bootstrapperBase, err error) {
	bb = new(bootstrapperBase)
	bb.params = params
	bb.Parameters = btpParams

	bb.logdslots = btpParams.LogMaxDimensions().Cols
	bb.dslots = 1 << bb.logdslots
	if maxLogSlots := params.LogMaxDimensions().Cols; bb.dslots < maxLogSlots {
		bb.dslots <<= 1
		bb.logdslots++
	}

	if bb.mod1Parameters, err = float.NewMod1ParametersFromLiteral(params, btpParams.Mod1ParametersLiteral); err != nil {
		return nil, err
	}

	scFac := bb.mod1Parameters.ScFac()
	K := bb.mod1Parameters.K() / scFac

	// Correcting factor for approximate division by Q
	// The second correcting factor for approximate multiplication by Q is included in the coefficients of the EvalMod polynomials
	qDiff := bb.mod1Parameters.QDiff()

	Q0 := params.Q()[0]

	// Q0/|m|
	bb.q0OverMessageRatio = math.Exp2(math.Round(math.Log2(float64(Q0) / bb.mod1Parameters.MessageRatio())))

	// If the scale used during the EvalMod step is smaller than Q0, then we cannot increase the scale during
	// the EvalMod step to get a free division by MessageRatio, and we need to do this division (totally or partly)
	// during the CoeffstoSlots step
	qDiv := bb.mod1Parameters.ScalingFactor().Float64() / math.Exp2(math.Round(math.Log2(float64(Q0))))

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

	if bb.ctsMatrices, err = float.NewDFTMatrixFromLiteral(params, bb.CoeffsToSlotsParameters, encoder); err != nil {
		return
	}

	// SlotsToCoeffs vectors
	// Rescaling factor to set the final ciphertext to the desired scale

	if bb.SlotsToCoeffsParameters.Scaling == nil {
		bb.SlotsToCoeffsParameters.Scaling = new(big.Float).SetFloat64(bb.params.DefaultScale().Float64() / (bb.mod1Parameters.ScalingFactor().Float64() / bb.mod1Parameters.MessageRatio()) * qDiff)
	} else {
		bb.SlotsToCoeffsParameters.Scaling.Mul(bb.SlotsToCoeffsParameters.Scaling, new(big.Float).SetFloat64(bb.params.DefaultScale().Float64()/(bb.mod1Parameters.ScalingFactor().Float64()/bb.mod1Parameters.MessageRatio())*qDiff))
	}

	if bb.stcMatrices, err = float.NewDFTMatrixFromLiteral(params, bb.SlotsToCoeffsParameters, encoder); err != nil {
		return
	}

	encoder = nil

	return
}
