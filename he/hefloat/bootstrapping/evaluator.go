package bootstrapping

import (
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/schemes/ckks"
)

// Evaluator is a struct to store a memory buffer with the plaintext matrices,
// the polynomial approximation, and the keys for the bootstrapping.
// It is used to evaluate the bootstrapping circuit on single ciphertexts.
type Evaluator struct {
	Parameters
	*hefloat.Evaluator
	*hefloat.DFTEvaluator
	*hefloat.Mod1Evaluator
	*EvaluationKeys

	ckks.DomainSwitcher

	// [1, x, x^2, x^4, ..., x^N1/2] / (X^N1 +1)
	xPow2N1 []ring.Poly
	// [1, x, x^2, x^4, ..., x^N2/2] / (X^N2 +1)
	xPow2N2 []ring.Poly
	// [1, x^-1, x^-2, x^-4, ..., x^-N2/2] / (X^N2 +1)
	xPow2InvN2 []ring.Poly

	Mod1Parameters hefloat.Mod1Parameters
	S2CDFTMatrix   hefloat.DFTMatrix
	C2SDFTMatrix   hefloat.DFTMatrix

	SkDebug *rlwe.SecretKey
}

// NewEvaluator creates a new Evaluator.
func NewEvaluator(btpParams Parameters, evk *EvaluationKeys) (eval *Evaluator, err error) {

	eval = &Evaluator{}

	paramsN1 := btpParams.ResidualParameters
	paramsN2 := btpParams.BootstrappingParameters

	switch paramsN1.RingType() {
	case ring.Standard:
		if paramsN1.N() != paramsN2.N() && (evk.EvkN1ToN2 == nil || evk.EvkN2ToN1 == nil) {
			return nil, fmt.Errorf("cannot NewBootstrapper: evk.(BootstrappingKeys) is missing EvkN1ToN2 and EvkN2ToN1")
		}
	case ring.ConjugateInvariant:
		if evk.EvkCmplxToReal == nil || evk.EvkRealToCmplx == nil {
			return nil, fmt.Errorf("cannot NewBootstrapper: evk.(BootstrappingKeys) is missing EvkN1ToN2 and EvkN2ToN1")
		}

		var err error
		if eval.DomainSwitcher, err = ckks.NewDomainSwitcher(paramsN2.Parameters, evk.EvkCmplxToReal, evk.EvkRealToCmplx); err != nil {
			return nil, fmt.Errorf("cannot NewBootstrapper: ckks.NewDomainSwitcher: %w", err)
		}

		// The switch to standard to conjugate invariant multiplies the scale by 2
		btpParams.SlotsToCoeffsParameters.Scaling = new(big.Float).SetFloat64(0.5)
	}

	eval.Parameters = btpParams

	if paramsN1.N() != paramsN2.N() {
		eval.xPow2N1 = rlwe.GenXPow2(paramsN1.RingQ().AtLevel(0), paramsN2.LogN(), false)
		eval.xPow2N2 = rlwe.GenXPow2(paramsN2.RingQ().AtLevel(0), paramsN2.LogN(), false)
		eval.xPow2InvN2 = rlwe.GenXPow2(paramsN2.RingQ(), paramsN2.LogN(), true)
	}

	if btpParams.Mod1ParametersLiteral.Mod1Type == hefloat.SinContinuous && btpParams.Mod1ParametersLiteral.DoubleAngle != 0 {
		return nil, fmt.Errorf("cannot use double angle formula for Mod1Type = Sin -> must use Mod1Type = Cos")
	}

	if btpParams.Mod1ParametersLiteral.Mod1Type == hefloat.CosDiscrete && btpParams.Mod1ParametersLiteral.Mod1Degree < 2*(btpParams.Mod1ParametersLiteral.K-1) {
		return nil, fmt.Errorf("Mod1Type 'hefloat.CosDiscrete' uses a minimum degree of 2*(K-1) but EvalMod degree is smaller")
	}

	switch btpParams.CircuitOrder {
	case ModUpThenEncode:
		if btpParams.CoeffsToSlotsParameters.LevelStart-btpParams.CoeffsToSlotsParameters.Depth(true) != btpParams.Mod1ParametersLiteral.LevelStart {
			return nil, fmt.Errorf("starting level and depth of CoeffsToSlotsParameters inconsistent starting level of Mod1ParametersLiteral")
		}

		if btpParams.Mod1ParametersLiteral.LevelStart-btpParams.Mod1ParametersLiteral.Depth() != btpParams.SlotsToCoeffsParameters.LevelStart {
			return nil, fmt.Errorf("starting level and depth of Mod1ParametersLiteral inconsistent starting level of CoeffsToSlotsParameters")
		}
	case DecodeThenModUp:
		if btpParams.BootstrappingParameters.MaxLevel()-btpParams.CoeffsToSlotsParameters.Depth(true) != btpParams.Mod1ParametersLiteral.LevelStart {
			return nil, fmt.Errorf("starting level and depth of Mod1ParametersLiteral inconsistent starting level of CoeffsToSlotsParameters")
		}
	case Custom:
	default:
		return nil, fmt.Errorf("invalid CircuitOrder value")
	}

	if err = eval.initialize(btpParams); err != nil {
		return
	}

	if err = eval.checkKeys(evk); err != nil {
		return
	}

	params := btpParams.BootstrappingParameters

	eval.EvaluationKeys = evk

	eval.Evaluator = hefloat.NewEvaluator(params, evk)

	eval.DFTEvaluator = hefloat.NewDFTEvaluator(params, eval.Evaluator)

	eval.Mod1Evaluator = hefloat.NewMod1Evaluator(eval.Evaluator, hefloat.NewPolynomialEvaluator(params, eval.Evaluator), eval.Mod1Parameters)

	return
}

// ShallowCopy creates a shallow copy of this Evaluator in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Evaluator can be used concurrently.
func (eval Evaluator) ShallowCopy() *Evaluator {
	heEvaluator := eval.Evaluator.ShallowCopy()

	paramsN1 := eval.ResidualParameters
	paramsN2 := eval.BootstrappingParameters

	var DomainSwitcher ckks.DomainSwitcher
	if paramsN1.RingType() == ring.ConjugateInvariant {
		var err error
		if DomainSwitcher, err = ckks.NewDomainSwitcher(paramsN2.Parameters, eval.EvkCmplxToReal, eval.EvkRealToCmplx); err != nil {
			panic(fmt.Errorf("cannot NewBootstrapper: ckks.NewDomainSwitcher: %w", err))
		}
	}

	return &Evaluator{
		Parameters:     eval.Parameters,
		EvaluationKeys: eval.EvaluationKeys,
		Mod1Parameters: eval.Mod1Parameters,
		S2CDFTMatrix:   eval.S2CDFTMatrix,
		C2SDFTMatrix:   eval.C2SDFTMatrix,
		Evaluator:      heEvaluator,
		xPow2N1:        eval.xPow2N1,
		xPow2N2:        eval.xPow2N2,
		xPow2InvN2:     eval.xPow2InvN2,
		DomainSwitcher: DomainSwitcher,
		DFTEvaluator:   hefloat.NewDFTEvaluator(paramsN2, heEvaluator),
		Mod1Evaluator:  hefloat.NewMod1Evaluator(heEvaluator, hefloat.NewPolynomialEvaluator(paramsN2, heEvaluator), eval.Mod1Parameters),
		SkDebug:        eval.SkDebug,
	}
}

// CheckKeys checks if all the necessary keys are present in the instantiated Evaluator
func (eval Evaluator) checkKeys(evk *EvaluationKeys) (err error) {

	if _, err = evk.GetRelinearizationKey(); err != nil {
		return
	}

	for _, galEl := range eval.GaloisElements(eval.BootstrappingParameters) {
		if _, err = evk.GetGaloisKey(galEl); err != nil {
			return
		}
	}

	if evk.EvkDenseToSparse == nil && eval.EphemeralSecretWeight != 0 {
		return fmt.Errorf("rlwe.EvaluationKey key dense to sparse is nil")
	}

	if evk.EvkSparseToDense == nil && eval.EphemeralSecretWeight != 0 {
		return fmt.Errorf("rlwe.EvaluationKey key sparse to dense is nil")
	}

	return
}

func (eval *Evaluator) initialize(btpParams Parameters) (err error) {

	eval.Parameters = btpParams
	params := btpParams.BootstrappingParameters

	if eval.Mod1Parameters, err = hefloat.NewMod1ParametersFromLiteral(params, btpParams.Mod1ParametersLiteral); err != nil {
		return
	}

	scFac := eval.Mod1Parameters.ScFac()
	K := eval.Mod1Parameters.K() / scFac

	// Correcting factor for approximate division by Q
	// The second correcting factor for approximate multiplication by Q is included in the coefficients of the EvalMod polynomials
	qDiff := eval.Mod1Parameters.QDiff()

	// If the scale used during the EvalMod step is smaller than Q0, then we cannot increase the scale during
	// the EvalMod step to get a free division by MessageRatio, and we need to do this division (totally or partly)
	// during the CoeffstoSlots step
	qDiv := eval.Mod1Parameters.ScalingFactor().Float64() / math.Exp2(math.Round(math.Log2(float64(params.Q()[0]))))

	// Sets qDiv to 1 if there is enough room for the division to happen using scale manipulation.
	if qDiv > 1 {
		qDiv = 1
	}

	encoder := hefloat.NewEncoder(params)

	// CoeffsToSlots vectors
	// Change of variable for the evaluation of the Chebyshev polynomial + cancelling factor for the DFT and SubSum + eventual scaling factor for the double angle formula

	scale := eval.BootstrappingParameters.DefaultScale().Float64()
	offset := eval.Mod1Parameters.ScalingFactor().Float64() / eval.Mod1Parameters.MessageRatio()

	C2SScaling := new(big.Float).SetFloat64(qDiv / (K * scFac * qDiff))
	StCScaling := new(big.Float).SetFloat64(scale / offset)

	if eval.CoeffsToSlotsParameters.Scaling == nil {
		eval.CoeffsToSlotsParameters.Scaling = C2SScaling
	} else {
		eval.CoeffsToSlotsParameters.Scaling.Mul(eval.CoeffsToSlotsParameters.Scaling, C2SScaling)
	}

	if eval.SlotsToCoeffsParameters.Scaling == nil {
		eval.SlotsToCoeffsParameters.Scaling = StCScaling
	} else {
		eval.SlotsToCoeffsParameters.Scaling.Mul(eval.SlotsToCoeffsParameters.Scaling, StCScaling)
	}

	if eval.C2SDFTMatrix, err = hefloat.NewDFTMatrixFromLiteral(params, eval.CoeffsToSlotsParameters, encoder); err != nil {
		return
	}

	if eval.S2CDFTMatrix, err = hefloat.NewDFTMatrixFromLiteral(params, eval.SlotsToCoeffsParameters, encoder); err != nil {
		return
	}

	encoder = nil // For the GC

	return
}
