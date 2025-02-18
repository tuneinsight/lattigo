package bootstrapping

import (
	"fmt"
	"math"
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// Evaluator is a struct to store a memory buffer with the plaintext matrices,
// the polynomial approximation, and the keys for the bootstrapping.
// It is used to evaluate the bootstrapping circuit on single ciphertexts.
type Evaluator struct {
	Parameters
	*ckks.Evaluator
	DFTEvaluator  *dft.Evaluator
	Mod1Evaluator *mod1.Evaluator
	*EvaluationKeys

	ckks.DomainSwitcher

	// [1, x, x^2, x^4, ..., x^N1/2] / (X^N1 +1)
	xPow2N1 []ring.Poly
	// [1, x, x^2, x^4, ..., x^N2/2] / (X^N2 +1)
	xPow2N2 []ring.Poly
	// [1, x^-1, x^-2, x^-4, ..., x^-N1/2] / (X^N1 +1)
	xPow2InvN1 []ring.Poly
	// [1, x^-1, x^-2, x^-4, ..., x^-N2/2] / (X^N2 +1)
	xPow2InvN2 []ring.Poly

	Mod1Parameters mod1.Parameters
	S2CDFTMatrix   dft.Matrix
	C2SDFTMatrix   dft.Matrix

	SkDebug *rlwe.SecretKey
}

// NewEvaluator creates a new [Evaluator].
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
		if eval.DomainSwitcher, err = ckks.NewDomainSwitcher(paramsN2, evk.EvkCmplxToReal, evk.EvkRealToCmplx); err != nil {
			return nil, fmt.Errorf("cannot NewBootstrapper: ckks.NewDomainSwitcher: %w", err)
		}

		// The switch to standard to conjugate invariant multiplies the scale by 2
		btpParams.SlotsToCoeffsParameters.Scaling = new(big.Float).SetFloat64(0.5)
	}

	eval.Parameters = btpParams

	if paramsN1.N() != paramsN2.N() {
		eval.xPow2N1 = rlwe.GenXPow2NTT(paramsN1.RingQ().AtLevel(0), paramsN2.LogN(), false)
		eval.xPow2InvN1 = rlwe.GenXPow2NTT(paramsN1.RingQ(), paramsN1.LogN(), true)
	}
	eval.xPow2N2 = rlwe.GenXPow2NTT(paramsN2.RingQ().AtLevel(0), paramsN2.LogN(), false)
	eval.xPow2InvN2 = rlwe.GenXPow2NTT(paramsN2.RingQ(), paramsN2.LogN(), true)

	if btpParams.Mod1ParametersLiteral.Mod1Type == mod1.SinContinuous && btpParams.Mod1ParametersLiteral.DoubleAngle != 0 {
		return nil, fmt.Errorf("cannot use double angle formula for Mod1Type = Sin -> must use Mod1Type = Cos")
	}

	if btpParams.Mod1ParametersLiteral.Mod1Type == mod1.CosDiscrete && btpParams.Mod1ParametersLiteral.Mod1Degree < 2*(btpParams.Mod1ParametersLiteral.K-1) {
		return nil, fmt.Errorf("Mod1Type 'mod1.CosDiscrete' uses a minimum degree of 2*(K-1) but EvalMod degree is smaller")
	}

	switch btpParams.CircuitOrder {
	case ModUpThenEncode:
		if btpParams.CoeffsToSlotsParameters.LevelQ-btpParams.CoeffsToSlotsParameters.Depth(true) != btpParams.Mod1ParametersLiteral.LevelQ {
			return nil, fmt.Errorf("starting level and depth of CoeffsToSlotsParameters inconsistent starting level of Mod1ParametersLiteral")
		}

		if btpParams.Mod1ParametersLiteral.LevelQ-btpParams.Mod1ParametersLiteral.Depth() != btpParams.SlotsToCoeffsParameters.LevelQ {
			return nil, fmt.Errorf("starting level and depth of Mod1ParametersLiteral inconsistent starting level of CoeffsToSlotsParameters")
		}
	case DecodeThenModUp:
		if btpParams.BootstrappingParameters.MaxLevel()-btpParams.CoeffsToSlotsParameters.Depth(true) != btpParams.Mod1ParametersLiteral.LevelQ {
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

	eval.Evaluator = ckks.NewEvaluator(params, evk)

	eval.DFTEvaluator = dft.NewEvaluator(params, eval.Evaluator)

	eval.Mod1Evaluator = mod1.NewEvaluator(eval.Evaluator, polynomial.NewEvaluator(params, eval.Evaluator), eval.Mod1Parameters)

	return
}

// ShallowCopy creates a shallow copy of this [Evaluator] in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Evaluator can be used concurrently.
func (eval Evaluator) ShallowCopy() *Evaluator {
	heEvaluator := eval.Evaluator.ShallowCopy()

	paramsN1 := eval.ResidualParameters
	paramsN2 := eval.BootstrappingParameters

	var DomainSwitcher ckks.DomainSwitcher
	if paramsN1.RingType() == ring.ConjugateInvariant {
		var err error
		if DomainSwitcher, err = ckks.NewDomainSwitcher(paramsN2, eval.EvkCmplxToReal, eval.EvkRealToCmplx); err != nil {
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
		DFTEvaluator:   dft.NewEvaluator(paramsN2, heEvaluator),
		Mod1Evaluator:  mod1.NewEvaluator(heEvaluator, polynomial.NewEvaluator(paramsN2, heEvaluator), eval.Mod1Parameters),
		SkDebug:        eval.SkDebug,
	}
}

// CheckKeys checks if all the necessary keys are present in the instantiated [Evaluator]
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

	if eval.Mod1Parameters, err = mod1.NewParametersFromLiteral(params, btpParams.Mod1ParametersLiteral); err != nil {
		return
	}

	// [-K, K]
	K := eval.Mod1Parameters.K

	// Correcting factor for approximate division by Q
	// The second correcting factor for approximate multiplication by Q is included in the coefficients of the EvalMod polynomials
	qDiff := eval.Mod1Parameters.QDiff

	// If the scale used during the EvalMod step is smaller than Q0, then we cannot increase the scale during
	// the EvalMod step to get a free division by MessageRatio, and we need to do this division (totally or partly)
	// during the CoeffstoSlots step
	qDiv := eval.Mod1Parameters.ScalingFactor().Float64() / math.Exp2(math.Round(math.Log2(float64(params.Q()[0]))))

	// Sets qDiv to 1 if there is enough room for the division to happen using scale manipulation.
	if qDiv > 1 {
		qDiv = 1
	}

	encoder := ckks.NewEncoder(params)

	// CoeffsToSlots vectors
	// Change of variable for the evaluation of the Chebyshev polynomial + cancelling factor for the DFT and SubSum + eventual scaling factor for the double angle formula

	scale := eval.BootstrappingParameters.DefaultScale().Float64()
	offset := eval.Mod1Parameters.ScalingFactor().Float64() / eval.Mod1Parameters.MessageRatio()

	C2SScaling := new(big.Float).SetFloat64(qDiv / (K * qDiff))
	StCScaling := new(big.Float).SetFloat64(scale / offset)

	if btpParams.CoeffsToSlotsParameters.Scaling == nil {
		eval.CoeffsToSlotsParameters.Scaling = C2SScaling
	} else {
		eval.CoeffsToSlotsParameters.Scaling = new(big.Float).Mul(btpParams.CoeffsToSlotsParameters.Scaling, C2SScaling)
	}

	if btpParams.SlotsToCoeffsParameters.Scaling == nil {
		eval.SlotsToCoeffsParameters.Scaling = StCScaling
	} else {
		eval.SlotsToCoeffsParameters.Scaling = new(big.Float).Mul(btpParams.SlotsToCoeffsParameters.Scaling, StCScaling)
	}

	if eval.C2SDFTMatrix, err = dft.NewMatrixFromLiteral(params, eval.CoeffsToSlotsParameters, encoder); err != nil {
		return
	}

	if eval.S2CDFTMatrix, err = dft.NewMatrixFromLiteral(params, eval.SlotsToCoeffsParameters, encoder); err != nil {
		return
	}

	encoder = nil // For the GC

	return
}

// Bootstrap bootstraps a single ciphertext and returns the bootstrapped ciphertext.
func (eval Evaluator) Bootstrap(ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	cts := []rlwe.Ciphertext{*ct}
	cts, err := eval.BootstrapMany(cts)
	if err != nil {
		return nil, err
	}
	return &cts[0], nil
}

// BootstrapMany bootstraps a list of ciphertexts and returns the list of bootstrapped ciphertexts.
func (eval Evaluator) BootstrapMany(cts []rlwe.Ciphertext) ([]rlwe.Ciphertext, error) {

	var err error

	switch eval.ResidualParameters.RingType() {
	case ring.ConjugateInvariant:

		for i := 0; i < len(cts); i = i + 2 {

			even, odd := i, i+1

			ct0 := &cts[even]

			var ct1 *rlwe.Ciphertext
			if odd < len(cts) {
				ct1 = &cts[odd]
			}

			if ct0, ct1, err = eval.EvaluateConjugateInvariant(ct0, ct1); err != nil {
				return nil, fmt.Errorf("cannot BootstrapMany: %w", err)
			}

			cts[even] = *ct0

			if ct1 != nil {
				cts[odd] = *ct1
			}
		}

	default:

		ctsPacked, ctxt1, ctxt2, err := eval.PackAndSwitchN1ToN2(cts)
		if err != nil {
			return nil, fmt.Errorf("cannot BootstrapMany: %w", err)
		}

		cts = ctsPacked
		for i := range cts {
			var ct *rlwe.Ciphertext
			if ct, err = eval.Evaluate(&cts[i]); err != nil {
				return nil, fmt.Errorf("cannot BootstrapMany: %w", err)
			}
			cts[i] = *ct
		}

		if cts, err = eval.UnpackAndSwitchN2ToN1(cts, ctxt1, ctxt2); err != nil {
			return nil, fmt.Errorf("cannot BootstrapMany: %w", err)
		}
	}

	for i := range cts {
		cts[i].Scale = eval.ResidualParameters.DefaultScale()
	}

	return cts, err
}

// Depth returns the multiplicative depth (number of levels consumed) of the bootstrapping circuit.
func (eval Evaluator) Depth() int {
	return eval.BootstrappingParameters.MaxLevel() - eval.ResidualParameters.MaxLevel()
}

// OutputLevel returns the output level after the evaluation of the bootstrapping circuit.
func (eval Evaluator) OutputLevel() int {
	return eval.ResidualParameters.MaxLevel()
}

// MinimumInputLevel returns the minimum level at which a ciphertext must be to be bootstrapped.
func (eval Evaluator) MinimumInputLevel() int {
	return eval.BootstrappingParameters.LevelsConsumedPerRescaling()
}

// Evaluate re-encrypts a ciphertext to a ciphertext at MaxLevel - k where k is the depth of the bootstrapping circuit.
// If the input ciphertext level is zero, the input scale must be an exact power of two smaller than Q[0]/MessageRatio
// (it can't be equal since Q[0] is not a power of two).
// The message ratio is an optional field in the bootstrapping parameters, by default it set to 2^{LogMessageRatio = 8}.
// See the bootstrapping parameters for more information about the message ratio or other parameters related to the bootstrapping.
// If the input ciphertext is at level one or more, the input scale does not need to be an exact power of two as one level
// can be used to do a scale matching.
//
// The circuit consists in 5 steps.
//  1. ScaleDown: scales the ciphertext to q/|m| and bringing it down to q
//  2. ModUp: brings the modulus from q to Q
//  3. CoeffsToSlots: homomorphic encoding
//  4. EvalMod: homomorphic modular reduction
//  5. SlotsToCoeffs: homomorphic decoding
func (eval Evaluator) Evaluate(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext, err error) {

	if eval.IterationsParameters == nil && eval.ResidualParameters.PrecisionMode() != ckks.PREC128 {
		ctOut, _, err = eval.bootstrap(ctIn)
		return

	} else {

		var errScale *rlwe.Scale
		// [M^{d}/q1 + e^{d-logprec}]
		if ctOut, errScale, err = eval.bootstrap(ctIn.CopyNew()); err != nil {
			return nil, err
		}

		// Stores by how much a ciphertext must be scaled to get back
		// to the input scale
		// Error correcting factor of the approximate division by q1
		// diffScale = ctIn.Scale / (ctOut.Scale * errScale)
		diffScale := ctIn.Scale.Div(ctOut.Scale)
		diffScale = diffScale.Div(*errScale)

		// [M^{d} + e^{d-logprec}]
		if err = eval.Evaluator.Mul(ctOut, diffScale.BigInt(), ctOut); err != nil {
			return nil, err
		}
		ctOut.Scale = ctIn.Scale

		if eval.IterationsParameters != nil {

			QiReserved := eval.BootstrappingParameters.Q()[eval.ResidualParameters.MaxLevel()+1]

			var totLogPrec float64

			for i := 0; i < len(eval.IterationsParameters.BootstrappingPrecision); i++ {

				logPrec := eval.IterationsParameters.BootstrappingPrecision[i]

				totLogPrec += logPrec

				// prec = round(2^{logprec})
				log2 := bignum.Log(new(big.Float).SetPrec(256).SetUint64(2))
				log2TimesLogPrec := log2.Mul(log2, new(big.Float).SetFloat64(totLogPrec))
				prec := new(big.Int)
				log2TimesLogPrec.Add(bignum.Exp(log2TimesLogPrec), new(big.Float).SetFloat64(0.5)).Int(prec)

				// Corrects the last iteration 2^{logprec} such that diffScale / prec * QReserved is as close to an integer as possible.
				// This is necessary to not lose bits of precision during the last iteration is a reserved prime is used.
				// If this correct is not done, what can happen is that there is a loss of up to 2^{logprec/2} bits from the last iteration.
				if eval.IterationsParameters.ReservedPrimeBitSize != 0 && i == len(eval.IterationsParameters.BootstrappingPrecision)-1 {

					// 1) Computes the scale = diffScale / prec * QReserved
					scale := new(big.Float).Quo(&diffScale.Value, new(big.Float).SetInt(prec))
					scale.Mul(scale, new(big.Float).SetUint64(QiReserved))

					// 2) Finds the closest integer to scale with scale = round(scale)
					scale.Add(scale, new(big.Float).SetFloat64(0.5))
					tmp := new(big.Int)
					scale.Int(tmp)
					scale.SetInt(tmp)

					// 3) Computes the corrected precision = diffScale * QReserved / round(scale)
					preccorrected := new(big.Float).Quo(&diffScale.Value, scale)
					preccorrected.Mul(preccorrected, new(big.Float).SetUint64(QiReserved))
					preccorrected.Add(preccorrected, new(big.Float).SetFloat64(0.5))

					// 4) Updates with the corrected precision
					preccorrected.Int(prec)
				}

				// round(q1/logprec)
				scale := new(big.Int).Set(diffScale.BigInt())
				bignum.DivRound(scale, prec, scale)

				// Checks that round(q1/logprec) >= 2^{logprec}
				requiresReservedPrime := scale.Cmp(new(big.Int).SetUint64(1)) < 0

				if requiresReservedPrime && eval.IterationsParameters.ReservedPrimeBitSize == 0 {
					return ctOut, fmt.Errorf("warning: early stopping at iteration k=%d: reason: round(q1/2^{logprec}) < 1 and no reserverd prime was provided", i+1)
				}

				// [M^{d} + e^{d-logprec}] - [M^{d}] -> [e^{d-logprec}]
				tmp, err := eval.Evaluator.SubNew(ctOut, ctIn)

				if err != nil {
					return nil, err
				}

				// prec * [e^{d-logprec}] -> [e^{d}]
				if err = eval.Evaluator.Mul(tmp, prec, tmp); err != nil {
					return nil, err
				}

				tmp.Scale = ctOut.Scale

				// [e^{d}] -> [e^{d}/q1] -> [e^{d}/q1 + e'^{d-logprec}]
				if tmp, errScale, err = eval.bootstrap(tmp); err != nil {
					return nil, err
				}

				tmp.Scale = tmp.Scale.Mul(*errScale)

				// [[e^{d}/q1 + e'^{d-logprec}] * q1/logprec -> [e^{d-logprec} + e'^{d-2logprec}*q1]
				if eval.IterationsParameters.ReservedPrimeBitSize == 0 {
					if err = eval.Evaluator.Mul(tmp, scale, tmp); err != nil {
						return nil, err
					}
				} else {

					// Else we compute the floating point ratio
					scale := new(big.Float).SetInt(diffScale.BigInt())
					scale.Quo(scale, new(big.Float).SetInt(prec))

					if new(big.Float).Mul(scale, new(big.Float).SetUint64(QiReserved)).Cmp(new(big.Float).SetUint64(1)) == -1 {
						return ctOut, fmt.Errorf("warning: early stopping at iteration k=%d: reason: maximum precision achieved", i+1)
					}

					// Do a scaled multiplication by the last prime
					if err = eval.Evaluator.Mul(tmp, scale, tmp); err != nil {
						return nil, err
					}

					// And rescale
					if err = eval.Evaluator.Rescale(tmp, tmp); err != nil {
						return nil, err
					}
				}

				// This is a given
				tmp.Scale = ctOut.Scale

				// [M^{d} + e^{d-logprec}] - [e^{d-logprec} + e'^{d-2logprec}*q1] -> [M^{d} + e'^{d-2logprec}*q1]
				if err = eval.Evaluator.Sub(ctOut, tmp, ctOut); err != nil {
					return nil, err
				}
			}
		}

		for ctOut.Level() > eval.ResidualParameters.MaxLevel() {
			eval.Evaluator.DropLevel(ctOut, 1)
		}
	}

	return
}

// EvaluateConjugateInvariant takes two ciphertext in the Conjugate Invariant ring, repacks them in a single ciphertext in the standard ring
// using the real and imaginary part, bootstrap both ciphertext, and then extract back the real and imaginary part before repacking them
// individually in two new ciphertexts in the Conjugate Invariant ring.
func (eval Evaluator) EvaluateConjugateInvariant(ctLeftN1Q0, ctRightN1Q0 *rlwe.Ciphertext) (ctLeftN1QL, ctRightN1QL *rlwe.Ciphertext, err error) {

	if ctLeftN1Q0 == nil {
		return nil, nil, fmt.Errorf("ctLeftN1Q0 cannot be nil")
	}

	// Switches ring from ring.ConjugateInvariant to ring.Standard
	ctLeftN2Q0 := eval.RealToComplexNew(ctLeftN1Q0)

	// Repacks ctRightN1Q0 into the imaginary part of ctLeftN1Q0
	// which is zero since it comes from the Conjugate Invariant ring)
	if ctRightN1Q0 != nil {
		ctRightN2Q0 := eval.RealToComplexNew(ctRightN1Q0)

		if err = eval.Evaluator.Mul(ctRightN2Q0, 1i, ctRightN2Q0); err != nil {
			return nil, nil, fmt.Errorf("cannot BootstrapMany: %w", err)
		}

		if err = eval.Evaluator.Add(ctLeftN2Q0, ctRightN2Q0, ctLeftN2Q0); err != nil {
			return nil, nil, fmt.Errorf("cannot BootstrapMany: %w", err)
		}
	}

	// Bootstraps in the ring.Standard
	var ctLeftAndRightN2QL *rlwe.Ciphertext
	if ctLeftAndRightN2QL, err = eval.Evaluate(ctLeftN2Q0); err != nil {
		return nil, nil, fmt.Errorf("cannot BootstrapMany: %w", err)
	}

	// The SlotsToCoeffs transformation scales the ciphertext by 0.5
	// This is done to compensate for the 2x factor introduced by ringStandardToConjugate(*).
	ctLeftAndRightN2QL.Scale = ctLeftAndRightN2QL.Scale.Mul(rlwe.NewScale(1 / 2.0))

	// Switches ring from ring.Standard to ring.ConjugateInvariant
	ctLeftN1QL = eval.ComplexToRealNew(ctLeftAndRightN2QL)

	// Extracts the imaginary part
	if ctRightN1Q0 != nil {
		if err = eval.Evaluator.Mul(ctLeftAndRightN2QL, -1i, ctLeftAndRightN2QL); err != nil {
			return nil, nil, fmt.Errorf("cannot BootstrapMany: %w", err)
		}
		ctRightN1QL = eval.ComplexToRealNew(ctLeftAndRightN2QL)
	}

	return
}

// checks if the current message ratio is greater or equal to the last prime times the target message ratio.
func checkMessageRatio(ct *rlwe.Ciphertext, msgRatio float64, r *ring.Ring) bool {
	level := ct.Level()
	currentMessageRatio := rlwe.NewScale(r.ModulusAtLevel[level])
	currentMessageRatio = currentMessageRatio.Div(ct.Scale)
	return currentMessageRatio.Cmp(rlwe.NewScale(r.SubRings[level].Modulus).Mul(rlwe.NewScale(msgRatio))) > -1
}

func (eval Evaluator) bootstrap(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext, errScale *rlwe.Scale, err error) {

	// Step 1: scale to q/|m|
	if ctOut, errScale, err = eval.ScaleDown(ctIn); err != nil {
		return
	}

	// Step 2 : Extend the basis from q to Q
	if ctOut, err = eval.ModUp(ctOut); err != nil {
		return
	}

	// Step 3 : CoeffsToSlots (Homomorphic encoding)
	// ctReal = Ecd(real)
	// ctImag = Ecd(imag)
	// If n < N/2 then ctReal = Ecd(real||imag)
	var ctReal, ctImag *rlwe.Ciphertext
	if ctReal, ctImag, err = eval.CoeffsToSlots(ctOut); err != nil {
		return
	}

	// Step 4 : EvalMod (Homomorphic modular reduction)
	if ctReal, err = eval.EvalMod(ctReal); err != nil {
		return
	}

	// Step 4 : EvalMod (Homomorphic modular reduction)
	if ctImag != nil {
		if ctImag, err = eval.EvalMod(ctImag); err != nil {
			return
		}
	}

	// Step 5 : SlotsToCoeffs (Homomorphic decoding)
	if ctOut, err = eval.SlotsToCoeffs(ctReal, ctImag); err != nil {
		return
	}

	return
}

// ScaleDown brings the ciphertext level to zero and scaling factor to Q[0]/MessageRatio
// It multiplies the ciphertexts by round(currentMessageRatio / targetMessageRatio) where:
//   - currentMessageRatio = Q/ctIn.Scale
//   - targetMessageRatio = q/|m|
//
// and updates the scale of ctIn accordingly
// It then rescales the ciphertext down to q if necessary and also returns the rescaling error from this process
func (eval Evaluator) ScaleDown(ctIn *rlwe.Ciphertext) (*rlwe.Ciphertext, *rlwe.Scale, error) {

	params := &eval.BootstrappingParameters

	r := params.RingQ()

	// Removes unecessary primes
	for ctIn.Level() != 0 && checkMessageRatio(ctIn, eval.Mod1Parameters.MessageRatio(), r) {
		ctIn.Resize(ctIn.Degree(), ctIn.Level()-1)
	}

	// Current Message Ratio
	currentMessageRatio := rlwe.NewScale(r.ModulusAtLevel[ctIn.Level()])
	currentMessageRatio = currentMessageRatio.Div(ctIn.Scale)

	// Desired Message Ratio
	targetMessageRatio := rlwe.NewScale(eval.Mod1Parameters.MessageRatio())

	// (Current Message Ratio) / (Desired Message Ratio)
	scaleUp := currentMessageRatio.Div(targetMessageRatio)

	if scaleUp.Cmp(rlwe.NewScale(0.5)) == -1 {
		return nil, nil, fmt.Errorf("initial Q/Scale = %f < 0.5*Q[0]/MessageRatio = %f", currentMessageRatio.Float64(), targetMessageRatio.Float64())
	}

	scaleUpBigint := scaleUp.BigInt()

	if err := eval.Evaluator.Mul(ctIn, scaleUpBigint, ctIn); err != nil {
		return nil, nil, err
	}

	ctIn.Scale = ctIn.Scale.Mul(rlwe.NewScale(scaleUpBigint))

	// errScale = CtIn.Scale/(Q[0]/MessageRatio)
	targetScale := new(big.Float).SetPrec(256).SetInt(r.ModulusAtLevel[0])
	targetScale.Quo(targetScale, new(big.Float).SetFloat64(eval.Mod1Parameters.MessageRatio()))

	if ctIn.Level() != 0 {
		if err := eval.RescaleTo(ctIn, rlwe.NewScale(targetScale), ctIn); err != nil {
			return nil, nil, err
		}
	}

	// Rescaling error (if any)
	errScale := ctIn.Scale.Div(rlwe.NewScale(targetScale))

	return ctIn, &errScale, nil
}

// ModUp raise the modulus from q to Q, scales the message  and applies the Trace if the ciphertext is sparsely packed.
func (eval Evaluator) ModUp(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext, err error) {

	// Switch to the sparse key
	if eval.EvkDenseToSparse != nil {
		if err := eval.ApplyEvaluationKey(ctIn, eval.EvkDenseToSparse, ctIn); err != nil {
			return nil, err
		}
	}

	params := eval.BootstrappingParameters

	ringQ := params.RingQ().AtLevel(ctIn.Level())
	ringP := params.RingP()

	for i := range ctIn.Value {
		ringQ.INTT(ctIn.Value[i], ctIn.Value[i])
	}

	// Extend the ciphertext from q to Q with zero values.
	ctIn.Resize(ctIn.Degree(), params.MaxLevel())

	levelQ := params.QCount() - 1
	levelP := params.PCount() - 1

	ringQ = ringQ.AtLevel(levelQ)

	Q := ringQ.ModuliChain()
	P := ringP.ModuliChain()
	q := Q[0]
	BRCQ := ringQ.BRedConstants()
	BRCP := ringP.BRedConstants()

	var coeff, tmp, pos, neg uint64

	N := ringQ.N()

	// ModUp q->Q for ctIn[0] centered around q
	for j := 0; j < N; j++ {

		coeff = ctIn.Value[0].Coeffs[0][j]
		pos, neg = 1, 0
		if coeff >= (q >> 1) {
			coeff = q - coeff
			pos, neg = 0, 1
		}

		for i := 1; i < levelQ+1; i++ {
			tmp = ring.BRedAdd(coeff, Q[i], BRCQ[i])
			ctIn.Value[0].Coeffs[i][j] = tmp*pos + (Q[i]-tmp)*neg
		}
	}

	if eval.EvkSparseToDense != nil {

		ks := eval.Evaluator.Evaluator

		// ModUp q->QP for ctIn[1] centered around q
		for j := 0; j < N; j++ {

			coeff = ctIn.Value[1].Coeffs[0][j]
			pos, neg = 1, 0
			if coeff > (q >> 1) {
				coeff = q - coeff
				pos, neg = 0, 1
			}

			for i := 0; i < levelQ+1; i++ {
				tmp = ring.BRedAdd(coeff, Q[i], BRCQ[i])
				ks.BuffDecompQP[0].Q.Coeffs[i][j] = tmp*pos + (Q[i]-tmp)*neg

			}

			for i := 0; i < levelP+1; i++ {
				tmp = ring.BRedAdd(coeff, P[i], BRCP[i])
				ks.BuffDecompQP[0].P.Coeffs[i][j] = tmp*pos + (P[i]-tmp)*neg
			}
		}

		for i := len(ks.BuffDecompQP) - 1; i >= 0; i-- {
			ringQ.NTT(ks.BuffDecompQP[0].Q, ks.BuffDecompQP[i].Q)
		}

		for i := len(ks.BuffDecompQP) - 1; i >= 0; i-- {
			ringP.NTT(ks.BuffDecompQP[0].P, ks.BuffDecompQP[i].P)
		}

		ringQ.NTT(ctIn.Value[0], ctIn.Value[0])

		// Scale the message from Q0/|m| to QL/|m|, where QL is the largest modulus used during the bootstrapping.
		if scale := (eval.Mod1Parameters.ScalingFactor().Float64() / eval.Mod1Parameters.MessageRatio()) / ctIn.Scale.Float64(); scale > 1 {

			scalar := uint64(math.Round(scale))

			for i := len(ks.BuffDecompQP) - 1; i >= 0; i-- {
				ringQ.MulScalar(ks.BuffDecompQP[0].Q, scalar, ks.BuffDecompQP[i].Q)
			}

			for i := len(ks.BuffDecompQP) - 1; i >= 0; i-- {
				ringP.MulScalar(ks.BuffDecompQP[0].P, scalar, ks.BuffDecompQP[i].P)
			}

			ringQ.MulScalar(ctIn.Value[0], scalar, ctIn.Value[0])

			ctIn.Scale = ctIn.Scale.Mul(rlwe.NewScale(scale))
		}

		ctTmp := &rlwe.Ciphertext{}
		ctTmp.Value = []ring.Poly{ks.BuffQP[1].Q, ctIn.Value[1]}
		ctTmp.MetaData = ctIn.MetaData

		// Switch back to the dense key
		ks.GadgetProductHoisted(levelQ, ks.BuffDecompQP, &eval.EvkSparseToDense.GadgetCiphertext, ctTmp)
		ringQ.Add(ctIn.Value[0], ctTmp.Value[0], ctIn.Value[0])

	} else {

		for j := 0; j < N; j++ {

			coeff = ctIn.Value[1].Coeffs[0][j]
			pos, neg = 1, 0
			if coeff >= (q >> 1) {
				coeff = q - coeff
				pos, neg = 0, 1
			}

			for i := 1; i < levelQ+1; i++ {
				tmp = ring.BRedAdd(coeff, Q[i], BRCQ[i])
				ctIn.Value[1].Coeffs[i][j] = tmp*pos + (Q[i]-tmp)*neg
			}
		}

		ringQ.NTT(ctIn.Value[0], ctIn.Value[0])
		ringQ.NTT(ctIn.Value[1], ctIn.Value[1])

		// Scale the message from Q0/|m| to QL/|m|, where QL is the largest modulus used during the bootstrapping.
		if scale := (eval.Mod1Parameters.ScalingFactor().Float64() / eval.Mod1Parameters.MessageRatio()) / ctIn.Scale.Float64(); scale > 1 {

			scalar := uint64(math.Round(scale))

			ringQ.MulScalar(ctIn.Value[0], scalar, ctIn.Value[0])
			ringQ.MulScalar(ctIn.Value[1], scalar, ctIn.Value[1])

			ctIn.Scale = ctIn.Scale.Mul(rlwe.NewScale(scale))
		}
	}

	//SubSum X -> (N/dslots) * Y^dslots
	return ctIn, eval.Trace(ctIn, eval.CoeffsToSlotsParameters.LogSlots, ctIn)
}

// CoeffsToSlots applies the homomorphic decoding
func (eval Evaluator) CoeffsToSlots(ctIn *rlwe.Ciphertext) (ctReal, ctImag *rlwe.Ciphertext, err error) {
	return eval.DFTEvaluator.CoeffsToSlotsNew(ctIn, eval.C2SDFTMatrix)
}

// EvalMod applies the homomorphic modular reduction by q.
func (eval Evaluator) EvalMod(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext, err error) {
	if ctOut, err = eval.Mod1Evaluator.EvaluateNew(ctIn); err != nil {
		return nil, err
	}

	ctOut.Scale = eval.BootstrappingParameters.DefaultScale()
	return
}

// EvalModAndScale applies the homomorphic modular reduction by q and scales the output value (without
// consuming an additional level).
func (eval Evaluator) EvalModAndScale(ctIn *rlwe.Ciphertext, scaling complex128) (ctOut *rlwe.Ciphertext, err error) {
	if ctOut, err = eval.Mod1Evaluator.EvaluateAndScaleNew(ctIn, scaling); err != nil {
		return nil, err
	}

	ctOut.Scale = eval.BootstrappingParameters.DefaultScale()
	return
}

func (eval Evaluator) SlotsToCoeffs(ctReal, ctImag *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext, err error) {
	return eval.DFTEvaluator.SlotsToCoeffsNew(ctReal, ctImag, eval.S2CDFTMatrix)
}

func (eval Evaluator) switchRingDegreeN1ToN2New(ctN1 *rlwe.Ciphertext) (ctN2 *rlwe.Ciphertext) {
	ctN2 = ckks.NewCiphertext(eval.BootstrappingParameters, 1, ctN1.Level())

	// Sanity check, this error should never happen unless this algorithm has been improperly
	// modified to pass invalid inputs.
	if err := eval.Evaluator.ApplyEvaluationKey(ctN1, eval.EvkN1ToN2, ctN2); err != nil {
		panic(err)
	}
	return
}

func (eval Evaluator) switchRingDegreeN2ToN1New(ctN2 *rlwe.Ciphertext) (ctN1 *rlwe.Ciphertext) {
	ctN1 = ckks.NewCiphertext(eval.ResidualParameters, 1, ctN2.Level())

	// Sanity check, this error should never happen unless this algorithm has been improperly
	// modified to pass invalid inputs.
	if err := eval.Evaluator.ApplyEvaluationKey(ctN2, eval.EvkN2ToN1, ctN1); err != nil {
		panic(err)
	}
	return
}

func (eval Evaluator) ComplexToRealNew(ctCmplx *rlwe.Ciphertext) (ctReal *rlwe.Ciphertext) {
	ctReal = ckks.NewCiphertext(eval.ResidualParameters, 1, ctCmplx.Level())

	// Sanity check, this error should never happen unless this algorithm has been improperly
	// modified to pass invalid inputs.
	if err := eval.DomainSwitcher.ComplexToReal(eval.Evaluator, ctCmplx, ctReal); err != nil {
		panic(err)
	}
	return
}

func (eval Evaluator) RealToComplexNew(ctReal *rlwe.Ciphertext) (ctCmplx *rlwe.Ciphertext) {
	ctCmplx = ckks.NewCiphertext(eval.BootstrappingParameters, 1, ctReal.Level())

	// Sanity check, this error should never happen unless this algorithm has been improperly
	// modified to pass invalid inputs.
	if err := eval.DomainSwitcher.RealToComplex(eval.Evaluator, ctReal, ctCmplx); err != nil {
		panic(err)
	}
	return
}

// packingContext contains the parameters used when packing (with Pack())
type packingContext struct {
	Params           *ckks.Parameters // Parameters of the ring we are packing to or unpacking from
	LogMaxDimensions ring.Dimensions  // maximum dimension of a packed ciphertext (logMaxDimensions <= params.LogMaxDimensions())
	LogSlots         int              // number of slots in a ct before packing (resp. after unpacking)
	NbPackedCTs      int              // number of cts to be packed (resp. to be unpacked into)
}

// PackAndSwitchN1ToN2 packs the ciphertexts into N1 and switch to N2 if N1 < N2
// then it packs the ciphertexts into N2.
func (eval Evaluator) PackAndSwitchN1ToN2(cts []rlwe.Ciphertext) ([]rlwe.Ciphertext, *packingContext, *packingContext, error) {

	var err error
	var packN1, packN2 *packingContext

	// If N1 < N2, we pack ciphertexts into N1 and then switch to N2
	if eval.ResidualParameters.N() != eval.BootstrappingParameters.N() {

		packN1 = &packingContext{&eval.ResidualParameters, eval.ResidualParameters.LogMaxDimensions(), cts[0].LogSlots(), len(cts)}

		// If the bootstrapping max slots are smaller than the max slots of N1, we only pack up to the former
		if eval.Parameters.LogMaxSlots() < eval.ResidualParameters.LogMaxSlots() {
			packN1.LogMaxDimensions = eval.Parameters.LogMaxDimensions()
		}
		if cts, err = eval.pack(cts, *packN1, eval.xPow2N1); err != nil {
			return nil, nil, nil, fmt.Errorf("cannot PackAndSwitchN1ToN2: PackN1: %w", err)
		}

		for i := range cts {
			cts[i] = *eval.switchRingDegreeN1ToN2New(&cts[i])
		}
	}

	// Packing ciphertexts into N2 (up to eval.Parameters.LogMaxDimensions())
	packN2 = &packingContext{&eval.BootstrappingParameters, eval.Parameters.LogMaxDimensions(), cts[0].LogSlots(), len(cts)}

	if cts, err = eval.pack(cts, *packN2, eval.xPow2N2); err != nil {
		return nil, nil, nil, fmt.Errorf("cannot PackAndSwitchN1ToN2: PackN1: %w", err)
	}

	return cts, packN1, packN2, nil
}

// UnpackAndSwitchN2ToN1 unpacks the ciphertexts into N2 and, if N1 < N2, it switches the ciphertexts
// to N1 and unpacks further into N1
func (eval Evaluator) UnpackAndSwitchN2ToN1(cts []rlwe.Ciphertext, ctxtN1, ctxtN2 *packingContext) ([]rlwe.Ciphertext, error) {

	var ctsOut []rlwe.Ciphertext

	logSlots := ctxtN2.LogSlots

	// Unpack ciphertexts in N2
	for i := range cts {
		ctsUnpack, err := eval.unpack(&cts[i], *ctxtN2, eval.xPow2InvN2)

		if err != nil {
			return nil, fmt.Errorf("cannot UnpackAndSwitchN2Tn1: UnpackN2: %w", err)
		}

		ctsOut = append(ctsOut, ctsUnpack...)
		ctxtN2.NbPackedCTs -= len(ctsUnpack)
	}

	// If N1 != N2 (i.e. ctxtN1 != nil): 1) switch cts to N1 2) unpack the cts in N1
	if ctxtN1 != nil {
		var ctsN1 []rlwe.Ciphertext
		logSlots = ctxtN1.LogSlots

		for i := range ctsOut {
			ctsOut[i] = *eval.switchRingDegreeN2ToN1New(&ctsOut[i])
		}

		for i := range ctsOut {
			ctsUnpack, err := eval.unpack(&ctsOut[i], *ctxtN1, eval.xPow2InvN1)
			if err != nil {
				return nil, fmt.Errorf("cannot UnpackAndSwitchN2Tn1: UnpackN1: %w", err)
			}

			ctsN1 = append(ctsN1, ctsUnpack...)
			ctxtN1.NbPackedCTs -= len(ctsUnpack)
		}

		ctsOut = ctsN1
	}

	// Set back the dimension of cts to its original value
	for i := range ctsOut {
		ctsOut[i].LogDimensions.Cols = logSlots
	}

	return ctsOut, nil
}

// unpack unpacks one sparse ciphertext of (log) dimension ctxt.logMaxDimensions
// into ctxt.NbPackedCTs ciphertexts of (log) dimension {0, ctxt.LogSlots}
func (eval Evaluator) unpack(ct *rlwe.Ciphertext, ctxt packingContext, xPow2Inv []ring.Poly) ([]rlwe.Ciphertext, error) {
	logPackCTs := ctxt.LogMaxDimensions.Cols - ctxt.LogSlots // log of number of CTs that can be packed in one ct

	cts := []rlwe.Ciphertext{*ct}
	if logPackCTs == 0 {
		return cts, nil
	}

	n := utils.Min(ctxt.NbPackedCTs, 1<<logPackCTs) // #cts to unpack from ct
	cts = append(cts, make([]rlwe.Ciphertext, n-1)...)

	for i := 1; i < len(cts); i++ {
		cts[i] = *ct.CopyNew()
	}

	r := ctxt.Params.RingQ().AtLevel(cts[0].Level())

	logGap := (ctxt.Params.LogMaxSlots() - ctxt.LogSlots) - 1 // log gap of CTs with params.N (minus one)

	/* #nosec G115 -- n-1 cannot be negative */
	for i := 0; i < utils.Min(bits.Len64(uint64(n-1)), logPackCTs); i++ {

		step := 1 << (i + 1)

		for j := 0; j < n; j += step {

			for k := step >> 1; k < step; k++ {

				if (j + k) >= n {
					break
				}

				r.MulCoeffsMontgomery(cts[j+k].Value[0], xPow2Inv[logGap-i], cts[j+k].Value[0])
				r.MulCoeffsMontgomery(cts[j+k].Value[1], xPow2Inv[logGap-i], cts[j+k].Value[1])
			}
		}
	}

	return cts, nil
}

// pack packs ctxt.NbPackedCTs sparse ciphertexts of (log) dimension {0, ctxt.LogSlots}
// into one ciphertext of (log) dimension ctxt.logMaxDimensions
func (eval Evaluator) pack(cts []rlwe.Ciphertext, ctxt packingContext, xPow2 []ring.Poly) ([]rlwe.Ciphertext, error) {

	var logSlots = ctxt.LogSlots
	var logMaxSlots = ctxt.LogMaxDimensions.Cols
	ringDegree := ctxt.Params.N()

	if logSlots > logMaxSlots {
		return nil, fmt.Errorf("cannot Pack: cts[0].LogSlots()=%d > logMaxSlots=%d", logSlots, logMaxSlots)
	}

	for i, ct := range cts {
		if s := ct.LogSlots(); s != logSlots {
			return nil, fmt.Errorf("cannot Pack: cts[%d].PlaintextLogSlots()=%d != cts[0].PlaintextLogSlots=%d", i, s, logSlots)
		}

		if N := ct.Value[0].N(); N != ringDegree {
			return nil, fmt.Errorf("cannot Pack: cts[%d].Value[0].N()=%d != params.N()=%d", i, N, ringDegree)
		}
	}

	logPackCTs := logMaxSlots - logSlots // log of number of CTs that can be packed in one ct
	logGap := (ctxt.Params.LogMaxSlots() - logSlots - 1)

	if logPackCTs == 0 {
		return cts, nil
	}

	for i := 0; i < logPackCTs; i++ {

		for j := 0; j < len(cts)>>1; j++ {

			eve := cts[j*2+0]
			odd := cts[j*2+1]

			level := utils.Min(eve.Level(), odd.Level())

			r := ctxt.Params.RingQ().AtLevel(level)

			r.MulCoeffsMontgomeryThenAdd(odd.Value[0], xPow2[logGap-i], eve.Value[0])
			r.MulCoeffsMontgomeryThenAdd(odd.Value[1], xPow2[logGap-i], eve.Value[1])

			cts[j] = eve
		}

		if len(cts)&1 == 1 {
			cts[len(cts)>>1] = cts[len(cts)-1]
			cts = cts[:len(cts)>>1+1]
		} else {
			cts = cts[:len(cts)>>1]
		}
	}

	for i := range cts {
		cts[i].LogDimensions = ctxt.LogMaxDimensions
	}

	return cts, nil
}
