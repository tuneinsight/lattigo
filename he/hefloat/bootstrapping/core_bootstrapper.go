package bootstrapping

import (
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/core/rlwe"
	"github.com/tuneinsight/lattigo/v4/he/hefloat"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// CoreBootstrapper is a struct to store a memory buffer with the plaintext matrices,
// the polynomial approximation, and the keys for the bootstrapping.
type CoreBootstrapper struct {
	*hefloat.Evaluator
	*hefloat.DFTEvaluator
	*hefloat.Mod1Evaluator
	*bootstrapperBase
	SkDebug *rlwe.SecretKey
}

type bootstrapperBase struct {
	Parameters
	*EvaluationKeys
	params hefloat.Parameters

	dslots    int // Number of plaintext slots after the re-encoding: min(2*slots, N/2)
	logdslots int // log2(dslots)

	mod1Parameters hefloat.Mod1Parameters
	stcMatrices    hefloat.DFTMatrix
	ctsMatrices    hefloat.DFTMatrix

	q0OverMessageRatio float64
}

// NewCoreBootstrapper creates a new CoreBootstrapper.
func NewCoreBootstrapper(btpParams Parameters, evk *EvaluationKeys) (btp *CoreBootstrapper, err error) {

	if btpParams.Mod1ParametersLiteral.Mod1Type == hefloat.SinContinuous && btpParams.Mod1ParametersLiteral.DoubleAngle != 0 {
		return nil, fmt.Errorf("cannot use double angle formula for Mod1Type = Sin -> must use Mod1Type = Cos")
	}

	if btpParams.Mod1ParametersLiteral.Mod1Type == hefloat.CosDiscrete && btpParams.Mod1ParametersLiteral.Mod1Degree < 2*(btpParams.Mod1ParametersLiteral.K-1) {
		return nil, fmt.Errorf("Mod1Type 'hefloat.CosDiscrete' uses a minimum degree of 2*(K-1) but EvalMod degree is smaller")
	}

	if btpParams.CoeffsToSlotsParameters.LevelStart-btpParams.CoeffsToSlotsParameters.Depth(true) != btpParams.Mod1ParametersLiteral.LevelStart {
		return nil, fmt.Errorf("starting level and depth of CoeffsToSlotsParameters inconsistent starting level of SineEvalParameters")
	}

	if btpParams.Mod1ParametersLiteral.LevelStart-btpParams.Mod1ParametersLiteral.Depth() != btpParams.SlotsToCoeffsParameters.LevelStart {
		return nil, fmt.Errorf("starting level and depth of SineEvalParameters inconsistent starting level of CoeffsToSlotsParameters")
	}

	params := btpParams.BootstrappingParameters

	btp = new(CoreBootstrapper)
	if btp.bootstrapperBase, err = newBootstrapperBase(params, btpParams, evk); err != nil {
		return
	}

	if err = btp.bootstrapperBase.CheckKeys(evk); err != nil {
		return nil, fmt.Errorf("invalid bootstrapping key: %w", err)
	}

	btp.EvaluationKeys = evk

	btp.Evaluator = hefloat.NewEvaluator(params, evk)

	btp.DFTEvaluator = hefloat.NewDFTEvaluator(params, btp.Evaluator)

	btp.Mod1Evaluator = hefloat.NewMod1Evaluator(btp.Evaluator, hefloat.NewPolynomialEvaluator(params, btp.Evaluator), btp.bootstrapperBase.mod1Parameters)

	return
}

// ShallowCopy creates a shallow copy of this CoreBootstrapper in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// CoreBootstrapper can be used concurrently.
func (btp CoreBootstrapper) ShallowCopy() *CoreBootstrapper {
	Evaluator := btp.Evaluator.ShallowCopy()
	params := btp.BootstrappingParameters
	return &CoreBootstrapper{
		Evaluator:        Evaluator,
		bootstrapperBase: btp.bootstrapperBase,
		DFTEvaluator:     hefloat.NewDFTEvaluator(params, Evaluator),
		Mod1Evaluator:    hefloat.NewMod1Evaluator(Evaluator, hefloat.NewPolynomialEvaluator(params, Evaluator), btp.bootstrapperBase.mod1Parameters),
	}
}

// CheckKeys checks if all the necessary keys are present in the instantiated CoreBootstrapper
func (bb *bootstrapperBase) CheckKeys(evk *EvaluationKeys) (err error) {

	if _, err = evk.GetRelinearizationKey(); err != nil {
		return
	}

	for _, galEl := range bb.GaloisElements(bb.params) {
		if _, err = evk.GetGaloisKey(galEl); err != nil {
			return
		}
	}

	if evk.EvkDenseToSparse == nil && bb.EphemeralSecretWeight != 0 {
		return fmt.Errorf("rlwe.EvaluationKey key dense to sparse is nil")
	}

	if evk.EvkSparseToDense == nil && bb.EphemeralSecretWeight != 0 {
		return fmt.Errorf("rlwe.EvaluationKey key sparse to dense is nil")
	}

	return
}

func newBootstrapperBase(params hefloat.Parameters, btpParams Parameters, evk *EvaluationKeys) (bb *bootstrapperBase, err error) {
	bb = new(bootstrapperBase)
	bb.params = params
	bb.Parameters = btpParams

	bb.logdslots = btpParams.LogMaxDimensions().Cols
	bb.dslots = 1 << bb.logdslots
	if maxLogSlots := params.LogMaxDimensions().Cols; bb.dslots < maxLogSlots {
		bb.dslots <<= 1
		bb.logdslots++
	}

	if bb.mod1Parameters, err = hefloat.NewMod1ParametersFromLiteral(params, btpParams.Mod1ParametersLiteral); err != nil {
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

	encoder := hefloat.NewEncoder(bb.params)

	// CoeffsToSlots vectors
	// Change of variable for the evaluation of the Chebyshev polynomial + cancelling factor for the DFT and SubSum + eventual scaling factor for the double angle formula

	if bb.CoeffsToSlotsParameters.Scaling == nil {
		bb.CoeffsToSlotsParameters.Scaling = new(big.Float).SetFloat64(qDiv / (K * scFac * qDiff))
	} else {
		bb.CoeffsToSlotsParameters.Scaling.Mul(bb.CoeffsToSlotsParameters.Scaling, new(big.Float).SetFloat64(qDiv/(K*scFac*qDiff)))
	}

	if bb.ctsMatrices, err = hefloat.NewDFTMatrixFromLiteral(params, bb.CoeffsToSlotsParameters, encoder); err != nil {
		return
	}

	// SlotsToCoeffs vectors
	// Rescaling factor to set the final ciphertext to the desired scale

	if bb.SlotsToCoeffsParameters.Scaling == nil {
		bb.SlotsToCoeffsParameters.Scaling = new(big.Float).SetFloat64(bb.params.DefaultScale().Float64() / (bb.mod1Parameters.ScalingFactor().Float64() / bb.mod1Parameters.MessageRatio()))
	} else {
		bb.SlotsToCoeffsParameters.Scaling.Mul(bb.SlotsToCoeffsParameters.Scaling, new(big.Float).SetFloat64(bb.params.DefaultScale().Float64()/(bb.mod1Parameters.ScalingFactor().Float64()/bb.mod1Parameters.MessageRatio())))
	}

	if bb.stcMatrices, err = hefloat.NewDFTMatrixFromLiteral(params, bb.SlotsToCoeffsParameters, encoder); err != nil {
		return
	}

	encoder = nil // For the GC

	return
}

func (btp CoreBootstrapper) MinimumInputLevel() int {
	return btp.params.LevelsConsumedPerRescaling()
}

func (btp CoreBootstrapper) OutputLevel() int {
	return btp.params.MaxLevel() - btp.Depth()
}

// Bootstrap re-encrypts a ciphertext to a ciphertext at MaxLevel - k where k is the depth of the bootstrapping circuit.
// If the input ciphertext level is zero, the input scale must be an exact power of two smaller than Q[0]/MessageRatio
// (it can't be equal since Q[0] is not a power of two).
// The message ratio is an optional field in the bootstrapping parameters, by default it set to 2^{LogMessageRatio = 8}.
// See the bootstrapping parameters for more information about the message ratio or other parameters related to the bootstrapping.
// If the input ciphertext is at level one or more, the input scale does not need to be an exact power of two as one level
// can be used to do a scale matching.
//
// The circuit has two variants, each consisting in 5 steps.
// Variant I:
// 1) ScaleDown: scales the ciphertext to q/|m| and bringing it down to q
// 2) ModUp: brings the modulus from q to Q
// 3) CoeffsToSlots: homomorphic encoding
// 4) EvalMod: homomorphic modular reduction
// 5) SlotsToCoeffs: homomorphic decoding
//
// Variant II:
// 1) SlotsToCoeffs: homomorphic decoding
// 2) ScaleDown: scales the ciphertext to q/|m| and bringing it down to q
// 3) ModUp: brings the modulus from q to Q
// 4) CoeffsToSlots: homomorphic encoding
// 5) EvalMod: homomorphic modular reduction
func (btp CoreBootstrapper) Bootstrap(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext, err error) {

	if btp.IterationsParameters == nil {
		ctOut, _, err = btp.bootstrap(ctIn)
		return

	} else {

		var errScale *rlwe.Scale
		// [M^{d}/q1 + e^{d-logprec}]
		if ctOut, errScale, err = btp.bootstrap(ctIn.CopyNew()); err != nil {
			return nil, err
		}

		// Stores by how much a ciphertext must be scaled to get back
		// to the input scale
		// Error correcting factor of the approximate division by q1
		// diffScale = ctIn.Scale / (ctOut.Scale * errScale)
		diffScale := ctIn.Scale.Div(ctOut.Scale)
		diffScale = diffScale.Div(*errScale)

		// [M^{d} + e^{d-logprec}]
		if err = btp.Evaluator.Mul(ctOut, diffScale.BigInt(), ctOut); err != nil {
			return nil, err
		}
		ctOut.Scale = ctIn.Scale

		var totLogPrec float64

		for i := 0; i < len(btp.IterationsParameters.BootstrappingPrecision); i++ {

			logPrec := btp.IterationsParameters.BootstrappingPrecision[i]

			totLogPrec += logPrec

			// prec = round(2^{logprec})
			log2 := bignum.Log(new(big.Float).SetPrec(256).SetUint64(2))
			log2TimesLogPrec := log2.Mul(log2, new(big.Float).SetFloat64(totLogPrec))
			prec := new(big.Int)
			log2TimesLogPrec.Add(bignum.Exp(log2TimesLogPrec), new(big.Float).SetFloat64(0.5)).Int(prec)

			// round(q1/logprec)
			scale := new(big.Int).Set(diffScale.BigInt())
			bignum.DivRound(scale, prec, scale)

			// Checks that round(q1/logprec) >= 2^{logprec}
			requiresReservedPrime := scale.Cmp(new(big.Int).SetUint64(1)) < 0

			if requiresReservedPrime && btp.IterationsParameters.ReservedPrimeBitSize == 0 {
				return ctOut, fmt.Errorf("warning: early stopping at iteration k=%d: reason: round(q1/2^{logprec}) < 1 and no reserverd prime was provided", i+1)
			}

			// [M^{d} + e^{d-logprec}] - [M^{d}] -> [e^{d-logprec}]
			tmp, err := btp.Evaluator.SubNew(ctOut, ctIn)

			if err != nil {
				return nil, err
			}

			// prec * [e^{d-logprec}] -> [e^{d}]
			if err = btp.Evaluator.Mul(tmp, prec, tmp); err != nil {
				return nil, err
			}

			tmp.Scale = ctOut.Scale

			// [e^{d}] -> [e^{d}/q1] -> [e^{d}/q1 + e'^{d-logprec}]
			if tmp, errScale, err = btp.bootstrap(tmp); err != nil {
				return nil, err
			}

			tmp.Scale = tmp.Scale.Mul(*errScale)

			// [[e^{d}/q1 + e'^{d-logprec}] * q1/logprec -> [e^{d-logprec} + e'^{d-2logprec}*q1]
			// If scale > 2^{logprec}, then we ensure a precision of at least 2^{logprec} even with a rounding of the scale
			if !requiresReservedPrime {
				if err = btp.Evaluator.Mul(tmp, scale, tmp); err != nil {
					return nil, err
				}
			} else {

				// Else we compute the floating point ratio
				ss := new(big.Float).SetInt(diffScale.BigInt())
				ss.Quo(ss, new(big.Float).SetInt(prec))

				// Do a scaled multiplication by the last prime
				if err = btp.Evaluator.Mul(tmp, ss, tmp); err != nil {
					return nil, err
				}

				// And rescale
				if err = btp.Evaluator.Rescale(tmp, tmp); err != nil {
					return nil, err
				}
			}

			// This is a given
			tmp.Scale = ctOut.Scale

			// [M^{d} + e^{d-logprec}] - [e^{d-logprec} + e'^{d-2logprec}*q1] -> [M^{d} + e'^{d-2logprec}*q1]
			if err = btp.Evaluator.Sub(ctOut, tmp, ctOut); err != nil {
				return nil, err
			}
		}
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

// ScaleDown brings the ciphertext level to zero and scaling factor to Q[0]/MessageRatio
// It multiplies the ciphertexts by round(currentMessageRatio / targetMessageRatio) where:
// - currentMessageRatio = Q/ctIn.Scale
// - targetMessageRatio = q/|m|
// and updates the scale of ctIn accordingly
// It then rescales the ciphertext down to q if necessary and also returns the rescaling error from this process
func (btp CoreBootstrapper) ScaleDown(ctIn *rlwe.Ciphertext) (*rlwe.Ciphertext, *rlwe.Scale, error) {

	params := &btp.params

	r := params.RingQ()

	// Removes unecessary primes
	for ctIn.Level() != 0 && checkMessageRatio(ctIn, btp.Mod1Parameters.MessageRatio(), r) {
		ctIn.Resize(ctIn.Degree(), ctIn.Level()-1)
	}

	// Current Message Ratio
	currentMessageRatio := rlwe.NewScale(r.ModulusAtLevel[ctIn.Level()])
	currentMessageRatio = currentMessageRatio.Div(ctIn.Scale)

	// Desired Message Ratio
	targetMessageRatio := rlwe.NewScale(btp.Mod1Parameters.MessageRatio())

	// (Current Message Ratio) / (Desired Message Ratio)
	scaleUp := currentMessageRatio.Div(targetMessageRatio)

	if scaleUp.Cmp(rlwe.NewScale(0.5)) == -1 {
		return nil, nil, fmt.Errorf("initial Q/Scale < 0.5*Q[0]/MessageRatio")
	}

	scaleUpBigint := scaleUp.BigInt()

	if err := btp.Evaluator.Mul(ctIn, scaleUpBigint, ctIn); err != nil {
		return nil, nil, err
	}

	ctIn.Scale = ctIn.Scale.Mul(rlwe.NewScale(scaleUpBigint))

	// errScale = CtIn.Scale/(Q[0]/MessageRatio)
	targetScale := new(big.Float).SetPrec(256).SetInt(r.ModulusAtLevel[0])
	targetScale.Quo(targetScale, new(big.Float).SetFloat64(btp.Mod1Parameters.MessageRatio()))

	if ctIn.Level() != 0 {
		if err := btp.RescaleTo(ctIn, rlwe.NewScale(targetScale), ctIn); err != nil {
			return nil, nil, err
		}
	}

	// Rescaling error (if any)
	errScale := ctIn.Scale.Div(rlwe.NewScale(targetScale))

	return ctIn, &errScale, nil
}

// ModUp raise the modulus from q to Q, scales the message  and applies the Trace if the ciphertext is sparsely packed.
func (btp *CoreBootstrapper) ModUp(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext, err error) {

	// Switch to the sparse key
	if btp.EvkDenseToSparse != nil {
		if err := btp.ApplyEvaluationKey(ctIn, btp.EvkDenseToSparse, ctIn); err != nil {
			return nil, err
		}
	}

	ringQ := btp.params.RingQ().AtLevel(ctIn.Level())
	ringP := btp.params.RingP()

	for i := range ctIn.Value {
		ringQ.INTT(ctIn.Value[i], ctIn.Value[i])
	}

	// Extend the ciphertext from q to Q with zero values.
	ctIn.Resize(ctIn.Degree(), btp.params.MaxLevel())

	levelQ := btp.params.QCount() - 1
	levelP := btp.params.PCount() - 1

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

	if btp.EvkSparseToDense != nil {

		ks := btp.Evaluator.Evaluator

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

		ctTmp := &rlwe.Ciphertext{}
		ctTmp.Value = []ring.Poly{ks.BuffQP[1].Q, ctIn.Value[1]}
		ctTmp.MetaData = ctIn.MetaData

		// Switch back to the dense key
		ks.GadgetProductHoisted(levelQ, ks.BuffDecompQP, &btp.EvkSparseToDense.GadgetCiphertext, ctTmp)
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
	}

	// Scale the message from Q0/|m| to QL/|m|, where QL is the largest modulus used during the bootstrapping.
	if scale := (btp.Mod1Parameters.ScalingFactor().Float64() / btp.Mod1Parameters.MessageRatio()) / ctIn.Scale.Float64(); scale > 1 {
		if err = btp.ScaleUp(ctIn, rlwe.NewScale(scale), ctIn); err != nil {
			return nil, err
		}
	}

	//SubSum X -> (N/dslots) * Y^dslots
	return ctIn, btp.Trace(ctIn, ctIn.LogDimensions.Cols, ctIn)
}

// CoeffsToSlots applies the homomorphic decoding
func (btp *CoreBootstrapper) CoeffsToSlots(ctIn *rlwe.Ciphertext) (ctReal, ctImag *rlwe.Ciphertext, err error) {
	return btp.CoeffsToSlotsNew(ctIn, btp.ctsMatrices)
}

// EvalMod applies the homomorphic modular reduction by q.
func (btp *CoreBootstrapper) EvalMod(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext, err error) {

	if ctOut, err = btp.Mod1Evaluator.EvaluateNew(ctIn); err != nil {
		return nil, err
	}
	ctOut.Scale = btp.params.DefaultScale()
	return
}

func (btp *CoreBootstrapper) SlotsToCoeffs(ctReal, ctImag *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext, err error) {
	return btp.SlotsToCoeffsNew(ctReal, ctImag, btp.stcMatrices)
}

func (btp *CoreBootstrapper) bootstrap(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext, errScale *rlwe.Scale, err error) {

	// Step 1: scale to q/|m|
	if ctOut, errScale, err = btp.ScaleDown(ctIn); err != nil {
		return
	}

	// Step 2 : Extend the basis from q to Q
	if ctOut, err = btp.ModUp(ctOut); err != nil {
		return
	}

	// Step 3 : CoeffsToSlots (Homomorphic encoding)
	// ctReal = Ecd(real)
	// ctImag = Ecd(imag)
	// If n < N/2 then ctReal = Ecd(real||imag)
	var ctReal, ctImag *rlwe.Ciphertext
	if ctReal, ctImag, err = btp.CoeffsToSlots(ctOut); err != nil {
		return
	}

	// Step 4 : EvalMod (Homomorphic modular reduction)
	if ctReal, err = btp.EvalMod(ctReal); err != nil {
		return
	}

	// Step 4 : EvalMod (Homomorphic modular reduction)
	if ctImag != nil {
		if ctImag, err = btp.EvalMod(ctImag); err != nil {
			return
		}
	}

	// Step 5 : SlotsToCoeffs (Homomorphic decoding)
	if ctOut, err = btp.SlotsToCoeffs(ctReal, ctImag); err != nil {
		return
	}

	return
}
