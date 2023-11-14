// Package bootstrapping implements bootstrapping for fixed-point encrypted
// approximate homomorphic encryption over the complex/real numbers.
package bootstrapping

import (
	"fmt"
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/schemes/ckks"
	"github.com/tuneinsight/lattigo/v5/utils"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
)

// Evaluate re-encrypts a ciphertext to a ciphertext at MaxLevel - k where k is the depth of the bootstrapping circuit.
// If the input ciphertext level is zero, the input scale must be an exact power of two smaller than Q[0]/MessageRatio
// (it can't be equal since Q[0] is not a power of two).
// The message ratio is an optional field in the bootstrapping parameters, by default it set to 2^{LogMessageRatio = 8}.
// See the bootstrapping parameters for more information about the message ratio or other parameters related to the bootstrapping.
// If the input ciphertext is at level one or more, the input scale does not need to be an exact power of two as one level
// can be used to do a scale matching.
//
// The circuit consists in 5 steps.
// 1) ScaleDown: scales the ciphertext to q/|m| and bringing it down to q
// 2) ModUp: brings the modulus from q to Q
// 3) CoeffsToSlots: homomorphic encoding
// 4) EvalMod: homomorphic modular reduction
// 5) SlotsToCoeffs: homomorphic decoding
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
// - currentMessageRatio = Q/ctIn.Scale
// - targetMessageRatio = q/|m|
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
	}

	// Scale the message from Q0/|m| to QL/|m|, where QL is the largest modulus used during the bootstrapping.
	if scale := (eval.Mod1Parameters.ScalingFactor().Float64() / eval.Mod1Parameters.MessageRatio()) / ctIn.Scale.Float64(); scale > 1 {
		if err = eval.ScaleUp(ctIn, rlwe.NewScale(scale), ctIn); err != nil {
			return nil, err
		}
	}

	//SubSum X -> (N/dslots) * Y^dslots
	return ctIn, eval.Trace(ctIn, ctIn.LogDimensions.Cols, ctIn)
}

// CoeffsToSlots applies the homomorphic decoding
func (eval Evaluator) CoeffsToSlots(ctIn *rlwe.Ciphertext) (ctReal, ctImag *rlwe.Ciphertext, err error) {
	return eval.CoeffsToSlotsNew(ctIn, eval.C2SDFTMatrix)
}

// EvalMod applies the homomorphic modular reduction by q.
func (eval Evaluator) EvalMod(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext, err error) {

	if ctOut, err = eval.Mod1Evaluator.EvaluateNew(ctIn); err != nil {
		return nil, err
	}
	ctOut.Scale = eval.BootstrappingParameters.DefaultScale()
	return
}

func (eval Evaluator) SlotsToCoeffs(ctReal, ctImag *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext, err error) {
	return eval.SlotsToCoeffsNew(ctReal, ctImag, eval.S2CDFTMatrix)
}

func (eval Evaluator) SwitchRingDegreeN1ToN2New(ctN1 *rlwe.Ciphertext) (ctN2 *rlwe.Ciphertext) {
	ctN2 = hefloat.NewCiphertext(eval.BootstrappingParameters, 1, ctN1.Level())

	// Sanity check, this error should never happen unless this algorithm has been improperly
	// modified to pass invalid inputs.
	if err := eval.Evaluator.ApplyEvaluationKey(ctN1, eval.EvkN1ToN2, ctN2); err != nil {
		panic(err)
	}
	return
}

func (eval Evaluator) SwitchRingDegreeN2ToN1New(ctN2 *rlwe.Ciphertext) (ctN1 *rlwe.Ciphertext) {
	ctN1 = hefloat.NewCiphertext(eval.ResidualParameters, 1, ctN2.Level())

	// Sanity check, this error should never happen unless this algorithm has been improperly
	// modified to pass invalid inputs.
	if err := eval.Evaluator.ApplyEvaluationKey(ctN2, eval.EvkN2ToN1, ctN1); err != nil {
		panic(err)
	}
	return
}

func (eval Evaluator) ComplexToRealNew(ctCmplx *rlwe.Ciphertext) (ctReal *rlwe.Ciphertext) {
	ctReal = hefloat.NewCiphertext(eval.ResidualParameters, 1, ctCmplx.Level())

	// Sanity check, this error should never happen unless this algorithm has been improperly
	// modified to pass invalid inputs.
	if err := eval.DomainSwitcher.ComplexToReal(&eval.Evaluator.Evaluator, ctCmplx, ctReal); err != nil {
		panic(err)
	}
	return
}

func (eval Evaluator) RealToComplexNew(ctReal *rlwe.Ciphertext) (ctCmplx *rlwe.Ciphertext) {
	ctCmplx = hefloat.NewCiphertext(eval.BootstrappingParameters, 1, ctReal.Level())

	// Sanity check, this error should never happen unless this algorithm has been improperly
	// modified to pass invalid inputs.
	if err := eval.DomainSwitcher.RealToComplex(&eval.Evaluator.Evaluator, ctReal, ctCmplx); err != nil {
		panic(err)
	}
	return
}

func (eval Evaluator) PackAndSwitchN1ToN2(cts []rlwe.Ciphertext) ([]rlwe.Ciphertext, error) {

	var err error

	if eval.ResidualParameters.N() != eval.BootstrappingParameters.N() {
		if cts, err = eval.Pack(cts, eval.ResidualParameters, eval.xPow2N1); err != nil {
			return nil, fmt.Errorf("cannot PackAndSwitchN1ToN2: PackN1: %w", err)
		}

		for i := range cts {
			cts[i] = *eval.SwitchRingDegreeN1ToN2New(&cts[i])
		}
	}

	if cts, err = eval.Pack(cts, eval.BootstrappingParameters, eval.xPow2N2); err != nil {
		return nil, fmt.Errorf("cannot PackAndSwitchN1ToN2: PackN2: %w", err)
	}

	return cts, nil
}

func (eval Evaluator) UnpackAndSwitchN2Tn1(cts []rlwe.Ciphertext, LogSlots, Nb int) ([]rlwe.Ciphertext, error) {

	var err error

	if eval.ResidualParameters.N() != eval.BootstrappingParameters.N() {
		if cts, err = eval.UnPack(cts, eval.BootstrappingParameters, LogSlots, Nb, eval.xPow2InvN2); err != nil {
			return nil, fmt.Errorf("cannot UnpackAndSwitchN2Tn1: UnpackN2: %w", err)
		}

		for i := range cts {
			cts[i] = *eval.SwitchRingDegreeN2ToN1New(&cts[i])
		}
	}

	for i := range cts {
		cts[i].LogDimensions.Cols = LogSlots
	}

	return cts, nil
}

func (eval Evaluator) UnPack(cts []rlwe.Ciphertext, params hefloat.Parameters, LogSlots, Nb int, xPow2Inv []ring.Poly) ([]rlwe.Ciphertext, error) {
	LogGap := params.LogMaxSlots() - LogSlots

	if LogGap == 0 {
		return cts, nil
	}

	cts = append(cts, make([]rlwe.Ciphertext, Nb-1)...)

	for i := 1; i < len(cts); i++ {
		cts[i] = *cts[0].CopyNew()
	}

	r := params.RingQ().AtLevel(cts[0].Level())

	N := len(cts)

	for i := 0; i < utils.Min(bits.Len64(uint64(N-1)), LogGap); i++ {

		step := 1 << (i + 1)

		for j := 0; j < N; j += step {

			for k := step >> 1; k < step; k++ {

				if (j + k) >= N {
					break
				}

				r.MulCoeffsMontgomery(cts[j+k].Value[0], xPow2Inv[i], cts[j+k].Value[0])
				r.MulCoeffsMontgomery(cts[j+k].Value[1], xPow2Inv[i], cts[j+k].Value[1])
			}
		}
	}

	return cts, nil
}

func (eval Evaluator) Pack(cts []rlwe.Ciphertext, params hefloat.Parameters, xPow2 []ring.Poly) ([]rlwe.Ciphertext, error) {

	var LogSlots = cts[0].LogSlots()
	RingDegree := params.N()

	for i, ct := range cts {
		if N := ct.LogSlots(); N != LogSlots {
			return nil, fmt.Errorf("cannot Pack: cts[%d].PlaintextLogSlots()=%d != cts[0].PlaintextLogSlots=%d", i, N, LogSlots)
		}

		if N := ct.Value[0].N(); N != RingDegree {
			return nil, fmt.Errorf("cannot Pack: cts[%d].Value[0].N()=%d != params.N()=%d", i, N, RingDegree)
		}
	}

	LogGap := params.LogMaxSlots() - LogSlots

	if LogGap == 0 {
		return cts, nil
	}

	for i := 0; i < LogGap; i++ {

		for j := 0; j < len(cts)>>1; j++ {

			eve := cts[j*2+0]
			odd := cts[j*2+1]

			level := utils.Min(eve.Level(), odd.Level())

			r := params.RingQ().AtLevel(level)

			r.MulCoeffsMontgomeryThenAdd(odd.Value[0], xPow2[i], eve.Value[0])
			r.MulCoeffsMontgomeryThenAdd(odd.Value[1], xPow2[i], eve.Value[1])

			cts[j] = eve
		}

		if len(cts)&1 == 1 {
			cts[len(cts)>>1] = cts[len(cts)-1]
			cts = cts[:len(cts)>>1+1]
		} else {
			cts = cts[:len(cts)>>1]
		}
	}

	LogMaxDimensions := params.LogMaxDimensions()
	for i := range cts {
		cts[i].LogDimensions = LogMaxDimensions
	}

	return cts, nil
}
