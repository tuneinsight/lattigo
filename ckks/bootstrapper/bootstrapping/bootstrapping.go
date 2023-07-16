// Package bootstrapping implement the bootstrapping for the CKKS scheme.
package bootstrapping

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

func (btp Bootstrapper) StartingLevel() int {
	return btp.params.PlaintextScaleToModuliRatio() - 1
}

// Bootstrap re-encrypts a ciphertext to a ciphertext at MaxLevel - k where k is the depth of the bootstrapping circuit.
// If the input ciphertext level is zero, the input scale must be an exact power of two smaller than Q[0]/MessageRatio
// (it can't be equal since Q[0] is not a power of two).
// The message ratio is an optional field in the bootstrapping parameters, by default it set to 2^{LogMessageRatio = 8}.
// See the bootstrapping parameters for more information about the message ratio or other parameters related to the bootstrapping.
// If the input ciphertext is at level one or more, the input scale does not need to be an exact power of two as one level
// can be used to do a scale matching.
func (btp Bootstrapper) Bootstrap(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext, err error) {

	// Pre-processing
	ctDiff := ctIn.CopyNew()

	var errScale *rlwe.Scale

	// [M^{d}/q1]
	if ctDiff, errScale, err = btp.scaleDownToQ0OverMessageRatio(ctDiff); err != nil {
		return nil, err
	}

	// [M^{d}/q1 + e^{d-logprec}]
	if ctOut, err = btp.bootstrap(ctDiff.CopyNew()); err != nil {
		return nil, err
	}

	// Error correcting factor of the approximate division by q1
	ctOut.PlaintextScale = ctOut.PlaintextScale.Mul(*errScale)

	// Stores by how much a ciphertext must be scaled to get back
	// to the input scale
	diffScale := ctIn.PlaintextScale.Div(ctOut.PlaintextScale).Bigint()

	// [M^{d} + e^{d-logprec}]
	if err = btp.Mul(ctOut, diffScale, ctOut); err != nil {
		return nil, err
	}
	ctOut.PlaintextScale = ctIn.PlaintextScale

	if btp.IterationsParameters != nil {

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
			scale := new(big.Int).Set(diffScale)
			bignum.DivRound(scale, prec, scale)

			// Checks that round(q1/logprec) >= 2^{logprec}
			requiresReservedPrime := scale.Cmp(new(big.Int).SetUint64(1)) < 0

			if requiresReservedPrime && btp.IterationsParameters.ReservedPrimeBitSize == 0 {
				return ctOut, fmt.Errorf("warning: early stopping at iteration k=%d: reason: round(q1/2^{logprec}) < 1 and no reserverd prime was provided", i+1)
			}

			// [M^{d} + e^{d-logprec}] - [M^{d}] -> [e^{d-logprec}]
			tmp, err := btp.SubNew(ctOut, ctIn)

			if err != nil {
				return nil, err
			}

			// prec * [e^{d-logprec}] -> [e^{d}]
			if err = btp.Mul(tmp, prec, tmp); err != nil {
				return nil, err
			}

			tmp.PlaintextScale = ctOut.PlaintextScale

			// [e^{d}] / q1 -> [e^{d}/q1]
			if tmp, errScale, err = btp.scaleDownToQ0OverMessageRatio(tmp); err != nil {
				return nil, err
			}

			// [e^{d}/q1] -> [e^{d}/q1 + e'^{d-logprec}]
			if tmp, err = btp.bootstrap(tmp); err != nil {
				return nil, err
			}

			tmp.PlaintextScale = tmp.PlaintextScale.Mul(*errScale)

			// [[e^{d}/q1 + e'^{d-logprec}] * q1/logprec -> [e^{d-logprec} + e'^{d-2logprec}*q1]
			// If scale > 2^{logprec}, then we ensure a precision of at least 2^{logprec} even with a rounding of the scale
			if !requiresReservedPrime {
				if err = btp.Mul(tmp, scale, tmp); err != nil {
					return nil, err
				}
			} else {

				// Else we compute the floating point ratio
				ss := new(big.Float).SetInt(diffScale)
				ss.Quo(ss, new(big.Float).SetInt(prec))

				// Do a scaled multiplication by the last prime
				if err = btp.Mul(tmp, ss, tmp); err != nil {
					return nil, err
				}

				// And rescale
				if err = btp.Rescale(tmp, btp.params.PlaintextScale(), tmp); err != nil {
					return nil, err
				}
			}

			// This is a given
			tmp.PlaintextScale = ctOut.PlaintextScale

			// [M^{d} + e^{d-logprec}] - [e^{d-logprec} + e'^{d-2logprec}*q1] -> [M^{d} + e'^{d-2logprec}*q1]
			if err = btp.Sub(ctOut, tmp, ctOut); err != nil {
				return nil, err
			}
		}
	}

	return
}

func currentMessageRatioIsGreaterOrEqualToLastPrimeTimesTargetMessageRatio(ct *rlwe.Ciphertext, msgRatio float64, r *ring.Ring) bool {
	level := ct.Level()
	currentMessageRatio := rlwe.NewScale(r.ModulusAtLevel[level])
	currentMessageRatio = currentMessageRatio.Div(ct.PlaintextScale)
	return currentMessageRatio.Cmp(rlwe.NewScale(r.SubRings[level].Modulus).Mul(rlwe.NewScale(msgRatio))) > -1
}

// The purpose of this pre-processing step is to bring the ciphertext level to zero and scaling factor to Q[0]/MessageRatio
func (btp Bootstrapper) scaleDownToQ0OverMessageRatio(ctIn *rlwe.Ciphertext) (*rlwe.Ciphertext, *rlwe.Scale, error) {

	params := &btp.params

	r := params.RingQ()

	// Removes unecessary primes
	for ctIn.Level() != 0 && currentMessageRatioIsGreaterOrEqualToLastPrimeTimesTargetMessageRatio(ctIn, btp.evalModPoly.MessageRatio(), r) {
		ctIn.Resize(ctIn.Degree(), ctIn.Level()-1)
	}

	// Current Message Ratio
	currentMessageRatio := rlwe.NewScale(r.ModulusAtLevel[ctIn.Level()])
	currentMessageRatio = currentMessageRatio.Div(ctIn.PlaintextScale)

	// Desired Message Ratio
	targetMessageRatio := rlwe.NewScale(btp.evalModPoly.MessageRatio())

	// (Current Message Ratio) / (Desired Message Ratio)
	scaleUp := currentMessageRatio.Div(targetMessageRatio)

	if scaleUp.Cmp(rlwe.NewScale(0.5)) == -1 {
		return nil, nil, fmt.Errorf("cannot scaleDownToQ0OverMessageRatio: initial Q/PlaintextScale < 0.5*Q[0]/MessageRatio")
	}

	scaleUpBigint := scaleUp.Bigint()

	if err := btp.Mul(ctIn, scaleUpBigint, ctIn); err != nil {
		return nil, nil, fmt.Errorf("cannot scaleDownToQ0OverMessageRatio: %w", err)
	}

	ctIn.PlaintextScale = ctIn.PlaintextScale.Mul(rlwe.NewScale(scaleUpBigint))

	// errScale = CtIn.Scale/(Q[0]/MessageRatio)
	targetScale := new(big.Float).SetPrec(256).SetInt(r.ModulusAtLevel[0])
	targetScale.Quo(targetScale, new(big.Float).SetFloat64(btp.evalModPoly.MessageRatio()))

	if ctIn.Level() != 0 {
		if err := btp.Rescale(ctIn, rlwe.NewScale(targetScale), ctIn); err != nil {
			return nil, nil, fmt.Errorf("cannot scaleDownToQ0OverMessageRatio: %w", err)
		}
	}

	errScale := ctIn.PlaintextScale.Div(rlwe.NewScale(targetScale))

	return ctIn, &errScale, nil
}

func (btp *Bootstrapper) bootstrap(ctIn *rlwe.Ciphertext) (opOut *rlwe.Ciphertext, err error) {

	// Step 1 : Extend the basis from q to Q
	if opOut, err = btp.modUpFromQ0(ctIn); err != nil {
		return
	}

	// Scale the message from Q0/|m| to QL/|m|, where QL is the largest modulus used during the bootstrapping.
	if scale := (btp.evalModPoly.ScalingFactor().Float64() / btp.evalModPoly.MessageRatio()) / opOut.PlaintextScale.Float64(); scale > 1 {
		if err = btp.ScaleUp(opOut, rlwe.NewScale(scale), opOut); err != nil {
			return nil, err
		}
	}

	//SubSum X -> (N/dslots) * Y^dslots
	if err = btp.Trace(opOut, opOut.PlaintextLogDimensions[1], opOut); err != nil {
		return nil, err
	}

	// Step 2 : CoeffsToSlots (Homomorphic encoding)
	ctReal, ctImag, err := btp.CoeffsToSlotsNew(opOut, btp.ctsMatrices)
	if err != nil {
		return nil, err
	}

	// Step 3 : EvalMod (Homomorphic modular reduction)
	// ctReal = Ecd(real)
	// ctImag = Ecd(imag)
	// If n < N/2 then ctReal = Ecd(real|imag)
	if ctReal, err = btp.EvalModNew(ctReal, btp.evalModPoly); err != nil {
		return nil, err
	}
	ctReal.PlaintextScale = btp.params.PlaintextScale()

	if ctImag != nil {
		if ctImag, err = btp.EvalModNew(ctImag, btp.evalModPoly); err != nil {
			return nil, err
		}
		ctImag.PlaintextScale = btp.params.PlaintextScale()
	}

	// Step 4 : SlotsToCoeffs (Homomorphic decoding)
	opOut, err = btp.SlotsToCoeffsNew(ctReal, ctImag, btp.stcMatrices)

	return
}

func (btp *Bootstrapper) modUpFromQ0(ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {

	if btp.EvkDtS != nil {
		if err := btp.ApplyEvaluationKey(ct, btp.EvkDtS, ct); err != nil {
			return nil, err
		}
	}

	ringQ := btp.params.RingQ().AtLevel(ct.Level())
	ringP := btp.params.RingP()

	for i := range ct.Value {
		ringQ.INTT(ct.Value[i], ct.Value[i])
	}

	// Extend the ciphertext with zero polynomials.
	ct.Resize(ct.Degree(), btp.params.MaxLevel())

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

	// ModUp q->Q for ct[0] centered around q
	for j := 0; j < N; j++ {

		coeff = ct.Value[0].Coeffs[0][j]
		pos, neg = 1, 0
		if coeff >= (q >> 1) {
			coeff = q - coeff
			pos, neg = 0, 1
		}

		for i := 1; i < levelQ+1; i++ {
			tmp = ring.BRedAdd(coeff, Q[i], BRCQ[i])
			ct.Value[0].Coeffs[i][j] = tmp*pos + (Q[i]-tmp)*neg
		}
	}

	if btp.EvkStD != nil {

		ks := btp.Evaluator.Evaluator

		// ModUp q->QP for ct[1] centered around q
		for j := 0; j < N; j++ {

			coeff = ct.Value[1].Coeffs[0][j]
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

		ringQ.NTT(ct.Value[0], ct.Value[0])

		ctTmp := &rlwe.Ciphertext{}
		ctTmp.Value = []ring.Poly{ks.BuffQP[1].Q, ct.Value[1]}
		ctTmp.MetaData = ct.MetaData

		ks.GadgetProductHoisted(levelQ, ks.BuffDecompQP, &btp.EvkStD.GadgetCiphertext, ctTmp)
		ringQ.Add(ct.Value[0], ctTmp.Value[0], ct.Value[0])

	} else {

		for j := 0; j < N; j++ {

			coeff = ct.Value[1].Coeffs[0][j]
			pos, neg = 1, 0
			if coeff >= (q >> 1) {
				coeff = q - coeff
				pos, neg = 0, 1
			}

			for i := 1; i < levelQ+1; i++ {
				tmp = ring.BRedAdd(coeff, Q[i], BRCQ[i])
				ct.Value[1].Coeffs[i][j] = tmp*pos + (Q[i]-tmp)*neg
			}
		}

		ringQ.NTT(ct.Value[0], ct.Value[0])
		ringQ.NTT(ct.Value[1], ct.Value[1])
	}

	return ct, nil
}
