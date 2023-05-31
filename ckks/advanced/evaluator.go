// Package advanced implements advanced operations for the CKKS scheme.
package advanced

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

type Evaluator struct {
	*ckks.Evaluator
}

// NewEvaluator creates a new Evaluator.
func NewEvaluator(params ckks.Parameters, evk rlwe.EvaluationKeySetInterface) *Evaluator {
	return &Evaluator{ckks.NewEvaluator(params, evk)}
}

// ShallowCopy creates a shallow copy of this Evaluator in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Evaluators can be used concurrently.
func (eval *Evaluator) ShallowCopy() *Evaluator {
	return &Evaluator{eval.Evaluator.ShallowCopy()}
}

// WithKey creates a shallow copy of the receiver Evaluator for which the new EvaluationKey is evaluationKey
// and where the temporary buffers are shared. The receiver and the returned Evaluators cannot be used concurrently.
func (eval *Evaluator) WithKey(evk rlwe.EvaluationKeySetInterface) *Evaluator {
	return &Evaluator{eval.Evaluator.WithKey(evk)}
}

// CoeffsToSlotsNew applies the homomorphic encoding and returns the result on new ciphertexts.
// Homomorphically encodes a complex vector vReal + i*vImag.
// If the packing is sparse (n < N/2), then returns ctReal = Ecd(vReal || vImag) and ctImag = nil.
// If the packing is dense (n == N/2), then returns ctReal = Ecd(vReal) and ctImag = Ecd(vImag).
func (eval *Evaluator) CoeffsToSlotsNew(ctIn *rlwe.Ciphertext, ctsMatrices HomomorphicDFTMatrix) (ctReal, ctImag *rlwe.Ciphertext) {
	ctReal = ckks.NewCiphertext(eval.Parameters(), 1, ctsMatrices.LevelStart)

	if maxLogSlots := eval.Parameters().MaxLogSlots()[1]; ctsMatrices.LogSlots == maxLogSlots {
		ctImag = ckks.NewCiphertext(eval.Parameters(), 1, ctsMatrices.LevelStart)
	}

	eval.CoeffsToSlots(ctIn, ctsMatrices, ctReal, ctImag)
	return
}

// CoeffsToSlots applies the homomorphic encoding and returns the results on the provided ciphertexts.
// Homomorphically encodes a complex vector vReal + i*vImag of size n on a real vector of size 2n.
// If the packing is sparse (n < N/2), then returns ctReal = Ecd(vReal || vImag) and ctImag = nil.
// If the packing is dense (n == N/2), then returns ctReal = Ecd(vReal) and ctImag = Ecd(vImag).
func (eval *Evaluator) CoeffsToSlots(ctIn *rlwe.Ciphertext, ctsMatrices HomomorphicDFTMatrix, ctReal, ctImag *rlwe.Ciphertext) {

	if ctsMatrices.RepackImag2Real {

		zV := ctIn.CopyNew()

		eval.dft(ctIn, ctsMatrices.Matrices, zV)

		eval.Conjugate(zV, ctReal)

		var tmp *rlwe.Ciphertext
		if ctImag != nil {
			tmp = ctImag
		} else {
			tmp = rlwe.NewCiphertextAtLevelFromPoly(ctReal.Level(), eval.BuffCt.Value[:2])
			tmp.IsNTT = true
		}

		// Imag part
		eval.Sub(zV, ctReal, tmp)
		eval.Mul(tmp, -1i, tmp)

		// Real part
		eval.Add(ctReal, zV, ctReal)

		// If repacking, then ct0 and ct1 right n/2 slots are zero.
		if maxLogSlots := eval.Parameters().MaxLogSlots()[1]; ctsMatrices.LogSlots < maxLogSlots {
			eval.Rotate(tmp, ctIn.Slots()[1], tmp)
			eval.Add(ctReal, tmp, ctReal)
		}

		zV = nil

	} else {
		eval.dft(ctIn, ctsMatrices.Matrices, ctReal)
	}
}

// SlotsToCoeffsNew applies the homomorphic decoding and returns the result on a new ciphertext.
// Homomorphically decodes a real vector of size 2n on a complex vector vReal + i*vImag of size n.
// If the packing is sparse (n < N/2) then ctReal = Ecd(vReal || vImag) and ctImag = nil.
// If the packing is dense (n == N/2), then ctReal = Ecd(vReal) and ctImag = Ecd(vImag).
func (eval *Evaluator) SlotsToCoeffsNew(ctReal, ctImag *rlwe.Ciphertext, stcMatrices HomomorphicDFTMatrix) (ctOut *rlwe.Ciphertext) {

	if ctReal.Level() < stcMatrices.LevelStart || (ctImag != nil && ctImag.Level() < stcMatrices.LevelStart) {
		panic("ctReal.Level() or ctImag.Level() < HomomorphicDFTMatrix.LevelStart")
	}

	ctOut = ckks.NewCiphertext(eval.Parameters(), 1, stcMatrices.LevelStart)
	eval.SlotsToCoeffs(ctReal, ctImag, stcMatrices, ctOut)
	return

}

// SlotsToCoeffs applies the homomorphic decoding and returns the result on the provided ciphertext.
// Homomorphically decodes a real vector of size 2n on a complex vector vReal + i*vImag of size n.
// If the packing is sparse (n < N/2) then ctReal = Ecd(vReal || vImag) and ctImag = nil.
// If the packing is dense (n == N/2), then ctReal = Ecd(vReal) and ctImag = Ecd(vImag).
func (eval *Evaluator) SlotsToCoeffs(ctReal, ctImag *rlwe.Ciphertext, stcMatrices HomomorphicDFTMatrix, ctOut *rlwe.Ciphertext) {
	// If full packing, the repacking can be done directly using ct0 and ct1.
	if ctImag != nil {
		eval.Mul(ctImag, 1i, ctOut)
		eval.Add(ctOut, ctReal, ctOut)
		eval.dft(ctOut, stcMatrices.Matrices, ctOut)
	} else {
		eval.dft(ctReal, stcMatrices.Matrices, ctOut)
	}
}

func (eval *Evaluator) dft(ctIn *rlwe.Ciphertext, plainVectors []rlwe.LinearTransform, ctOut *rlwe.Ciphertext) {

	inputLogSlots := ctIn.LogSlots

	// Sequentially multiplies w with the provided dft matrices.
	scale := ctIn.Scale
	var in, out *rlwe.Ciphertext
	for i, plainVector := range plainVectors {
		in, out = ctOut, ctOut
		if i == 0 {
			in, out = ctIn, ctOut
		}

		eval.LinearTransform(in, plainVector, []*rlwe.Ciphertext{out})

		if err := eval.Rescale(out, scale, out); err != nil {
			panic(err)
		}
	}

	// Encoding matrices are a special case of `fractal` linear transform
	// that doesn't change the underlying plaintext polynomial Y = X^{N/n}
	// of the input ciphertext.
	ctOut.LogSlots = inputLogSlots
}

// EvalModNew applies a homomorphic mod Q on a vector scaled by Delta, scaled down to mod 1 :
//
//  1. Delta * (Q/Delta * I(X) + m(X)) (Delta = scaling factor, I(X) integer poly, m(X) message)
//  2. Delta * (I(X) + Delta/Q * m(X)) (divide by Q/Delta)
//  3. Delta * (Delta/Q * m(X)) (x mod 1)
//  4. Delta * (m(X)) (multiply back by Q/Delta)
//
// Since Q is not a power of two, but Delta is, then does an approximate division by the closest
// power of two to Q instead. Hence, it assumes that the input plaintext is already scaled by
// the correcting factor Q/2^{round(log(Q))}.
//
// !! Assumes that the input is normalized by 1/K for K the range of the approximation.
//
// Scaling back error correction by 2^{round(log(Q))}/Q afterward is included in the polynomial
func (eval *Evaluator) EvalModNew(ct *rlwe.Ciphertext, evalModPoly EvalModPoly) *rlwe.Ciphertext {

	if ct.Level() < evalModPoly.LevelStart() {
		panic("ct.Level() < evalModPoly.LevelStart")
	}

	if ct.Level() > evalModPoly.LevelStart() {
		eval.DropLevel(ct, ct.Level()-evalModPoly.LevelStart())
	}

	// Stores default scales
	prevScaleCt := ct.Scale

	// Normalize the modular reduction to mod by 1 (division by Q)
	ct.Scale = evalModPoly.ScalingFactor()

	var err error

	// Compute the scales that the ciphertext should have before the double angle
	// formula such that after it it has the scale it had before the polynomial
	// evaluation

	Qi := eval.Parameters().Q()

	targetScale := ct.Scale
	for i := 0; i < evalModPoly.doubleAngle; i++ {
		targetScale = targetScale.Mul(rlwe.NewScale(Qi[evalModPoly.levelStart-evalModPoly.sinePoly.Depth()-evalModPoly.doubleAngle+i+1]))
		targetScale.Value.Sqrt(&targetScale.Value)
	}

	// Division by 1/2^r and change of variable for the Chebyshev evaluation
	if evalModPoly.sineType == CosDiscrete || evalModPoly.sineType == CosContinuous {
		offset := new(big.Float).Sub(&evalModPoly.sinePoly.B, &evalModPoly.sinePoly.A)
		offset.Mul(offset, new(big.Float).SetFloat64(evalModPoly.scFac))
		offset.Quo(new(big.Float).SetFloat64(-0.5), offset)
		eval.Add(ct, offset, ct)
	}

	// Chebyshev evaluation
	if ct, err = eval.Polynomial(ct, evalModPoly.sinePoly, rlwe.NewScale(targetScale)); err != nil {
		panic(err)
	}

	// Double angle
	sqrt2pi := evalModPoly.sqrt2Pi
	for i := 0; i < evalModPoly.doubleAngle; i++ {
		sqrt2pi *= sqrt2pi
		eval.MulRelin(ct, ct, ct)
		eval.Add(ct, ct, ct)
		eval.Add(ct, -sqrt2pi, ct)
		if err := eval.Rescale(ct, rlwe.NewScale(targetScale), ct); err != nil {
			panic(err)
		}
	}

	// ArcSine
	if evalModPoly.arcSinePoly != nil {
		if ct, err = eval.Polynomial(ct, evalModPoly.arcSinePoly, ct.Scale); err != nil {
			panic(err)
		}
	}

	// Multiplies back by q
	ct.Scale = prevScaleCt
	return ct
}
