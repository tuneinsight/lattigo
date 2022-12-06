// Package bootstrapping implement the bootstrapping for the CKKS scheme.
package bootstrapping

import (
	"math"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// Bootstrap re-encrypts a ciphertext at lvl Q0 to a ciphertext at MaxLevel-k where k is the depth of the bootstrapping circuit.
// If the input ciphertext level is zero, the input scale must be an exact power of two smaller or equal to round(Q0/2^{10}).
// If the input ciphertext is at level one or more, the input scale does not need to be an exact power of two as one level
// can be used to do a scale matching.
func (btp *Bootstrapper) Bootstrap(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {

	ctOut = ctIn.CopyNew()

	// Drops the level to 1
	for ctOut.Level() > 1 {
		btp.DropLevel(ctOut, 1)
	}

	// Brings the ciphertext scale to Q0/MessageRatio
	if ctOut.Level() == 1 {

		// If one level is available, then uses it to match the scale
		btp.SetScale(ctOut, rlwe.NewScale(btp.q0OverMessageRatio))

		// Then drops to level 0
		for ctOut.Level() != 0 {
			btp.DropLevel(ctOut, 1)
		}

	} else {

		// Does an integer constant mult by round((Q0/Delta_m)/ctscle)
		if btp.q0OverMessageRatio < ctOut.Scale.Float64() {
			panic("Cannot bootstrap: ciphetext scale > q/||m||)")
		}

		btp.ScaleUp(ctOut, rlwe.NewScale(math.Round(btp.q0OverMessageRatio/ctOut.Scale.Float64())), ctOut)
	}

	// Scales the message to Q0/|m|, which is the maximum possible before ModRaise to avoid plaintext overflow.
	if scale := math.Round((btp.params.QiFloat64(0) / btp.evalModPoly.MessageRatio()) / ctOut.Scale.Float64()); scale > 1 {
		btp.ScaleUp(ctOut, rlwe.NewScale(scale), ctOut)
	}

	// Step 1 : Extend the basis from q to Q
	ctOut = btp.modUpFromQ0(ctOut)

	// Scale the message from Q0/|m| to QL/|m|, where QL is the largest modulus used during the bootstrapping.
	if scale := (btp.evalModPoly.ScalingFactor().Float64() / btp.evalModPoly.MessageRatio()) / ctOut.Scale.Float64(); scale > 1 {
		btp.ScaleUp(ctOut, rlwe.NewScale(scale), ctOut)
	}

	//SubSum X -> (N/dslots) * Y^dslots
	btp.Trace(ctOut, btp.params.LogSlots(), ctOut)

	// Step 2 : CoeffsToSlots (Homomorphic encoding)
	ctReal, ctImag := btp.CoeffsToSlotsNew(ctOut, btp.ctsMatrices)

	// Step 3 : EvalMod (Homomorphic modular reduction)
	// ctReal = Ecd(real)
	// ctImag = Ecd(imag)
	// If n < N/2 then ctReal = Ecd(real|imag)
	ctReal = btp.EvalModNew(ctReal, btp.evalModPoly)
	ctReal.Scale = btp.params.DefaultScale()

	if ctImag != nil {
		ctImag = btp.EvalModNew(ctImag, btp.evalModPoly)
		ctImag.Scale = btp.params.DefaultScale()
	}

	// Step 4 : SlotsToCoeffs (Homomorphic decoding)
	ctOut = btp.SlotsToCoeffsNew(ctReal, ctImag, btp.stcMatrices)

	return
}

func (btp *Bootstrapper) modUpFromQ0(ct *rlwe.Ciphertext) *rlwe.Ciphertext {

	if btp.swkDtS != nil {
		btp.SwitchKeys(ct, btp.swkDtS, ct)
	}

	ringQ := btp.params.RingQ()
	ringP := btp.params.RingP()

	for i := range ct.Value {
		ringQ.InvNTTLvl(ct.Level(), ct.Value[i], ct.Value[i])
	}

	// Extend the ciphertext with zero polynomials.
	ct.Resize(ct.Degree(), btp.params.MaxLevel())

	levelQ := btp.params.QCount() - 1
	levelP := btp.params.PCount() - 1

	Q := ringQ.Moduli()
	P := ringP.Moduli()
	q := Q[0]
	bredparamsQ := ringQ.BRedParams()
	bredparamsP := ringP.BRedParams()

	var coeff, tmp, pos, neg uint64

	// ModUp q->Q for ct[0] centered around q
	for j := 0; j < btp.params.N(); j++ {

		coeff = ct.Value[0].Coeffs[0][j]
		pos, neg = 1, 0
		if coeff >= (q >> 1) {
			coeff = q - coeff
			pos, neg = 0, 1
		}

		for i := 1; i < levelQ+1; i++ {
			tmp = ring.BRedAdd(coeff, Q[i], bredparamsQ[i])
			ct.Value[0].Coeffs[i][j] = tmp*pos + (Q[i]-tmp)*neg
		}
	}

	if btp.swkStD != nil {

		ks := btp.GetRLWEEvaluator()

		// ModUp q->QP for ct[1] centered around q
		for j := 0; j < btp.params.N(); j++ {

			coeff = ct.Value[1].Coeffs[0][j]
			pos, neg = 1, 0
			if coeff > (q >> 1) {
				coeff = q - coeff
				pos, neg = 0, 1
			}

			for i := 0; i < levelQ+1; i++ {
				tmp = ring.BRedAdd(coeff, Q[i], bredparamsQ[i])
				ks.BuffDecompQP[0].Q.Coeffs[i][j] = tmp*pos + (Q[i]-tmp)*neg

			}

			for i := 0; i < levelP+1; i++ {
				tmp = ring.BRedAdd(coeff, P[i], bredparamsP[i])
				ks.BuffDecompQP[0].P.Coeffs[i][j] = tmp*pos + (P[i]-tmp)*neg
			}
		}

		for i := len(ks.BuffDecompQP) - 1; i >= 0; i-- {
			ringQ.NTTLvl(levelQ, ks.BuffDecompQP[0].Q, ks.BuffDecompQP[i].Q)
		}

		for i := len(ks.BuffDecompQP) - 1; i >= 0; i-- {
			ringP.NTTLvl(levelP, ks.BuffDecompQP[0].P, ks.BuffDecompQP[i].P)
		}

		ringQ.NTTLvl(levelQ, ct.Value[0], ct.Value[0])

		ks.KeyswitchHoisted(levelQ, ks.BuffDecompQP, btp.swkStD, ks.BuffQP[1].Q, ct.Value[1], ks.BuffQP[1].P, ks.BuffQP[2].P)
		ringQ.AddLvl(levelQ, ct.Value[0], ks.BuffQP[1].Q, ct.Value[0])

	} else {

		for j := 0; j < btp.params.N(); j++ {

			coeff = ct.Value[1].Coeffs[0][j]
			pos, neg = 1, 0
			if coeff >= (q >> 1) {
				coeff = q - coeff
				pos, neg = 0, 1
			}

			for i := 1; i < levelQ+1; i++ {
				tmp = ring.BRedAdd(coeff, Q[i], bredparamsQ[i])
				ct.Value[1].Coeffs[i][j] = tmp*pos + (Q[i]-tmp)*neg
			}
		}

		ringQ.NTTLvl(levelQ, ct.Value[0], ct.Value[0])
		ringQ.NTTLvl(levelQ, ct.Value[1], ct.Value[1])
	}

	return ct
}
