package bootstrapping

import (
	"math"

	"github.com/tuneinsight/lattigo/v3/ckks"
	"github.com/tuneinsight/lattigo/v3/ring"
)

// Bootstrapp re-encrypt a ciphertext at lvl Q0 to a ciphertext at MaxLevel-k where k is the depth of the bootstrapping circuit.
// If the input ciphertext level is zero, the input scale must be an exact power of two smaller or equal to round(Q0/2^{10}).
// If the input ciphertext is at level one or more, the input scale does not need to be an exact power of two as one level
// can be used to do a scale matching.
func (btp *Bootstrapper) Bootstrapp(ctIn *ckks.Ciphertext) (ctOut *ckks.Ciphertext) {

	ctOut = ctIn.CopyNew()

	// Drops the level to 1
	for ctOut.Level() > 1 {
		btp.DropLevel(ctOut, 1)
	}

	// Brings the ciphertext scale to Q0/MessageRatio
	if ctOut.Level() == 1 {

		// If one level is available, then uses it to match the scale
		btp.SetScale(ctOut, btp.q0OverMessageRatio)

		// Then drops to level 0
		for ctOut.Level() != 0 {
			btp.DropLevel(ctOut, 1)
		}

	} else {

		// Does an integer constant mult by round((Q0/Delta_m)/ctscle)
		if btp.q0OverMessageRatio < ctOut.Scale {
			panic("ciphetext scale > q/||m||)")
		}

		btp.ScaleUp(ctOut, math.Round(btp.q0OverMessageRatio/ctOut.Scale), ctOut)
	}

	// Scales the message to Q0/|m|, which is the maximum possible before ModRaise to avoid plaintext overflow.
	if math.Round((btp.params.QiFloat64(0)/btp.evalModPoly.MessageRatio())/ctOut.Scale) > 1 {
		btp.ScaleUp(ctOut, math.Round((btp.params.QiFloat64(0)/btp.evalModPoly.MessageRatio())/ctOut.Scale), ctOut)
	}

	// Step 1 : Extend the basis from q to Q
	ctOut = btp.modUpFromQ0(ctOut)

	// Scale the message from Q0/|m| to QL/|m|, where QL is the largest modulus used during the bootstrapping.
	if (btp.evalModPoly.ScalingFactor()/btp.evalModPoly.MessageRatio())/ctOut.Scale > 1 {
		btp.ScaleUp(ctOut, math.Round((btp.evalModPoly.ScalingFactor()/btp.evalModPoly.MessageRatio())/ctOut.Scale), ctOut)
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

func (btp *Bootstrapper) modUpFromQ0(ct *ckks.Ciphertext) *ckks.Ciphertext {

<<<<<<< btp_eprint
	if btp.swkDtS != nil {
		btp.SwitchKeys(ct, btp.swkDtS, ct)
	}

=======
<<<<<<< 83ae36f5f9908381fe0d957ce0daa4f037d38e6f
=======
	btp.SwitchKeys(ct, btp.swkDtS, ct)

	ks := btp.GetRLWEEvaluator()

>>>>>>> rebased onto btp_eprint
>>>>>>> rebased onto btp_eprint
	ringQ := btp.params.RingQ()
	ringP := btp.params.RingP()

	for i := range ct.Value {
		ringQ.InvNTTLvl(ct.Level(), ct.Value[i], ct.Value[i])
	}

	// Extend the ciphertext with zero polynomials.
	for u := range ct.Value {
		ct.Value[u].Coeffs = append(ct.Value[u].Coeffs, make([][]uint64, btp.params.MaxLevel())...)
		for i := 1; i < btp.params.MaxLevel()+1; i++ {
			ct.Value[u].Coeffs[i] = make([]uint64, btp.params.N())
		}
	}

	levelQ := btp.params.QCount() - 1
	levelP := btp.params.PCount() - 1

	Q := ringQ.Modulus
	P := ringP.Modulus
	q := Q[0]
	bredparamsQ := ringQ.BredParams
	bredparamsP := ringP.BredParams

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

		ks := btp.GetKeySwitcher()

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

<<<<<<< btp_eprint
		ks.KeyswitchHoisted(levelQ, ks.BuffDecompQP, btp.swkStD, ks.BuffQP[1].Q, ct.Value[1], ks.BuffQP[1].P, ks.BuffQP[2].P)
		ringQ.AddLvl(levelQ, ct.Value[0], ks.BuffQP[1].Q, ct.Value[0])
=======
<<<<<<< 83ae36f5f9908381fe0d957ce0daa4f037d38e6f
			for i := 1; i < btp.params.MaxLevel()+1; i++ {
=======
		for i := 0; i < levelQ+1; i++ {
			tmp = ring.BRedAdd(coeff, Q[i], bredparamsQ[i])
			ks.BuffDecompQP[0].Q.Coeffs[i][j] = tmp*pos + (Q[i]-tmp)*neg
>>>>>>> rebased onto btp_eprint
>>>>>>> rebased onto btp_eprint

	} else {

		for j := 0; j < btp.params.N(); j++ {

<<<<<<< btp_eprint
			coeff = ct.Value[1].Coeffs[0][j]
			pos, neg = 1, 0
			if coeff >= (q >> 1) {
				coeff = q - coeff
				pos, neg = 0, 1
			}

			for i := 1; i < levelQ+1; i++ {
				tmp = ring.BRedAdd(coeff, Q[i], bredparamsQ[i])
				ct.Value[1].Coeffs[i][j] = tmp*pos + (Q[i]-tmp)*neg
=======
<<<<<<< 83ae36f5f9908381fe0d957ce0daa4f037d38e6f
				if coeff > (Q >> 1) {
					ct.Value[u].Coeffs[i][j] = qi - ring.BRedAdd(Q-coeff, qi, bredparams[i])
				} else {
					ct.Value[u].Coeffs[i][j] = ring.BRedAdd(coeff, qi, bredparams[i])
				}
>>>>>>> rebased onto btp_eprint
			}
		}

		ringQ.NTTLvl(levelQ, ct.Value[0], ct.Value[0])
		ringQ.NTTLvl(levelQ, ct.Value[1], ct.Value[1])
	}

=======
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

>>>>>>> rebased onto btp_eprint
	return ct
}
