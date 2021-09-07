package bootstrapping

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
	"math"
)

// Bootstrapp re-encrypt a ciphertext at lvl Q0 to a ciphertext at MaxLevel-k where k is the depth of the bootstrapping circuit.
// If the input ciphertext level is zero, the input scale must be an exact power of two smaller or equal to round(Q0/2^{10}).
// If the input ciphertext is at level one or more, the input scale does not need to be an exact power of two as one level
// can be used to do a scale matching.
func (btp *Bootstrapper) Bootstrapp(ctIn *ckks.Ciphertext) (ctOut *ckks.Ciphertext) {

	ctOut = ctIn.CopyNew()

	bootstrappingScale := math.Exp2(math.Round(math.Log2(btp.params.QiFloat64(0) / btp.evalModPoly.MessageRatio())))

	// Drops the level to 1
	for ctOut.Level() > 1 {
		btp.DropLevel(ctOut, 1)
	}

	// Brings the ciphertext scale to Q0/MessageRatio
	if ctOut.Level() == 1 {

		// If one level is available, then uses it to match the scale
		btp.SetScale(ctOut, bootstrappingScale)

		// Then drops to level 0
		for ctOut.Level() != 0 {
			btp.DropLevel(ctOut, 1)
		}

	} else {

		// Does an integer constant mult by round((Q0/Delta_m)/ctscle)
		if bootstrappingScale < ctOut.Scale {
			panic("ciphetext scale > q/||m||)")
		}

		btp.ScaleUp(ctOut, math.Round(bootstrappingScale/ctOut.Scale), ctOut)
	}

	// Step 1 : Extend the basis from q to Q
	ctOut = btp.modUpFromQ0(ctOut)

	// Brings the ciphertext scale to sineQi/(Q0/scale) if Q0 < sineQi
	// Does it after modUp to avoid plaintext overflow
	// Reduces the additive error of the next steps
	btp.ScaleUp(ctOut, math.Round((btp.evalModPoly.ScalingFactor()/btp.evalModPoly.MessageRatio())/ctOut.Scale), ctOut)

	//SubSum X -> (N/dslots) * Y^dslots
	btp.Trace(ctOut, btp.params.LogSlots(), btp.params.LogN()-1, ctOut)

	// Step 2 : CoeffsToSlots (Homomorphic encoding)
	ctReal, ctImag := btp.CoeffsToSlotsNew(ctOut, btp.ctsMatrices)

	// Step 3 : EvalMod (Homomorphic modular reduction)
	// ctReal = Ecd(real)
	// ctImag = Ecd(imag)
	// If n < N/2 then ctReal = Ecd(real|imag)
	ctReal = btp.EvalModNew(ctReal, btp.evalModPoly)
	ctReal.Scale = btp.params.Scale()

	if ctImag != nil {
		ctImag = btp.EvalModNew(ctImag, btp.evalModPoly)
		ctImag.Scale = btp.params.Scale()
	}

	// Step 4 : SlotsToCoeffs (Homomorphic decoding)
	ctOut = btp.SlotsToCoeffsNew(ctReal, ctImag, btp.stcMatrices)

	return
}

func (btp *Bootstrapper) modUpFromQ0(ct *ckks.Ciphertext) *ckks.Ciphertext {

	ringQ := btp.params.RingQ()

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

	//Centers the values around Q0 and extends the basis from Q0 to QL
	Q := ringQ.Modulus[0]
	bredparams := ringQ.BredParams

	var coeff, qi uint64
	for u := range ct.Value {

		for j := 0; j < btp.params.N(); j++ {

			coeff = ct.Value[u].Coeffs[0][j]

			for i := 1; i < btp.params.MaxLevel()+1; i++ {

				qi = ringQ.Modulus[i]

				if coeff > (Q >> 1) {
					ct.Value[u].Coeffs[i][j] = qi - ring.BRedAdd(Q-coeff, qi, bredparams[i])
				} else {
					ct.Value[u].Coeffs[i][j] = ring.BRedAdd(coeff, qi, bredparams[i])
				}
			}
		}
	}

	for i := range ct.Value {
		ringQ.NTTLvl(ct.Level(), ct.Value[i], ct.Value[i])
	}

	return ct
}
