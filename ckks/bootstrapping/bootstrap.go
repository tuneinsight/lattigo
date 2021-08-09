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
func (btp *Bootstrapper) Bootstrapp(ct *ckks.Ciphertext) *ckks.Ciphertext {

	bootstrappingScale := math.Exp2(math.Round(math.Log2(btp.params.QiFloat64(0) / btp.evalModPoly.MessageRatio)))

	var ct0, ct1 *ckks.Ciphertext

	// Drops the level to 1
	for ct.Level() > 1 {
		btp.DropLevel(ct, 1)
	}

	// Brings the ciphertext scale to Q0/MessageRatio
	if ct.Level() == 1 {

		// If one level is available, then uses it to match the scale
		btp.SetScale(ct, bootstrappingScale)

		// Then drops to level 0
		for ct.Level() != 0 {
			btp.DropLevel(ct, 1)
		}

	} else {

		// Does an integer constant mult by round((Q0/Delta_m)/ctscle)
		if bootstrappingScale < ct.Scale {
			panic("ciphetext scale > q/||m||)")
		}

		btp.ScaleUp(ct, math.Round(bootstrappingScale/ct.Scale), ct)
	}

	// Step 1 : Extend the basis from q to Q
	ct = btp.modUpFromQ0(ct)

	// Brings the ciphertext scale to sineQi/(Q0/scale) if Q0 < sineQi
	// Does it after modUp to avoid plaintext overflow
	// Reduces the additive error of the next steps
	btp.ScaleUp(ct, math.Round((btp.evalModPoly.ScalingFactor/btp.evalModPoly.MessageRatio)/ct.Scale), ct)

	//SubSum X -> (N/dslots) * Y^dslots
	ct = btp.Trace(ct, btp.params.LogSlots(), btp.params.LogN()-1)

	// Step 2 : CoeffsToSlots (Homomorphic encoding)
	ct0, ct1 = btp.CoeffsToSlots(ct, btp.ctsMatrices)

	// Step 3 : EvalMod (Homomorphic modular reduction)
	// ct0 = Ecd(real)
	// ct1 = Ecd(imag)
	// If n < N/2 then ct0 = Ecd(real|imag)
	ct0 = btp.EvalMod(ct0, btp.evalModPoly)
	ct0.Scale = btp.params.Scale()

	if ct1 != nil {
		ct1 = btp.EvalMod(ct1, btp.evalModPoly)
		ct1.Scale = btp.params.Scale()
	}

	// Step 4 : SlotsToCoeffs (Homomorphic decoding)
	ct0 = btp.SlotsToCoeffs(ct0, ct1, btp.stcMatrices)

	return ct0
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
