package ckks

import (
	"github.com/ldsec/lattigo/v2/ring"
	"math"
)

// Bootstrapp re-encrypt a ciphertext at lvl Q0 to a ciphertext at MaxLevel-k where k is the depth of the bootstrapping circuit.
// If the input ciphertext level is zero, the input scale must be an exact power of two smaller or equal to round(Q0/2^{10}).
// If the input ciphertext is at level one or more, the input scale does not need to be an exact power of two as one level
// can be used to do a scale matching.
func (btp *Bootstrapper) Bootstrapp(ct *Ciphertext) *Ciphertext {

	//var t time.Time
	var ct0, ct1 *Ciphertext

	// Drops the level to 1
	for ct.Level() > 1 {
		btp.evaluator.DropLevel(ct, 1)
	}

	// Brings the ciphertext scale to Q0/2^{10}
	if ct.Level() == 1 {

		// if one level is available, then uses it to match the scale
		btp.evaluator.SetScale(ct, btp.prescale)

		// then drops to level 0
		for ct.Level() != 0 {
			btp.evaluator.DropLevel(ct, 1)
		}

	} else {

		// else drop to level 0
		for ct.Level() != 0 {
			btp.evaluator.DropLevel(ct, 1)
		}

		// and does an integer constant mult by round((Q0/Delta_m)/ctscle)

		if btp.prescale < ct.Scale() {
			panic("ciphetext scale > Q[0]/(Q[0]/Delta_m)")
		}
		btp.evaluator.ScaleUp(ct, math.Round(btp.prescale/ct.Scale()), ct)
	}

	// ModUp ct_{Q_0} -> ct_{Q_L}
	//t = time.Now()
	ct = btp.modUp(ct)
	//log.Println("After ModUp  :", time.Now().Sub(t), ct.Level(), ct.Scale())

	// Brings the ciphertext scale to sineQi/(Q0/scale) if its under
	btp.evaluator.ScaleUp(ct, math.Round(btp.postscale/ct.Scale()), ct)

	//SubSum X -> (N/dslots) * Y^dslots
	//t = time.Now()
	ct = btp.subSum(ct)
	//log.Println("After SubSum :", time.Now().Sub(t), ct.Level(), ct.Scale())
	// Part 1 : Coeffs to slots

	//t = time.Now()
	ct0, ct1 = btp.coeffsToSlots(ct)
	//log.Println("After CtS    :", time.Now().Sub(t), ct0.Level(), ct0.Scale())

	// Part 2 : SineEval
	//t = time.Now()
	ct0, ct1 = btp.evaluateSine(ct0, ct1)
	//log.Println("After Sine   :", time.Now().Sub(t), ct0.Level(), ct0.Scale())

	// Part 3 : Slots to coeffs
	//t = time.Now()
	ct0 = btp.slotsToCoeffs(ct0, ct1)

	ct0.SetScale(math.Exp2(math.Round(math.Log2(ct0.Scale())))) // rounds to the nearest power of two
	//log.Println("After StC    :", time.Now().Sub(t), ct0.Level(), ct0.Scale())
	return ct0
}

func (btp *Bootstrapper) subSum(ct *Ciphertext) *Ciphertext {

	for i := btp.params.logSlots; i < btp.params.MaxLogSlots(); i++ {

		btp.evaluator.Rotate(ct, 1<<i, btp.evaluator.ctxpool)

		btp.evaluator.Add(ct, btp.evaluator.ctxpool, ct)
	}

	return ct
}

func (btp *Bootstrapper) modUp(ct *Ciphertext) *Ciphertext {

	ringQ := btp.evaluator.ringQ

	ct.InvNTT(ringQ, ct.El())

	// Extend the ciphertext with zero polynomials.
	for u := range ct.Value() {
		ct.Value()[u].Coeffs = append(ct.Value()[u].Coeffs, make([][]uint64, btp.params.MaxLevel())...)
		for i := 1; i < btp.params.MaxLevel()+1; i++ {
			ct.Value()[u].Coeffs[i] = make([]uint64, btp.params.N())
		}
	}

	//Centers the values around Q0 and extends the basis from Q0 to QL
	Q := ringQ.Modulus[0]
	bredparams := ringQ.BredParams

	var coeff, qi uint64
	for u := range ct.Value() {

		for j := 0; j < btp.params.N(); j++ {

			coeff = ct.Value()[u].Coeffs[0][j]

			for i := 1; i < btp.params.MaxLevel()+1; i++ {

				qi = ringQ.Modulus[i]

				if coeff > (Q >> 1) {
					ct.Value()[u].Coeffs[i][j] = qi - ring.BRedAdd(Q-coeff, qi, bredparams[i])
				} else {
					ct.Value()[u].Coeffs[i][j] = ring.BRedAdd(coeff, qi, bredparams[i])
				}
			}
		}
	}

	ct.NTT(ringQ, ct.El())

	return ct
}

func (btp *Bootstrapper) coeffsToSlots(vec *Ciphertext) (ct0, ct1 *Ciphertext) {

	var zV, zVconj *Ciphertext

	zV = btp.dft(vec, btp.pDFTInv, true)

	zVconj = btp.ConjugateNew(zV)

	// The real part is stored in ct0
	ct0 = btp.AddNew(zV, zVconj)

	// The imaginary part is stored in ct1
	ct1 = btp.SubNew(zV, zVconj)

	btp.DivByi(ct1, ct1)

	// If repacking, then ct0 and ct1 right n/2 slots are zero.
	if btp.repack {
		btp.Rotate(ct1, int(btp.params.Slots()), ct1)
		btp.Add(ct0, ct1, ct0)
		return ct0, nil
	}

	zV = nil
	zVconj = nil

	return ct0, ct1
}

func (btp *Bootstrapper) slotsToCoeffs(ct0, ct1 *Ciphertext) (ct *Ciphertext) {

	// If full packing, the repacking can be done directly using ct0 and ct1.
	if !btp.repack {
		btp.MultByi(ct1, ct1)
		btp.Add(ct0, ct1, ct0)
	}

	ct1 = nil

	return btp.dft(ct0, btp.pDFT, false)
}

func (btp *Bootstrapper) dft(vec *Ciphertext, plainVectors []*PtDiagMatrix, forward bool) *Ciphertext {

	// Sequentially multiplies w with the provided dft matrices.
	for _, plainVector := range plainVectors {
		vec = btp.LinearTransform(vec, plainVector)[0]
		if err := btp.Rescale(vec, btp.scale, vec); err != nil {
			panic(err)
		}
	}

	return vec
}

// Sine Evaluation ct0 = Q/(2pi) * sin((2pi/Q) * ct0)
func (btp *Bootstrapper) evaluateSine(ct0, ct1 *Ciphertext) (*Ciphertext, *Ciphertext) {

	ct0.MulScale(btp.MessageRatio)
	btp.evaluator.scale = btp.sinescale // Reference scale is changed to the Qi used for the SineEval (which is also close to the new ciphetext scale)

	ct0 = btp.evaluateCheby(ct0)

	ct0.DivScale(btp.MessageRatio * btp.postscale / btp.params.scale)

	if ct1 != nil {
		ct1.MulScale(btp.MessageRatio)
		ct1 = btp.evaluateCheby(ct1)
		ct1.DivScale(btp.MessageRatio * btp.postscale / btp.params.scale)
	}

	// Reference scale is changed back to the current ciphertext's scale.
	btp.evaluator.scale = ct0.Scale()

	return ct0, ct1
}

func (btp *Bootstrapper) evaluateCheby(ct *Ciphertext) *Ciphertext {

	var err error

	cheby := btp.sineEvalPoly

	targetScale := btp.sinescale

	// Compute the scales that the ciphertext should have before the double angle
	// formula such that after it it has the scale it had before the polynomial
	// evaluation
	for i := 0; i < btp.SinRescal; i++ {
		targetScale = math.Sqrt(targetScale * float64(btp.SineEvalModuli.Qi[i]))
	}

	// Division by 1/2^r and change of variable for the Chebysehev evaluation
	if btp.SinType == Cos1 || btp.SinType == Cos2 {
		btp.AddConst(ct, -0.5/(complex(btp.scFac, 0)*(cheby.b-cheby.a)), ct)
	}

	// Chebyshev evaluation
	if ct, err = btp.EvaluateCheby(ct, cheby, targetScale); err != nil {
		panic(err)
	}

	// Double angle
	sqrt2pi := btp.sqrt2pi
	for i := 0; i < btp.SinRescal; i++ {
		sqrt2pi *= sqrt2pi
		btp.MulRelin(ct, ct, ct)
		btp.Add(ct, ct, ct)
		btp.AddConst(ct, -sqrt2pi, ct)
		if err := btp.Rescale(ct, btp.evaluator.scale, ct); err != nil {
			panic(err)
		}
	}

	// ArcSine
	if btp.ArcSineDeg > 0 {
		if ct, err = btp.EvaluatePoly(ct, btp.arcSinePoly, ct.Scale()); err != nil {
			panic(err)
		}
	}

	return ct
}
