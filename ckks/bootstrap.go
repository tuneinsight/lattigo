package ckks

import (
	"github.com/ldsec/lattigo/v2/ring"
	//"log"
	"math"
	//"time"
	"fmt"
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

		// and does an integer constant mult by round((Q0/2^{10})/ctscle)
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

		btp.evaluator.Rotate(ct, 1<<i, btp.rotkeys, btp.ctxpool)

		btp.evaluator.Add(ct, btp.ctxpool, ct)
	}

	return ct
}

func (btp *Bootstrapper) modUp(ct *Ciphertext) *Ciphertext {

	ringQ := btp.evaluator.(*evaluator).ringQ

	ct.InvNTT(ringQ, ct.El())

	// Extend the ciphertext with zero polynomials.
	for u := range ct.Value() {
		ct.Value()[u].Coeffs = append(ct.Value()[u].Coeffs, make([][]uint64, btp.params.MaxLevel())...)
		for i := uint64(1); i < btp.params.MaxLevel()+1; i++ {
			ct.Value()[u].Coeffs[i] = make([]uint64, btp.params.N())
		}
	}

	//Centers the values around Q0 and extends the basis from Q0 to QL
	Q := ringQ.Modulus[0]
	bredparams := ringQ.BredParams

	var coeff, qi uint64
	for u := range ct.Value() {

		for j := uint64(0); j < btp.params.N(); j++ {

			coeff = ct.Value()[u].Coeffs[0][j]

			for i := uint64(1); i < btp.params.MaxLevel()+1; i++ {

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

	eval := btp.evaluator

	var zV, zVconj *Ciphertext

	zV = btp.dft(vec, btp.pDFTInv, true)

	// Extraction of real and imaginary parts.
	zVconj = eval.ConjugateNew(zV, btp.rotkeys)

	// The real part is stored in ct0
	ct0 = eval.AddNew(zV, zVconj)

	// The imaginary part is stored in ct1
	ct1 = eval.SubNew(zV, zVconj)

	// Switches imaginary and real and negates
	eval.DivByi(ct1, ct1)

	// If repacking, then ct0 and ct1 right n/2 slots are zero.
	if btp.repack {
		// The imaginary part is put in the right n/2 slots of ct0.
		eval.Rotate(ct1, btp.params.Slots(), btp.rotkeys, ct1)
		eval.Add(ct0, ct1, ct0)
		return ct0, nil
	}

	zV = nil
	zVconj = nil

	return ct0, ct1
}

func (btp *Bootstrapper) slotsToCoeffs(ct0, ct1 *Ciphertext) (ct *Ciphertext) {

	// If full packing, the repacking can be done directly using ct0 and ct1.
	if !btp.repack {
		btp.evaluator.MultByi(ct1, ct1)
		btp.evaluator.Add(ct0, ct1, ct0)
	}

	ct1 = nil

	return btp.dft(ct0, btp.pDFT, false)
}

func (btp *Bootstrapper) dft(vec *Ciphertext, plainVectors []*PtDiagMatrix, forward bool) *Ciphertext {

	eval := btp.evaluator.(*evaluator)
	// Sequencially multiplies vec with the provided dft matrices.
	for _, plainVector := range plainVectors {
		vec = eval.LinearTransform(vec, plainVector, btp.rotkeys)[0]
		eval.Rescale(vec, eval.scale, vec)
	}

	return vec
}

// Sine Evaluation ct0 = Q/(2pi) * sin((2pi/Q) * ct0)
func (btp *Bootstrapper) evaluateSine(ct0, ct1 *Ciphertext) (*Ciphertext, *Ciphertext) {

	fmt.Println(ct0.Level(), ct0.Scale())

	evaluator := btp.evaluator.(*evaluator)

	ct0.MulScale(btp.MessageRatio)
	evaluator.scale = btp.sinescale // Reference scale is changed to the Qi used for the SineEval (which is also close to the new ciphetext scale)

	ct0 = btp.evaluateCheby(ct0)

	ct0.DivScale(btp.MessageRatio * btp.postscale / btp.params.scale)

	if ct1 != nil {
		ct1.MulScale(btp.MessageRatio)
		ct1 = btp.evaluateCheby(ct1)
		ct1.DivScale(btp.MessageRatio * btp.postscale / btp.params.scale)
	}

	// Reference scale is changed back to the current ciphertext's scale.
	evaluator.scale = ct0.Scale()

	return ct0, ct1
}

func (btp *Bootstrapper) evaluateCheby(ct *Ciphertext) *Ciphertext {

	fmt.Println(ct.Level(), ct.Scale())

	eval := btp.evaluator.(*evaluator)

	cheby := btp.chebycoeffs

	targetScale := btp.sinescale

	// Compute the scales that the ciphertext should have before the double angle
	// formula such that after it it has the scale it had before the polynomial 
	// evaluation
	for i := uint64(0); i < btp.SinRescal; i++ {
		targetScale = math.Sqrt(targetScale * float64(btp.SineEvalModuli.Qi[i]))
	}

	// If arcsine is evaluated, then the multiplication by 1/2pi is merged with it
	var sqrt2pi float64
	if btp.ArcSineDeg > 0 {
		sqrt2pi = math.Pow(1, 1.0/float64(int(1<<btp.SinRescal)))
	} else {
		sqrt2pi = math.Pow(0.15915494309189535, 1.0/float64(int(1<<btp.SinRescal)))
	}

	// Division by 1/2^r and change of variable for the Chebysehev evaluation
	if btp.SinType == Cos1 || btp.SinType == Cos2 {
		scfac := complex(float64(int(1<<btp.SinRescal)), 0)
		eval.AddConst(ct, -0.5/(scfac*(cheby.b-cheby.a)), ct)
	}

	// Chebyshev evaluation
	ct, _ = eval.EvaluateCheby(ct, cheby, targetScale, btp.relinkey)

	// Double angle
	for i := uint64(0); i < btp.SinRescal; i++ {
		sqrt2pi *= sqrt2pi
		eval.MulRelin(ct, ct, btp.relinkey, ct)
		eval.Add(ct, ct, ct)
		eval.AddConst(ct, -sqrt2pi, ct)
		if err := eval.Rescale(ct, eval.scale, ct); err != nil {
			panic(err)
		}
	}

	// ArcSine
	if btp.ArcSineDeg > 0 {
		poly := NewPoly([]complex128{0, 0.15915494309189535, 0, 1.0 / 6.0 * 0.15915494309189535, 0, 3.0 / 40.0 * 0.15915494309189535, 0, 5.0 / 112.0 * 0.15915494309189535}[:btp.ArcSineDeg+1])
		ct, _ = eval.EvaluatePoly(ct, poly, ct.Scale(), btp.relinkey)
	}

	return ct
}
