package ckks

import (
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
	"log"
	"math"
	"time"
)

// Bootstrapp re-encrypt a ciphertext at lvl Q0 to a ciphertext at MaxLevel-k where k is the depth of the bootstrapping circuit.
func (bootcontext *BootContext) Bootstrapp(ct *Ciphertext) *Ciphertext {
	var t time.Time
	var ct0, ct1 *Ciphertext

	for ct.Level() != 0 {
		bootcontext.evaluator.DropLevel(ct, 1)
	}

	// TODO : better management of the initial scale
	bootcontext.evaluator.ScaleUp(ct, math.Round(bootcontext.sinScale/ct.Scale()), ct)

	// ModUp ct_{Q_0} -> ct_{Q_L}

	t = time.Now()
	ct = bootcontext.modUp(ct)
	log.Println("After ModUp  :", time.Now().Sub(t), ct.Level(), ct.Scale())

	//SubSum X -> (N/dslots) * Y^dslots
	t = time.Now()
	ct = bootcontext.subSum(ct)
	log.Println("After SubSum :", time.Now().Sub(t), ct.Level(), ct.Scale())
	// Part 1 : Coeffs to slots

	t = time.Now()
	ct0, ct1 = bootcontext.coeffsToSlots(ct)
	log.Println("After CtS    :", time.Now().Sub(t), ct0.Level(), ct0.Scale())

	// Part 2 : SineEval
	t = time.Now()
	ct0, ct1 = bootcontext.evaluateSine(ct0, ct1)
	log.Println("After Sine   :", time.Now().Sub(t), ct0.Level(), ct0.Scale())

	// Part 3 : Slots to coeffs
	t = time.Now()
	ct0 = bootcontext.slotsToCoeffs(ct0, ct1)
	log.Println("After StC    :", time.Now().Sub(t), ct0.Level(), ct0.Scale())
	return ct0
}

func (bootcontext *BootContext) subSum(ct *Ciphertext) *Ciphertext {

	for i := bootcontext.LogSlots; i < bootcontext.LogN-1; i++ {

		bootcontext.evaluator.RotateColumns(ct, 1<<i, bootcontext.rotkeys, bootcontext.ctxpool[0])

		bootcontext.evaluator.Add(ct, bootcontext.ctxpool[0], ct)
	}

	return ct
}

func (bootcontext *BootContext) modUp(ct *Ciphertext) *Ciphertext {

	contextQ := bootcontext.evaluator.(*evaluator).ckksContext.contextQ

	ct.InvNTT(contextQ, ct.Element())

	// Extend the ciphertext with zero polynomials.
	for u := range ct.Value() {
		ct.Value()[u].Coeffs = append(ct.Value()[u].Coeffs, make([][]uint64, bootcontext.MaxLevel)...)
		for i := uint64(1); i < bootcontext.MaxLevel+1; i++ {
			ct.Value()[u].Coeffs[i] = make([]uint64, bootcontext.n)
		}
	}

	//Centers the values around Q0 and extends the basis from Q0 to QL
	Q := contextQ.Modulus[0]
	bredparams := contextQ.GetBredParams()

	var coeff, qi uint64
	for u := range ct.Value() {

		for j := uint64(0); j < bootcontext.n; j++ {

			coeff = ct.Value()[u].Coeffs[0][j]

			for i := uint64(1); i < bootcontext.MaxLevel+1; i++ {

				qi = contextQ.Modulus[i]

				if coeff > (Q >> 1) {
					ct.Value()[u].Coeffs[i][j] = qi - ring.BRedAdd(Q-coeff, qi, bredparams[i])
				} else {
					ct.Value()[u].Coeffs[i][j] = ring.BRedAdd(coeff, qi, bredparams[i])
				}
			}
		}
	}

	ct.NTT(contextQ, ct.Element())

	return ct
}

func (bootcontext *BootContext) coeffsToSlots(vec *Ciphertext) (ct0, ct1 *Ciphertext) {

	evaluator := bootcontext.evaluator

	var zV, zVconj *Ciphertext

	zV = bootcontext.dft(vec, bootcontext.pDFTInv, true)

	// Extraction of real and imaginary parts.
	zVconj = evaluator.ConjugateNew(zV, bootcontext.rotkeys)

	// The real part is stored in ct0
	ct0 = evaluator.AddNew(zV, zVconj)

	// The imaginary part is stored in ct1
	ct1 = evaluator.SubNew(zV, zVconj)

	evaluator.DivByi(ct1, ct1)

	// If repacking, then ct0 and ct1 right n/2 slots are zero.
	if bootcontext.repack {

		// The imaginary part is put in the right n/2 slots of ct0.
		evaluator.RotateColumns(ct1, bootcontext.slots, bootcontext.rotkeys, ct1)

		evaluator.Add(ct0, ct1, ct0)

		return ct0, nil
	}

	return ct0, ct1
}

func (bootcontext *BootContext) slotsToCoeffs(ct0, ct1 *Ciphertext) (ct *Ciphertext) {

	// If full packing, the repacking can be done directly using ct0 and ct1.
	if !bootcontext.repack {

		bootcontext.evaluator.MultByi(ct1, ct1)

		bootcontext.evaluator.Add(ct0, ct1, ct0)
	}

	return bootcontext.dft(ct0, bootcontext.pDFT, false)
}

func (bootcontext *BootContext) dft(vec *Ciphertext, plainVectors []*dftvectors, forward bool) *Ciphertext {

	evaluator := bootcontext.evaluator.(*evaluator)

	// Sequencially multiplies w with the provided dft matrices.
	for _, plainVector := range plainVectors {
		vec = bootcontext.multiplyByDiagMatrice(vec, plainVector)
		evaluator.Rescale(vec, evaluator.ckksContext.scale, vec)
	}

	return vec
}

func (bootcontext *BootContext) multiplyByDiagMatrice(vec *Ciphertext, plainVectors *dftvectors) (res *Ciphertext) {

	eval := bootcontext.evaluator.(*evaluator)
	contextQ := eval.ckksContext.contextQ
	contextP := eval.ckksContext.contextP

	levelQ := vec.Level()
	levelP := uint64(len(contextP.Modulus) - 1)

	var N1 uint64

	res = NewCiphertext(&bootcontext.Parameters, 1, vec.Level(), vec.Scale())

	// N1*N2 = N
	N1 = plainVectors.N1

	// Computes the rotations indexes of the non-zero rows of the diagonalized DFT matrix for the baby-step giang-step algorithm
	index := make(map[uint64][]uint64)
	rotations := []uint64{}
	for key := range plainVectors.Vec {

		idx1 := key / N1
		idx2 := key & (N1 - 1)

		if index[idx1] == nil {
			index[idx1] = []uint64{idx2}
		} else {
			index[idx1] = append(index[idx1], idx2)
		}

		if !utils.IsInSliceUint64(idx2, rotations) {
			rotations = append(rotations, idx2)
		}
	}

	// Pre-rotates ciphertext for the baby-step giant-step algorithm, does not divide by P yet
	vecRotQ, vecRotP := eval.RotateHoistedNoModDown(vec, rotations, bootcontext.rotkeys)

	tmpQ0 := bootcontext.poolQ[0]
	tmpQ1 := bootcontext.poolQ[1]
	tmpQ2 := bootcontext.poolQ[2]
	tmpQ3 := bootcontext.poolQ[3]

	tmpP0 := bootcontext.poolP[0]
	tmpP1 := bootcontext.poolP[1]
	tmpP2 := bootcontext.poolP[2]
	tmpP3 := bootcontext.poolP[3]

	tmpResQ0 := bootcontext.poolQ[4]
	tmpResQ1 := bootcontext.poolQ[5]

	pool2Q := eval.poolQ[1]
	pool3Q := eval.poolQ[2]
	pool2P := eval.poolP[1]
	pool3P := eval.poolP[2]

	N1Rot := 0
	N2Rot := 0

	// OUTER LOOP
	cnt0 := 0
	for j := range index {

		if j != 0 {

			// INNER LOOP
			state := false
			cnt1 := 0
			for _, i := range index[j] {

				if i == 0 {
					state = true
				} else {

					N1Rot++

					plaintextQ := plainVectors.Vec[N1*j+uint64(i)][0]
					plaintextP := plainVectors.Vec[N1*j+uint64(i)][1]

					if cnt1 == 0 {
						contextQ.MulCoeffsMontgomeryLvl(levelQ, plaintextQ, vecRotQ[i].value[0], tmpQ2) // phi(d0_Q) * plaintext
						contextQ.MulCoeffsMontgomeryLvl(levelQ, plaintextQ, vecRotQ[i].value[1], tmpQ3) // phi(d1_Q) * plaintext
						contextP.MulCoeffsMontgomery(plaintextP, vecRotP[i].value[0], tmpP2)            // phi(d0_P) * plaintext
						contextP.MulCoeffsMontgomery(plaintextP, vecRotP[i].value[1], tmpP3)            // phi(d1_P) * plaintext
					} else {
						contextQ.MulCoeffsMontgomeryAndAddLvl(levelQ, plaintextQ, vecRotQ[i].value[0], tmpQ2) // phi(d0_Q) * plaintext
						contextQ.MulCoeffsMontgomeryAndAddLvl(levelQ, plaintextQ, vecRotQ[i].value[1], tmpQ3) // phi(d1_Q) * plaintext
						contextP.MulCoeffsMontgomeryAndAdd(plaintextP, vecRotP[i].value[0], tmpP2)            // phi(d0_P) * plaintext
						contextP.MulCoeffsMontgomeryAndAdd(plaintextP, vecRotP[i].value[1], tmpP3)            // phi(d1_P) * plaintext
					}

					cnt1++
				}
			}

			// Hoisting of the ModDown of sum(sum(phi(d0) * plaintext)) and sum(sum(phi(d1) * plaintext))
			eval.baseconverter.ModDownSplitedNTTPQ(levelQ, tmpQ2, tmpP2, tmpResQ0) // sum(phi(d0) * plaintext)/P
			eval.baseconverter.ModDownSplitedNTTPQ(levelQ, tmpQ3, tmpP3, tmpResQ1) // sum(phi(d1) * plaintext)/P

			// Inner loop rotation of c0 (in bas Q only)
			for _, i := range index[j] {
				if i != 0 {
					ring.PermuteNTTWithIndexLvl(levelQ, vec.value[0], bootcontext.rotkeys.permuteNTTLeftIndex[i], tmpQ2) // phi(c0)
					contextQ.MulCoeffsMontgomeryAndAddLvl(levelQ, plainVectors.Vec[N1*j+uint64(i)][0], tmpQ2, tmpResQ0)  // phi(c0) * plaintext mod Q
				}
			}

			// If a rotation by zero happens in the inner loop
			if state {
				N1Rot++
				contextQ.MulCoeffsMontgomeryAndAddLvl(levelQ, plainVectors.Vec[N1*j][0], vec.value[0], tmpResQ0) // c0 * plaintext + sum(phi(d0) * plaintext)/P + phi(c0) * plaintext mod Q
				contextQ.MulCoeffsMontgomeryAndAddLvl(levelQ, plainVectors.Vec[N1*j][0], vec.value[1], tmpResQ1) // c1 * plaintext + sum(phi(d1) * plaintext)/P + phi(c1) * plaintext mod Q
			}

			// Outer loop rotations
			ring.PermuteNTTWithIndexLvl(levelQ, tmpResQ0, bootcontext.rotkeys.permuteNTTLeftIndex[N1*j], tmpQ2) // phi(tmpRes_0)
			contextQ.AddLvl(levelQ, res.value[0], tmpQ2, res.value[0])                                          // res += phi(tmpRes)

			eval.switchKeysInPlaceNoModDown(levelQ, tmpResQ1, bootcontext.rotkeys.evakeyRotColLeft[N1*j], pool2Q, pool2P, pool3Q, pool3P) // Switchkey(phi(tmpRes_1)) = (d0, d1) in base QP

			rot := bootcontext.rotkeys.permuteNTTLeftIndex[N1*j]

			N2Rot++

			if cnt0 == 0 {
				ring.PermuteNTTWithIndexLvl(levelQ, pool2Q, rot, tmpQ0) // sum(phi(d0_Q))
				ring.PermuteNTTWithIndexLvl(levelQ, pool3Q, rot, tmpQ1) // sum(phi(d1_Q))
				ring.PermuteNTTWithIndexLvl(levelP, pool2P, rot, tmpP0) // sum(phi(d0_P))
				ring.PermuteNTTWithIndexLvl(levelP, pool3P, rot, tmpP1) // sum(phi(d1_P))
			} else {
				ring.PermuteNTTWithIndexAndAddNoModLvl(levelQ, pool2Q, rot, tmpQ0) // sum(phi(d0_Q))
				ring.PermuteNTTWithIndexAndAddNoModLvl(levelQ, pool3Q, rot, tmpQ1) // sum(phi(d1_Q))
				ring.PermuteNTTWithIndexAndAddNoModLvl(levelP, pool2P, rot, tmpP0) // sum(phi(d0_P))
				ring.PermuteNTTWithIndexAndAddNoModLvl(levelP, pool3P, rot, tmpP1) // sum(phi(d1_P))
			}

			if cnt0 == 7 {
				contextQ.ReduceLvl(levelQ, tmpQ0, tmpQ0)
				contextQ.ReduceLvl(levelQ, tmpQ1, tmpQ1)
				contextP.Reduce(tmpP0, tmpP0)
				contextP.Reduce(tmpP1, tmpP1)
			}

			cnt0++
		}
	}

	if cnt0 != 7 {
		contextQ.ReduceLvl(levelQ, tmpQ0, tmpQ0)
		contextQ.ReduceLvl(levelQ, tmpQ1, tmpQ1)
		contextP.Reduce(tmpP0, tmpP0)
		contextP.Reduce(tmpP1, tmpP1)
	}

	eval.baseconverter.ModDownSplitedNTTPQ(levelQ, tmpQ0, tmpP0, tmpQ0) // sum(phi(d0_QP))/P
	eval.baseconverter.ModDownSplitedNTTPQ(levelQ, tmpQ1, tmpP1, tmpQ1) // sum(phi(d1_QP))/P

	contextQ.AddLvl(levelQ, res.value[0], tmpQ0, res.value[0]) // res += sum(phi(d0_QP))/P
	contextQ.AddLvl(levelQ, res.value[1], tmpQ1, res.value[1]) // res += sum(phi(d1_QP))/P

	// INNER LOOP ROTATION NON ZERO WITH OUTER LOOP ROTATION ZERO
	state := false
	cnt0 = 0
	for _, i := range index[0] {

		if i == 0 {
			state = true
		} else {

			plaintextQ := plainVectors.Vec[uint64(i)][0]
			plaintextP := plainVectors.Vec[uint64(i)][1]

			N1Rot++

			// keyswitch(c1_Q) = (d0_QP, d1_QP)

			if cnt0 == 0 {
				contextQ.MulCoeffsMontgomeryLvl(levelQ, plaintextQ, vecRotQ[i].value[0], tmpQ2) // phi(d0_Q) * plaintext
				contextQ.MulCoeffsMontgomeryLvl(levelQ, plaintextQ, vecRotQ[i].value[1], tmpQ3) // phi(d1_Q) * plaintext
				contextP.MulCoeffsMontgomery(plaintextP, vecRotP[i].value[0], tmpP2)            // phi(d0_P) * plaintext
				contextP.MulCoeffsMontgomery(plaintextP, vecRotP[i].value[1], tmpP3)            // phi(d1_P) * plaintext
			} else {
				contextQ.MulCoeffsMontgomeryAndAddLvl(levelQ, plaintextQ, vecRotQ[i].value[0], tmpQ2) // phi(d0_Q) * plaintext
				contextQ.MulCoeffsMontgomeryAndAddLvl(levelQ, plaintextQ, vecRotQ[i].value[1], tmpQ3) // phi(d1_Q) * plaintext
				contextP.MulCoeffsMontgomeryAndAdd(plaintextP, vecRotP[i].value[0], tmpP2)            // phi(d0_P) * plaintext
				contextP.MulCoeffsMontgomeryAndAdd(plaintextP, vecRotP[i].value[1], tmpP3)            // phi(d1_P) * plaintext
			}
			cnt0++
		}
	}

	eval.baseconverter.ModDownSplitedNTTPQ(levelQ, tmpQ2, tmpP2, tmpResQ0) // sum(phi(d0_QP) * plaintext)/P
	eval.baseconverter.ModDownSplitedNTTPQ(levelQ, tmpQ3, tmpP3, tmpResQ1) // sum(phi(d1_QP) * plaintext)/P

	contextQ.AddLvl(levelQ, res.value[0], tmpResQ0, res.value[0]) // res += sum(phi(d0_QP) * plaintext)/P
	contextQ.AddLvl(levelQ, res.value[1], tmpResQ1, res.value[1]) // res += sum(phi(d1_QP) * plaintext)/P

	// INNER LOOP ROTATION ZERO WITH OUTER LOOP ROTATION ZERO
	for _, i := range index[0] { // sum(phi(d0_Q) * plaintext)/P + phi(c0_Q) * plaintext
		if i != 0 {
			ring.PermuteNTTWithIndexLvl(levelQ, vec.value[0], bootcontext.rotkeys.permuteNTTLeftIndex[i], tmpQ2)
			contextQ.MulCoeffsMontgomeryAndAddLvl(levelQ, plainVectors.Vec[uint64(i)][0], tmpQ2, res.value[0]) // res += phi(c0_Q) * plaintext
		}
	}

	if state { // If a rotation by zero should have happened in the inner loop
		N1Rot++
		contextQ.MulCoeffsMontgomeryAndAddLvl(levelQ, plainVectors.Vec[0][0], vec.value[0], res.value[0]) // res += c0_Q * plaintext
		contextQ.MulCoeffsMontgomeryAndAddLvl(levelQ, plainVectors.Vec[0][0], vec.value[1], res.value[1]) // res += c1_Q * plaintext
	}

	res.SetScale(plainVectors.Scale * vec.Scale())

	//log.Println(N1Rot, N2Rot)

	return
}

// Sine Evaluation ct0 = Q/(2pi) * sin((2pi/Q) * ct0)
func (bootcontext *BootContext) evaluateSine(ct0, ct1 *Ciphertext) (*Ciphertext, *Ciphertext) {

	evaluator := bootcontext.evaluator.(*evaluator)

	ct0.MulScale(1024)
	evaluator.ckksContext.scale = ct0.Scale() // Reference scale is changed to the new ciphertext's scale.

	// pre-computes the target scale for the output of the polynomial evaluation such that
	// the output scale after the polynomial evaluation followed by the double angle formula
	// does not change the scale of the ciphertext.
	for i := uint64(0); i < bootcontext.SinRescal; i++ {
		evaluator.ckksContext.scale *= float64(evaluator.params.Qi[bootcontext.StCLevel[0]+i+1])
		evaluator.ckksContext.scale = math.Sqrt(evaluator.ckksContext.scale)
	}

	ct0 = bootcontext.evaluateCheby(ct0)

	ct0.DivScale(1024 * bootcontext.sinScale / bootcontext.Scale)

	if ct1 != nil {
		ct1.MulScale(1024)
		ct1 = bootcontext.evaluateCheby(ct1)
		ct1.DivScale(1024 * bootcontext.sinScale / bootcontext.Scale)
	}

	// Reference scale is changed back to the current ciphertext's scale.
	evaluator.ckksContext.scale = ct0.Scale()

	return ct0, ct1
}

func (bootcontext *BootContext) evaluateCheby(ct *Ciphertext) (res *Ciphertext) {

	eval := bootcontext.evaluator.(*evaluator)

	C := make(map[uint64]*Ciphertext)
	C[1] = ct.CopyNew().Ciphertext()

	cheby := bootcontext.chebycoeffs
	if bootcontext.SinType == Sin {
		eval.AddConst(C[1], (-cheby.a-cheby.b)/(cheby.b-cheby.a), C[1])
	} else {
		sc_fac := complex(float64(int(1<<bootcontext.SinRescal)), 0)
		eval.AddConst(C[1], (-cheby.a-cheby.b-(2*0.25/sc_fac))/(cheby.b-cheby.a), C[1])
	}

	res = eval.evalCheby(cheby, C, bootcontext.relinkey)

	sqrt2pi := math.Pow(0.15915494309189535, 1.0/float64(int(1<<bootcontext.SinRescal)))

	for i := uint64(0); i < bootcontext.SinRescal; i++ {
		sqrt2pi *= sqrt2pi
		eval.MulRelin(res, res, bootcontext.relinkey, res)
		eval.Add(res, res, res)
		eval.AddConst(res, -sqrt2pi, res)
		eval.Rescale(res, eval.ckksContext.scale, res)
	}

	return
}
