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
					ct.Value()[u].Coeffs[i][j] = qi - ring.BRedAdd(Q-coeff+1, qi, bredparams[i])
				} else {
					ct.Value()[u].Coeffs[i][j] = ring.BRedAdd(coeff, qi, bredparams[i])
				}
			}

			qi = contextQ.Modulus[0]

			if coeff > (Q >> 1) {
				ct.Value()[u].Coeffs[0][j] = qi - ring.BRedAdd(Q-coeff+1, qi, bredparams[0])
			} else {
				ct.Value()[u].Coeffs[0][j] = ring.BRedAdd(coeff, qi, bredparams[0])
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

	evaluator := bootcontext.evaluator

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

	// Pre-rotates ciphertext for the baby-step giant-step algorithm
	vecRot := evaluator.RotateHoisted(vec, rotations, bootcontext.rotkeys)

	var tmpVec, tmp *Ciphertext

	tmpVec = NewCiphertext(&bootcontext.Parameters, 1, bootcontext.MaxLevel, vec.Scale())
	tmp = NewCiphertext(&bootcontext.Parameters, 1, bootcontext.MaxLevel, vec.Scale())

	for j := range index {

		tmpVec.Value()[0].Zero()
		tmpVec.Value()[1].Zero()

		for _, i := range index[j] {
			evaluator.MulRelin(vecRot[uint64(i)], plainVectors.Vec[N1*j+uint64(i)], nil, tmp)
			evaluator.Add(tmpVec, tmp, tmpVec)
		}

		evaluator.RotateColumns(tmpVec, N1*j, bootcontext.rotkeys, tmp)

		evaluator.Add(res, tmp, res)
	}

	return res
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
