package ckks

import (
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
	"math"
	"math/bits"
)

// Bootstrapp re-encrypt a ciphertext at lvl Q0 to a ciphertext at MaxLevel-k where k is the depth of the bootstrapping circuit.
func (bootcontext *BootContext) Bootstrapp(ct *Ciphertext) *Ciphertext {

	var ct0, ct1 *Ciphertext

	for ct.Level() != 0 {
		bootcontext.evaluator.DropLevel(ct, 1)
	}

	// TODO : better management of the initial scale
	bootcontext.evaluator.ScaleUp(ct, math.Round(bootcontext.sinScale/ct.Scale()), ct)

	// ModUp ct_{Q_0} -> ct_{Q_L}
	ct = bootcontext.modUp(ct)

	//SubSum X -> (N/dslots) * Y^dslots
	ct = bootcontext.subSum(ct)

	// Part 1 : Coeffs to slots
	ct0, ct1 = bootcontext.coeffsToSlots(ct)

	// Part 2 : SineEval
	ct0, ct1 = bootcontext.evaluateSine(ct0, ct1)

	// Part 3 : Slots to coeffs
	return bootcontext.slotsToCoeffs(ct0, ct1)
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

func (bootcontext *BootContext) evaluateSine(ct0, ct1 *Ciphertext) (*Ciphertext, *Ciphertext) {

	evaluator := bootcontext.evaluator.(*evaluator)

	// Reference scale is changed to the new ciphertext's scale.
	evaluator.ckksContext.scale = float64(bootcontext.Qi[ct0.Level()-1])

	ct0.MulScale(1024)

	// Sine Evaluation ct0 = Q/(2pi) * sin((2pi/Q) * ct0)
	if bootcontext.SinType == Sin {
		ct0 = bootcontext.evaluateChebySin(ct0)
	} else if bootcontext.SinType == Cos {
		ct0 = bootcontext.evaluateChebyCos(ct0)
	} else {
		panic("bootstrapp -> evaluate sine -> invalid sineType")
	}

	ct0.SetScale(bootcontext.Scale)

	if ct1 != nil {

		ct1.MulScale(1024)

		// Sine Evaluation ct1 = Q/(2pi) * sin((2pi/Q) * ct1)
		if bootcontext.SinType == Sin {
			ct1 = bootcontext.evaluateChebySin(ct1)
		} else if bootcontext.SinType == Cos {
			ct1 = bootcontext.evaluateChebyCos(ct1)
		} else {
			panic("bootstrapp -> evaluate sine -> invalid sineType")
		}

		ct1.SetScale(bootcontext.Scale)
	}

	// Reference scale is changed back to the current ciphertext's scale.
	evaluator.ckksContext.scale = ct0.Scale()

	return ct0, ct1
}

func (bootcontext *BootContext) evaluateChebySin(ct *Ciphertext) (res *Ciphertext) {

	evaluator := bootcontext.evaluator.(*evaluator)

	// Chebyshev params
	a := bootcontext.chebycoeffs.a
	b := bootcontext.chebycoeffs.b
	degree := bootcontext.chebycoeffs.degree()

	C := make(map[uint64]*Ciphertext)
	C[1] = ct.CopyNew().Ciphertext()

	evaluator.AddConst(C[1], (-a-b)/(b-a), C[1])

	M := uint64(bits.Len64(degree - 1))
	L := bootcontext.BabySplit

	for i := uint64(2); i < (1<<L)+1; i++ {
		computePowerBasisCheby(i, C, evaluator, bootcontext.relinkey)
	}

	for i := L + 1; i < M; i++ {
		computePowerBasisCheby(1<<i, C, evaluator, bootcontext.relinkey)
	}

	res = recurseCheby(degree, L, M, bootcontext.chebycoeffs.Poly(), C, evaluator, bootcontext.relinkey)

	return
}

func (bootcontext *BootContext) evaluateChebyCos(ct *Ciphertext) (res *Ciphertext) {

	evaluator := bootcontext.evaluator.(*evaluator)

	cheby := bootcontext.chebycoeffs
	sc_fac := complex(float64(int(1<<bootcontext.SinRescal)), 0)

	// Chebyshev params
	a := cheby.a
	b := cheby.b
	degree := cheby.maxDeg - 1

	C := make(map[uint64]*Ciphertext)
	C[1] = ct.CopyNew().Ciphertext()

	evaluator.AddConst(C[1], (-a-b-(2*0.25/sc_fac))/(b-a), C[1])

	M := uint64(bits.Len64(degree - 1))
	L := bootcontext.BabySplit

	for i := uint64(2); i < (1<<L)+1; i++ {
		computePowerBasisCheby(i, C, evaluator, bootcontext.relinkey)
	}

	for i := L + 1; i < M; i++ {
		computePowerBasisCheby(1<<i, C, evaluator, bootcontext.relinkey)
	}

	res = recurseCheby(degree, L, M, cheby.Poly(), C, evaluator, bootcontext.relinkey)

	/*
		for i := uint64(0); i < bootcontext.SinRescal; i++ {
			evaluator.MulRelin(res, res, bootcontext.relinkey, res)
			evaluator.MultByConst(res, 2, res)
			evaluator.AddConst(res, -1, res)
			evaluator.Rescale(res, evaluator.ckksContext.scale, res)
		}
	*/

	if bootcontext.SinRescal == 1 {
		// r = 2*y2 - a
		a := -1.0 / 6.283185307179586

		evaluator.MulRelin(res, res, bootcontext.relinkey, res)
		evaluator.Rescale(res, evaluator.ckksContext.scale, res)
		evaluator.AddConst(res, a, res)
	}

	if bootcontext.SinRescal == 2 {

		// r = 4 * y2 * (y2 - a) + b

		a := -0.5641895835477563
		b := 1.0 / 6.283185307179586

		evaluator.MulRelin(res, res, bootcontext.relinkey, res)
		evaluator.Rescale(res, evaluator.ckksContext.scale, res)

		y := evaluator.AddConstNew(res, a)

		evaluator.MulRelin(res, y, bootcontext.relinkey, res)

		evaluator.MultByConst(res, 4, res)
		evaluator.AddConst(res, b, res)

		evaluator.Rescale(res, evaluator.ckksContext.scale, res)
	}

	if bootcontext.SinRescal == 3 {

		// r = e*(y4 * (a*y4 - b*y2 + c) - d*y2) + f

		a := 4.0
		b := -6.00900435571954
		c := 2.8209479177387813
		d := -0.42377720812375763
		e := 16.0
		f := 0.15915494309189535

		// y2 (10, 16)
		y2 := evaluator.MulRelinNew(res, res, bootcontext.relinkey)
		evaluator.Rescale(y2, evaluator.ckksContext.scale, y2)

		// tmp1 (10, 33)
		tmp1 := y2.CopyNew().Ciphertext()
		evaluator.MultByConst(tmp1, b, tmp1)

		// tmp2 (10, 33)
		tmp2 := y2.CopyNew().Ciphertext()
		evaluator.MultByConst(tmp2, d, tmp2)

		// y4 (10, 33)
		y4 := evaluator.MulRelinNew(y2, y2, bootcontext.relinkey)

		// res (10, 33)
		res = y4.CopyNew().Ciphertext()

		// y4 (9, 16)
		evaluator.Rescale(y4, evaluator.ckksContext.scale, y4)

		// res (10, 33)
		evaluator.MultByConst(res, a, res)

		// res (10, 33) + tmp1 (10, 33)
		evaluator.Add(res, tmp1, res)
		evaluator.AddConst(res, c, res)

		// res (9, 16)
		evaluator.Rescale(res, evaluator.ckksContext.scale, res)

		// res (9, 16) * y4 (9, 16) = res (9, 33)
		evaluator.MulRelin(res, y4, bootcontext.relinkey, res)

		// res (9, 33) + tmp2 (10, 33)
		evaluator.Add(res, tmp2, res)

		evaluator.MultByConst(res, e, res)
		evaluator.AddConst(res, f, res)

		// res (8, 16)
		evaluator.Rescale(res, evaluator.ckksContext.scale, res)

	}

	return
}
