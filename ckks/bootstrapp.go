package ckks

import (
	"fmt"
	"github.com/ldsec/lattigo/ckks/bettersine"
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
	"math"
	"math/bits"
	"math/cmplx"
)

// BootContext stores the parameters for the bootstrapping.
type bootParams struct {
	params *Parameters

	n    uint64
	logn uint64

	slots    uint64
	logSlots uint64

	dslots    uint64
	logdslots uint64

	encoder   Encoder
	evaluator Evaluator

	plaintextSize uint64

	ctsRescale bool
	stcRescale bool

	ctsDepth uint64
	stcDepth uint64
	sinDepth uint64

	repack    bool
	sineScale float64

	ctxpool [3]*Ciphertext
}

type bootSine struct {
	chebycoeffs *ChebyshevInterpolation
}

type bootDFT struct {
	coeffsToSlotsDiffScale complex128
	slotsToCoeffsDiffScale complex128
	pDFT                   []*dftvectors
	pDFTInv                []*dftvectors
}

type bootKeys struct {
	relinkey *EvaluationKey
	rotkeys  *RotationKeys
}

type BootContext struct {
	*bootParams
	*bootSine
	*bootDFT
	*bootKeys
}

type BootContextBetterSine struct {
	BootContext
	sc_num int
}

type dftvectors struct {
	N1  uint64
	Vec map[uint64]*Plaintext
}

func sin2pi2pi(x complex128) complex128 {
	return cmplx.Sin(6.283185307179586*x) / 6.283185307179586
}

func showcoeffs(decryptor Decryptor, encoder Encoder, slots uint64, ciphertext *Ciphertext, message string) (coeffs []complex128) {

	coeffs = encoder.Decode(decryptor.DecryptNew(ciphertext), slots)

	if slots == 2 {
		fmt.Printf(message+"%22.10f %22.10f...\n", coeffs[0], coeffs[1])
	} else {
		fmt.Printf(message+"%22.10f %22.10f %22.10f %22.10f...\n", coeffs[0], coeffs[1], coeffs[2], coeffs[3])
	}

	return coeffs
}

// NewBootContext creates a new bootcontext.
func (bootcontext *BootContext) newBootParams(bootparams *BootParams) {

	params := &bootparams.Parameters

	bootcontext.bootParams = new(bootParams)

	bootcontext.params = params.Copy()

	bootcontext.n = uint64(1 << params.LogN)
	bootcontext.slots = uint64(1 << params.LogSlots)

	if params.LogSlots < params.LogN-1 {
		bootcontext.repack = true
	}

	bootcontext.dslots = bootcontext.slots
	if params.LogSlots < params.LogN-1 {
		bootcontext.dslots <<= 1
	}

	bootcontext.sineScale = 1 << 45

	bootcontext.ctsDepth = bootparams.ctsDepth
	bootcontext.stcDepth = bootparams.stcDepth

	bootcontext.ctsRescale = bootparams.ctsRescale
	bootcontext.stcRescale = bootparams.stcRescale

	bootcontext.encoder = NewEncoder(params)
	bootcontext.evaluator = NewEvaluator(params)

	bootcontext.ctxpool[0] = NewCiphertext(params, 1, params.MaxLevel(), 0)
	bootcontext.ctxpool[1] = NewCiphertext(params, 1, params.MaxLevel(), 0)
	bootcontext.ctxpool[2] = NewCiphertext(params, 1, params.MaxLevel(), 0)
}

func (bootcontext *BootContext) newBootKeys(sk *SecretKey) {

	bootcontext.bootKeys = new(bootKeys)

	slots := bootcontext.slots

	// List of the rotation key values to needed for the bootstrapp
	rotations := []uint64{}

	//SubSum rotation needed X -> Y^slots rotations
	for i := bootcontext.params.LogSlots; i < bootcontext.params.LogN-1; i++ {
		if !utils.IsInSliceUint64(1<<i, rotations) {
			rotations = append(rotations, 1<<i)
		}
	}

	var index uint64
	// Coeffs to Slots rotations
	for i := range bootcontext.pDFTInv {
		for j := range bootcontext.pDFTInv[i].Vec {

			index = ((j / bootcontext.pDFTInv[i].N1) * bootcontext.pDFTInv[i].N1) & (slots - 1)

			if !utils.IsInSliceUint64(index, rotations) {
				rotations = append(rotations, index)
			}

			index = j & (bootcontext.pDFTInv[i].N1 - 1)

			if !utils.IsInSliceUint64(index, rotations) {
				rotations = append(rotations, index)
			}
		}
	}

	// Slots to Coeffs rotations
	for i := range bootcontext.pDFT {
		for j := range bootcontext.pDFT[i].Vec {

			if bootcontext.repack && i == 0 {
				// Sparse repacking, occuring during the first DFT matrix of the CoeffsToSlots.
				index = ((j / bootcontext.pDFT[i].N1) * bootcontext.pDFT[i].N1) & (2*slots - 1)
			} else {
				// Other cases
				index = ((j / bootcontext.pDFT[i].N1) * bootcontext.pDFT[i].N1) & (slots - 1)
			}

			if !utils.IsInSliceUint64(index, rotations) {
				rotations = append(rotations, index)
			}

			index = j & (bootcontext.pDFT[i].N1 - 1)

			if !utils.IsInSliceUint64(index, rotations) {
				rotations = append(rotations, index)
			}
		}
	}

	fmt.Println("DFT vector size (GB) :", float64(bootcontext.plaintextSize)/float64(1000000000))
	fmt.Println("Switching-Keys size (GB) :", float64(bootcontext.n*2*uint64(len(rotations))*bootcontext.params.Beta()*uint64(len(bootcontext.params.Qi)+len(bootcontext.params.Pi))*8)/float64(1000000000), "(", len(rotations), "keys)")

	kgen := NewKeyGenerator(bootcontext.params)

	bootcontext.rotkeys = NewRotationKeys()

	kgen.GenRot(Conjugate, sk, 0, bootcontext.rotkeys)

	for _, i := range rotations {
		kgen.GenRot(RotationLeft, sk, uint64(i), bootcontext.rotkeys)
	}

	bootcontext.relinkey = kgen.GenRelinKey(sk)

	return
}

func (bootcontext *BootContext) newBootSine() {

	bootcontext.bootSine = new(bootSine)

	sineDeg := 127

	K := complex(15, 0)

	bootcontext.chebycoeffs = Approximate(sin2pi2pi, -K, K, sineDeg)

	bootcontext.sinDepth = uint64(math.Ceil(math.Log2(float64(sineDeg))) + 1)

	if !bootcontext.ctsRescale {
		bootcontext.sinDepth++
	}
}

func (bootcontext *BootContextBetterSine) newBootBetterSine() {

	bootcontext.bootSine = new(bootSine)

	K := 16
	deg := 40
	dev := 10
	sc_num := 2

	bootcontext.sc_num = sc_num

	sc_fac := complex(float64(int(1<<sc_num)), 0)

	cheby := new(ChebyshevInterpolation)
	cheby.coeffs = bettersine.Approximate(K, deg, dev, sc_num)

	if sc_num == 1 {
		for i := range cheby.coeffs {
			cheby.coeffs[i] *= 0.5641895835477563
		}
	}

	if sc_num == 2 {
		for i := range cheby.coeffs {
			cheby.coeffs[i] *= 0.7511255444649425
		}
	}

	cheby.maxDeg = uint64(deg) + 1
	cheby.a = complex(float64(-K), 0) / sc_fac
	cheby.b = complex(float64(K), 0) / sc_fac

	bootcontext.chebycoeffs = cheby

	bootcontext.sinDepth = uint64(math.Ceil(math.Log2(float64(deg))) + float64(sc_num) + 1)

	if !bootcontext.ctsRescale {
		bootcontext.sinDepth++
	}
}

func (bootcontext *BootContext) newBootDFT() {

	bootcontext.bootDFT = new(bootDFT)

	// ==========================================================
	// Constant multiplication to be merged with the DFT vectors
	// ==========================================================

	if bootcontext.ctsRescale {
		// Change of variable + SubSum + CoeffsToSlots cancelling factor
		a := real(bootcontext.chebycoeffs.a)
		b := real(bootcontext.chebycoeffs.b)
		n := float64(bootcontext.n)

		// Change of variable for the evaluation of the Chebyshev polynomial + cancelling factor for the DFT and SubSum
		bootcontext.coeffsToSlotsDiffScale = complex(math.Pow(2/((b-a)*n), 1.0/float64(bootcontext.ctsDepth)), 0)

	} else {
		bootcontext.coeffsToSlotsDiffScale = complex(1, 0)

	}

	if bootcontext.stcRescale {
		// Rescaling factor to set the final ciphertext to the desired scale
		bootcontext.slotsToCoeffsDiffScale = complex(math.Pow(bootcontext.sineScale/bootcontext.params.Scale, 1.0/float64(bootcontext.stcDepth)), 0)
	} else {
		bootcontext.slotsToCoeffsDiffScale = complex(1, 0)
	}

	// Computation and encoding of the matrices for CoeffsToSlots and SlotsToCoeffs.
	bootcontext.computePlaintextVectors()
}

func (bootcontext *BootContextBetterSine) newBootDFTBetterSine() {

	bootcontext.bootDFT = new(bootDFT)

	// ==========================================================
	// Constant multiplication to be merged with the DFT vectors
	// ==========================================================

	if bootcontext.ctsRescale {

		// Change of variable + SubSum + CoeffsToSlots cancelling factor
		a := real(bootcontext.chebycoeffs.a)
		b := real(bootcontext.chebycoeffs.b)
		n := float64(bootcontext.n)
		sc_fac := float64(int(1 << bootcontext.sc_num))

		// Change of variable for the evaluation of the Chebyshev polynomial + cancelling factor for the DFT and SubSum
		bootcontext.coeffsToSlotsDiffScale = complex(math.Pow(2/((b-a)*n*sc_fac), 1.0/float64(bootcontext.ctsDepth)), 0)

	} else {
		bootcontext.coeffsToSlotsDiffScale = complex(1, 0)
	}

	if bootcontext.stcRescale {
		// Rescaling factor to set the final ciphertext to the desired scale
		bootcontext.slotsToCoeffsDiffScale = complex(math.Pow(bootcontext.sineScale/bootcontext.params.Scale, 1.0/float64(bootcontext.stcDepth)), 0)
	} else {
		bootcontext.slotsToCoeffsDiffScale = complex(1, 0)
	}

	// Computation and encoding of the matrices for CoeffsToSlots and SlotsToCoeffs.
	bootcontext.computePlaintextVectors()

	return
}

// NewBootContext creates a new bootcontext.
func NewBootContext(bootparams *BootParams, sk *SecretKey) (bootcontext *BootContext) {

	bootcontext = new(BootContext)
	bootcontext.newBootParams(bootparams)
	bootcontext.newBootSine()
	bootcontext.newBootDFT()
	bootcontext.newBootKeys(sk)

	return bootcontext
}

// NewBootContext creates a new bootcontext.
func NewBootContextBetterSine(bootparams *BootParams, sk *SecretKey) (bootcontext *BootContextBetterSine) {

	bootcontext = new(BootContextBetterSine)

	bootcontext.newBootParams(bootparams)
	bootcontext.newBootBetterSine()
	bootcontext.newBootDFTBetterSine()
	bootcontext.newBootKeys(sk)

	return bootcontext
}

// Bootstrapp re-encrypt a ciphertext at lvl Q0 to a ciphertext at MaxLevel.
func (bootcontext *BootContext) Bootstrapp(ct *Ciphertext) *Ciphertext {

	var ct0, ct1 *Ciphertext

	for ct.Level() != 0 {
		bootcontext.evaluator.DropLevel(ct, 1)
	}

	// TODO : better management of the initial scale
	bootcontext.evaluator.ScaleUp(ct, math.Round(bootcontext.sineScale/ct.Scale()), ct)

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

// Bootstrapp re-encrypt a ciphertext at lvl Q0 to a ciphertext at MaxLevel.
func (bootcontext *BootContextBetterSine) Bootstrapp(ct *Ciphertext) *Ciphertext {

	var ct0, ct1 *Ciphertext

	for ct.Level() != 0 {
		bootcontext.evaluator.DropLevel(ct, 1)
	}

	// TODO : better management of the initial scale
	bootcontext.evaluator.ScaleUp(ct, math.Round(bootcontext.sineScale/ct.Scale()), ct)

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

	for i := bootcontext.params.LogSlots; i < bootcontext.params.LogN-1; i++ {

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
		ct.Value()[u].Coeffs = append(ct.Value()[u].Coeffs, make([][]uint64, bootcontext.params.MaxLevel())...)
		for i := uint64(1); i < bootcontext.params.MaxLevel()+1; i++ {
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

			for i := uint64(1); i < bootcontext.params.MaxLevel()+1; i++ {

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

func (bootcontext *BootContext) evaluateSine(ct0, ct1 *Ciphertext) (*Ciphertext, *Ciphertext) {

	evaluator := bootcontext.evaluator.(*evaluator)

	// Reference scale is changed to the new ciphertext's scale.
	evaluator.ckksContext.scale = float64(bootcontext.params.Qi[ct0.Level()-1])

	// TODO : manage scale dynamicly depending on Q_0, the Qi of the SineEval and the ciphertext's scale.
	ct0.MulScale(1024)
	// Sine Evaluation ct0 = Q/(2pi) * sin((2pi/Q) * ct0)
	ct0 = bootcontext.evaluateChebyBoot(ct0)

	if bootcontext.stcRescale {
		ct0.SetScale(bootcontext.params.Scale)
	} else {
		ct0.DivScale(1024)
	}

	if ct1 != nil {

		ct1.MulScale(1024)
		// Sine Evaluation ct0 = Q/(2pi) * sin((2pi/Q) * ct0)
		ct1 = bootcontext.evaluateChebyBoot(ct1)

		if bootcontext.stcRescale {
			ct1.SetScale(bootcontext.params.Scale)
		} else {
			ct1.DivScale(1024)
		}
	}

	// Reference scale is changed back to the current ciphertext's scale.
	evaluator.ckksContext.scale = ct0.Scale()

	return ct0, ct1
}

func (bootcontext *BootContextBetterSine) evaluateSine(ct0, ct1 *Ciphertext) (*Ciphertext, *Ciphertext) {

	evaluator := bootcontext.evaluator.(*evaluator)

	// Reference scale is changed to the new ciphertext's scale.
	evaluator.ckksContext.scale = float64(bootcontext.params.Qi[ct0.Level()-1])

	// TODO : manage scale dynamicly depending on Q_0, the Qi of the SineEval and the ciphertext's scale.
	ct0.MulScale(1024)

	// Sine Evaluation ct0 = Q/(2pi) * sin((2pi/Q) * ct0)
	ct0 = bootcontext.evaluateChebyBootBetterSine(ct0)

	if bootcontext.stcRescale {
		ct0.SetScale(bootcontext.params.Scale)
	} else {
		ct0.DivScale(1024)
	}

	if ct1 != nil {

		ct1.MulScale(1024)

		ct1 = bootcontext.evaluateChebyBootBetterSine(ct1)

		if bootcontext.stcRescale {
			ct1.SetScale(bootcontext.params.Scale)
		} else {
			ct1.DivScale(1024)
		}
	}

	// Reference scale is changed back to the current ciphertext's scale.
	evaluator.ckksContext.scale = ct0.Scale()

	return ct0, ct1
}

func (bootcontext *BootContext) evaluateChebyBoot(ct *Ciphertext) (res *Ciphertext) {

	evaluator := bootcontext.evaluator.(*evaluator)

	// Chebyshev params
	a := bootcontext.chebycoeffs.a
	b := bootcontext.chebycoeffs.b
	degree := bootcontext.chebycoeffs.degree()

	C := make(map[uint64]*Ciphertext)
	C[1] = ct.CopyNew().Ciphertext()

	if bootcontext.ctsRescale {
		evaluator.AddConst(C[1], (-a-b)/(b-a), C[1])
	} else {
		n := complex(float64(bootcontext.n), 0)
		evaluator.MultByConst(C[1], 2/((b-a)*n), C[1])
		evaluator.AddConst(C[1], (-a-b)/(b-a), C[1])
		evaluator.Rescale(C[1], evaluator.ckksContext.scale, C[1])
	}

	M := uint64(bits.Len64(degree - 1))
	L := uint64(M >> 1)

	for i := uint64(2); i < (1<<L)+1; i++ {
		computePowerBasisCheby(i, C, evaluator, bootcontext.relinkey)
	}

	for i := L + 1; i < M; i++ {
		computePowerBasisCheby(1<<i, C, evaluator, bootcontext.relinkey)
	}

	res = recurseCheby(degree, L, M, bootcontext.chebycoeffs.Poly(), C, evaluator, bootcontext.relinkey)

	return
}

func (bootcontext *BootContextBetterSine) evaluateChebyBootBetterSine(ct *Ciphertext) (res *Ciphertext) {

	evaluator := bootcontext.evaluator.(*evaluator)

	cheby := bootcontext.chebycoeffs
	sc_num := bootcontext.sc_num
	sc_fac := complex(float64(int(1<<sc_num)), 0)

	// Chebyshev params
	a := cheby.a
	b := cheby.b
	degree := cheby.maxDeg - 1

	C := make(map[uint64]*Ciphertext)
	C[1] = ct.CopyNew().Ciphertext()

	if bootcontext.ctsRescale {
		evaluator.AddConst(C[1], (-a-b-(2*0.25/sc_fac))/(b-a), C[1])
	} else {
		n := complex(float64(bootcontext.n), 0)
		evaluator.AddConst(C[1], -0.25*n, C[1])
		evaluator.MultByConst(C[1], 2/((b-a)*n*sc_fac), C[1])
		evaluator.AddConst(C[1], (-a-b)/(b-a), C[1])
		evaluator.Rescale(C[1], evaluator.ckksContext.scale, C[1])
	}

	M := uint64(bits.Len64(degree - 1))
	L := uint64(2)

	for i := uint64(2); i < (1<<L)+1; i++ {
		computePowerBasisCheby(i, C, evaluator, bootcontext.relinkey)
	}

	for i := L + 1; i < M; i++ {
		computePowerBasisCheby(1<<i, C, evaluator, bootcontext.relinkey)
	}

	res = recurseCheby(degree, L, M, cheby.Poly(), C, evaluator, bootcontext.relinkey)

	if sc_num == 1 {
		evaluator.MulRelin(res, res, bootcontext.relinkey, res)
		evaluator.Rescale(res, evaluator.ckksContext.scale, res)
		evaluator.AddConst(res, -1.0/6.283185307179586, res)
	}

	if sc_num == 2 {
		//fmt.Println()

		//fmt.Println("Mul", res.Level(), res.Level(), res.Scale())
		evaluator.MulRelin(res, res, bootcontext.relinkey, res)
		evaluator.Rescale(res, evaluator.ckksContext.scale, res)

		y := evaluator.AddConstNew(res, -0.5641895835477563)
		//fmt.Println("Mul", res.Level(), y.Level())
		evaluator.MulRelin(res, y, bootcontext.relinkey, res)

		evaluator.MultByConst(res, 4, res)
		evaluator.AddConst(res, 1.0/6.283185307179586, res)

		evaluator.Rescale(res, evaluator.ckksContext.scale, res)
	}

	return
}

func (bootcontext *BootContext) multiplyByDiagMatrice(vec *Ciphertext, plainVectors *dftvectors) (res *Ciphertext) {

	evaluator := bootcontext.evaluator

	var N1 uint64

	res = NewCiphertext(bootcontext.params, 1, vec.Level(), vec.Scale())

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

	tmpVec = NewCiphertext(bootcontext.params, 1, bootcontext.params.MaxLevel(), vec.Scale())
	tmp = NewCiphertext(bootcontext.params, 1, bootcontext.params.MaxLevel(), vec.Scale())

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

func computeRoots(N uint64) (roots []complex128) {

	var angle float64

	m := N << 1

	roots = make([]complex128, m)

	roots[0] = 1

	for i := uint64(1); i < m; i++ {
		angle = 6.283185307179586 * float64(i) / float64(m)
		roots[i] = complex(math.Cos(angle), math.Sin(angle))
	}

	return
}

func fftPlainVec(N uint64, roots []complex128, pow5 []uint64) (a, b, c [][]complex128) {

	var logN, m, index, tt, gap, k, mask, idx1, idx2 uint64

	logN = uint64(bits.Len64(N) - 1)

	a = make([][]complex128, logN)
	b = make([][]complex128, logN)
	c = make([][]complex128, logN)

	index = 0
	for m = 2; m <= N; m <<= 1 {

		a[index] = make([]complex128, 2*N)
		b[index] = make([]complex128, 2*N)
		c[index] = make([]complex128, 2*N)

		tt = m >> 1

		for i := uint64(0); i < N; i += m {

			gap = N / m
			mask = (m << 2) - 1

			for j := uint64(0); j < m>>1; j++ {

				k = (pow5[j] & mask) * gap

				idx1 = i + j
				idx2 = i + j + tt

				for u := uint64(0); u < 2; u++ {
					a[index][idx1+u*N] = 1
					a[index][idx2+u*N] = -roots[k]
					b[index][idx1+u*N] = roots[k]
					c[index][idx2+u*N] = 1
				}
			}
		}

		index++
	}

	return
}

func fftInvPlainVec(N uint64, roots []complex128, pow5 []uint64) (a, b, c [][]complex128) {

	var logN, m, index, tt, gap, k, mask, idx1, idx2 uint64

	logN = uint64(bits.Len64(N) - 1)

	a = make([][]complex128, logN)
	b = make([][]complex128, logN)
	c = make([][]complex128, logN)

	index = 0
	for m = N; m >= 2; m >>= 1 {

		a[index] = make([]complex128, 2*N)
		b[index] = make([]complex128, 2*N)
		c[index] = make([]complex128, 2*N)

		tt = m >> 1

		for i := uint64(0); i < N; i += m {

			gap = N / m
			mask = (m << 2) - 1

			for j := uint64(0); j < m>>1; j++ {

				k = ((m << 2) - (pow5[j] & mask)) * gap

				idx1 = i + j
				idx2 = i + j + tt

				for u := uint64(0); u < 2; u++ {

					a[index][idx1+u*N] = 1
					a[index][idx2+u*N] = -roots[k]
					b[index][idx1+u*N] = 1
					c[index][idx2+u*N] = roots[k]
				}
			}
		}

		index++
	}

	return
}

func (bootcontext *BootContext) computePlaintextVectors() {

	slots := bootcontext.slots
	dslots := bootcontext.dslots

	ctsDepth := bootcontext.ctsDepth
	stcDepth := bootcontext.stcDepth

	roots := computeRoots(slots << 1)
	pow5 := make([]uint64, (slots<<1)+1)
	pow5[0] = 1
	for i := uint64(1); i < (slots<<1)+1; i++ {
		pow5[i] = pow5[i-1] * 5
		pow5[i] &= (slots << 2) - 1
	}

	// CoeffsToSlots vectors
	bootcontext.pDFTInv = make([]*dftvectors, ctsDepth)

	pVecDFTInv := bootcontext.computeDFTPlaintextVectors(roots, pow5, bootcontext.coeffsToSlotsDiffScale, true)

	for i := uint64(0); i < bootcontext.ctsDepth; i++ {

		bootcontext.pDFTInv[i] = new(dftvectors)

		bootcontext.pDFTInv[i].N1 = findbestbabygiantstepsplit(pVecDFTInv[i], dslots)

		bootcontext.encodePVec(pVecDFTInv[i], bootcontext.pDFTInv[i], i, true)
	}

	// SlotsToCoeffs vectors
	bootcontext.pDFT = make([]*dftvectors, stcDepth)

	pVecDFT := bootcontext.computeDFTPlaintextVectors(roots, pow5, bootcontext.slotsToCoeffsDiffScale, false)

	for i := uint64(0); i < stcDepth; i++ {

		bootcontext.pDFT[i] = new(dftvectors)

		bootcontext.pDFT[i].N1 = findbestbabygiantstepsplit(pVecDFT[i], dslots)

		bootcontext.encodePVec(pVecDFT[i], bootcontext.pDFT[i], i, false)
	}
}

// Finds the best N1*N2 = N for the baby-step giant-step algorithm for matrix multiplication.
func findbestbabygiantstepsplit(vector map[uint64][]complex128, maxN uint64) (minN uint64) {

	var sum uint64

	sum = maxN

	for N1 := uint64(1); N1 < maxN; N1 <<= 1 {

		index := make(map[uint64][]uint64)

		for key := range vector {

			idx1 := key / N1
			idx2 := key & (N1 - 1)

			if index[idx1] == nil {
				index[idx1] = []uint64{idx2}
			} else {
				index[idx1] = append(index[idx1], idx2)
			}
		}

		if uint64(len(index)+len(index[0])) < sum {
			minN = N1
			sum = uint64(len(index) + len(index[0]))
		}
	}

	if minN == 0 {
		minN = 1
	}

	return
}

func (bootcontext *BootContext) encodePVec(pVec map[uint64][]complex128, plaintextVec *dftvectors, k uint64, forward bool) {
	var N, N1, level uint64

	// N1*N2 = N
	N = bootcontext.n
	N1 = plaintextVec.N1

	index := make(map[uint64][]uint64)

	for key := range pVec {
		idx1 := key / N1
		idx2 := key & (N1 - 1)
		if index[idx1] == nil {
			index[idx1] = []uint64{idx2}
		} else {
			index[idx1] = append(index[idx1], idx2)
		}
	}

	plaintextVec.Vec = make(map[uint64]*Plaintext)

	var scale float64

	for j := range index {

		for _, i := range index[j] {

			if forward {
				level = bootcontext.params.MaxLevel() - k
				scale = float64(bootcontext.params.Qi[level])
			} else {
				level = bootcontext.params.MaxLevel() - k/2 - bootcontext.ctsDepth - bootcontext.sinDepth

				if bootcontext.stcRescale {
					scale = 1 << 30
				} else {
					scale = 1 << 25
				}
			}

			plaintextVec.Vec[N1*j+uint64(i)] = NewPlaintext(bootcontext.params, level, scale)

			bootcontext.plaintextSize += (level + 1) * 8 * bootcontext.n

			bootcontext.encoder.Encode(plaintextVec.Vec[N1*j+uint64(i)], rotate(pVec[N1*j+uint64(i)], (N>>1)-(N1*j))[:bootcontext.dslots], bootcontext.dslots)
		}
	}
}

func (bootcontext *BootContext) computeDFTPlaintextVectors(roots []complex128, pow5 []uint64, diffscale complex128, forward bool) (plainVector []map[uint64][]complex128) {

	var level, depth, nextLevel, slots, logSlots uint64

	slots = bootcontext.slots

	logSlots = uint64(bits.Len64(bootcontext.slots) - 1)

	level = logSlots

	var a, b, c [][]complex128
	var maxDepth uint64

	if forward {
		maxDepth = bootcontext.ctsDepth
		a, b, c = fftInvPlainVec(slots, roots, pow5)
	} else {
		maxDepth = bootcontext.stcDepth
		a, b, c = fftPlainVec(slots, roots, pow5)
	}

	plainVector = make([]map[uint64][]complex128, maxDepth)

	// We compute the chain of merge in order or reverse order depending if its DFT or InvDFT because
	// the way the levels are collapsed has an inpact on the total number of rotations and keys to be
	// stored. Ex. instead of using 255 + 64 plaintext vectors, we can use 127 + 128 plaintext vectors
	// by reversing the order of the merging.
	merge := make([]uint64, maxDepth)
	for i := uint64(0); i < maxDepth; i++ {

		depth = uint64(math.Ceil(float64(level) / float64(maxDepth-i)))

		if forward {
			merge[i] = depth
		} else {
			merge[uint64(len(merge))-i-1] = depth

		}

		level -= depth
	}

	level = logSlots
	for i := uint64(0); i < maxDepth; i++ {

		if bootcontext.repack && !forward && i == 0 {

			// Special initial matrix for the repacking before SlotsToCoeffs
			plainVector[i] = genWfftRepack(logSlots, level)

			// Merges this special initial matrix with the first layer of SlotsToCoeffs DFT
			plainVector[i] = nextLevelfft(plainVector[i], logSlots, 2<<logSlots, level, a[logSlots-level], b[logSlots-level], c[logSlots-level], forward)

			// Continues the merging with the next layers if the total depth requires it.
			nextLevel = level - 1
			for j := uint64(0); j < merge[i]-1; j++ {
				plainVector[i] = nextLevelfft(plainVector[i], logSlots, 2<<logSlots, nextLevel, a[logSlots-nextLevel], b[logSlots-nextLevel], c[logSlots-nextLevel], forward)
				nextLevel--
			}

		} else {
			// First layer of the i-th level of the DFT
			plainVector[i] = genWfft(logSlots, level, a[logSlots-level], b[logSlots-level], c[logSlots-level], forward)

			// Merges the layer with the next levels of the DFT if the total depth requires it.
			nextLevel = level - 1
			for j := uint64(0); j < merge[i]-1; j++ {
				plainVector[i] = nextLevelfft(plainVector[i], logSlots, 1<<logSlots, nextLevel, a[logSlots-nextLevel], b[logSlots-nextLevel], c[logSlots-nextLevel], forward)
				nextLevel--
			}
		}

		level -= merge[i]
	}

	// Repacking after the CoeffsToSlots (we multiply the last DFT matrix with the vector [1, 1, ..., 1, 1, 0, 0, ..., 0, 0]).
	if bootcontext.repack && forward {
		for j := range plainVector[maxDepth-1] {
			for x := uint64(0); x < slots; x++ {
				plainVector[maxDepth-1][j][x+bootcontext.slots] = complex(0, 0)
			}
		}
	}

	if !forward {
		// Rescaling of the DFT matrix of the SlotsToCoeffs to match the desired output scale
		for j := range plainVector {
			for x := range plainVector[j] {
				for i := range plainVector[j][x] {
					plainVector[j][x][i] /= diffscale
				}
			}
		}
	} else {
		// Rescaling of the DFT matrix of the SlotsToCoeffs to operate the change of variable for
		// the evaluation of the Chebyshev polynomial and the DFT + SubSum cancellation factor.
		for j := range plainVector {
			for x := range plainVector[j] {
				for i := range plainVector[j][x] {
					plainVector[j][x][i] *= diffscale
				}
			}
		}
	}

	return
}

func genWfft(logL, level uint64, a, b, c []complex128, forward bool) (vectors map[uint64][]complex128) {

	var rot uint64

	if forward {
		rot = 1 << (level - 1)
	} else {
		rot = 1 << (logL - level)
	}

	vectors = make(map[uint64][]complex128)

	addToDicVector(vectors, 0, a)
	addToDicVector(vectors, rot, b)
	addToDicVector(vectors, (1<<logL)-rot, c)

	return
}

func genWfftRepack(logL, level uint64) (vectors map[uint64][]complex128) {

	vectors = make(map[uint64][]complex128)

	a := make([]complex128, 2<<logL)
	b := make([]complex128, 2<<logL)

	for i := uint64(0); i < 1<<logL; i++ {
		a[i] = complex(1, 0)
		a[i+(1<<logL)] = complex(0, 1)

		b[i] = complex(0, 1)
		b[i+(1<<logL)] = complex(1, 0)
	}

	addToDicVector(vectors, 0, a)
	addToDicVector(vectors, (1 << logL), b)

	return
}

func nextLevelfft(vec map[uint64][]complex128, logL, N, nextLevel uint64, a, b, c []complex128, forward bool) (newVec map[uint64][]complex128) {

	var rot uint64

	newVec = make(map[uint64][]complex128)

	if forward {
		rot = (1 << (nextLevel - 1)) & (N - 1)
	} else {
		rot = (1 << (logL - nextLevel)) & (N - 1)
	}

	for i := range vec {
		addToDicVector(newVec, i, mul(vec[i], a))
		addToDicVector(newVec, (i+rot)&(N-1), mul(rotate(vec[i], rot), b))
		addToDicVector(newVec, (i+N-rot)&(N-1), mul(rotate(vec[i], N-rot), c))
	}

	return
}

func addToDicVector(dic map[uint64][]complex128, index uint64, vec []complex128) {
	if dic[index] == nil {
		dic[index] = vec
	} else {
		dic[index] = add(dic[index], vec)
	}
}

func rotate(x []complex128, n uint64) (y []complex128) {

	y = make([]complex128, len(x))

	mask := uint64(len(x) - 1)

	// Rotates to the left
	for i := uint64(0); i < uint64(len(x)); i++ {
		y[i] = x[(i+n)&mask]
	}

	return
}

func mul(a, b []complex128) (res []complex128) {

	res = make([]complex128, len(a))

	for i := 0; i < len(a); i++ {
		res[i] = a[i] * b[i]
	}

	return
}

func add(a, b []complex128) (res []complex128) {

	res = make([]complex128, len(a))

	for i := 0; i < len(a); i++ {
		res[i] = a[i] + b[i]
	}

	return
}
