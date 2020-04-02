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
type BootContext struct {
	BootParams

	n    uint64 // Ring degree
	logn uint64 // log of the Ring degree

	slots    uint64 // Number of plaintext slots
	logSlots uint64 // Log of the number of plaintext slots

	dslots    uint64 // Number of plaintext slots after the re-encoding
	logdslots uint64 // Log of the number of plaintext slots after the re-encoding

	encoder   Encoder   // Encoder
	evaluator Evaluator // Evaluator

	plaintextSize uint64 // Byte size of the plaintext DFT matrices

	repack      bool                    // If true then can repack the CoeffsToSlots into on ciphertext
	sinScale    float64                 // Input scale to the SineEval
	chebycoeffs *ChebyshevInterpolation // Coefficients of the Chebyshev Interpolation of sin(2*pi*x) or cos(2*pi*x/r)

	coeffsToSlotsDiffScale complex128    // Matrice rescaling
	slotsToCoeffsDiffScale complex128    // Matrice rescaling
	pDFT                   []*dftvectors // Matrice vectors
	pDFTInv                []*dftvectors // Matrice vectors

	relinkey *EvaluationKey // Relinearization key
	rotkeys  *RotationKeys  // Rotation and conjugation keys

	ctxpool [3]*Ciphertext // Memory pool

	decryptor Decryptor
}

type dftvectors struct {
	N1  uint64
	Vec map[uint64]*Plaintext
}

func sin2pi2pi(x complex128) complex128 {
	return cmplx.Sin(6.283185307179586*x) / 6.283185307179586
}

func (b *BootContext) printDebug(message string, ciphertext *Ciphertext) {

	coeffs := b.encoder.Decode(b.decryptor.DecryptNew(ciphertext), b.dslots)

	if b.dslots == 2 {
		fmt.Printf(message+"%.10f %.10f...\n", coeffs[0], coeffs[1])
	} else {
		fmt.Printf(message+"%.10f %.10f %.10f %.10f...\n", coeffs[0], coeffs[1], coeffs[2], coeffs[3])
	}
}

// NewBootContext creates a new bootcontext.
func NewBootContext(bootparams *BootParams, sk *SecretKey) (bootcontext *BootContext) {

	bootcontext = new(BootContext)

	bootcontext.BootParams = *bootparams

	bootcontext.n = uint64(1 << bootparams.Parameters.LogN)
	bootcontext.slots = uint64(1 << bootparams.Parameters.LogSlots)

	if bootparams.Parameters.LogSlots < bootparams.Parameters.LogN-1 {
		bootcontext.repack = true
	}

	bootcontext.dslots = bootcontext.slots
	if bootparams.Parameters.LogSlots < bootparams.Parameters.LogN-1 {
		bootcontext.dslots <<= 1
	}

	bootcontext.sinScale = 1 << 45

	bootcontext.encoder = NewEncoder(&bootparams.Parameters)

	bootcontext.newBootSine()
	bootcontext.newBootDFT()
	bootcontext.newBootKeys(sk)

	bootcontext.evaluator = NewEvaluator(&bootparams.Parameters)
	bootcontext.decryptor = NewDecryptor(&bootparams.Parameters, sk)

	bootcontext.ctxpool[0] = NewCiphertext(&bootparams.Parameters, 1, bootparams.Parameters.MaxLevel, 0)
	bootcontext.ctxpool[1] = NewCiphertext(&bootparams.Parameters, 1, bootparams.Parameters.MaxLevel, 0)
	bootcontext.ctxpool[2] = NewCiphertext(&bootparams.Parameters, 1, bootparams.Parameters.MaxLevel, 0)

	return bootcontext
}

func (bootcontext *BootContext) newBootDFT() {

	if bootcontext.SinType == Sin {

		if bootcontext.CtSRescale {
			// Change of variable + SubSum + CoeffsToSlots cancelling factor
			a := real(bootcontext.chebycoeffs.a)
			b := real(bootcontext.chebycoeffs.b)
			n := float64(bootcontext.n)

			// Change of variable for the evaluation of the Chebyshev polynomial + cancelling factor for the DFT and SubSum
			bootcontext.coeffsToSlotsDiffScale = complex(math.Pow(2/((b-a)*n), 1.0/float64(bootcontext.CtSDepth)), 0)

		} else {
			bootcontext.coeffsToSlotsDiffScale = complex(1, 0)
		}

		if bootcontext.StCRescale {
			// Rescaling factor to set the final ciphertext to the desired scale
			bootcontext.slotsToCoeffsDiffScale = complex(math.Pow(bootcontext.sinScale/bootcontext.Scale, 1.0/float64(bootcontext.StCDepth)), 0)
		} else {
			bootcontext.slotsToCoeffsDiffScale = complex(1, 0)
		}

	} else if bootcontext.SinType == Cos {

		if bootcontext.CtSRescale {

			// Change of variable + SubSum + CoeffsToSlots cancelling factor
			a := real(bootcontext.chebycoeffs.a)
			b := real(bootcontext.chebycoeffs.b)
			n := float64(bootcontext.n)
			sc_fac := float64(int(1 << bootcontext.SinRescal))

			// Change of variable for the evaluation of the Chebyshev polynomial + cancelling factor for the DFT and SubSum
			bootcontext.coeffsToSlotsDiffScale = complex(math.Pow(2/((b-a)*n*sc_fac), 1.0/float64(bootcontext.CtSDepth)), 0)

		} else {
			bootcontext.coeffsToSlotsDiffScale = complex(1, 0)
		}

		if bootcontext.StCRescale {
			// Rescaling factor to set the final ciphertext to the desired scale
			bootcontext.slotsToCoeffsDiffScale = complex(math.Pow(bootcontext.sinScale/bootcontext.Scale, 1.0/float64(bootcontext.StCDepth)), 0)
		} else {
			bootcontext.slotsToCoeffsDiffScale = complex(1, 0)
		}

	} else {
		panic("bootcontext -> invalid sineType")
	}

	// Computation and encoding of the matrices for CoeffsToSlots and SlotsToCoeffs.
	bootcontext.computePlaintextVectors()

	return
}

func (bootcontext *BootContext) newBootSine() {

	if bootcontext.SinType == Sin {

		K := complex(float64(bootcontext.SinRange), 0)

		bootcontext.chebycoeffs = Approximate(sin2pi2pi, -K, K, int(bootcontext.SinDeg))

	} else if bootcontext.SinType == Cos {

		K := int(bootcontext.SinRange)
		deg := int(bootcontext.SinDeg)
		dev := 10
		sc_fac := complex(float64(int(1<<bootcontext.SinRescal)), 0)

		cheby := new(ChebyshevInterpolation)

		cheby.coeffs = bettersine.Approximate(K, deg, dev, int(bootcontext.SinRescal))

		if int(bootcontext.SinRescal) == 1 {
			for i := range cheby.coeffs {
				cheby.coeffs[i] *= 0.5641895835477563
			}
		}

		if int(bootcontext.SinRescal) == 2 {
			for i := range cheby.coeffs {
				cheby.coeffs[i] *= 0.7511255444649425
			}
		}

		if int(bootcontext.SinRescal) == 3 {
			for i := range cheby.coeffs {
				cheby.coeffs[i] *= 0.8666749935615672
			}
		}

		cheby.maxDeg = uint64(deg) + 1
		cheby.a = complex(float64(-K), 0) / sc_fac
		cheby.b = complex(float64(K), 0) / sc_fac

		bootcontext.chebycoeffs = cheby

	} else {
		panic("bootcontext -> invalid sineType")
	}
}

func (bootcontext *BootContext) newBootKeys(sk *SecretKey) {

	slots := bootcontext.slots

	// List of the rotation key values to needed for the bootstrapp
	rotations := []uint64{}

	//SubSum rotation needed X -> Y^slots rotations
	for i := bootcontext.LogSlots; i < bootcontext.LogN-1; i++ {
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
	fmt.Println("Switching-Keys size (GB) :", float64(bootcontext.n*2*uint64(len(rotations))*bootcontext.Beta*uint64(len(bootcontext.Qi)+len(bootcontext.Pi))*8)/float64(1000000000), "(", len(rotations), "keys)")

	kgen := NewKeyGenerator(&bootcontext.Parameters)

	bootcontext.rotkeys = NewRotationKeys()

	kgen.GenRot(Conjugate, sk, 0, bootcontext.rotkeys)

	for _, i := range rotations {
		kgen.GenRot(RotationLeft, sk, uint64(i), bootcontext.rotkeys)
	}

	bootcontext.relinkey = kgen.GenRelinKey(sk)

	return
}

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

func (bootcontext *BootContext) evaluateSine(ct0, ct1 *Ciphertext) (*Ciphertext, *Ciphertext) {

	evaluator := bootcontext.evaluator.(*evaluator)

	// Reference scale is changed to the new ciphertext's scale.
	evaluator.ckksContext.scale = float64(bootcontext.Qi[ct0.Level()-1])

	// TODO : manage scale dynamicly depending on Q_0, the Qi of the SineEval and the ciphertext's scale.
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

	if bootcontext.CtSRescale {
		evaluator.AddConst(C[1], (-a-b)/(b-a), C[1])
	} else {
		n := complex(float64(bootcontext.n), 0)
		evaluator.MultByConst(C[1], 2/((b-a)*n), C[1])
		evaluator.AddConst(C[1], (-a-b)/(b-a), C[1])
		evaluator.Rescale(C[1], evaluator.ckksContext.scale, C[1])
	}

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

	if bootcontext.CtSRescale {
		evaluator.AddConst(C[1], (-a-b-(2*0.25/sc_fac))/(b-a), C[1])
	} else {
		n := complex(float64(bootcontext.n), 0)
		evaluator.AddConst(C[1], -0.25*n, C[1])
		evaluator.MultByConst(C[1], 2/((b-a)*n*sc_fac), C[1])
		evaluator.AddConst(C[1], (-a-b)/(b-a), C[1])
		evaluator.Rescale(C[1], evaluator.ckksContext.scale, C[1])
	}

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

	CtSDepth := bootcontext.CtSDepth
	StCDepth := bootcontext.StCDepth

	roots := computeRoots(slots << 1)
	pow5 := make([]uint64, (slots<<1)+1)
	pow5[0] = 1
	for i := uint64(1); i < (slots<<1)+1; i++ {
		pow5[i] = pow5[i-1] * 5
		pow5[i] &= (slots << 2) - 1
	}

	// CoeffsToSlots vectors
	bootcontext.pDFTInv = make([]*dftvectors, CtSDepth)

	pVecDFTInv := bootcontext.computeDFTPlaintextVectors(roots, pow5, bootcontext.coeffsToSlotsDiffScale, true)

	for i := uint64(0); i < bootcontext.CtSDepth; i++ {

		bootcontext.pDFTInv[i] = new(dftvectors)

		bootcontext.pDFTInv[i].N1 = findbestbabygiantstepsplit(pVecDFTInv[i], dslots)

		bootcontext.encodePVec(pVecDFTInv[i], bootcontext.pDFTInv[i], i, true)
	}

	// SlotsToCoeffs vectors
	bootcontext.pDFT = make([]*dftvectors, StCDepth)

	pVecDFT := bootcontext.computeDFTPlaintextVectors(roots, pow5, bootcontext.slotsToCoeffsDiffScale, false)

	for i := uint64(0); i < StCDepth; i++ {

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
				level = bootcontext.MaxLevel - k
				scale = float64(bootcontext.Qi[level])
			} else {
				level = bootcontext.MaxLevel - uint64((float64(k)/2.0)+0.5) - bootcontext.CtSDepth - bootcontext.SinDepth

				// If the first moduli
				if bootcontext.LogQi[level] > 30 {
					scale = float64(uint64(1 << (bootcontext.LogQi[level] >> 1)))
				} else {
					scale = float64(bootcontext.Qi[level])
				}
			}

			plaintextVec.Vec[N1*j+uint64(i)] = NewPlaintext(&bootcontext.Parameters, level, scale)

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
		maxDepth = bootcontext.CtSDepth
		a, b, c = fftInvPlainVec(slots, roots, pow5)
	} else {
		maxDepth = bootcontext.StCDepth
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
