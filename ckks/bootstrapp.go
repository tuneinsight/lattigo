package ckks

import (
	"fmt"
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
	"math"
	"math/bits"
	"math/cmplx"
)

type BootContext struct {
	ckkscontext *Context

	encoder *Encoder
	slots   uint64
	dslots  uint64

	//Sine evaluation
	chebycoeffs *ChebyshevInterpolation

	plaintextSize uint64

	//Coeffs to slots and slots to coeffs
	ctsDepth uint64
	stcDepth uint64
	sinDepth uint64
	repack   bool
	pDFT     []*dftvectors
	pDFTInv  []*dftvectors

	relinkey *EvaluationKey
	rotkeys  *RotationKeys

	ctxpool [3]*Ciphertext
}

type dftvectors struct {
	N1  uint64
	Vec map[uint64]*Plaintext
}

func sin2pi2pi(x complex128) complex128 {
	return cmplx.Sin(6.283185307179586*(1.0/1.0)*x) * (1.0 / 6.283185307179586)
}

func showcoeffs(decryptor *Decryptor, encoder *Encoder, slots uint64, ciphertext *Ciphertext, message string) (coeffs []complex128) {

	coeffs = encoder.Decode(decryptor.DecryptNew(ciphertext), slots)

	if slots == 2 {
		fmt.Printf(message+"%22.10f %22.10f...\n", coeffs[0], coeffs[1])
	} else {
		fmt.Printf(message+"%22.10f %22.10f %22.10f %22.10f...\n", coeffs[0], coeffs[1], coeffs[2], coeffs[3])
	}

	return coeffs
}

func (ckkscontext *Context) NewBootContext(slots uint64, sk *SecretKey, ctsDepth, stcDepth uint64) (bootcontext *BootContext, err error) {

	bootcontext = new(BootContext)

	bootcontext.ckkscontext = ckkscontext

	bootcontext.ctsDepth = ctsDepth
	bootcontext.stcDepth = stcDepth

	if slots < ckkscontext.maxSlots {
		bootcontext.repack = true
	}

	bootcontext.slots = slots

	if slots == ckkscontext.maxSlots {
		bootcontext.dslots = ckkscontext.maxSlots
	} else {
		bootcontext.dslots = slots << 1
	}

	bootcontext.encoder = ckkscontext.NewEncoder()

	sineDeg := 127

	bootcontext.chebycoeffs = Approximate(sin2pi2pi, -15, 15, sineDeg)

	bootcontext.sinDepth = uint64(math.Ceil(math.Log2(float64(sineDeg))) + 2)

	// List of the rotation key values to needed for the bootstrapp
	rotations := []uint64{}

	//SubSum rotation needed X -> Y^slots rotations
	logSlots := uint64(bits.Len64(bootcontext.slots) - 1)
	logMaxSlots := uint64(bits.Len64(bootcontext.ckkscontext.maxSlots) - 1)
	for i := logSlots; i < logMaxSlots; i++ {
		if !utils.IsInSliceUint64(1<<i, rotations) {
			rotations = append(rotations, 1<<i)
		}
	}

	// Computation and encoding of the matrices for CoeffsToSlots and SlotsToCoeffs.
	bootcontext.computePlaintextVectors()

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
	fmt.Println("Switching-Keys size (GB) :", float64(ckkscontext.n*2*uint64(len(rotations))*ckkscontext.Beta()*uint64(len(ckkscontext.contextKeys.Modulus))*8)/float64(1000000000), "(", len(rotations), "keys)")

	kgen := ckkscontext.NewKeyGenerator()

	bootcontext.rotkeys = ckkscontext.NewRotationKeys()

	kgen.GenRot(Conjugate, sk, 0, bootcontext.rotkeys)

	for _, i := range rotations {
		kgen.GenRot(RotationLeft, sk, uint64(i), bootcontext.rotkeys)
	}

	bootcontext.relinkey = kgen.NewRelinKey(sk)

	bootcontext.ctxpool[0] = ckkscontext.NewCiphertext(1, ckkscontext.levels-1, 0)
	bootcontext.ctxpool[1] = ckkscontext.NewCiphertext(1, ckkscontext.levels-1, 0)
	bootcontext.ctxpool[2] = ckkscontext.NewCiphertext(1, ckkscontext.levels-1, 0)

	return bootcontext, nil
}

func (evaluator *Evaluator) Bootstrapp(ct *Ciphertext, bootcontext *BootContext) *Ciphertext {

	var ct0, ct1 *Ciphertext

	for ct.Level() != 0 {
		evaluator.DropLevel(ct, 1)
	}

	// TODO : better management of the initial scale
	evaluator.ScaleUp(ct, math.Round(float64(1<<45)/ct.Scale()), ct)

	// ModUp ct_{Q_0} -> ct_{Q_L}
	ct = bootcontext.modUp(ct)

	//SubSum X -> (N/dslots) * Y^dslots
	logSlots := uint64(bits.Len64(bootcontext.slots) - 1)
	logMaxSlots := uint64(bits.Len64(bootcontext.ckkscontext.maxSlots) - 1)
	for i := logSlots; i < logMaxSlots; i++ {

		evaluator.RotateColumns(ct, 1<<i, bootcontext.rotkeys, bootcontext.ctxpool[0])

		evaluator.Add(ct, bootcontext.ctxpool[0], ct)
	}

	// Part 1 : Coeffs to slots
	ct0, ct1 = bootcontext.coeffsToSlots(evaluator, ct)

	// Part 2 : SineEval
	ct0, ct1 = bootcontext.evaluateSine(ct0, ct1, evaluator)

	// Part 3 : Slots to coeffs
	return bootcontext.slotsToCoeffs(evaluator, ct0, ct1)
}

func (bootcontext *BootContext) modUp(ct *Ciphertext) *Ciphertext {

	ct.InvNTT(bootcontext.ckkscontext, ct.Element())

	// Extend the ciphertext with zero polynomials.
	for u := range ct.Value() {
		ct.Value()[u].Coeffs = append(ct.Value()[u].Coeffs, make([][]uint64, bootcontext.ckkscontext.levels-1)...)
		for i := uint64(1); i < bootcontext.ckkscontext.levels; i++ {
			ct.Value()[u].Coeffs[i] = make([]uint64, bootcontext.ckkscontext.n)
		}
	}

	//Centers the values around Q0 and extends the basis from Q0 to QL
	Q := bootcontext.ckkscontext.moduli[0]
	bredparams := bootcontext.ckkscontext.contextQ.GetBredParams()

	var coeff, qi uint64
	for u := range ct.Value() {

		for j := uint64(0); j < bootcontext.ckkscontext.n; j++ {

			coeff = ct.Value()[u].Coeffs[0][j]

			for i := uint64(1); i < bootcontext.ckkscontext.levels; i++ {

				qi = bootcontext.ckkscontext.moduli[i]

				if coeff > (Q >> 1) {
					ct.Value()[u].Coeffs[i][j] = qi - ring.BRedAdd(Q-coeff+1, qi, bredparams[i])
				} else {
					ct.Value()[u].Coeffs[i][j] = ring.BRedAdd(coeff, qi, bredparams[i])
				}
			}

			qi = bootcontext.ckkscontext.moduli[0]

			if coeff > (Q >> 1) {
				ct.Value()[u].Coeffs[0][j] = qi - ring.BRedAdd(Q-coeff+1, qi, bredparams[0])
			} else {
				ct.Value()[u].Coeffs[0][j] = ring.BRedAdd(coeff, qi, bredparams[0])
			}
		}
	}

	ct.NTT(bootcontext.ckkscontext, ct.Element())

	return ct
}

func (bootcontext *BootContext) coeffsToSlots(evaluator *Evaluator, vec *Ciphertext) (ct0, ct1 *Ciphertext) {

	var zV, zVconj *Ciphertext

	zV = bootcontext.dft(evaluator, vec, bootcontext.pDFTInv, true)

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

func (bootcontext *BootContext) slotsToCoeffs(evaluator *Evaluator, ct0, ct1 *Ciphertext) (ct *Ciphertext) {

	// If full packing, the repacking can be done directly using ct0 and ct1.
	if !bootcontext.repack {

		evaluator.MultByi(ct1, ct1)

		evaluator.Add(ct0, ct1, ct0)
	}

	return bootcontext.dft(evaluator, ct0, bootcontext.pDFT, false)
}

func (bootcontext *BootContext) dft(evaluator *Evaluator, vec *Ciphertext, plainVectors []*dftvectors, forward bool) (w *Ciphertext) {

	levels := uint64(len(plainVectors))

	w = vec.CopyNew().Ciphertext()

	// Sequencially multiplies w with the provided dft matrices.
	for i := uint64(0); i < levels; i++ {

		w = bootcontext.multiplyByDiagMatrice(evaluator, w, plainVectors[i])

		evaluator.Rescale(w, evaluator.ckkscontext.scale, w)
	}

	return w
}

func (bootcontext *BootContext) evaluateSine(ct0, ct1 *Ciphertext, evaluator *Evaluator) (*Ciphertext, *Ciphertext) {

	// Reference scale is changed to the new ciphertext's scale.
	bootcontext.ckkscontext.scale = bootcontext.ckkscontext.scalechain[ct0.Level()]

	// TODO : manage scale dynamicly depending on Q_0, the Qi of the SineEval and the ciphertext's scale.
	ct0.MulScale(1024)
	// Sine Evaluation ct0 = Q/(2pi) * sin((2pi/Q) * ct0)
	ct0 = bootcontext.evaluateChebyBoot(evaluator, ct0)

	ct0.DivScale(1024)

	if ct1 != nil {

		ct1.MulScale(1024)
		// Sine Evaluation ct0 = Q/(2pi) * sin((2pi/Q) * ct0)
		ct1 = bootcontext.evaluateChebyBoot(evaluator, ct1)
		ct1.DivScale(1024)
	}

	// Reference scale is changed back to the current ciphertext's scale.
	bootcontext.ckkscontext.scale = ct0.Scale()

	return ct0, ct1
}

func (bootcontext *BootContext) evaluateChebyBoot(evaluator *Evaluator, ct *Ciphertext) (res *Ciphertext) {

	// Chebyshev params
	a := bootcontext.chebycoeffs.a
	b := bootcontext.chebycoeffs.b
	degree := bootcontext.chebycoeffs.degree
	coeffs := bootcontext.chebycoeffs.coeffs

	// SubSum + CoeffsToSlots cancelling factor
	n := complex(float64(evaluator.ckkscontext.n), 0)

	C := make(map[uint64]*Ciphertext)
	C[1] = ct.CopyNew().Ciphertext()

	evaluator.MultByConst(C[1], 2/((b-a)*n), C[1])
	evaluator.AddConst(C[1], (-a-b)/(b-a), C[1])
	evaluator.Rescale(C[1], evaluator.ckkscontext.scale, C[1])

	M := uint64(bits.Len64(degree - 1))
	L := uint64(M >> 1)

	for i := uint64(2); i < (1<<L)+1; i++ {
		computePowerBasisCheby(i, C, evaluator, bootcontext.relinkey)
	}

	for i := L + 1; i < M; i++ {
		computePowerBasisCheby(1<<i, C, evaluator, bootcontext.relinkey)
	}

	return recurseCheby(degree, L, M, coeffs, C, evaluator, bootcontext.relinkey)
}

func (bootcontext *BootContext) multiplyByDiagMatrice(evaluator *Evaluator, vec *Ciphertext, plainVectors *dftvectors) (res *Ciphertext) {

	var N1 uint64

	res = bootcontext.ckkscontext.NewCiphertext(1, vec.Level(), vec.Scale())

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
	vec_rot := evaluator.RotateHoisted(vec, rotations, bootcontext.rotkeys)

	var tmp_vec, tmp *Ciphertext

	tmp_vec = bootcontext.ckkscontext.NewCiphertext(1, bootcontext.ckkscontext.levels-1, vec.Scale())
	tmp = bootcontext.ckkscontext.NewCiphertext(1, bootcontext.ckkscontext.levels-1, vec.Scale())

	for j := range index {

		tmp_vec.Value()[0].Zero()
		tmp_vec.Value()[1].Zero()

		for _, i := range index[j] {
			evaluator.MulRelin(vec_rot[uint64(i)], plainVectors.Vec[N1*j+uint64(i)], nil, tmp)
			evaluator.Add(tmp_vec, tmp, tmp_vec)
		}

		evaluator.RotateColumns(tmp_vec, N1*j, bootcontext.rotkeys, tmp)

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

		index += 1
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

		index += 1
	}

	return
}

func (bootcontext *BootContext) computePlaintextVectors() {

	roots := computeRoots(bootcontext.slots << 1)
	pow5 := make([]uint64, (bootcontext.slots<<1)+1)
	pow5[0] = 1
	for i := uint64(1); i < (bootcontext.slots<<1)+1; i++ {
		pow5[i] = pow5[i-1] * 5
		pow5[i] &= (bootcontext.slots << 2) - 1
	}

	// CoeffsToSlots vectors
	bootcontext.pDFTInv = make([]*dftvectors, bootcontext.ctsDepth)

	pVecDFTInv := bootcontext.computeDFTPlaintextVectors(roots, pow5, true)

	for i := uint64(0); i < bootcontext.ctsDepth; i++ {

		bootcontext.pDFTInv[i] = new(dftvectors)

		bootcontext.pDFTInv[i].N1 = findbestbabygiantstepsplit(pVecDFTInv[i], bootcontext.dslots)

		bootcontext.encodePVec(pVecDFTInv[i], bootcontext.pDFTInv[i], i, true)
	}

	// SlotsToCoeffs vectors
	bootcontext.pDFT = make([]*dftvectors, bootcontext.stcDepth)

	pVecDFT := bootcontext.computeDFTPlaintextVectors(roots, pow5, false)

	for i := uint64(0); i < bootcontext.stcDepth; i++ {

		bootcontext.pDFT[i] = new(dftvectors)

		bootcontext.pDFT[i].N1 = findbestbabygiantstepsplit(pVecDFT[i], bootcontext.dslots)

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
	N = bootcontext.ckkscontext.n
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

	for j := range index {

		for _, i := range index[j] {

			if forward {
				level = bootcontext.ckkscontext.Levels() - 1 - k
			} else {
				level = bootcontext.ckkscontext.Levels() - 1 - k - bootcontext.ctsDepth - bootcontext.sinDepth
			}

			plaintextVec.Vec[N1*j+uint64(i)] = bootcontext.ckkscontext.NewPlaintext(level, bootcontext.ckkscontext.scalechain[level])

			bootcontext.plaintextSize += (level + 1) * 8 * bootcontext.ckkscontext.n

			bootcontext.encoder.Encode(plaintextVec.Vec[N1*j+uint64(i)], rotate(pVec[N1*j+uint64(i)], (N>>1)-(N1*j))[:bootcontext.dslots], bootcontext.dslots)
		}
	}
}

func (bootcontext *BootContext) computeDFTPlaintextVectors(roots []complex128, pow5 []uint64, forward bool) (plainVector []map[uint64][]complex128) {

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
				nextLevel -= 1
			}

		} else {
			// First layer of the i-th level of the DFT
			plainVector[i] = genWfft(logSlots, level, a[logSlots-level], b[logSlots-level], c[logSlots-level], forward)

			// Merges the layer with the next levels of the DFT if the total depth requires it.
			nextLevel = level - 1
			for j := uint64(0); j < merge[i]-1; j++ {
				plainVector[i] = nextLevelfft(plainVector[i], logSlots, 1<<logSlots, nextLevel, a[logSlots-nextLevel], b[logSlots-nextLevel], c[logSlots-nextLevel], forward)
				nextLevel -= 1
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
