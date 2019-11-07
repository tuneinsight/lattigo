package ckks

import (
	"errors"
	"fmt"
	"github.com/ldsec/lattigo/ring"
	"math"
	"math/bits"
	"math/cmplx"
)

type BootContext struct {
	bootcontext *CkksContext
	ckkscontext *CkksContext

	encoder *Encoder
	slots   uint64
	dslots  uint64

	//Sine evaluation
	chebycoeffs *ChebyshevInterpolation

	plaintextSize uint64

	//Coeffs to slots and slots to coeffs
	ctsDepth     uint64
	stcDepth     uint64
	sinDepth     uint64
	repack       bool
	repackVecStC *Plaintext
	repackVecCtS *Plaintext
	pDFT         []*dftvectors
	pDFTInv      []*dftvectors

	relinkey *EvaluationKey
	rotkeys  *RotationKey

	ctxpool [3]*Ciphertext

	ptxpool *Plaintext
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

func (ckkscontext *CkksContext) NewBootContext(slots uint64, sk *SecretKey, ctsDepth, stcDepth uint64, repack bool) (bootcontext *BootContext, err error) {

	bootcontext = new(BootContext)

	bootcontext.ckkscontext = ckkscontext

	bootcontext.ctsDepth = ctsDepth
	bootcontext.stcDepth = stcDepth
	if slots == ckkscontext.maxSlots && repack {
		return nil, errors.New("invalide bootcontext -> cannot repack if ciphertext is already fully packed")
	}
	bootcontext.repack = repack

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

	if bootcontext.repack {

		pVec := make([]complex128, bootcontext.dslots)

		vecLevel := ckkscontext.Levels() - 1 - ctsDepth

		bootcontext.repackVecCtS = bootcontext.ckkscontext.NewPlaintext(vecLevel, ckkscontext.scalechain[vecLevel])

		bootcontext.plaintextSize += (vecLevel + 1) * 8 * bootcontext.ckkscontext.n

		// Repacking from Y^(N/n) to X (we multiply with a vector [1, 1, ..., 1, 1, 0, 0, ..., 0, 0])
		// TODO : include in the DFT to save a level.
		for i := uint64(0); i < bootcontext.slots; i++ {
			pVec[i] = complex(1, 0)
			pVec[i+bootcontext.slots] = complex(0, 0)
		}

		bootcontext.encoder.Encode(bootcontext.repackVecCtS, pVec, bootcontext.dslots)

		vecLevel = ckkscontext.Levels() - 1 - ctsDepth - 1 - bootcontext.sinDepth

		bootcontext.repackVecStC = bootcontext.ckkscontext.NewPlaintext(vecLevel, ckkscontext.scalechain[vecLevel])

		bootcontext.plaintextSize += (vecLevel + 1) * 8 * bootcontext.ckkscontext.n

		// Repacking from X to Y^(N/n) : ct0 = ct0 + rotate(ct0 * [1, 1, ..., 1, 1, -i, -i, ..., -i, -i], slots).
		// TODO : include in the DFT to save a level.
		for i := uint64(0); i < bootcontext.slots; i++ {
			pVec[i] = complex(1, 0)
			pVec[i+bootcontext.slots] = complex(0, -1)
		}

		bootcontext.encoder.Encode(bootcontext.repackVecStC, pVec, bootcontext.dslots)
	}

	// List of the rotation key values to needed for the bootstrapp
	rotations := []uint64{}

	//SubSum rotation needed X -> Y^slots rotations
	logSlots := uint64(bits.Len64(bootcontext.slots) - 1)
	logMaxSlots := uint64(bits.Len64(bootcontext.ckkscontext.maxSlots) - 1)
	for i := logSlots; i < logMaxSlots; i++ {
		if !IsInSlice(1<<i, rotations) {
			rotations = append(rotations, 1<<i)
		}
	}

	bootcontext.computePlaintextVectors()

	var index uint64

	// Coeffs to Slots rotations
	for i := range bootcontext.pDFTInv {
		for j := range bootcontext.pDFTInv[i].Vec {

			index = ((j / bootcontext.pDFTInv[i].N1) * bootcontext.pDFTInv[i].N1) & (slots - 1)

			if !IsInSlice(index, rotations) {
				rotations = append(rotations, index)
			}

			index = j % bootcontext.pDFTInv[i].N1

			if !IsInSlice(j%bootcontext.pDFTInv[i].N1, rotations) {
				rotations = append(rotations, index)
			}
		}
	}

	// Slots to Coeffs rotations
	for i := range bootcontext.pDFT {
		for j := range bootcontext.pDFT[i].Vec {

			index = ((j / bootcontext.pDFT[i].N1) * bootcontext.pDFT[i].N1) & (slots - 1)

			if !IsInSlice(index, rotations) {
				rotations = append(rotations, ((j/bootcontext.pDFT[i].N1)*bootcontext.pDFT[i].N1)&(slots-1))
			}

			index = j % bootcontext.pDFT[i].N1

			if !IsInSlice(index, rotations) {
				rotations = append(rotations, index)
			}
		}
	}

	fmt.Println("DFT vector size (GB) :", float64(bootcontext.plaintextSize)/float64(1000000000))
	fmt.Println("Switching-Keys size (GB) :", float64(ckkscontext.n*2*uint64(len(rotations))*ckkscontext.Beta()*uint64(len(ckkscontext.contextKeys.Modulus))*8)/float64(1000000000), "(", len(rotations), "keys)")

	kgen := ckkscontext.NewKeyGenerator()

	if bootcontext.rotkeys, err = kgen.NewRotationKeys(sk, rotations, nil, true); err != nil {
		return nil, err
	}

	bootcontext.relinkey = kgen.NewRelinKey(sk)

	bootcontext.ctxpool[0] = ckkscontext.NewCiphertext(1, ckkscontext.levels-1, 0)
	bootcontext.ctxpool[1] = ckkscontext.NewCiphertext(1, ckkscontext.levels-1, 0)
	bootcontext.ctxpool[2] = ckkscontext.NewCiphertext(1, ckkscontext.levels-1, 0)

	bootcontext.ptxpool = ckkscontext.NewPlaintext(ckkscontext.levels-1, 0)

	return bootcontext, nil
}

func (evaluator *Evaluator) Bootstrapp(ct *Ciphertext, bootcontext *BootContext) (*Ciphertext, error) {

	var err error
	var ct0, ct1 *Ciphertext

	for ct.Level() != 0 {
		evaluator.DropLevel(ct.Element(), 1)
	}

	evaluator.ScaleUp(ct, math.Round(float64(1<<45)/ct.Scale()), ct)

	ct = bootcontext.modUp(ct)

	//SubSum X -> (N/n) * Y^n
	logSlots := uint64(bits.Len64(bootcontext.slots) - 1)
	logMaxSlots := uint64(bits.Len64(bootcontext.ckkscontext.maxSlots) - 1)
	for i := logSlots; i < logMaxSlots; i++ {

		if err = evaluator.RotateColumns(ct, 1<<i, bootcontext.rotkeys, bootcontext.ctxpool[0]); err != nil {
			return nil, err
		}

		if err = evaluator.Add(ct, bootcontext.ctxpool[0], ct); err != nil {
			return nil, err
		}
	}

	// Part 1 : Coeffs to slots
	if ct0, ct1, err = bootcontext.coeffsToSlots(evaluator, ct); err != nil {
		return nil, err
	}

	// Part 2 : homomorphic evaluation of the modulo function
	if ct0, ct1, err = bootcontext.evaluateSine(ct0, ct1, evaluator); err != nil {
		return nil, err
	}

	// Part 3 : Slots to coeffs
	return bootcontext.slotsToCoeffs(evaluator, ct0, ct1)
}

func (bootcontext *BootContext) modUp(ct *Ciphertext) *Ciphertext {

	ct.InvNTT(bootcontext.ckkscontext, ct.Element())

	for u := range ct.Value() {
		ct.Value()[u].Coeffs = append(ct.Value()[u].Coeffs, make([][]uint64, bootcontext.ckkscontext.levels-1)...)
		for i := uint64(1); i < bootcontext.ckkscontext.levels; i++ {
			ct.Value()[u].Coeffs[i] = make([]uint64, bootcontext.ckkscontext.n)
		}
	}

	ct.SetCurrentModulus(bootcontext.ckkscontext.bigintChain[bootcontext.ckkscontext.levels-1])

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

func (bootcontext *BootContext) coeffsToSlots(evaluator *Evaluator, vec *Ciphertext) (ct0, ct1 *Ciphertext, err error) {

	var zV, zVconj *Ciphertext

	if zV, err = bootcontext.dft(evaluator, vec, bootcontext.pDFTInv); err != nil {
		return nil, nil, err
	}

	if bootcontext.repack {

		if err = evaluator.MulRelin(zV, bootcontext.repackVecCtS, nil, zV); err != nil {
			return nil, nil, err
		}

		if err = evaluator.Rescale(zV, bootcontext.ckkscontext.scale, zV); err != nil {
			return nil, nil, err
		}
	}

	// Extraction of real and imaginary parts
	if zVconj, err = evaluator.ConjugateNew(zV, bootcontext.rotkeys); err != nil {
		return nil, nil, err
	}

	if ct0, err = evaluator.AddNew(zV, zVconj); err != nil {
		return nil, nil, err
	}

	if ct1, err = evaluator.SubNew(zV, zVconj); err != nil {
		return nil, nil, err
	}

	if err = evaluator.DivByi(ct1, ct1); err != nil {
		return nil, nil, err
	}

	if bootcontext.repack {

		// Real part is put in the first n/2 slots, imaginary part is put in the last n/2 slots
		if err = evaluator.RotateColumns(ct1, bootcontext.slots, bootcontext.rotkeys, ct1); err != nil {
			return nil, nil, err
		}

		if err = evaluator.Add(ct0, ct1, ct0); err != nil {
			return nil, nil, err
		}
	}

	return ct0, ct1, nil
}

func (bootcontext *BootContext) slotsToCoeffs(evaluator *Evaluator, ct0, ct1 *Ciphertext) (ct *Ciphertext, err error) {

	if bootcontext.repack {

		if evaluator.MulRelin(ct0, bootcontext.repackVecStC, nil, ct0); err != nil {
			return nil, err
		}

		if err = evaluator.RotateColumns(ct0, bootcontext.slots, bootcontext.rotkeys, ct1); err != nil {
			return nil, err
		}

		if err = evaluator.Add(ct0, ct1, ct0); err != nil {
			return nil, err
		}

		if err = evaluator.Rescale(ct0, bootcontext.ckkscontext.scale, ct0); err != nil {
			return nil, err
		}

	} else {

		if err = evaluator.DivByi(ct1, ct1); err != nil {
			return nil, err
		}

		if err = evaluator.Add(ct0, ct1, ct0); err != nil {
			return nil, err
		}
	}

	return bootcontext.dft(evaluator, ct0, bootcontext.pDFT)
}

func (bootcontext *BootContext) evaluateSine(ct0, ct1 *Ciphertext, evaluator *Evaluator) (*Ciphertext, *Ciphertext, error) {

	var err error

	// Sine Evaluation ct0 = Q/(2pi) * sin((2pi/Q) * ct0)
	bootcontext.ckkscontext.scale = bootcontext.ckkscontext.scalechain[ct0.Level()]

	ct0.MulScale(1024)
	if ct0, err = evaluator.EvaluateChebyBoot(ct0, bootcontext.chebycoeffs, bootcontext.relinkey); err != nil {
		return nil, nil, err
	}
	ct0.DivScale(1024)

	if !bootcontext.repack {

		ct1.MulScale(1024)
		if ct1, err = evaluator.EvaluateChebyBoot(ct1, bootcontext.chebycoeffs, bootcontext.relinkey); err != nil {
			return nil, nil, err
		}
		ct1.DivScale(1024)
	}

	bootcontext.ckkscontext.scale = ct0.Scale()

	return ct0, ct1, nil
}

func (evaluator *Evaluator) EvaluateChebyBoot(ct *Ciphertext, cheby *ChebyshevInterpolation, evakey *EvaluationKey) (res *Ciphertext, err error) {

	C := make(map[uint64]*Ciphertext)

	C[1] = ct.CopyNew().Ciphertext()

	evaluator.MultConst(C[1], 2/((cheby.b-cheby.a)*complex(float64(evaluator.ckkscontext.n), 0)), C[1])
	evaluator.AddConst(C[1], (-cheby.a-cheby.b)/(cheby.b-cheby.a), C[1])
	evaluator.Rescale(C[1], evaluator.ckkscontext.scale, C[1])

	C[2], _ = evaluator.MulRelinNew(C[1], C[1], evakey)
	evaluator.Rescale(C[2], evaluator.ckkscontext.scale, C[2])

	evaluator.Add(C[2], C[2], C[2])
	evaluator.AddConst(C[2], -1, C[2])

	M := uint64(bits.Len64(cheby.degree - 1))
	L := uint64(M >> 1)

	for i := uint64(3); i < (1<<L)+1; i++ {
		if err = computePowerBasis(i, C, evaluator, evakey); err != nil {
			return nil, err
		}
	}

	for i := L + 1; i < M; i++ {
		if err = computePowerBasis(1<<i, C, evaluator, evakey); err != nil {
			return nil, err
		}
	}

	return recurse(cheby.degree, L, M, cheby.coeffs, C, evaluator, evakey)
}

func (bootcontext *BootContext) dft(evaluator *Evaluator, vec *Ciphertext, plainVectors []*dftvectors) (w *Ciphertext, err error) {

	levels := uint64(len(plainVectors))

	w = vec.CopyNew().Ciphertext()

	for i := uint64(0); i < levels; i++ {

		if w, err = bootcontext.multiplyByDiagMatrice(evaluator, w, plainVectors[i]); err != nil {
			return nil, err
		}

		evaluator.Rescale(w, evaluator.ckkscontext.scale, w)
	}

	return w, nil
}

func (bootcontext *BootContext) multiplyByDiagMatrice(evaluator *Evaluator, vec *Ciphertext, plainVectors *dftvectors) (res *Ciphertext, err error) {

	var N1 uint64

	res = bootcontext.ckkscontext.NewCiphertext(1, vec.Level(), vec.Scale())

	// N1*N2 = N
	N1 = plainVectors.N1

	vec_rot := make(map[uint64]*Ciphertext, N1)

	index := make(map[uint64][]uint64)

	for key := range plainVectors.Vec {
		if index[key/N1] == nil {
			index[key/N1] = []uint64{key % N1}
		} else {
			index[key/N1] = append(index[key/N1], key%N1)
		}

		if vec_rot[key%N1] == nil {
			if vec_rot[key%N1], err = evaluator.RotateColumnsNew(vec, key%N1, bootcontext.rotkeys); err != nil {

				return nil, err
			}
		}
	}

	var tmp_vec, tmp *Ciphertext

	//plaintextVec := bootcontext.ckkscontext.NewPlaintext(vec.Level(), evaluator.ckkscontext.scalechain[vec.Level()])
	tmp_vec = bootcontext.ckkscontext.NewCiphertext(1, bootcontext.ckkscontext.levels-1, vec.Scale())
	tmp = bootcontext.ckkscontext.NewCiphertext(1, bootcontext.ckkscontext.levels-1, vec.Scale())

	for j := range index {

		tmp_vec.Value()[0].Zero()
		tmp_vec.Value()[1].Zero()

		for _, i := range index[j] {

			if err = evaluator.MulRelin(vec_rot[uint64(i)], plainVectors.Vec[N1*j+uint64(i)], nil, tmp); err != nil {
				return nil, err
			}

			if err = evaluator.Add(tmp_vec, tmp, tmp_vec); err != nil {
				return nil, err
			}
		}

		if err = evaluator.RotateColumns(tmp_vec, N1*j, bootcontext.rotkeys, tmp); err != nil {
			return nil, err
		}

		if err = evaluator.Add(res, tmp, res); err != nil {
			return nil, err
		}
	}

	return res, nil
}

func computeRoots(N uint64) (roots, rootsInv []complex128) {
	m := N << 1

	roots = make([]complex128, m)
	rootsInv = make([]complex128, m)

	angle := 6.283185307179586 / float64(m)

	psi := complex(math.Cos(angle), math.Sin(angle))

	roots[0] = 1
	rootsInv[0] = 1
	for i := uint64(1); i < m; i++ {
		roots[i] = roots[i-1] * psi
		rootsInv[i] = 1.0 / roots[i]
	}

	return
}

func fftPlainVec(N uint64, roots []complex128) (a, b, c [][]complex128) {

	var logN, m, index, tt, gap, k, mask uint64

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

				k = 1
				for u := uint64(0); u < j; u++ {
					k *= 5
					k &= mask
				}
				k *= gap

				a[index][i+j] = 1
				a[index][i+j+tt] = -roots[k]
				b[index][i+j] = roots[k]
				c[index][i+j+tt] = 1

				a[index][N+i+j] = 1
				a[index][N+i+j+tt] = -roots[k]
				b[index][N+i+j] = roots[k]
				c[index][N+i+j+tt] = 1

			}
		}

		index += 1
	}

	return
}

func fftInvPlainVec(N uint64, rootsInv []complex128) (a, b, c [][]complex128) {

	var logN, m, index, tt, gap, k, mask uint64

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

				k = 1
				for u := uint64(0); u < j; u++ {
					k *= 5
					k &= mask
				}
				k = (((m << 2) - k) * gap)

				a[index][i+j] = 1
				a[index][i+j+tt] = -rootsInv[k]
				b[index][i+j] = 1
				c[index][i+j+tt] = rootsInv[k]

				a[index][N+i+j] = 1
				a[index][N+i+j+tt] = -rootsInv[k]
				b[index][N+i+j] = 1
				c[index][N+i+j+tt] = rootsInv[k]

			}
		}

		index += 1
	}

	return
}

func (bootcontext *BootContext) computePlaintextVectors() (err error) {
	var N, logL, L uint64
	var a, b, c [][]complex128

	N = bootcontext.slots << 1
	logL = uint64(bits.Len64(bootcontext.slots) - 1)
	L = bootcontext.slots

	roots, rootsInv := computeRoots(N)

	bootcontext.pDFTInv = make([]*dftvectors, bootcontext.ctsDepth)

	a, b, c = fftInvPlainVec(L, roots)
	pVecDFTInv := computeDFTPlaintextVectors(logL, bootcontext.ctsDepth, a, b, c, true)

	for i := uint64(0); i < bootcontext.ctsDepth; i++ {

		bootcontext.pDFTInv[i] = new(dftvectors)

		bootcontext.pDFTInv[i].N1 = findbestbabygiantstepsplit(pVecDFTInv[i], bootcontext.dslots)

		if err = bootcontext.encodePVec(pVecDFTInv[i], bootcontext.pDFTInv[i], i, true); err != nil {
			return err
		}
	}

	bootcontext.pDFT = make([]*dftvectors, bootcontext.stcDepth)

	a, b, c = fftPlainVec(L, rootsInv)
	pVecDFT := computeDFTPlaintextVectors(logL, bootcontext.stcDepth, a, b, c, false)

	for i := uint64(0); i < bootcontext.stcDepth; i++ {

		bootcontext.pDFT[i] = new(dftvectors)

		bootcontext.pDFT[i].N1 = findbestbabygiantstepsplit(pVecDFT[i], bootcontext.dslots)

		if err = bootcontext.encodePVec(pVecDFT[i], bootcontext.pDFT[i], i, false); err != nil {
			return err
		}
	}

	return
}

func findbestbabygiantstepsplit(vector map[uint64][]complex128, maxN uint64) (minN uint64) {

	var sum uint64

	sum = maxN

	for N1 := uint64(1); N1 < maxN; N1 <<= 1 {

		index := make(map[uint64][]uint64)

		for key := range vector {
			if index[key/N1] == nil {
				index[key/N1] = []uint64{key % N1}
			} else {
				index[key/N1] = append(index[key/N1], key%N1)
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

func (bootcontext *BootContext) encodePVec(pVec map[uint64][]complex128, plaintextVec *dftvectors, k uint64, forward bool) (err error) {
	var N, N1, level uint64

	// N1*N2 = N
	N = bootcontext.ckkscontext.n
	N1 = plaintextVec.N1

	index := make(map[uint64][]uint64)

	for key := range pVec {
		if index[key/N1] == nil {
			index[key/N1] = []uint64{key % N1}
		} else {
			index[key/N1] = append(index[key/N1], key%N1)
		}
	}

	plaintextVec.Vec = make(map[uint64]*Plaintext)

	for j := range index {

		for _, i := range index[j] {

			if forward {
				level = bootcontext.ckkscontext.Levels() - 1 - k
			} else {
				level = bootcontext.ckkscontext.Levels() - 1 - k - bootcontext.ctsDepth - bootcontext.sinDepth

				// We take into account the levels needed for repacking
				if bootcontext.repack {
					level -= 2
				}
			}

			plaintextVec.Vec[N1*j+uint64(i)] = bootcontext.ckkscontext.NewPlaintext(level, bootcontext.ckkscontext.scalechain[level])

			bootcontext.plaintextSize += (level + 1) * 8 * bootcontext.ckkscontext.n

			if err = bootcontext.encoder.Encode(plaintextVec.Vec[N1*j+uint64(i)], rotate(pVec[N1*j+uint64(i)], (N>>1)-(N1*j))[:bootcontext.dslots], bootcontext.dslots); err != nil {
				return err
			}
		}
	}

	return
}

func computeDFTPlaintextVectors(logL, maxDepth uint64, a, b, c [][]complex128, forward bool) (plainVector []map[uint64][]complex128) {

	var level, depth, nextLevel uint64

	plainVector = make([]map[uint64][]complex128, maxDepth)

	level = logL

	// We compute the chain of merge in order or reverse order depending if its DFT or InvDFT because
	// the way the levels are collapsed has an inpact on the total number of rotations and keys to be
	// stored. Ex. instead of using 255 + 64 plaintext vectors, we can use 127 + 128 plaintext vectors.
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

	level = logL
	for i := uint64(0); i < maxDepth; i++ {

		plainVector[i] = genWfft(logL, level, a, b, c, forward)

		nextLevel = level - 1
		for j := uint64(0); j < merge[i]-1; j++ {
			plainVector[i] = nextLevelfft(plainVector[i], logL, nextLevel, a, b, c, forward)
			nextLevel -= 1
		}

		level -= merge[i]
	}

	return
}

func genWfft(logL, level uint64, a, b, c [][]complex128, forward bool) (vectors map[uint64][]complex128) {

	var rot uint64

	if forward {
		rot = 1 << (level - 1)
	} else {
		rot = 1 << (logL - level)
	}

	vectors = make(map[uint64][]complex128)

	addToDicVector(vectors, 0, a[logL-level])
	addToDicVector(vectors, rot, b[logL-level])
	addToDicVector(vectors, (1<<logL)-rot, c[logL-level])

	return
}

func nextLevelfft(vec map[uint64][]complex128, logL, nextLevel uint64, a, b, c [][]complex128, forward bool) (newVec map[uint64][]complex128) {

	var N, rot uint64

	N = 1 << logL

	newVec = make(map[uint64][]complex128)

	if forward {
		rot = (1 << (nextLevel - 1)) & (N - 1)
	} else {
		rot = (1 << (logL - nextLevel)) & (N - 1)
	}

	for i := range vec {
		addToDicVector(newVec, i, mul(vec[i], a[logL-nextLevel]))
		addToDicVector(newVec, (i+rot)&(N-1), mul(rotate(vec[i], rot), b[logL-nextLevel]))
		addToDicVector(newVec, (i+N-rot)&(N-1), mul(rotate(vec[i], N-rot), c[logL-nextLevel]))
	}

	return
}

// Rotates to the left
func rotate(x []complex128, n uint64) (y []complex128) {

	y = make([]complex128, len(x))

	mask := uint64(len(x) - 1)

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

func addToDicVector(dic map[uint64][]complex128, index uint64, vec []complex128) {
	if dic[index] == nil {
		dic[index] = vec
	} else {
		dic[index] = add(dic[index], vec)
	}
}
