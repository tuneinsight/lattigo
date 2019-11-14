package ckks

import (
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
	ctsDepth uint64
	stcDepth uint64
	sinDepth uint64
	repack   bool
	pDFT     []*dftvectors
	pDFTInv  []*dftvectors

	relinkey *EvaluationKey
	rotkeys  *RotationKey

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

func (ckkscontext *CkksContext) NewBootContext(slots uint64, sk *SecretKey, ctsDepth, stcDepth uint64) (bootcontext *BootContext, err error) {

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
		if !IsInSlice(1<<i, rotations) {
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

			if !IsInSlice(index, rotations) {
				rotations = append(rotations, index)
			}

			index = j % bootcontext.pDFTInv[i].N1

			if !IsInSlice(index, rotations) {
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

			if !IsInSlice(index, rotations) {
				rotations = append(rotations, index)
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

	return bootcontext, nil
}

func (evaluator *Evaluator) Bootstrapp(ct *Ciphertext, bootcontext *BootContext) (*Ciphertext, error) {

	var err error
	var ct0, ct1 *Ciphertext

	for ct.Level() != 0 {
		evaluator.DropLevel(ct.Element(), 1)
	}

	// TODO : better management of the initial scale
	evaluator.ScaleUp(ct, math.Round(float64(1<<45)/ct.Scale()), ct)

	// ModUp ct_{Q_0} -> ct_{Q_L}
	ct = bootcontext.modUp(ct)

	//SubSum X -> (N/dslots) * Y^dslots
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

	// Part 2 : SineEval
	if ct0, ct1, err = bootcontext.evaluateSine(ct0, ct1, evaluator); err != nil {
		return nil, err
	}

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

	if zV, err = bootcontext.dft(evaluator, vec, bootcontext.pDFTInv, true); err != nil {
		return nil, nil, err
	}

	// Extraction of real and imaginary parts.
	if zVconj, err = evaluator.ConjugateNew(zV, bootcontext.rotkeys); err != nil {
		return nil, nil, err
	}

	// The real part is stored in ct0
	if ct0, err = evaluator.AddNew(zV, zVconj); err != nil {
		return nil, nil, err
	}

	// The imaginary part is stored in ct1
	if ct1, err = evaluator.SubNew(zV, zVconj); err != nil {
		return nil, nil, err
	}

	if err = evaluator.DivByi(ct1, ct1); err != nil {
		return nil, nil, err
	}

	// If repacking, then ct0 and ct1 right n/2 slots are zero.
	if bootcontext.repack {

		// The imaginary part is put in the right n/2 slots of ct0.
		if err = evaluator.RotateColumns(ct1, bootcontext.slots, bootcontext.rotkeys, ct1); err != nil {
			return nil, nil, err
		}

		if err = evaluator.Add(ct0, ct1, ct0); err != nil {
			return nil, nil, err
		}

		return ct0, nil, nil
	}

	return ct0, ct1, nil
}

func (bootcontext *BootContext) slotsToCoeffs(evaluator *Evaluator, ct0, ct1 *Ciphertext) (ct *Ciphertext, err error) {

	// If full packing, the repacking can be done directly using ct0 and ct1.
	if !bootcontext.repack {

		if err = evaluator.MultByi(ct1, ct1); err != nil {
			return nil, err
		}

		if err = evaluator.Add(ct0, ct1, ct0); err != nil {
			return nil, err
		}
	}

	return bootcontext.dft(evaluator, ct0, bootcontext.pDFT, false)
}

func (bootcontext *BootContext) dft(evaluator *Evaluator, vec *Ciphertext, plainVectors []*dftvectors, forward bool) (w *Ciphertext, err error) {

	levels := uint64(len(plainVectors))

	w = vec.CopyNew().Ciphertext()

	// Sequencially multiplies w with the provided dft matrices.
	for i := uint64(0); i < levels; i++ {

		if w, err = bootcontext.multiplyByDiagMatrice(evaluator, w, plainVectors[i]); err != nil {
			return nil, err
		}

		evaluator.Rescale(w, evaluator.ckkscontext.scale, w)
	}

	return w, nil
}

func (bootcontext *BootContext) evaluateSine(ct0, ct1 *Ciphertext, evaluator *Evaluator) (*Ciphertext, *Ciphertext, error) {

	var err error

	// Reference scale is changed to the new ciphertext's scale.
	bootcontext.ckkscontext.scale = bootcontext.ckkscontext.scalechain[ct0.Level()]

	// TODO : manage scale dynamicly depending on Q_0, the Qi of the SineEval and the ciphertext's scale.
	ct0.MulScale(1024)
	// Sine Evaluation ct0 = Q/(2pi) * sin((2pi/Q) * ct0)
	if ct0, err = bootcontext.evaluateChebyBoot(evaluator, ct0); err != nil {
		return nil, nil, err
	}
	ct0.DivScale(1024)

	if ct1 != nil {

		ct1.MulScale(1024)
		// Sine Evaluation ct0 = Q/(2pi) * sin((2pi/Q) * ct0)
		if ct1, err = bootcontext.evaluateChebyBoot(evaluator, ct1); err != nil {
			return nil, nil, err
		}
		ct1.DivScale(1024)
	}

	// Reference scale is changed back to the current ciphertext's scale.
	bootcontext.ckkscontext.scale = ct0.Scale()

	return ct0, ct1, nil
}

func (bootcontext *BootContext) evaluateChebyBoot(evaluator *Evaluator, ct *Ciphertext) (res *Ciphertext, err error) {

	// Chebyshev params
	a := bootcontext.chebycoeffs.a
	b := bootcontext.chebycoeffs.b
	degree := bootcontext.chebycoeffs.degree
	coeffs := bootcontext.chebycoeffs.coeffs

	// SubSum + CoeffsToSlots cancelling factor
	n := complex(float64(evaluator.ckkscontext.n), 0)

	C := make(map[uint64]*Ciphertext)
	C[1] = ct.CopyNew().Ciphertext()

	evaluator.MultConst(C[1], 2/((b-a)*n), C[1])
	evaluator.AddConst(C[1], (-a-b)/(b-a), C[1])
	evaluator.Rescale(C[1], evaluator.ckkscontext.scale, C[1])

	C[2], _ = evaluator.MulRelinNew(C[1], C[1], bootcontext.relinkey)
	evaluator.Rescale(C[2], evaluator.ckkscontext.scale, C[2])

	evaluator.Add(C[2], C[2], C[2])
	evaluator.AddConst(C[2], -1, C[2])

	M := uint64(bits.Len64(degree - 1))
	L := uint64(M >> 1)

	for i := uint64(3); i < (1<<L)+1; i++ {
		if err = computePowerBasis(i, C, evaluator, bootcontext.relinkey); err != nil {
			return nil, err
		}
	}

	for i := L + 1; i < M; i++ {
		if err = computePowerBasis(1<<i, C, evaluator, bootcontext.relinkey); err != nil {
			return nil, err
		}
	}

	return recurse(degree, L, M, coeffs, C, evaluator, bootcontext.relinkey)
}

func (bootcontext *BootContext) multiplyByDiagMatrice(evaluator *Evaluator, vec *Ciphertext, plainVectors *dftvectors) (res *Ciphertext, err error) {

	var N1 uint64

	res = bootcontext.ckkscontext.NewCiphertext(1, vec.Level(), vec.Scale())

	// N1*N2 = N
	N1 = plainVectors.N1

	vec_rot := make(map[uint64]*Ciphertext, N1)

	// Pre-computation for rotations using hoisting
	contextQ := bootcontext.ckkscontext.contextQ
	contextP := bootcontext.ckkscontext.contextP

	c2NTT := vec.value[1]
	c2InvNTT := contextQ.NewPoly()
	contextQ.InvNTTLvl(vec.Level(), c2NTT, c2InvNTT)

	c2_qiQDecomp := make([]*ring.Poly, bootcontext.ckkscontext.beta)
	c2_qiPDecomp := make([]*ring.Poly, bootcontext.ckkscontext.beta)

	alpha := evaluator.ckkscontext.alpha
	beta := uint64(math.Ceil(float64(vec.Level()+1) / float64(alpha)))

	for i := uint64(0); i < beta; i++ {
		c2_qiQDecomp[i] = contextQ.NewPoly()
		c2_qiPDecomp[i] = contextP.NewPoly()
		evaluator.decomposeAndSplitNTT(vec.Level(), i, c2NTT, c2InvNTT, c2_qiQDecomp[i], c2_qiPDecomp[i])
	}

	// Computes the rotations indexes of the non-zero rows of the diagonalized DFT matrix for the baby-step giang-step algorithm
	index := make(map[uint64][]uint64)
	for key := range plainVectors.Vec {
		if index[key/N1] == nil {
			index[key/N1] = []uint64{key % N1}
		} else {
			index[key/N1] = append(index[key/N1], key%N1)
		}

		// Pre-rotates ciphertext for the baby-step giant-step algorithm
		if vec_rot[key%N1] == nil {

			// Rotation using hoisting
			if key%N1 == 0 {
				vec_rot[key%N1] = vec.CopyNew().Ciphertext()
			} else {
				vec_rot[key%N1] = bootcontext.ckkscontext.NewCiphertext(1, vec.Level(), vec.Scale())
				evaluator.rotateLeftHoisted(vec, c2_qiQDecomp, c2_qiPDecomp, key%N1, bootcontext.rotkeys, vec_rot[key%N1])
			}
		}
	}

	var tmp_vec, tmp *Ciphertext

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

func (evaluator *Evaluator) rotateLeftHoisted(ctIn *Ciphertext, c2_qiQDecomp, c2_qiPDecomp []*ring.Poly, k uint64, evakey *RotationKey, ctOut *Ciphertext) {

	var level, reduce uint64

	level = ctOut.Level()

	contextQ := evaluator.ckkscontext.contextQ
	contextP := evaluator.ckkscontext.contextP

	if ctIn != ctOut {
		ring.PermuteNTT(ctIn.value[0], evaluator.ckkscontext.galElRotColLeft[k], evaluator.ringpool[0])
		contextQ.CopyLvl(level, evaluator.ringpool[0], ctOut.value[0])
		ring.PermuteNTT(ctIn.value[1], evaluator.ckkscontext.galElRotColLeft[k], evaluator.ringpool[0])
		contextQ.CopyLvl(level, evaluator.ringpool[0], ctOut.value[1])
	} else {
		ring.PermuteNTT(ctIn.value[0], evaluator.ckkscontext.galElRotColLeft[k], ctOut.value[0])
		ring.PermuteNTT(ctIn.value[1], evaluator.ckkscontext.galElRotColLeft[k], ctOut.value[1])
	}

	for i := range evaluator.poolQ {
		evaluator.poolQ[i].Zero()
	}

	for i := range evaluator.poolP {
		evaluator.poolP[i].Zero()
	}

	c2_qiQPermute := evaluator.poolQ[0]
	c2_qiPPermute := evaluator.poolP[0]

	pool2Q := evaluator.poolQ[1]
	pool2P := evaluator.poolP[1]

	pool3Q := evaluator.poolQ[2]
	pool3P := evaluator.poolP[2]

	reduce = 0

	alpha := evaluator.ckkscontext.alpha
	beta := uint64(math.Ceil(float64(level+1) / float64(alpha)))

	// Key switching with crt decomposition for the Qi
	for i := uint64(0); i < beta; i++ {

		ring.PermuteNTT(c2_qiQDecomp[i], evaluator.ckkscontext.galElRotColLeft[k], c2_qiQPermute)
		ring.PermuteNTT(c2_qiPDecomp[i], evaluator.ckkscontext.galElRotColLeft[k], c2_qiPPermute)

		contextQ.MulCoeffsMontgomeryAndAddNoModLvl(level, evakey.evakey_rot_col_L[k].evakey[i][0], c2_qiQPermute, pool2Q)
		contextQ.MulCoeffsMontgomeryAndAddNoModLvl(level, evakey.evakey_rot_col_L[k].evakey[i][1], c2_qiQPermute, pool3Q)

		// We continue with the keyswitch primes.
		for j, keysindex := uint64(0), evaluator.ckkscontext.levels; j < uint64(len(evaluator.ckkscontext.specialprimes)); j, keysindex = j+1, keysindex+1 {

			pj := contextP.Modulus[j]
			mredParams := contextP.GetMredParams()[j]

			key0 := evakey.evakey_rot_col_L[k].evakey[i][0].Coeffs[keysindex]
			key1 := evakey.evakey_rot_col_L[k].evakey[i][1].Coeffs[keysindex]
			p2tmp := pool2P.Coeffs[j]
			p3tmp := pool3P.Coeffs[j]
			c2tmp := c2_qiPPermute.Coeffs[j]

			for y := uint64(0); y < contextP.N; y++ {
				p2tmp[y] += ring.MRed(key0[y], c2tmp[y], pj, mredParams)
				p3tmp[y] += ring.MRed(key1[y], c2tmp[y], pj, mredParams)
			}
		}

		if reduce&7 == 1 {
			contextQ.ReduceLvl(level, pool2Q, pool2Q)
			contextQ.ReduceLvl(level, pool3Q, pool3Q)
			contextP.Reduce(pool2P, pool2P)
			contextP.Reduce(pool3P, pool3P)
		}

		reduce++
	}

	if (reduce-1)&7 != 1 {
		contextQ.ReduceLvl(level, pool2Q, pool2Q)
		contextQ.ReduceLvl(level, pool3Q, pool3Q)
		contextP.Reduce(pool2P, pool2P)
		contextP.Reduce(pool3P, pool3P)
	}

	//Independant of context (parameter : level)
	// Computes pool2Q = pool2Q/pool2P and pool3Q = pool3Q/pool3P
	evaluator.baseconverter.ModDownSplitedNTT(contextQ, contextP, evaluator.ckkscontext.rescaleParamsKeys, level, pool2Q, pool2P, pool2Q, evaluator.keyswitchpool[0])
	evaluator.baseconverter.ModDownSplitedNTT(contextQ, contextP, evaluator.ckkscontext.rescaleParamsKeys, level, pool3Q, pool3P, pool3Q, evaluator.keyswitchpool[0])

	//Independant of context (parameter : level)
	contextQ.AddLvl(level, ctOut.value[0], pool2Q, ctOut.value[0])
	contextQ.AddLvl(level, ctOut.value[1], pool3Q, ctOut.value[1])
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

func fftInvPlainVec(N uint64, roots []complex128) (a, b, c [][]complex128) {

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
				a[index][i+j+tt] = -roots[k]
				b[index][i+j] = 1
				c[index][i+j+tt] = roots[k]

				a[index][N+i+j] = 1
				a[index][N+i+j+tt] = -roots[k]
				b[index][N+i+j] = 1
				c[index][N+i+j+tt] = roots[k]

			}
		}

		index += 1
	}

	return
}

func (bootcontext *BootContext) computePlaintextVectors() (err error) {

	roots := computeRoots(bootcontext.slots << 1)

	// CoeffsToSlots vectors
	bootcontext.pDFTInv = make([]*dftvectors, bootcontext.ctsDepth)

	pVecDFTInv := bootcontext.computeDFTPlaintextVectors(roots, true)

	for i := uint64(0); i < bootcontext.ctsDepth; i++ {

		bootcontext.pDFTInv[i] = new(dftvectors)

		bootcontext.pDFTInv[i].N1 = findbestbabygiantstepsplit(pVecDFTInv[i], bootcontext.dslots)

		if err = bootcontext.encodePVec(pVecDFTInv[i], bootcontext.pDFTInv[i], i, true); err != nil {
			return err
		}
	}

	// SlotsToCoeffs vectors
	bootcontext.pDFT = make([]*dftvectors, bootcontext.stcDepth)

	pVecDFT := bootcontext.computeDFTPlaintextVectors(roots, false)

	for i := uint64(0); i < bootcontext.stcDepth; i++ {

		bootcontext.pDFT[i] = new(dftvectors)

		bootcontext.pDFT[i].N1 = findbestbabygiantstepsplit(pVecDFT[i], bootcontext.dslots)

		if err = bootcontext.encodePVec(pVecDFT[i], bootcontext.pDFT[i], i, false); err != nil {
			return err
		}
	}

	return
}

// Finds the best N1*N2 = N for the baby-step giant-step algorithm for matrix multiplication.
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

func (bootcontext *BootContext) computeDFTPlaintextVectors(roots []complex128, forward bool) (plainVector []map[uint64][]complex128) {

	var level, depth, nextLevel, slots, logSlots uint64

	slots = bootcontext.slots

	logSlots = uint64(bits.Len64(bootcontext.slots) - 1)

	level = logSlots

	var a, b, c [][]complex128
	var maxDepth uint64

	if forward {
		maxDepth = bootcontext.ctsDepth
		a, b, c = fftInvPlainVec(slots, roots)
	} else {
		maxDepth = bootcontext.stcDepth
		a, b, c = fftPlainVec(slots, roots)
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
