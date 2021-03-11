package ckks

import (
	"fmt"
	"math"
	"math/cmplx"

	"github.com/ldsec/lattigo/v2/ckks/bettersine"
	"github.com/ldsec/lattigo/v2/utils"
)

// Bootstrapper is a struct to stores a memory pool the plaintext matrices
// the polynomial approximation and the keys for the bootstrapping.
type Bootstrapper struct {
	*evaluator
	BootstrappingParameters
	*BootstrappingKey
	params *Parameters

	dslots    int // Number of plaintext slots after the re-encoding
	logdslots int

	encoder Encoder // Encoder

	plaintextSize int // Byte size of the plaintext DFT matrices

	repack       bool                    // If true then can repack the CoeffsToSlots into on ciphertext
	prescale     float64                 // Q[0]/(Q[0]/|m|)
	postscale    float64                 // Qi sineeval/(Q[0]/|m|)
	sinescale    float64                 // Qi sineeval
	sqrt2pi      float64                 // (1/2pi)^{-2^r}
	scFac        float64                 // 2^{r}
	sineEvalPoly *ChebyshevInterpolation // Coefficients of the Chebyshev Interpolation of sin(2*pi*x) or cos(2*pi*x/r)
	arcSinePoly  *Poly                   // Coefficients of the Taylor series of arcsine(x)

	coeffsToSlotsDiffScale complex128      // Matrice rescaling
	slotsToCoeffsDiffScale complex128      // Matrice rescaling
	pDFT                   []*PtDiagMatrix // Matrice vectors
	pDFTInv                []*PtDiagMatrix // Matrice vectors
	ctsLevel               []int           // index of the Qi used by CoeffsToSlots
	stcLevel               []int           // index of the Qi used by SlotsToCoeffs

	rotKeyIndex []int // a list of the required rotation keys
}

func sin2pi2pi(x complex128) complex128 {
	return cmplx.Sin(6.283185307179586*x) / 6.283185307179586
}

func cos2pi(x complex128) complex128 {
	return cmplx.Cos(6.283185307179586 * x)
}

// NewBootstrapper creates a new Bootstrapper.
func NewBootstrapper(params *Parameters, btpParams *BootstrappingParameters, btpKey BootstrappingKey) (btp *Bootstrapper, err error) {

	if btpParams.SinType == SinType(Sin) && btpParams.SinRescal != 0 {
		return nil, fmt.Errorf("cannot use double angle formul for SinType = Sin -> must use SinType = Cos")
	}

	btp = newBootstrapper(params, btpParams)

	btp.BootstrappingKey = &BootstrappingKey{btpKey.Rlk, btpKey.Rtks}
	if err = btp.CheckKeys(); err != nil {
		return nil, fmt.Errorf("invalid bootstrapping key: %w", err)
	}

	btp.evaluator = btp.evaluator.WithKey(EvaluationKey{btpKey.Rlk, btpKey.Rtks}).(*evaluator)

	return btp, nil
}

// newBootstrapper is a constructor of "dummy" bootstrapper to enable the generation of bootstrapping-related constants
// without providing a bootstrapping key. To be replaced by a propper factorization of the bootstrapping pre-computations.
func newBootstrapper(params *Parameters, btpParams *BootstrappingParameters) (btp *Bootstrapper) {
	btp = new(Bootstrapper)

	btp.params = params.Copy()
	btp.BootstrappingParameters = *btpParams.Copy()

	btp.dslots = params.Slots()
	btp.logdslots = params.LogSlots()
	if params.logSlots < params.MaxLogSlots() {
		btp.repack = true
		btp.dslots <<= 1
		btp.logdslots++
	}

	btp.prescale = math.Exp2(math.Round(math.Log2(float64(params.qi[0]) / btp.MessageRatio)))
	btp.sinescale = math.Exp2(math.Round(math.Log2(btp.SineEvalModuli.ScalingFactor)))
	btp.postscale = btp.sinescale / btp.MessageRatio

	btp.ctsLevel = []int{}
	for i := range btp.CoeffsToSlotsModuli.Qi {
		for range btp.CoeffsToSlotsModuli.ScalingFactor[btp.CtSDepth(true)-1-i] {
			btp.ctsLevel = append(btp.ctsLevel, btp.params.MaxLevel()-i)
		}
	}

	btp.stcLevel = []int{}
	for i := range btp.SlotsToCoeffsModuli.Qi {
		for range btp.SlotsToCoeffsModuli.ScalingFactor[btp.StCDepth(true)-1-i] {
			btp.stcLevel = append(btp.stcLevel, btp.params.MaxLevel()-btp.CtSDepth(true)-btp.SineEvalDepth(true)-btp.ArcSineDepth()-i)
		}
	}

	btp.encoder = NewEncoder(params)
	btp.evaluator = NewEvaluator(params, EvaluationKey{}).(*evaluator) // creates an evaluator without keys for genDFTMatrices

	btp.genSinePoly()
	btp.genDFTMatrices()

	btp.ctxpool = NewCiphertext(params, 1, params.MaxLevel(), 0)

	return btp
}

// CheckKeys checks if all the necessary keys are present
func (btp *Bootstrapper) CheckKeys() (err error) {

	if btp.Rlk == nil {
		return fmt.Errorf("relinearization key is nil")
	}

	if btp.Rtks == nil {
		return fmt.Errorf("rotation key is nil")
	}

	rotMissing := []int{}
	for _, i := range btp.rotKeyIndex {
		galEl := btp.params.GaloisElementForColumnRotationBy(int(i))
		if _, generated := btp.Rtks.Keys[galEl]; !generated {
			rotMissing = append(rotMissing, i)
		}
	}

	if len(rotMissing) != 0 {
		return fmt.Errorf("rotation key(s) missing: %d", rotMissing)
	}

	return nil
}

func (btp *Bootstrapper) addMatrixRotToList(pVec *PtDiagMatrix, rotations []int, slots int, repack bool) {

	var index int
	for j := range pVec.Vec {

		N1 := pVec.N1

		index = ((j / N1) * N1)

		if repack {
			// Sparse repacking, occurring during the first DFT matrix of the CoeffsToSlots.
			index &= 2*slots - 1
		} else {
			// Other cases
			index &= slots - 1
		}

		if index != 0 && !utils.IsInSliceInt(index, rotations) {
			rotations = append(rotations, index)
		}

		index = j & (N1 - 1)

		if index != 0 && !utils.IsInSliceInt(index, rotations) {
			rotations = append(rotations, index)
		}
	}
}

func (btp *Bootstrapper) genDFTMatrices() {

	a := real(btp.sineEvalPoly.a)
	b := real(btp.sineEvalPoly.b)
	n := float64(btp.params.N())
	qDiff := float64(btp.params.qi[0]) / math.Exp2(math.Round(math.Log2(float64(btp.params.qi[0]))))

	// Change of variable for the evaluation of the Chebyshev polynomial + cancelling factor for the DFT and SubSum + evantual scaling factor for the double angle formula
	btp.coeffsToSlotsDiffScale = complex(math.Pow(2.0/((b-a)*n*btp.scFac*qDiff), 1.0/float64(btp.CtSDepth(false))), 0)

	// Rescaling factor to set the final ciphertext to the desired scale
	btp.slotsToCoeffsDiffScale = complex(math.Pow((qDiff*btp.params.scale)/btp.postscale, 1.0/float64(btp.StCDepth(false))), 0)

	// Computation and encoding of the matrices for CoeffsToSlots and SlotsToCoeffs.
	btp.computePlaintextVectors()

	// List of the rotation key values to needed for the bootstrapp
	btp.rotKeyIndex = []int{}

	//SubSum rotation needed X -> Y^slots rotations
	for i := btp.params.logSlots; i < btp.params.MaxLogSlots(); i++ {
		if !utils.IsInSliceInt(1<<i, btp.rotKeyIndex) {
			btp.rotKeyIndex = append(btp.rotKeyIndex, 1<<i)
		}
	}

	// Coeffs to Slots rotations
	for _, pVec := range btp.pDFTInv {
		btp.addMatrixRotToList(pVec, btp.rotKeyIndex, btp.params.Slots(), false)
	}

	// Slots to Coeffs rotations
	for i, pVec := range btp.pDFT {
		if i == 0 {
			btp.addMatrixRotToList(pVec, btp.rotKeyIndex, btp.params.Slots(), btp.repack)
		} else {
			btp.addMatrixRotToList(pVec, btp.rotKeyIndex, btp.params.Slots(), false)
		}
	}
}

func (btp *Bootstrapper) genSinePoly() {

	K := int(btp.SinRange)
	deg := int(btp.SinDeg)
	btp.scFac = float64(int(1 << btp.SinRescal))

	if btp.ArcSineDeg > 0 {
		btp.sqrt2pi = 1.0

		coeffs := make([]complex128, btp.ArcSineDeg+1)

		coeffs[1] = 0.15915494309189535

		for i := 3; i < btp.ArcSineDeg+1; i += 2 {

			coeffs[i] = coeffs[i-2] * complex(float64(i*i-4*i+4)/float64(i*i-i), 0)

		}

		btp.arcSinePoly = NewPoly(coeffs)

	} else {
		btp.sqrt2pi = math.Pow(0.15915494309189535, 1.0/btp.scFac)
	}

	if btp.SinType == Sin {

		btp.sineEvalPoly = Approximate(sin2pi2pi, -complex(float64(K)/btp.scFac, 0), complex(float64(K)/btp.scFac, 0), deg)

	} else if btp.SinType == Cos1 {

		btp.sineEvalPoly = new(ChebyshevInterpolation)

		btp.sineEvalPoly.coeffs = bettersine.Approximate(K, deg, btp.MessageRatio, int(btp.SinRescal))

		btp.sineEvalPoly.maxDeg = btp.sineEvalPoly.Degree()
		btp.sineEvalPoly.a = complex(float64(-K)/btp.scFac, 0)
		btp.sineEvalPoly.b = complex(float64(K)/btp.scFac, 0)
		btp.sineEvalPoly.lead = true

	} else if btp.SinType == Cos2 {

		btp.sineEvalPoly = Approximate(cos2pi, -complex(float64(K)/btp.scFac, 0), complex(float64(K)/btp.scFac, 0), deg)

	} else {
		panic("Bootstrapper -> invalid sineType")
	}

	for i := range btp.sineEvalPoly.coeffs {
		btp.sineEvalPoly.coeffs[i] *= complex(btp.sqrt2pi, 0)
	}
}

func computeRoots(N int) (roots []complex128) {

	var angle float64

	m := N << 1

	roots = make([]complex128, m)

	roots[0] = 1

	for i := 1; i < m; i++ {
		angle = 6.283185307179586 * float64(i) / float64(m)
		roots[i] = complex(math.Cos(angle), math.Sin(angle))
	}

	return
}

func fftPlainVec(logN, dslots int, roots []complex128, pow5 []int) (a, b, c [][]complex128) {

	var N, m, index, tt, gap, k, mask, idx1, idx2 int

	N = 1 << logN

	a = make([][]complex128, logN)
	b = make([][]complex128, logN)
	c = make([][]complex128, logN)

	var size int
	if 2*N == dslots {
		size = 2
	} else {
		size = 1
	}

	index = 0
	for m = 2; m <= N; m <<= 1 {

		a[index] = make([]complex128, dslots)
		b[index] = make([]complex128, dslots)
		c[index] = make([]complex128, dslots)

		tt = m >> 1

		for i := 0; i < N; i += m {

			gap = N / m
			mask = (m << 2) - 1

			for j := 0; j < m>>1; j++ {

				k = (pow5[j] & mask) * gap

				idx1 = i + j
				idx2 = i + j + tt

				for u := 0; u < size; u++ {
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

func fftInvPlainVec(logN, dslots int, roots []complex128, pow5 []int) (a, b, c [][]complex128) {

	var N, m, index, tt, gap, k, mask, idx1, idx2 int

	N = 1 << logN

	a = make([][]complex128, logN)
	b = make([][]complex128, logN)
	c = make([][]complex128, logN)

	var size int
	if 2*N == dslots {
		size = 2
	} else {
		size = 1
	}

	index = 0
	for m = N; m >= 2; m >>= 1 {

		a[index] = make([]complex128, dslots)
		b[index] = make([]complex128, dslots)
		c[index] = make([]complex128, dslots)

		tt = m >> 1

		for i := 0; i < N; i += m {

			gap = N / m
			mask = (m << 2) - 1

			for j := 0; j < m>>1; j++ {

				k = ((m << 2) - (pow5[j] & mask)) * gap

				idx1 = i + j
				idx2 = i + j + tt

				for u := 0; u < size; u++ {

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

func (btp *Bootstrapper) computePlaintextVectors() {

	slots := btp.params.Slots()

	ctsLevel := btp.ctsLevel
	stcLevel := btp.stcLevel

	roots := computeRoots(slots << 1)
	pow5 := make([]int, (slots<<1)+1)
	pow5[0] = 1
	for i := 1; i < (slots<<1)+1; i++ {
		pow5[i] = pow5[i-1] * 5
		pow5[i] &= (slots << 2) - 1
	}

	// CoeffsToSlots vectors
	btp.pDFTInv = make([]*PtDiagMatrix, len(ctsLevel))
	pVecDFTInv := btp.computeDFTMatrices(roots, pow5, btp.coeffsToSlotsDiffScale, true)
	cnt := 0
	for i := range btp.CoeffsToSlotsModuli.ScalingFactor {
		for j := range btp.CoeffsToSlotsModuli.ScalingFactor[btp.CtSDepth(true)-i-1] {
			btp.pDFTInv[cnt] = btp.encoder.EncodeDiagMatrixAtLvl(ctsLevel[cnt], pVecDFTInv[cnt], btp.CoeffsToSlotsModuli.ScalingFactor[btp.CtSDepth(true)-i-1][j], btp.MaxN1N2Ratio, btp.logdslots)
			cnt++
		}
	}

	// SlotsToCoeffs vectors
	btp.pDFT = make([]*PtDiagMatrix, len(stcLevel))
	pVecDFT := btp.computeDFTMatrices(roots, pow5, btp.slotsToCoeffsDiffScale, false)
	cnt = 0
	for i := range btp.SlotsToCoeffsModuli.ScalingFactor {
		for j := range btp.SlotsToCoeffsModuli.ScalingFactor[btp.StCDepth(true)-i-1] {
			btp.pDFT[cnt] = btp.encoder.EncodeDiagMatrixAtLvl(stcLevel[cnt], pVecDFT[cnt], btp.SlotsToCoeffsModuli.ScalingFactor[btp.StCDepth(true)-i-1][j], btp.MaxN1N2Ratio, btp.logdslots)
			cnt++
		}
	}
}

func (btp *Bootstrapper) computeDFTMatrices(roots []complex128, pow5 []int, diffscale complex128, forward bool) (plainVector []map[int][]complex128) {

	var level, depth, nextLevel, logSlots int

	logSlots = btp.params.logSlots
	level = logSlots

	var a, b, c [][]complex128
	var maxDepth int

	if forward {
		maxDepth = btp.CtSDepth(false)
		a, b, c = fftInvPlainVec(btp.params.logSlots, btp.dslots, roots, pow5)
	} else {
		maxDepth = btp.StCDepth(false)
		a, b, c = fftPlainVec(btp.params.logSlots, btp.dslots, roots, pow5)
	}

	plainVector = make([]map[int][]complex128, maxDepth)

	// We compute the chain of merge in order or reverse order depending if its DFT or InvDFT because
	// the way the levels are collapsed has an inpact on the total number of rotations and keys to be
	// stored. Ex. instead of using 255 + 64 plaintext vectors, we can use 127 + 128 plaintext vectors
	// by reversing the order of the merging.
	merge := make([]int, maxDepth)
	for i := 0; i < maxDepth; i++ {

		depth = int(math.Ceil(float64(level) / float64(maxDepth-i)))

		if forward {
			merge[i] = depth
		} else {
			merge[len(merge)-i-1] = depth

		}

		level -= depth
	}

	level = logSlots
	for i := 0; i < maxDepth; i++ {

		if btp.repack && !forward && i == 0 {

			// Special initial matrix for the repacking before SlotsToCoeffs
			plainVector[i] = genWfftRepack(logSlots, level)

			// Merges this special initial matrix with the first layer of SlotsToCoeffs DFT
			plainVector[i] = nextLevelfft(plainVector[i], logSlots, 2<<logSlots, level, a[logSlots-level], b[logSlots-level], c[logSlots-level], forward)

			// Continues the merging with the next layers if the total depth requires it.
			nextLevel = level - 1
			for j := 0; j < merge[i]-1; j++ {
				plainVector[i] = nextLevelfft(plainVector[i], logSlots, 2<<logSlots, nextLevel, a[logSlots-nextLevel], b[logSlots-nextLevel], c[logSlots-nextLevel], forward)
				nextLevel--
			}

		} else {
			// First layer of the i-th level of the DFT
			plainVector[i] = genWfft(logSlots, level, a[logSlots-level], b[logSlots-level], c[logSlots-level], forward)

			// Merges the layer with the next levels of the DFT if the total depth requires it.
			nextLevel = level - 1
			for j := 0; j < merge[i]-1; j++ {
				plainVector[i] = nextLevelfft(plainVector[i], logSlots, 1<<logSlots, nextLevel, a[logSlots-nextLevel], b[logSlots-nextLevel], c[logSlots-nextLevel], forward)
				nextLevel--
			}
		}

		level -= merge[i]
	}

	// Repacking after the CoeffsToSlots (we multiply the last DFT matrix with the vector [1, 1, ..., 1, 1, 0, 0, ..., 0, 0]).
	if btp.repack && forward {
		for j := range plainVector[maxDepth-1] {
			for x := 0; x < btp.params.Slots(); x++ {
				plainVector[maxDepth-1][j][x+btp.params.Slots()] = complex(0, 0)
			}
		}
	}

	// Rescaling of the DFT matrix of the SlotsToCoeffs/CoeffsToSlots
	for j := range plainVector {
		for x := range plainVector[j] {
			for i := range plainVector[j][x] {
				plainVector[j][x][i] *= diffscale
			}
		}
	}

	return
}

func genWfft(logL, level int, a, b, c []complex128, forward bool) (vectors map[int][]complex128) {

	var rot int

	if forward {
		rot = 1<<(level - 1)
	} else {
		rot = 1<<(logL - level)
	}

	vectors = make(map[int][]complex128)

	addToDicVector(vectors, 0, a)
	addToDicVector(vectors, rot, b)
	addToDicVector(vectors, (1<<logL)-rot, c)

	return
}

func genWfftRepack(logL, level int) (vectors map[int][]complex128) {

	vectors = make(map[int][]complex128)

	a := make([]complex128, 2<<logL)
	b := make([]complex128, 2<<logL)

	for i := 0; i < 1<<logL; i++ {
		a[i] = complex(1, 0)
		a[i+(1<<logL)] = complex(0, 1)

		b[i] = complex(0, 1)
		b[i+(1<<logL)] = complex(1, 0)
	}

	addToDicVector(vectors, 0, a)
	addToDicVector(vectors, (1 << logL), b)

	return
}

func nextLevelfft(vec map[int][]complex128, logL, N, nextLevel int, a, b, c []complex128, forward bool) (newVec map[int][]complex128) {

	var rot int

	newVec = make(map[int][]complex128)

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

func addToDicVector(dic map[int][]complex128, index int, vec []complex128) {
	if dic[index] == nil {
		dic[index] = vec
	} else {
		dic[index] = add(dic[index], vec)
	}
}

func rotate(x []complex128, n int) (y []complex128) {

	y = make([]complex128, len(x))

	mask := int(len(x) - 1)

	// Rotates to the left
	for i := 0; i < len(x); i++ {
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
