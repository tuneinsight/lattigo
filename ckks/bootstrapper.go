package ckks

import (
	"fmt"
	"log"
	"math"
	"math/cmplx"

	"github.com/ldsec/lattigo/v2/ckks/bettersine"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// Bootstrapper is a struct to stores a memory pool the plaintext matrices
// the polynomial approximation and the keys for the bootstrapping.
type Bootstrapper struct {
	BootstrappParams
	*BootstrappingKey
	params *Parameters

	dslots    uint64 // Number of plaintext slots after the re-encoding
	logdslots uint64

	encoder   Encoder   // Encoder
	evaluator Evaluator // Evaluator

	plaintextSize uint64 // Byte size of the plaintext DFT matrices

	repack      bool                    // If true then can repack the CoeffsToSlots into on ciphertext
	deviation   float64                 // Q[0]/Scale
	prescale    float64                 // Q[0]/1024
	postscale   float64                 // Qi sineeval/2^{10}
	chebycoeffs *ChebyshevInterpolation // Coefficients of the Chebyshev Interpolation of sin(2*pi*x) or cos(2*pi*x/r)

	coeffsToSlotsDiffScale complex128    // Matrice rescaling
	slotsToCoeffsDiffScale complex128    // Matrice rescaling
	pDFT                   []*dftvectors // Matrice vectors
	pDFTInv                []*dftvectors // Matrice vectors

	rotKeyIndex []uint64 // a list of the required rotation keys

	ctxpool *Ciphertext // Memory pool

	decryptor Decryptor

	poolQ [1]*ring.Poly // Memory pool for the matrix evaluation
	poolP [2]*ring.Poly // Memory pool for the matrix evaluation
}

type dftvectors struct {
	N1    uint64
	Level uint64
	Scale float64
	Vec   map[uint64][2]*ring.Poly
}

func sin2pi2pi(x complex128) complex128 {
	return cmplx.Sin(6.283185307179586*x) / 6.283185307179586
}

func cos2pi(x complex128) complex128 {
	return cmplx.Cos(6.283185307179586 * x)
}

func (btp *Bootstrapper) printDebug(message string, ciphertext *Ciphertext) {

	coeffs := btp.encoder.Decode(btp.decryptor.DecryptNew(ciphertext), btp.dslots)

	if btp.dslots == 2 {
		log.Printf(message+"%.10f %.10f...\n", coeffs[0], coeffs[1])
	} else {
		log.Printf(message+"%.10f %.10f %.10f %.10f...\n", coeffs[0], coeffs[1], coeffs[2], coeffs[3])
	}
}

// NewBootstrapper creates a new Bootstrapper.
func NewBootstrapper(params *Parameters, btpParams *BootstrappParams, btpKey *BootstrappingKey) (btp *Bootstrapper, err error) {

	if btpParams.SinType == SinType(Sin) && btpParams.SinRescal != 0 {
		return nil, fmt.Errorf("BootstrappParams: cannot use double angle formul for SinType = Sin -> must use SinType = Cos")
	}

	if btpParams.CtSLevel[0] != params.MaxLevel() {
		return nil, fmt.Errorf("BootstrappParams: CtSLevel start not consistent with MaxLevel")
	}

	btp = newBootstrapper(params, btpParams)

	btp.BootstrappingKey = btpKey
	if err = btp.CheckKeys(); err != nil {
		return nil, err
	}

	return btp, nil
}

// newBootstrapper is a constructor of "dummy" bootstrapper to enable the generation of bootstrapping-related constants
// without providing a bootstrapping key. To be replaced by a propper factorization of the bootstrapping pre-computations.
func newBootstrapper(params *Parameters, btpParams *BootstrappParams) (btp *Bootstrapper) {
	btp = new(Bootstrapper)

	btp.params = params.Copy()
	btp.BootstrappParams = *btpParams.Copy()

	btp.dslots = params.Slots()
	btp.logdslots = params.LogSlots()
	if params.logSlots < params.MaxLogSlots() {
		btp.repack = true
		btp.dslots <<= 1
		btp.logdslots++
	}

	btp.deviation = 1024.0
	btp.prescale = math.Exp2(math.Round(math.Log2(float64(params.qi[0]) / btp.deviation)))
	btp.postscale = math.Exp2(math.Round(math.Log2(float64(params.qi[len(params.qi)-1-len(btpParams.CtSLevel)])))) / btp.deviation

	btp.encoder = NewEncoder(params)
	btp.evaluator = NewEvaluator(params)

	btp.genSinePoly()
	btp.genDFTMatrices()

	btp.ctxpool = NewCiphertext(params, 1, params.MaxLevel(), 0)

	for i := range btp.poolQ {
		btp.poolQ[i] = params.NewPolyQ()
	}

	for i := range btp.poolP {
		btp.poolP[i] = params.NewPolyP()
	}
	return btp
}

// CheckKeys checks if all the necessary keys are present
func (btp *Bootstrapper) CheckKeys() (err error) {

	if btp.relinkey == nil || btp.rotkeys == nil {
		return fmt.Errorf("empty relinkkey and/or rotkeys")
	}

	if btp.rotkeys.evakeyConjugate == nil {
		return fmt.Errorf("missing conjugate key")
	}

	rotMissing := []uint64{}
	for _, i := range btp.rotKeyIndex {
		if btp.rotkeys.evakeyRotColLeft[i] == nil || btp.rotkeys.permuteNTTLeftIndex[i] == nil {
			rotMissing = append(rotMissing, i)
		}
	}

	if len(rotMissing) != 0 {
		return fmt.Errorf("missing rotation keys : %d", rotMissing)
	}

	return nil
}

func (btp *Bootstrapper) genDFTMatrices() {

	a := real(btp.chebycoeffs.a)
	b := real(btp.chebycoeffs.b)
	n := float64(btp.params.N())
	scFac := float64(int(1 << btp.SinRescal))
	qDiff := float64(btp.params.qi[0]) / math.Exp2(math.Round(math.Log2(float64(btp.params.qi[0]))))

	// Change of variable for the evaluation of the Chebyshev polynomial + cancelling factor for the DFT and SubSum + evantual scaling factor for the double angle formula
	btp.coeffsToSlotsDiffScale = complex(math.Pow(2.0/((b-a)*n*scFac*qDiff), 1.0/float64(len(btp.CtSLevel))), 0)

	// Rescaling factor to set the final ciphertext to the desired scale
	btp.slotsToCoeffsDiffScale = complex(math.Pow((qDiff*btp.params.scale)/btp.postscale, 1.0/float64(len(btp.StCLevel))), 0)

	// Computation and encoding of the matrices for CoeffsToSlots and SlotsToCoeffs.
	btp.computePlaintextVectors()

	// List of the rotation key values to needed for the bootstrapp
	btp.rotKeyIndex = []uint64{}

	//SubSum rotation needed X -> Y^slots rotations
	for i := btp.params.logSlots; i < btp.params.MaxLogSlots(); i++ {
		if !utils.IsInSliceUint64(1<<i, btp.rotKeyIndex) {
			btp.rotKeyIndex = append(btp.rotKeyIndex, 1<<i)
		}
	}

	var index uint64
	// Coeffs to Slots rotations
	for i := range btp.pDFTInv {
		for j := range btp.pDFTInv[i].Vec {

			index = ((j / btp.pDFTInv[i].N1) * btp.pDFTInv[i].N1) & (btp.params.Slots() - 1)

			if index != 0 && !utils.IsInSliceUint64(index, btp.rotKeyIndex) {
				btp.rotKeyIndex = append(btp.rotKeyIndex, index)
			}

			index = j & (btp.pDFTInv[i].N1 - 1)

			if index != 0 && !utils.IsInSliceUint64(index, btp.rotKeyIndex) {
				btp.rotKeyIndex = append(btp.rotKeyIndex, index)
			}
		}
	}

	// Slots to Coeffs rotations
	for i := range btp.pDFT {
		for j := range btp.pDFT[i].Vec {

			if btp.repack && i == 0 {
				// Sparse repacking, occuring during the first DFT matrix of the CoeffsToSlots.
				index = ((j / btp.pDFT[i].N1) * btp.pDFT[i].N1) & (2*btp.params.Slots() - 1)
			} else {
				// Other cases
				index = ((j / btp.pDFT[i].N1) * btp.pDFT[i].N1) & (btp.params.Slots() - 1)
			}

			if index != 0 && !utils.IsInSliceUint64(index, btp.rotKeyIndex) {
				btp.rotKeyIndex = append(btp.rotKeyIndex, index)
			}

			index = j & (btp.pDFT[i].N1 - 1)

			if index != 0 && !utils.IsInSliceUint64(index, btp.rotKeyIndex) {
				btp.rotKeyIndex = append(btp.rotKeyIndex, index)
			}
		}
	}

	/*
		log.Println("DFT vector size (GB) :", float64(btp.plaintextSize)/float64(1000000000))

		nbKeys := uint64(len(btp.rotKeyIndex)) + 2 //rot keys + conj key + relin key
		nbPoly := btp.params.Beta()
		nbCoefficients := 2 * btp.params.N() * btp.params.QPiCount()
		bytesPerCoeff := uint64(8)

		log.Println("Switching-Keys size (GB) :", float64(nbKeys*nbPoly*nbCoefficients*bytesPerCoeff)/float64(1000000000), "(", nbKeys, "keys)")
	*/

	return
}

func (btp *Bootstrapper) genSinePoly() {

	if btp.SinType == Sin {

		K := complex(float64(btp.SinRange), 0)
		btp.chebycoeffs = Approximate(sin2pi2pi, -K, K, int(btp.SinDeg))

	} else if btp.SinType == Cos1 {

		K := int(btp.SinRange)
		deg := int(btp.SinDeg)
		scFac := complex(float64(int(1<<btp.SinRescal)), 0)

		cheby := new(ChebyshevInterpolation)

		cheby.coeffs = bettersine.Approximate(K, deg, btp.deviation, int(btp.SinRescal))

		sqrt2pi := math.Pow(0.15915494309189535, 1.0/real(scFac))

		for i := range cheby.coeffs {
			cheby.coeffs[i] *= complex(sqrt2pi, 0)
		}

		cheby.maxDeg = cheby.Degree()
		cheby.a = complex(float64(-K), 0) / scFac
		cheby.b = complex(float64(K), 0) / scFac
		cheby.lead = true

		btp.chebycoeffs = cheby

	} else if btp.SinType == Cos2 {

		K := int(btp.SinRange)
		deg := int(btp.SinDeg)
		scFac := complex(float64(int(1<<btp.SinRescal)), 0)

		cheby := Approximate(cos2pi, -complex(float64(K), 0)/scFac, complex(float64(K), 0)/scFac, deg)
		sqrt2pi := math.Pow(0.15915494309189535, 1.0/real(scFac))

		for i := range cheby.coeffs {
			cheby.coeffs[i] *= complex(sqrt2pi, 0)
		}

		btp.chebycoeffs = cheby

	} else {
		panic("Bootstrapper -> invalid sineType")
	}
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

func fftPlainVec(logN uint64, roots []complex128, pow5 []uint64) (a, b, c [][]complex128) {

	var N, m, index, tt, gap, k, mask, idx1, idx2 uint64

	N = 1 << logN

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

func fftInvPlainVec(logN uint64, roots []complex128, pow5 []uint64) (a, b, c [][]complex128) {

	var N, m, index, tt, gap, k, mask, idx1, idx2 uint64

	N = 1 << logN

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

func (btp *Bootstrapper) computePlaintextVectors() {

	slots := btp.params.Slots()
	dslots := btp.dslots

	CtSLevel := btp.CtSLevel
	StCLevel := btp.StCLevel

	roots := computeRoots(slots << 1)
	pow5 := make([]uint64, (slots<<1)+1)
	pow5[0] = 1
	for i := uint64(1); i < (slots<<1)+1; i++ {
		pow5[i] = pow5[i-1] * 5
		pow5[i] &= (slots << 2) - 1
	}

	// CoeffsToSlots vectors
	btp.pDFTInv = make([]*dftvectors, len(CtSLevel))
	pVecDFTInv := btp.computeDFTPlaintextVectors(roots, pow5, btp.coeffsToSlotsDiffScale, true)
	for i, lvl := range CtSLevel {
		btp.pDFTInv[i] = new(dftvectors)
		btp.pDFTInv[i].N1 = findbestbabygiantstepsplit(pVecDFTInv[i], dslots, btp.MaxN1N2Ratio)
		btp.encodePVec(pVecDFTInv[i], btp.pDFTInv[i], lvl, true)
	}

	// SlotsToCoeffs vectors
	btp.pDFT = make([]*dftvectors, len(StCLevel))
	pVecDFT := btp.computeDFTPlaintextVectors(roots, pow5, btp.slotsToCoeffsDiffScale, false)
	for i, lvl := range StCLevel {
		btp.pDFT[i] = new(dftvectors)
		btp.pDFT[i].N1 = findbestbabygiantstepsplit(pVecDFT[i], dslots, btp.MaxN1N2Ratio)
		btp.encodePVec(pVecDFT[i], btp.pDFT[i], lvl, false)
	}
}

// Finds the best N1*N2 = N for the baby-step giant-step algorithm for matrix multiplication.
func findbestbabygiantstepsplit(vector map[uint64][]complex128, maxN uint64, maxRatio float64) (minN uint64) {

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

		if len(index[0]) > 0 {

			hoisted := len(index[0]) - 1
			normal := len(index) - 1

			// The matrice is very sparse already
			if normal == 0 {
				return N1 / 2
			}

			if hoisted > normal {
				// Finds the next split that has a ratio hoisted/normal greater or equal to maxRatio
				for float64(hoisted)/float64(normal) < maxRatio {

					if normal/2 == 0 {
						break
					}
					N1 *= 2
					hoisted = hoisted*2 + 1
					normal = normal / 2
				}
				return N1
			}
		}
	}

	return 1
}

func (btp *Bootstrapper) encodePVec(pVec map[uint64][]complex128, plaintextVec *dftvectors, level uint64, forward bool) {
	var N, N1 uint64
	var scale float64

	// N1*N2 = N
	N = btp.params.N()
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

	plaintextVec.Vec = make(map[uint64][2]*ring.Poly)

	if forward {
		scale = float64(btp.params.qi[level])
	} else {
		// If the first moduli
		logQi := math.Round(math.Log2(float64(btp.params.qi[level])))
		if logQi >= 56.0 {
			scale = math.Exp2(logQi / 2)
		} else {
			scale = float64(btp.params.qi[level])
		}
	}

	plaintextVec.Level = level
	plaintextVec.Scale = scale
	ringQ := btp.evaluator.(*evaluator).ringQ
	ringP := btp.evaluator.(*evaluator).ringP
	encoder := btp.encoder.(*encoderComplex128)

	for j := range index {

		for _, i := range index[j] {

			//  levels * n coefficients of 8 bytes each
			btp.plaintextSize += 8 * N * (level + 1 + btp.params.PiCount())

			encoder.Embed(rotate(pVec[N1*j+uint64(i)], (N>>1)-(N1*j))[:btp.dslots], btp.logdslots)

			plaintextQ := ring.NewPoly(N, level+1)
			encoder.ScaleUp(plaintextQ, scale, ringQ.Modulus[:level+1])
			ringQ.NTTLvl(level, plaintextQ, plaintextQ)
			ringQ.MFormLvl(level, plaintextQ, plaintextQ)

			plaintextP := ring.NewPoly(N, level+1)
			encoder.ScaleUp(plaintextP, scale, ringP.Modulus)
			ringP.NTT(plaintextP, plaintextP)
			ringP.MForm(plaintextP, plaintextP)

			plaintextVec.Vec[N1*j+uint64(i)] = [2]*ring.Poly{plaintextQ, plaintextP}

			encoder.WipeInternalMemory()

		}
	}
}

func (btp *Bootstrapper) computeDFTPlaintextVectors(roots []complex128, pow5 []uint64, diffscale complex128, forward bool) (plainVector []map[uint64][]complex128) {

	var level, depth, nextLevel, logSlots uint64

	logSlots = btp.params.logSlots
	level = logSlots

	var a, b, c [][]complex128
	var maxDepth uint64

	if forward {
		maxDepth = uint64(len(btp.CtSLevel))
		a, b, c = fftInvPlainVec(btp.params.logSlots, roots, pow5)
	} else {
		maxDepth = uint64(len(btp.StCLevel))
		a, b, c = fftPlainVec(btp.params.logSlots, roots, pow5)
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

		if btp.repack && !forward && i == 0 {

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
	if btp.repack && forward {
		for j := range plainVector[maxDepth-1] {
			for x := uint64(0); x < btp.params.Slots(); x++ {
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
