package ckks

import (
	"fmt"
	"github.com/ldsec/lattigo/ckks/bettersine"
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
	"log"
	"math"
	"math/bits"
	"math/cmplx"
)

// BootContext stores the parameters for the bootstrapping.
type BootContext struct { // TODO: change to "Bootstrapper" ?
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

	rotKeyIndex []uint64       // a list of the required rotation keys
	relinkey    *EvaluationKey // Relinearization key
	rotkeys     *RotationKeys  // Rotation and conjugation keys

	ctxpool [3]*Ciphertext // Memory pool

	decryptor Decryptor

	poolQ [6]*ring.Poly
	poolP [4]*ring.Poly
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

func (b *BootContext) printDebug(message string, ciphertext *Ciphertext) {

	coeffs := b.encoder.Decode(b.decryptor.DecryptNew(ciphertext), b.dslots)

	if b.dslots == 2 {
		fmt.Printf(message+"%.10f %.10f...\n", coeffs[0], coeffs[1])
	} else {
		fmt.Printf(message+"%.10f %.10f %.10f %.10f...\n", coeffs[0], coeffs[1], coeffs[2], coeffs[3])
	}
}

// NewBootContext creates a new bootcontext.
func NewBootContext(bootparams *BootParams) (bootcontext *BootContext) {

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
	bootcontext.evaluator = NewEvaluator(&bootparams.Parameters)

	bootcontext.newBootSine()
	bootcontext.newBootDFT()

	bootcontext.ctxpool[0] = NewCiphertext(&bootparams.Parameters, 1, bootparams.Parameters.MaxLevel, 0)
	bootcontext.ctxpool[1] = NewCiphertext(&bootparams.Parameters, 1, bootparams.Parameters.MaxLevel, 0)
	bootcontext.ctxpool[2] = NewCiphertext(&bootparams.Parameters, 1, bootparams.Parameters.MaxLevel, 0)

	eval := bootcontext.evaluator.(*evaluator)
	contextQ := eval.ckksContext.contextQ
	contextP := eval.ckksContext.contextP

	for i := range bootcontext.poolQ {
		bootcontext.poolQ[i] = contextQ.NewPoly()
	}

	for i := range bootcontext.poolP {
		bootcontext.poolP[i] = contextP.NewPoly()
	}

	return bootcontext
}

func (bootcontext *BootContext) GenBootKeys(sk *SecretKey) {

	log.Println("DFT vector size (GB) :", float64(bootcontext.plaintextSize)/float64(1000000000))

	nbKeys := uint64(len(bootcontext.rotKeyIndex)) + 2 //rot keys + conj key + relin key
	nbPoly := bootcontext.Beta
	nbCoefficients := 2 * bootcontext.n * uint64(len(bootcontext.Qi)+len(bootcontext.Pi))
	bytesPerCoeff := uint64(8)

	log.Println("Switching-Keys size (GB) :", float64(nbKeys*nbPoly*nbCoefficients*bytesPerCoeff)/float64(1000000000), "(", nbKeys, "keys)")

	kgen := NewKeyGenerator(&bootcontext.Parameters)

	bootcontext.rotkeys = NewRotationKeys()

	kgen.GenRot(Conjugate, sk, 0, bootcontext.rotkeys)

	for _, i := range bootcontext.rotKeyIndex {
		kgen.GenRot(RotationLeft, sk, uint64(i), bootcontext.rotkeys)
	}

	bootcontext.relinkey = kgen.GenRelinKey(sk)

	return
}

func (bootcontext *BootContext) ExportKeys() (rlk *EvaluationKey, rotkeys *RotationKeys) {
	return bootcontext.relinkey, bootcontext.rotkeys
}

func (bootcontext *BootContext) ImportKeys(rlk *EvaluationKey, rotkeys *RotationKeys) {
	bootcontext.relinkey = rlk
	bootcontext.rotkeys = rotkeys
}

func (bootcontext *BootContext) CheckKeys() (err error) {

	if bootcontext.relinkey == nil || bootcontext.rotkeys == nil {
		return fmt.Errorf("empty relinkkey and/or rotkeys")
	}

	if bootcontext.rotkeys.evakeyConjugate == nil {
		return fmt.Errorf("missing conjugate key")
	}

	rotMissing := []uint64{}
	for _, i := range bootcontext.rotKeyIndex {
		if bootcontext.rotkeys.evakeyRotColLeft[i] == nil || bootcontext.rotkeys.permuteNTTLeftIndex[i] == nil {
			rotMissing = append(rotMissing, i)
		}
	}

	if len(rotMissing) != 0 {
		return fmt.Errorf("missing rotation keys : %d", rotMissing)
	}

	return nil
}

func (bootcontext *BootContext) newBootDFT() {

	a := real(bootcontext.chebycoeffs.a)
	b := real(bootcontext.chebycoeffs.b)
	n := float64(bootcontext.n)
	sc_fac := float64(int(1 << bootcontext.SinRescal))
	qDiff := float64(bootcontext.Qi[0]) / float64(1<<55)

	// Change of variable for the evaluation of the Chebyshev polynomial + cancelling factor for the DFT and SubSum + evantual scaling factor for the double angle formula
	bootcontext.coeffsToSlotsDiffScale = complex(math.Pow(2.0/((b-a)*n*sc_fac*qDiff), 1.0/float64(len(bootcontext.CtSLevel))), 0)

	// Rescaling factor to set the final ciphertext to the desired scale
	bootcontext.slotsToCoeffsDiffScale = complex(math.Pow((qDiff*bootcontext.Scale)/bootcontext.sinScale, 1.0/float64(len(bootcontext.StCLevel))), 0)

	// Computation and encoding of the matrices for CoeffsToSlots and SlotsToCoeffs.
	bootcontext.computePlaintextVectors()

	slots := bootcontext.slots

	// List of the rotation key values to needed for the bootstrapp
	bootcontext.rotKeyIndex = []uint64{}

	//SubSum rotation needed X -> Y^slots rotations
	for i := bootcontext.LogSlots; i < bootcontext.LogN-1; i++ {
		if !utils.IsInSliceUint64(1<<i, bootcontext.rotKeyIndex) {
			bootcontext.rotKeyIndex = append(bootcontext.rotKeyIndex, 1<<i)
		}
	}

	var index uint64
	// Coeffs to Slots rotations
	for i := range bootcontext.pDFTInv {
		for j := range bootcontext.pDFTInv[i].Vec {

			index = ((j / bootcontext.pDFTInv[i].N1) * bootcontext.pDFTInv[i].N1) & (slots - 1)

			if index != 0 && !utils.IsInSliceUint64(index, bootcontext.rotKeyIndex) {
				bootcontext.rotKeyIndex = append(bootcontext.rotKeyIndex, index)
			}

			index = j & (bootcontext.pDFTInv[i].N1 - 1)

			if index != 0 && !utils.IsInSliceUint64(index, bootcontext.rotKeyIndex) {
				bootcontext.rotKeyIndex = append(bootcontext.rotKeyIndex, index)
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

			if index != 0 && !utils.IsInSliceUint64(index, bootcontext.rotKeyIndex) {
				bootcontext.rotKeyIndex = append(bootcontext.rotKeyIndex, index)
			}

			index = j & (bootcontext.pDFT[i].N1 - 1)

			if index != 0 && !utils.IsInSliceUint64(index, bootcontext.rotKeyIndex) {
				bootcontext.rotKeyIndex = append(bootcontext.rotKeyIndex, index)
			}
		}
	}

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

		sqrt2pi := math.Pow(0.15915494309189535, 1.0/real(sc_fac))

		for i := range cheby.coeffs {
			cheby.coeffs[i] *= complex(sqrt2pi, 0)
		}

		cheby.maxDeg = cheby.degree()
		cheby.a = complex(float64(-K), 0) / sc_fac
		cheby.b = complex(float64(K), 0) / sc_fac

		bootcontext.chebycoeffs = cheby

	} else {
		panic("bootcontext -> invalid sineType")
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

	CtSLevel := bootcontext.CtSLevel
	StCLevel := bootcontext.StCLevel

	roots := computeRoots(slots << 1)
	pow5 := make([]uint64, (slots<<1)+1)
	pow5[0] = 1
	for i := uint64(1); i < (slots<<1)+1; i++ {
		pow5[i] = pow5[i-1] * 5
		pow5[i] &= (slots << 2) - 1
	}

	// CoeffsToSlots vectors
	bootcontext.pDFTInv = make([]*dftvectors, len(CtSLevel))
	pVecDFTInv := bootcontext.computeDFTPlaintextVectors(roots, pow5, bootcontext.coeffsToSlotsDiffScale, true)
	for i, lvl := range CtSLevel {
		bootcontext.pDFTInv[i] = new(dftvectors)
		bootcontext.pDFTInv[i].N1 = findbestbabygiantstepsplit(pVecDFTInv[i], dslots, bootcontext.MaxN1N2Ratio)
		bootcontext.encodePVec(pVecDFTInv[i], bootcontext.pDFTInv[i], lvl, true)
	}

	// SlotsToCoeffs vectors
	bootcontext.pDFT = make([]*dftvectors, len(StCLevel))
	pVecDFT := bootcontext.computeDFTPlaintextVectors(roots, pow5, bootcontext.slotsToCoeffsDiffScale, false)
	for i, lvl := range StCLevel {
		bootcontext.pDFT[i] = new(dftvectors)
		bootcontext.pDFT[i].N1 = findbestbabygiantstepsplit(pVecDFT[i], dslots, bootcontext.MaxN1N2Ratio)
		bootcontext.encodePVec(pVecDFT[i], bootcontext.pDFT[i], lvl, false)
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

			if hoisted > normal {

				// Finds the next split that has a ratio hoisted/normal greater or equal to 3
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

func (bootcontext *BootContext) encodePVec(pVec map[uint64][]complex128, plaintextVec *dftvectors, level uint64, forward bool) {
	var N, N1 uint64
	var scale float64

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

	plaintextVec.Vec = make(map[uint64][2]*ring.Poly)

	if forward {
		scale = float64(bootcontext.Qi[level])
	} else {
		// If the first moduli
		if bootcontext.LogQi[level] > 30 {
			scale = float64(uint64(1 << (bootcontext.LogQi[level] >> 1)))
		} else {
			scale = float64(bootcontext.Qi[level])
		}
	}

	plaintextVec.Level = level
	plaintextVec.Scale = scale
	contextQ := bootcontext.evaluator.(*evaluator).ckksContext.contextQ
	contextP := bootcontext.evaluator.(*evaluator).ckksContext.contextP
	encoder := bootcontext.encoder.(*encoderComplex128)

	for j := range index {

		for _, i := range index[j] {

			//  levels * n coefficients of 8 bytes each
			bootcontext.plaintextSize += (level + 1) * 8 * bootcontext.n

			encoder.embed(rotate(pVec[N1*j+uint64(i)], (N>>1)-(N1*j))[:bootcontext.dslots], bootcontext.dslots)

			plaintextQ := ring.NewPoly(bootcontext.Parameters.N, level+1)
			encoder.scaleUp(plaintextQ, scale, contextQ.Modulus[:level+1])
			contextQ.NTTLvl(level, plaintextQ, plaintextQ)
			contextQ.MFormLvl(level, plaintextQ, plaintextQ)

			plaintextP := ring.NewPoly(bootcontext.Parameters.N, level+1)
			encoder.scaleUp(plaintextP, scale, contextP.Modulus)
			contextP.NTT(plaintextP, plaintextP)
			contextP.MForm(plaintextP, plaintextP)

			plaintextVec.Vec[N1*j+uint64(i)] = [2]*ring.Poly{plaintextQ, plaintextP}

			encoder.wipeInternalMemory()

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
		maxDepth = uint64(len(bootcontext.CtSLevel))
		a, b, c = fftInvPlainVec(slots, roots, pow5)
	} else {
		maxDepth = uint64(len(bootcontext.StCLevel))
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
