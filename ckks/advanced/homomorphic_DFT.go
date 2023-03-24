package advanced

import (
	"math"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// DFTType is a type used to distinguish different linear transformations.
type DFTType int

// Encode (IDFT) and Decode (DFT) are two available linear transformations for homomorphic encoding and decoding.
const (
	Encode = DFTType(0) // Homomorphic Encoding (IDFT)
	Decode = DFTType(1) // Homomorphic Decoding (DFT)
)

// HomomorphicDFTMatrix is a struct storing the factorized IDFT, DFT matrices, which are
// used to hommorphically encode and decode a ciphertext respectively.
type HomomorphicDFTMatrix struct {
	HomomorphicDFTMatrixLiteral
	Matrices []ckks.LinearTransform
}

// HomomorphicDFTMatrixLiteral is a struct storing the parameters to generate the factorized DFT/IDFT matrices.
// This struct has mandatory and optional fields.
//
// Mandatory:
//   - DFTType: Encode (a.k.a. CoeffsToSlots) or Decode (a.k.a. SlotsToCoeffs)
//   - LogN: log2(RingDegree)
//   - LogSlots: log2(slots)
//   - LevelStart: starting level of the linear transformation
//   - Levels: depth of the linear transform (i.e. the degree of factorization of the encoding matrix)
//
// Optional:
//   - RepackImag2Real: if true, the imaginary part is repacked into the right n slots of the real part
//   - Scaling: constant by which the matrix is multiplied
//   - BitReversed: if true, then applies the transformation bit-reversed and expects bit-reversed inputs
//   - LogBSGSRatio: log2 of the ratio between the inner and outer loop of the baby-step giant-step algorithm
type HomomorphicDFTMatrixLiteral struct {
	// Mandatory
	Type       DFTType
	LogN       int
	LogSlots   int
	LevelStart int
	Levels     []int
	// Optional
	RepackImag2Real bool    // Default: False.
	Scaling         float64 // Default 1.0.
	BitReversed     bool    // Default: False.
	LogBSGSRatio    int     // Default: 0.
}

// Depth returns the number of levels allocated to the linear transform.
// If actual == true then returns the number of moduli consumed, else
// returns the factorization depth.
func (d *HomomorphicDFTMatrixLiteral) Depth(actual bool) (depth int) {
	if actual {
		depth = len(d.Levels)
	} else {
		for _, d := range d.Levels {
			depth += d
		}
	}
	return
}

// GaloisElements returns the list of rotations performed during the CoeffsToSlot operation.
func (d *HomomorphicDFTMatrixLiteral) GaloisElements(params ckks.Parameters) (galEls []uint64) {
	rotations := []int{}

	logSlots := d.LogSlots
	logN := d.LogN
	slots := 1 << logSlots
	dslots := slots
	if logSlots < logN-1 && d.RepackImag2Real {
		dslots <<= 1
		if d.Type == Encode {
			rotations = append(rotations, slots)
		}
	}

	indexCtS := d.computeBootstrappingDFTIndexMap()

	// Coeffs to Slots rotations
	for i, pVec := range indexCtS {
		N1 := ckks.FindBestBSGSRatio(pVec, dslots, d.LogBSGSRatio)
		rotations = addMatrixRotToList(pVec, rotations, N1, slots, d.Type == Decode && logSlots < logN-1 && i == 0 && d.RepackImag2Real)
	}

	return params.GaloisElementsForRotations(rotations)
}

// NewHomomorphicDFTMatrixFromLiteral generates the factorized DFT/IDFT matrices for the homomorphic encoding/decoding.
func NewHomomorphicDFTMatrixFromLiteral(d HomomorphicDFTMatrixLiteral, encoder ckks.Encoder) HomomorphicDFTMatrix {

	logSlots := d.LogSlots
	logdSlots := logSlots
	if logdSlots < d.LogN-1 && d.RepackImag2Real {
		logdSlots++
	}

	params := encoder.Parameters()

	// DFT vectors
	matrices := []ckks.LinearTransform{}
	pVecDFT := d.GenMatrices()

	level := d.LevelStart
	var idx int
	for i := range d.Levels {

		scale := rlwe.NewScale(math.Pow(params.QiFloat64(level), 1.0/float64(d.Levels[i])))

		for j := 0; j < d.Levels[i]; j++ {
			matrices = append(matrices, ckks.GenLinearTransformBSGS(encoder, pVecDFT[idx], level, scale, d.LogBSGSRatio, logdSlots))
			idx++
		}

		level--
	}

	return HomomorphicDFTMatrix{HomomorphicDFTMatrixLiteral: d, Matrices: matrices}
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

func ifftPlainVec(logN, dslots int, roots []complex128, pow5 []int) (a, b, c [][]complex128) {

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

func addMatrixRotToList(pVec map[int]bool, rotations []int, N1, slots int, repack bool) []int {

	if len(pVec) < 3 {
		for j := range pVec {
			if !utils.IsInSlice(j, rotations) {
				rotations = append(rotations, j)
			}
		}
	} else {
		var index int
		for j := range pVec {

			index = (j / N1) * N1

			if repack {
				// Sparse repacking, occurring during the first IDFT matrix.
				index &= (2*slots - 1)
			} else {
				// Other cases
				index &= (slots - 1)
			}

			if index != 0 && !utils.IsInSlice(index, rotations) {
				rotations = append(rotations, index)
			}

			index = j & (N1 - 1)

			if index != 0 && !utils.IsInSlice(index, rotations) {
				rotations = append(rotations, index)
			}
		}
	}

	return rotations
}

func (d *HomomorphicDFTMatrixLiteral) computeBootstrappingDFTIndexMap() (rotationMap []map[int]bool) {

	logN := d.LogN
	logSlots := d.LogSlots
	ltType := d.Type
	repacki2r := d.RepackImag2Real
	bitreversed := d.BitReversed
	maxDepth := d.Depth(false)

	var level, depth, nextLevel int

	level = logSlots

	rotationMap = make([]map[int]bool, maxDepth)

	// We compute the chain of merge in order or reverse order depending if its DFT or InvDFT because
	// the way the levels are collapsed has an impact on the total number of rotations and keys to be
	// stored. Ex. instead of using 255 + 64 plaintext vectors, we can use 127 + 128 plaintext vectors
	// by reversing the order of the merging.
	merge := make([]int, maxDepth)
	for i := 0; i < maxDepth; i++ {

		depth = int(math.Ceil(float64(level) / float64(maxDepth-i)))

		if ltType == Encode {
			merge[i] = depth
		} else {
			merge[len(merge)-i-1] = depth

		}

		level -= depth
	}

	level = logSlots
	for i := 0; i < maxDepth; i++ {

		if logSlots < logN-1 && ltType == Decode && i == 0 && repacki2r {

			// Special initial matrix for the repacking before DFT
			rotationMap[i] = genWfftRepackIndexMap(logSlots, level)

			// Merges this special initial matrix with the first layer of DFT
			rotationMap[i] = nextLevelfftIndexMap(rotationMap[i], logSlots, 2<<logSlots, level, ltType, bitreversed)

			// Continues the merging with the next layers if the total depth requires it.
			nextLevel = level - 1
			for j := 0; j < merge[i]-1; j++ {
				rotationMap[i] = nextLevelfftIndexMap(rotationMap[i], logSlots, 2<<logSlots, nextLevel, ltType, bitreversed)
				nextLevel--
			}

		} else {

			// First layer of the i-th level of the DFT
			rotationMap[i] = genWfftIndexMap(logSlots, level, ltType, bitreversed)

			// Merges the layer with the next levels of the DFT if the total depth requires it.
			nextLevel = level - 1
			for j := 0; j < merge[i]-1; j++ {
				rotationMap[i] = nextLevelfftIndexMap(rotationMap[i], logSlots, 1<<logSlots, nextLevel, ltType, bitreversed)
				nextLevel--
			}
		}

		level -= merge[i]
	}

	return
}

func genWfftIndexMap(logL, level int, ltType DFTType, bitreversed bool) (vectors map[int]bool) {

	var rot int

	if ltType == Encode && !bitreversed || ltType == Decode && bitreversed {
		rot = 1 << (level - 1)
	} else {
		rot = 1 << (logL - level)
	}

	vectors = make(map[int]bool)
	vectors[0] = true
	vectors[rot] = true
	vectors[(1<<logL)-rot] = true
	return
}

func genWfftRepackIndexMap(logL, level int) (vectors map[int]bool) {
	vectors = make(map[int]bool)
	vectors[0] = true
	vectors[(1 << logL)] = true
	return
}

func nextLevelfftIndexMap(vec map[int]bool, logL, N, nextLevel int, ltType DFTType, bitreversed bool) (newVec map[int]bool) {

	var rot int

	newVec = make(map[int]bool)

	if ltType == Encode && !bitreversed || ltType == Decode && bitreversed {
		rot = (1 << (nextLevel - 1)) & (N - 1)
	} else {
		rot = (1 << (logL - nextLevel)) & (N - 1)
	}

	for i := range vec {
		newVec[i] = true
		newVec[(i+rot)&(N-1)] = true
		newVec[(i-rot)&(N-1)] = true
	}

	return
}

// GenMatrices returns the ordered list of factors of the non-zero diagonals of the IDFT (encoding) or DFT (decoding) matrix.
func (d *HomomorphicDFTMatrixLiteral) GenMatrices() (plainVector []map[int][]complex128) {

	logSlots := d.LogSlots
	slots := 1 << logSlots
	maxDepth := d.Depth(false)
	ltType := d.Type
	bitreversed := d.BitReversed

	logdSlots := logSlots
	if logdSlots < d.LogN-1 && d.RepackImag2Real {
		logdSlots++
	}

	roots := ckks.GetRootsFloat64(slots << 2)
	pow5 := make([]int, (slots<<1)+1)
	pow5[0] = 1
	for i := 1; i < (slots<<1)+1; i++ {
		pow5[i] = pow5[i-1] * 5
		pow5[i] &= (slots << 2) - 1
	}

	var fftLevel, depth, nextfftLevel int

	fftLevel = logSlots

	var a, b, c [][]complex128
	if ltType == Encode {
		a, b, c = ifftPlainVec(logSlots, 1<<logdSlots, roots, pow5)
	} else {
		a, b, c = fftPlainVec(logSlots, 1<<logdSlots, roots, pow5)
	}

	plainVector = make([]map[int][]complex128, maxDepth)

	// We compute the chain of merge in order or reverse order depending if its DFT or InvDFT because
	// the way the levels are collapsed has an impact on the total number of rotations and keys to be
	// stored. Ex. instead of using 255 + 64 plaintext vectors, we can use 127 + 128 plaintext vectors
	// by reversing the order of the merging.
	merge := make([]int, maxDepth)
	for i := 0; i < maxDepth; i++ {

		depth = int(math.Ceil(float64(fftLevel) / float64(maxDepth-i)))

		if ltType == Encode {
			merge[i] = depth
		} else {
			merge[len(merge)-i-1] = depth

		}

		fftLevel -= depth
	}

	fftLevel = logSlots
	for i := 0; i < maxDepth; i++ {

		if logSlots != logdSlots && ltType == Decode && i == 0 && d.RepackImag2Real {

			// Special initial matrix for the repacking before DFT
			plainVector[i] = genRepackMatrix(logSlots, bitreversed)

			// Merges this special initial matrix with the first layer of DFT
			plainVector[i] = multiplyFFTMatrixWithNextFFTLevel(plainVector[i], logSlots, 2*slots, fftLevel, a[logSlots-fftLevel], b[logSlots-fftLevel], c[logSlots-fftLevel], ltType, bitreversed)

			// Continues the merging with the next layers if the total depth requires it.
			nextfftLevel = fftLevel - 1
			for j := 0; j < merge[i]-1; j++ {
				plainVector[i] = multiplyFFTMatrixWithNextFFTLevel(plainVector[i], logSlots, 2*slots, nextfftLevel, a[logSlots-nextfftLevel], b[logSlots-nextfftLevel], c[logSlots-nextfftLevel], ltType, bitreversed)
				nextfftLevel--
			}

		} else {
			// First layer of the i-th level of the DFT
			plainVector[i] = genFFTDiagMatrix(logSlots, fftLevel, a[logSlots-fftLevel], b[logSlots-fftLevel], c[logSlots-fftLevel], ltType, bitreversed)

			// Merges the layer with the next levels of the DFT if the total depth requires it.
			nextfftLevel = fftLevel - 1
			for j := 0; j < merge[i]-1; j++ {
				plainVector[i] = multiplyFFTMatrixWithNextFFTLevel(plainVector[i], logSlots, slots, nextfftLevel, a[logSlots-nextfftLevel], b[logSlots-nextfftLevel], c[logSlots-nextfftLevel], ltType, bitreversed)
				nextfftLevel--
			}
		}

		fftLevel -= merge[i]
	}

	// Repacking after the IDFT (we multiply the last matrix with the vector [1, 1, ..., 1, 1, 0, 0, ..., 0, 0]).
	if logSlots != logdSlots && ltType == Encode && d.RepackImag2Real {
		for j := range plainVector[maxDepth-1] {
			v := plainVector[maxDepth-1][j]
			for x := 0; x < slots; x++ {
				v[x+slots] = complex(0, 0)
			}
		}
	}

	// Rescaling of the DFT matrices.
	scaling := complex(d.Scaling, 0)

	// If no scaling (Default); set to 1
	if scaling == 0 {
		scaling = 1.0
	}

	// If DFT matrix, rescale by 1/N
	if ltType == Encode {
		scaling /= complex(float64(slots), 0)

		// Real/Imag extraction 1/2 factor
		if d.RepackImag2Real {
			scaling /= 2
		}
	}

	// Spreads the scale across the matrices
	scaling = complex(math.Pow(real(scaling), 1.0/float64(d.Depth(false))), 0)

	for j := range plainVector {
		for x := range plainVector[j] {
			v := plainVector[j][x]
			for i := range v {
				v[i] *= scaling
			}
		}
	}

	return
}

func genFFTDiagMatrix(logL, fftLevel int, a, b, c []complex128, ltType DFTType, bitreversed bool) (vectors map[int][]complex128) {

	var rot int

	if ltType == Encode && !bitreversed || ltType == Decode && bitreversed {
		rot = 1 << (fftLevel - 1)
	} else {
		rot = 1 << (logL - fftLevel)
	}

	vectors = make(map[int][]complex128)

	if bitreversed {
		ckks.SliceBitReverseInPlaceComplex128(a, 1<<logL)
		ckks.SliceBitReverseInPlaceComplex128(b, 1<<logL)
		ckks.SliceBitReverseInPlaceComplex128(c, 1<<logL)

		if len(a) > 1<<logL {
			ckks.SliceBitReverseInPlaceComplex128(a[1<<logL:], 1<<logL)
			ckks.SliceBitReverseInPlaceComplex128(b[1<<logL:], 1<<logL)
			ckks.SliceBitReverseInPlaceComplex128(c[1<<logL:], 1<<logL)
		}
	}

	addToDiagMatrix(vectors, 0, a)
	addToDiagMatrix(vectors, rot, b)
	addToDiagMatrix(vectors, (1<<logL)-rot, c)

	return
}

func genRepackMatrix(logL int, bitreversed bool) (vectors map[int][]complex128) {

	vectors = make(map[int][]complex128)

	a := make([]complex128, 2<<logL)
	b := make([]complex128, 2<<logL)

	for i := 0; i < 1<<logL; i++ {
		a[i] = complex(1, 0)
		a[i+(1<<logL)] = complex(0, 1)

		b[i] = complex(0, 1)
		b[i+(1<<logL)] = complex(1, 0)
	}

	addToDiagMatrix(vectors, 0, a)
	addToDiagMatrix(vectors, (1 << logL), b)

	return
}

func multiplyFFTMatrixWithNextFFTLevel(vec map[int][]complex128, logL, N, nextLevel int, a, b, c []complex128, ltType DFTType, bitreversed bool) (newVec map[int][]complex128) {

	var rot int

	newVec = make(map[int][]complex128)

	if ltType == Encode && !bitreversed || ltType == Decode && bitreversed {
		rot = (1 << (nextLevel - 1)) & (N - 1)
	} else {
		rot = (1 << (logL - nextLevel)) & (N - 1)
	}

	if bitreversed {
		ckks.SliceBitReverseInPlaceComplex128(a, 1<<logL)
		ckks.SliceBitReverseInPlaceComplex128(b, 1<<logL)
		ckks.SliceBitReverseInPlaceComplex128(c, 1<<logL)

		if len(a) > 1<<logL {
			ckks.SliceBitReverseInPlaceComplex128(a[1<<logL:], 1<<logL)
			ckks.SliceBitReverseInPlaceComplex128(b[1<<logL:], 1<<logL)
			ckks.SliceBitReverseInPlaceComplex128(c[1<<logL:], 1<<logL)
		}
	}

	for i := range vec {
		addToDiagMatrix(newVec, i, mul(vec[i], a))
		addToDiagMatrix(newVec, (i+rot)&(N-1), mul(rotate(vec[i], rot), b))
		addToDiagMatrix(newVec, (i-rot)&(N-1), mul(rotate(vec[i], -rot), c))
	}

	return
}

func addToDiagMatrix(diagMat map[int][]complex128, index int, vec []complex128) {
	if diagMat[index] == nil {
		diagMat[index] = vec
	} else {
		diagMat[index] = add(diagMat[index], vec)
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
