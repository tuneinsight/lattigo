package ckks

import (
	"github.com/ldsec/lattigo/v2/ring"
	"math"
)

type MatrixMultiplier interface {
	GenPlaintextMatrices(params *Parameters, level, d uint64, encoder Encoder) (mmpt *MMPt)
	GenRotationKeys(mmpt *MMPt, kgen KeyGenerator, sk *SecretKey, rotKeys *RotationKeys)
	GenTransposeDiagMatrix(d, logSlots uint64) (diagMatrix map[uint64][]complex128)
	GenPermuteAMatrix(d, logSlots uint64) (diagMatrix map[uint64][]complex128)
	GenPermuteBMatrix(d, logSlots uint64) (diagMatrix map[uint64][]complex128)
	GenSubVectorRotationMatrix(lvl uint64, scale float64, vectorSize, k, logSlots uint64, encoder Encoder) (*PtDiagMatrix, map[uint64][]complex128)
}

type matrixMultiplier struct {
}

type MMPt struct {
	dimension uint64
	mPermuteA *PtDiagMatrix
	mPermuteB *PtDiagMatrix
	mRotRows  []*PtDiagMatrix
	mRotCols  []*PtDiagMatrix
}

func NewMatrixMultiplier(params *Parameters) MatrixMultiplier {

	return &matrixMultiplier{}
}

func (mm *matrixMultiplier) GenPlaintextMatrices(params *Parameters, level uint64, dimension uint64, encoder Encoder) (mmpt *MMPt) {

	mmpt = new(MMPt)

	mmpt.dimension = dimension

	var scale float64
	scale = float64(params.Qi()[level]) * math.Sqrt(float64(params.Qi()[level-2])/params.Scale())

	mmpt.mPermuteA = encoder.EncodeDiagMatrixAtLvl(level, mm.GenPermuteAMatrix(dimension, params.LogSlots()), scale, 16.0, params.LogSlots())
	mmpt.mPermuteB = encoder.EncodeDiagMatrixAtLvl(level, mm.GenPermuteBMatrix(dimension, params.LogSlots()), scale, 16.0, params.LogSlots())

	mmpt.mRotCols = make([]*PtDiagMatrix, dimension-1)
	mmpt.mRotRows = make([]*PtDiagMatrix, dimension-1)

	scale = float64(params.Qi()[level-1])

	for i := uint64(0); i < dimension-1; i++ {
		mmpt.mRotCols[i], _ = mm.GenSubVectorRotationMatrix(level-1, scale, dimension, i+1, params.LogSlots(), encoder)
		mmpt.mRotRows[i], _ = mm.GenSubVectorRotationMatrix(level-1, scale, dimension*dimension, (i+1)*dimension, params.LogSlots(), encoder)

	}
	return
}

func (mm *matrixMultiplier) GenRotationKeys(mmpt *MMPt, kgen KeyGenerator, sk *SecretKey, rotKeys *RotationKeys) {

	kgen.GenRotKeysForDiagMatrix(mmpt.mPermuteA, sk, rotKeys)
	kgen.GenRotKeysForDiagMatrix(mmpt.mPermuteB, sk, rotKeys)

	for i := range mmpt.mRotCols {
		kgen.GenRotKeysForDiagMatrix(mmpt.mRotCols[i], sk, rotKeys)
		kgen.GenRotKeysForDiagMatrix(mmpt.mRotRows[i], sk, rotKeys)
	}
}

func (eval *evaluator) MulMatrixAB(A, B *Ciphertext, mmpt *MMPt, rlk *EvaluationKey, rotKeys *RotationKeys) (ciphertextAB *Ciphertext) {

	ciphertextA := eval.LinearTransform(A, mmpt.mPermuteA, rotKeys)[0]
	eval.Rescale(ciphertextA, eval.params.Scale(), ciphertextA)

	ciphertextB := eval.LinearTransform(B, mmpt.mPermuteB, rotKeys)[0]
	eval.Rescale(ciphertextB, eval.params.Scale(), ciphertextB)

	ciphertextAB = eval.MulRelinNew(ciphertextA, ciphertextB, nil)

	alpha := eval.params.Alpha()
	beta := uint64(math.Ceil(float64(ciphertextA.Level()+1) / float64(alpha)))

	c2QiQDecompB := make([]*ring.Poly, beta)
	c2QiPDecompB := make([]*ring.Poly, beta)

	for i := uint64(0); i < beta; i++ {
		c2QiQDecompB[i] = eval.ringQ.NewPolyLvl(ciphertextA.Level())
		c2QiPDecompB[i] = eval.ringP.NewPoly()
	}

	eval.DecompInternal(ciphertextA.Level(), ciphertextA.value[1], eval.c2QiQDecomp, eval.c2QiPDecomp)
	eval.DecompInternal(ciphertextB.Level(), ciphertextB.value[1], c2QiQDecompB, c2QiPDecompB)

	tmpC := NewCiphertext(eval.params, 2, ciphertextA.Level()-1, ciphertextA.Scale())

	tmpA := NewCiphertext(eval.params, 1, ciphertextA.Level(), ciphertextA.Scale())
	tmpB := NewCiphertext(eval.params, 1, ciphertextB.Level(), ciphertextB.Scale())

	tmpARescale := NewCiphertext(eval.params, 1, ciphertextA.Level()-1, ciphertextA.Scale())
	tmpBRescale := NewCiphertext(eval.params, 1, ciphertextB.Level()-1, ciphertextB.Scale())

	for i := uint64(0); i < mmpt.dimension-1; i++ {

		eval.multiplyByDiabMatrix(ciphertextA, tmpA, mmpt.mRotCols[i], rotKeys, eval.c2QiQDecomp, eval.c2QiPDecomp)
		eval.multiplyByDiabMatrix(ciphertextB, tmpB, mmpt.mRotRows[i], rotKeys, c2QiQDecompB, c2QiPDecompB)

		if err := eval.Rescale(tmpA, eval.params.Scale(), tmpARescale); err != nil {
			panic(err)
		}

		if err := eval.Rescale(tmpB, eval.params.Scale(), tmpBRescale); err != nil {
			panic(err)
		}

		eval.MulRelin(tmpARescale, tmpBRescale, nil, tmpC)

		eval.Add(ciphertextAB, tmpC, ciphertextAB)
	}

	eval.Relinearize(ciphertextAB, rlk, ciphertextAB)
	eval.Rescale(ciphertextAB, eval.params.Scale(), ciphertextAB)

	return
}

func (mm *matrixMultiplier) GenPermuteAMatrix(dimension, logSlots uint64) (diagMatrix map[uint64][]complex128) {

	slots := uint64(1 << logSlots)

	diagMatrix = make(map[uint64][]complex128)

	d2 := int(dimension * dimension)

	for i := -int(dimension) + 1; i < int(dimension); i++ {

		m := make([]complex128, slots)

		for k := 0; k < d2; k++ {

			if i < 0 {
				for j := i; j < int(dimension); j++ {
					x := (d2 + k - (int(dimension)+i)*int(dimension)) % d2
					if x < int(dimension) && x >= -i {
						m[k] = 1
					}
				}
			} else {

				for j := i; j < int(dimension); j++ {
					if (d2+k-int(dimension)*i)%d2 < int(dimension)-i {
						m[k] = 1
					}
				}
			}
		}

		populateVector(m, d2, logSlots)

		diagMatrix[uint64((i+int(slots)))%slots] = m
	}

	return

}

func (mm *matrixMultiplier) GenPermuteBMatrix(dimension, logSlots uint64) (diagMatrix map[uint64][]complex128) {

	slots := uint64(1 << logSlots)

	diagMatrix = make(map[uint64][]complex128)

	d2 := int(dimension * dimension)

	if uint64(d2) < slots {

		for i := -int((dimension - 1) * dimension); i < d2; i = i + int(dimension) {

			m := make([]complex128, 1<<logSlots)

			if i >= 0 {
				for j := 0; j < d2-i; j = j + int(dimension) {
					m[i/int(dimension)+j] = 1
				}
			} else {
				for j := 0; j < d2+i; j = j + int(dimension) {
					m[-i+int(dimension)+(i/int(dimension))+j] = 1
				}
			}

			populateVector(m, d2, logSlots)

			diagMatrix[uint64((i+int(slots)))%slots] = m
		}
	} else {
		for i := 0; i < int(dimension); i++ {

			m := make([]complex128, 1<<logSlots)

			for j := 0; j < d2; j = j + int(dimension) {
				m[j+i] = 1
			}

			populateVector(m, d2, logSlots)

			/*
				fmt.Printf("%3d %2d", i, uint64(i)*d)
				for i := range m[:d2]{
					fmt.Printf("%2.f ", real(m[i]))
				}
				fmt.Println()
			*/

			diagMatrix[uint64(i)*dimension] = m
		}
	}

	return
}

// GenSubVectorRotationMatrix allows to generate a permutation matrix that roates subvectors independently.
// Given a vector of size N=2^"logSlots", partitionned into N/"vectorSize" subvectors each of size "vectorSize",
// rotates each subvector by "k" positions to the left.
//
// Example :
// Given v = [a_(0), a_(1), a_(2), ..., a_(N-3), a_(N-2), a_(N-1)],
// Then M x v = [rotate(a_(0), a_(1), ..., a_(vectorsize-1), k), ... , rotate(a_(N-vectorsize-1), a_(N-vectorsize), ..., a_(N-1), k)]
//
// If vectorSize does not divide N, then the last N%vectorSize slots are zero.
// If N = vectorSize, then no mask is generated and the evaluation is instead a single rotation.
//
// This is done by generating the two masks :
//       	 |     vectorsize     |, ..., |     vectorsize     |
// mask_0 = [{1, ..., 1, 0, ..., 0}, ..., {1, ..., 1, 0, ..., 0}]
// mask_1 = [{0, ..., 0, 1, ..., 1}, ..., {0, ..., 0, 1, ..., 1}]
//            0 ----- k                    0 ----- k
func (mm *matrixMultiplier) GenSubVectorRotationMatrix(level uint64, scale float64, vectorSize, k, logSlots uint64, encoder Encoder) (*PtDiagMatrix, map[uint64][]complex128) {

	k %= vectorSize

	diagMatrix := make(map[uint64][]complex128)

	slots := uint64(1 << logSlots)

	matrix := new(PtDiagMatrix)
	matrix.Vec = make(map[uint64][2]*ring.Poly)

	if vectorSize < slots {
		m0 := make([]complex128, slots)
		m1 := make([]complex128, slots)

		for i := uint64(0); i < slots/vectorSize; i++ {

			index := i * vectorSize

			for j := uint64(0); j < k; j++ {
				m0[j+index] = 1
			}

			for j := k; j < vectorSize; j++ {
				m1[j+index] = 1
			}
		}

		/*
			fmt.Printf("%4d", (slots) - (vectorSize - k))
			for i := range m0[:vectorSize]{
				fmt.Printf("%2.f ", real(m0[i]))
			}
			fmt.Println()

			fmt.Printf("%4d", k)
			for i := range m1[:vectorSize]{
				fmt.Printf("%2.f ", real(m1[i]))
			}
			fmt.Println()
		*/

		diagMatrix[slots-vectorSize+k] = m0
		diagMatrix[k] = m1

		// Encoding
		ringQ := encoder.(*encoderComplex128).ringQ
		ringP := encoder.(*encoderComplex128).ringP
		N := ringQ.N

		matrix.LogSlots = logSlots
		matrix.Level = level
		matrix.Scale = scale
		matrix.naive = true

		// Encode m0
		encoder.Embed(rotate(m0, slots-vectorSize+k), logSlots)

		m0Q := ring.NewPoly(N, level+1)
		encoder.ScaleUp(m0Q, scale, ringQ.Modulus[:level+1])
		ringQ.NTTLvl(level, m0Q, m0Q)
		ringQ.MFormLvl(level, m0Q, m0Q)

		m0P := ring.NewPoly(N, level+1)
		encoder.ScaleUp(m0P, scale, ringP.Modulus)
		ringP.NTT(m0P, m0P)
		ringP.MForm(m0P, m0P)

		matrix.Vec[slots-vectorSize+k] = [2]*ring.Poly{m0Q, m0P}

		encoder.WipeInternalMemory()

		// Encode m1
		encoder.Embed(rotate(m1, k), logSlots)

		m1Q := ring.NewPoly(N, level+1)
		encoder.ScaleUp(m1Q, scale, ringQ.Modulus[:level+1])
		ringQ.NTTLvl(level, m1Q, m1Q)
		ringQ.MFormLvl(level, m1Q, m1Q)

		m1P := ring.NewPoly(N, level+1)
		encoder.ScaleUp(m1P, scale, ringP.Modulus)
		ringP.NTT(m1P, m1P)
		ringP.MForm(m1P, m1P)

		matrix.Vec[k] = [2]*ring.Poly{m1Q, m1P}

		encoder.WipeInternalMemory()

	} else {

		matrix.rotOnly = true

		// If N = vectorSize, the we a single rotation without masking is sufficient
		matrix.Vec[k] = [2]*ring.Poly{nil, nil}

	}

	return matrix, diagMatrix
}

func (mm *matrixMultiplier) GenTransposeDiagMatrix(dimension, logSlots uint64) (diagMatrix map[uint64][]complex128) {

	slots := uint64(1 << logSlots)

	diagMatrix = make(map[uint64][]complex128)

	d2 := int(dimension * dimension)

	for i := -int(dimension) + 1; i < int(dimension); i++ {

		m := make([]complex128, slots)

		if i >= 0 {
			for j := 0; j < d2-i*int(dimension); j = j + int(dimension) + 1 {
				m[i+j] = 1
			}
		} else {
			for j := -i * int(dimension); j < d2; j = j + int(dimension) + 1 {
				m[j] = 1
			}
		}

		populateVector(m, d2, logSlots)

		diagMatrix[uint64(i*int(dimension-1)+int(slots))%slots] = m
	}

	return
}

func populateVector(m []complex128, d2 int, logSlots uint64) {

	slots := uint64(1 << logSlots)

	for k := d2; k < int(slots); k = k + d2 {

		if k+d2 > int(slots) {
			break
		}

		for j := 0; j < d2; j++ {
			m[k+j] = m[j]
		}
	}
}
