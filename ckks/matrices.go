package ckks

import (
	"fmt"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
	"math"
)

// MMPt is a struct holding all the linear transformation necessary for the homomorphic
// multiplication of square matrices.
type MMPt struct {
	dimension    int
	mPermuteRows *PtDiagMatrix
	mPermuteCols *PtDiagMatrix
	mRotRows     []*PtDiagMatrix
	mRotCols     []*PtDiagMatrix
}

// GenMatMulLinTrans generates the plaintext linear transformation necessary for the homomorphic
// multiplication of square matrices.
func GenMatMulLinTrans(params *Parameters, level, dimension int, encoder Encoder) (mmpt *MMPt) {

	mmpt = new(MMPt)

	mmpt.dimension = dimension

	var scale float64
	scale = float64(params.Qi()[level]) * math.Sqrt(float64(params.Qi()[level-2])/params.Scale())

	mmpt.mPermuteRows, _ = genPermuteRowsMatrix(level, scale, 16.0, dimension, params.LogSlots(), encoder)
	mmpt.mPermuteCols, _ = genPermuteColsMatrix(level, scale, 16.0, dimension, params.LogSlots(), encoder)

	mmpt.mRotCols = make([]*PtDiagMatrix, dimension-1)
	mmpt.mRotRows = make([]*PtDiagMatrix, dimension-1)

	scale = float64(params.Qi()[level-1])

	for i := 0; i < dimension-1; i++ {
		mmpt.mRotCols[i], _ = GenSubVectorRotationMatrix(level-1, scale, dimension, i+1, params.LogSlots(), encoder)
		mmpt.mRotRows[i], _ = GenSubVectorRotationMatrix(level-1, scale, dimension*dimension, (i+1)*dimension, params.LogSlots(), encoder)

	}
	return
}

func (eval *evaluator) MulMatrix(A, B *Ciphertext, mmpt *MMPt) (ciphertextAB *Ciphertext) {

	ciphertextA := eval.LinearTransform(A, mmpt.mPermuteRows)[0]
	if err := eval.Rescale(ciphertextA, eval.params.Scale(), ciphertextA); err != nil {
		panic(err)
	}

	ciphertextB := eval.LinearTransform(B, mmpt.mPermuteCols)[0]
	if err := eval.Rescale(ciphertextB, eval.params.Scale(), ciphertextB); err != nil {
		panic(err)
	}

	ciphertextAB = eval.MulNew(ciphertextA, ciphertextB)

	alpha := eval.params.Alpha()
	beta := int(math.Ceil(float64(ciphertextA.Level()+1) / float64(alpha)))

	c2QiQDecompB := make([]*ring.Poly, beta)
	c2QiPDecompB := make([]*ring.Poly, beta)

	for i := 0; i < beta; i++ {
		c2QiQDecompB[i] = eval.ringQ.NewPolyLvl(ciphertextA.Level())
		c2QiPDecompB[i] = eval.ringP.NewPoly()
	}

	eval.DecompInternal(ciphertextA.Level(), ciphertextA.value[1], eval.c2QiQDecomp, eval.c2QiPDecomp)
	eval.DecompInternal(ciphertextB.Level(), ciphertextB.value[1], c2QiQDecompB, c2QiPDecompB)

	tmpC := NewCiphertext(eval.params, 1, ciphertextA.Level()-1, ciphertextA.Scale())

	tmpA := NewCiphertext(eval.params, 1, ciphertextA.Level(), ciphertextA.Scale())
	tmpB := NewCiphertext(eval.params, 1, ciphertextB.Level(), ciphertextB.Scale())

	tmpARescale := NewCiphertext(eval.params, 1, ciphertextA.Level()-1, ciphertextA.Scale())
	tmpBRescale := NewCiphertext(eval.params, 1, ciphertextB.Level()-1, ciphertextB.Scale())

	for i := 0; i < mmpt.dimension-1; i++ {

		eval.multiplyByDiabMatrix(ciphertextA, tmpA, mmpt.mRotCols[i], eval.c2QiQDecomp, eval.c2QiPDecomp)
		eval.multiplyByDiabMatrix(ciphertextB, tmpB, mmpt.mRotRows[i], c2QiQDecompB, c2QiPDecompB)

		if err := eval.Rescale(tmpA, eval.params.Scale(), tmpARescale); err != nil {
			panic(err)
		}

		if err := eval.Rescale(tmpB, eval.params.Scale(), tmpBRescale); err != nil {
			panic(err)
		}

		eval.Mul(tmpARescale, tmpBRescale, tmpC)

		eval.Add(ciphertextAB, tmpC, ciphertextAB)
	}

	eval.Relinearize(ciphertextAB, ciphertextAB)
	eval.Rescale(ciphertextAB, eval.params.Scale(), ciphertextAB)

	return
}

// genPermuteRowsMatrix rotates each row of the matrix by k position, where k is the row index.
func genPermuteRowsMatrix(level int, scale, maxM1N2Ratio float64, dimension int, logSlots int, encoder Encoder) (*PtDiagMatrix, map[int][]complex128) {

	slots := 1 << logSlots

	diagMatrix := make(map[int][]complex128)

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

		diagMatrix[(i+int(slots))%slots] = m
	}

	return encoder.EncodeDiagMatrixAtLvl(level, diagMatrix, scale, maxM1N2Ratio, logSlots), diagMatrix

}

// genPermuteColsMatrix rotates each column of the matrix by k position, where k is the column index.
func genPermuteColsMatrix(level int, scale, maxM1N2Ratio float64, dimension int, logSlots int, encoder Encoder) (*PtDiagMatrix, map[int][]complex128) {

	slots := 1 << logSlots

	diagMatrix := make(map[int][]complex128)

	d2 := int(dimension * dimension)

	if d2 < slots {

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

			diagMatrix[(i+int(slots))%slots] = m

		}
	} else {
		for i := 0; i < int(dimension); i++ {

			m := make([]complex128, 1<<logSlots)

			for j := 0; j < d2; j = j + int(dimension) {
				m[j+i] = 1
			}

			populateVector(m, d2, logSlots)

			diagMatrix[i*int(dimension)] = m
		}
	}

	return encoder.EncodeDiagMatrixAtLvl(level, diagMatrix, scale, maxM1N2Ratio, logSlots), diagMatrix
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
func GenSubVectorRotationMatrix(level int, scale float64, vectorSize, k int, logSlots int, encoder Encoder) (*PtDiagMatrix, map[int][]complex128) {

	k %= vectorSize

	diagMatrix := make(map[int][]complex128)

	slots := 1 << logSlots

	matrix := new(PtDiagMatrix)
	matrix.Vec = make(map[int][2]*ring.Poly)

	if vectorSize < slots {
		m0 := make([]complex128, slots)
		m1 := make([]complex128, slots)

		for i := 0; i < slots/vectorSize; i++ {

			index := i * vectorSize

			for j := 0; j < k; j++ {
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
		matrix.LogSlots = logSlots
		matrix.Level = level
		matrix.Scale = scale
		matrix.naive = true

		// Encode m0
		matrix.Vec[slots-vectorSize+k] = encoder.(*encoderComplex128).encodeDiagonal(logSlots, level, scale, m0, slots-vectorSize+k)
		// Encode m1
		matrix.Vec[k] = encoder.(*encoderComplex128).encodeDiagonal(logSlots, level, scale, m1, k)

	} else {

		matrix.rotOnly = true

		// If N = vectorSize, the we a single rotation without masking is sufficient
		matrix.Vec[k] = [2]*ring.Poly{nil, nil}

	}

	return matrix, diagMatrix
}

// GenTransposeDiagMatrix generates the linear transform plaintext vectors for the transpose of a square matrix.
func GenTransposeDiagMatrix(level int, scale, maxM1N2Ratio float64, dimension int, logSlots int, encoder Encoder) (*PtDiagMatrix, map[int][]complex128) {

	slots := 1 << logSlots

	diagMatrix := make(map[int][]complex128)

	d2 := dimension * dimension

	for i := -dimension + 1; i < dimension; i++ {

		m := make([]complex128, slots)

		if i >= 0 {
			for j := 0; j < d2-i*dimension; j = j + dimension + 1 {
				m[i+j] = 1
			}
		} else {
			for j := -i * dimension; j < d2; j = j + dimension + 1 {
				m[j] = 1
			}
		}

		populateVector(m, d2, logSlots)

		diagMatrix[(i*(dimension-1)+slots)%slots] = m
	}

	return encoder.EncodeDiagMatrixAtLvl(level, diagMatrix, scale, maxM1N2Ratio, logSlots), diagMatrix
}

func populateVector(m []complex128, d2, logSlots int) {

	slots := 1 << logSlots

	for k := d2; k < int(slots); k = k + d2 {

		if k+d2 > int(slots) {
			break
		}

		for j := 0; j < d2; j++ {
			m[k+j] = m[j]
		}
	}
}

// Matrix is a struct holding a row flatened complex matrix.
type Matrix struct {
	rows, cols int
	M          []complex128
}

// NewMatrix creates a new matrix.
func NewMatrix(rows, cols int) (m *Matrix) {
	m = new(Matrix)
	m.M = make([]complex128, rows*cols)
	m.rows = rows
	m.cols = cols
	return
}

// Rows returns the number of rows of the matrix.
func (m *Matrix) Rows() int {
	return m.rows
}

// Cols returns the number of columns of the matrix.
func (m *Matrix) Cols() int {
	return m.cols
}

// Add adds matrix A and B and stores the result on the target.
func (m *Matrix) Add(A, B *Matrix) {

	if len(A.M) != len(B.M) {
		panic("input matrices are incompatible for addition")
	}

	if m.M == nil {
		m.M = make([]complex128, A.Rows()*B.Rows())
	} else if len(m.M) > len(A.M) {
		m.M = m.M[:len(A.M)]
	} else if len(m.M) < len(A.M) {
		m.M = append(m.M, make([]complex128, len(A.M)-len(m.M))...)
	}

	for i := range A.M {
		m.M[i] = A.M[i] + B.M[i]
	}
}

// MulMat multiplies A with B and returns the result on the target.
func (m *Matrix) MulMat(A, B *Matrix) {

	if A.Cols() != B.Rows() {
		panic("matrices are incompatible for multiplication")
	}

	rowsA := A.Rows()
	colsA := A.Cols()
	colsB := B.Cols()

	acc := make([]complex128, rowsA*colsB)

	for i := 0; i < rowsA; i++ {
		for j := 0; j < colsB; j++ {
			for k := 0; k < colsA; k++ {
				acc[i*colsA+j] += A.M[i*colsA+k] * B.M[j+k*colsB]
			}
		}
	}

	if len(m.M) < rowsA*colsB {
		m.M = append(m.M, make([]complex128, m.Rows()*m.Cols()-rowsA*colsB)...)
	} else {
		m.M = m.M[:rowsA*colsB]
	}

	for i := range m.M {
		m.M[i] = acc[i]
	}
}

// GenRandomComplexMatrices generates a list of complex matrices.
func GenRandomComplexMatrices(rows, cols, n int) (Matrices []*Matrix) {

	Matrices = make([]*Matrix, n)

	for k := range Matrices {
		m := NewMatrix(rows, cols)
		for i := 0; i < rows*cols; i++ {
			m.M[i] = complex(utils.RandFloat64(-1, 1), utils.RandFloat64(-1, 1))
		}
		Matrices[k] = m
	}

	return
}

// PermuteRows rotates each row by k where k is the row index.
// Equivalent to Transpoe(PermuteCols(Transpose(M)))
func (m *Matrix) PermuteRows() {
	var index int
	tmp := make([]complex128, m.Cols())
	for i := 0; i < m.Rows(); i++ {
		index = i * m.Cols()
		for j := range tmp {
			tmp[j] = m.M[index+j]
		}

		tmp = append(tmp[i:], tmp[:i]...)

		for j, c := range tmp {
			m.M[index+j] = c
		}
	}
}

// PermuteCols rotates each column by k, where k is the column index.
// Equivalent to Transpoe(PermuteRows(Transpose(M)))
func (m *Matrix) PermuteCols() {
	tmp := make([]complex128, m.Rows())
	for i := 0; i < m.Cols(); i++ {
		for j := range tmp {
			tmp[j] = m.M[i+j*m.Cols()]
		}

		tmp = append(tmp[i:], tmp[:i]...)

		for j, c := range tmp {
			m.M[i+j*m.Cols()] = c
		}
	}
}

// RotateCols rotates each column by k position to the left.
func (m *Matrix) RotateCols(k int) {

	k %= m.Cols()
	var index int
	tmp := make([]complex128, m.Cols())
	for i := 0; i < m.Rows(); i++ {
		index = i * m.Cols()
		for j := range tmp {
			tmp[j] = m.M[index+j]
		}

		tmp = append(tmp[k:], tmp[:k]...)

		for j, c := range tmp {
			m.M[index+j] = c
		}
	}
}

// RotateRows rotates each row by k positions to the left.
func (m *Matrix) RotateRows(k int) {
	k %= m.Rows()
	m.M = append(m.M[k*m.Cols():], m.M[:k*m.Cols()]...)
}

// Transpose transposes the matrix.
func (m *Matrix) Transpose() (mT *Matrix) {
	rows := m.Rows()
	cols := m.Cols()
	mT = NewMatrix(cols, rows)

	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			mT.M[cols*i+j] = m.M[j*cols+i]
		}
	}
	return
}

// Print prints the target matrix.
func (m *Matrix) Print() {

	for i := 0; i < m.Cols(); i++ {
		for j := 0; j < m.Rows(); j++ {
			fmt.Printf("%7.4f ", m.M[i*m.Cols()+j])
		}
		fmt.Printf("\n")
	}
	fmt.Printf("\n")
}
