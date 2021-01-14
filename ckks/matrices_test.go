package ckks

import (
	"fmt"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
	"github.com/stretchr/testify/require"
	"math"
	"math/rand"
	"testing"
	"time"
)

type Matrix struct {
	rows, cols uint64
	M          []complex128
}

func NewMatrix(rows, cols uint64) (m *Matrix) {
	m = new(Matrix)
	m.M = make([]complex128, rows*cols)
	m.rows = rows
	m.cols = cols
	return
}

func (m *Matrix) Rows() uint64 {
	return m.rows
}

func (m *Matrix) Cols() uint64 {
	return m.cols
}

func MulMat(A, B *Matrix) (AB *Matrix) {
	rowsA := A.Rows()
	colsA := A.Cols()
	colsB := B.Cols()

	AB = NewMatrix(rowsA, colsB)

	for i := uint64(0); i < rowsA; i++ {
		for j := uint64(0); j < colsB; j++ {
			for k := uint64(0); k < colsA; k++ {
				AB.M[i*colsA+j] += A.M[i*colsA+k] * B.M[j+k*colsB]
			}
		}
	}

	return
}

func GenMatrices(rows, cols, n uint64) (Matrices []*Matrix) {

	Matrices = make([]*Matrix, n)

	for k := range Matrices {
		m := NewMatrix(rows, cols)
		for i := uint64(0); i < rows*cols; i++ {
			m.M[i] = complex(utils.RandFloat64(-1, 1), 0) //float64(i*int(d) + j+1  + k*int(d*d)), 0) //randomComplex(-1, 1)float64(i*int(d) + j+1), 0))
		}
		Matrices[k] = m
	}

	return
}

// Rotates each row by k where k is the row index.
// Equivalent to Transpoe(PermuteB(Transpose(M)))
func (m *Matrix) PermuteA() {
	var index uint64
	tmp := make([]complex128, m.Cols())
	for i := uint64(0); i < m.Rows(); i++ {
		index = i * m.Cols()
		for j := range tmp {
			tmp[j] = m.M[index+uint64(j)]
		}

		tmp = append(tmp[i:], tmp[:i]...)

		for j, c := range tmp {
			m.M[index+uint64(j)] = c
		}
	}
}

// Rotates each column by k, where k is the column index.
// Equivalent to Transpoe(PermuteA(Transpose(M)))
func (m *Matrix) PermuteB() {
	tmp := make([]complex128, m.Rows())
	for i := uint64(0); i < m.Cols(); i++ {
		for j := range tmp {
			tmp[j] = m.M[i+uint64(j)*m.Cols()]
		}

		tmp = append(tmp[i:], tmp[:i]...)

		for j, c := range tmp {
			m.M[i+uint64(j)*m.Cols()] = c
		}
	}
}

// RotateCols rotates each column by k position to the left.
func (m *Matrix) RotateCols(k uint64) {

	k %= m.Cols()
	var index uint64
	tmp := make([]complex128, m.Cols())
	for i := uint64(0); i < m.Rows(); i++ {
		index = i * m.Cols()
		for j := range tmp {
			tmp[j] = m.M[index+uint64(j)]
		}

		tmp = append(tmp[k:], tmp[:k]...)

		for j, c := range tmp {
			m.M[index+uint64(j)] = c
		}
	}
}

func (m *Matrix) RotateRows(k uint64) {
	k %= m.Rows()
	m.M = append(m.M[k*m.Cols():], m.M[:k*m.Cols()]...)
}

func (m *Matrix) Transpose() (mT *Matrix) {
	rows := m.Rows()
	cols := m.Cols()
	mT = NewMatrix(cols, rows)

	for i := uint64(0); i < rows; i++ {
		for j := uint64(0); j < cols; j++ {
			mT.M[cols*i+j] = m.M[j*cols+i]
		}
	}
	return
}

func (m *Matrix) Print() {

	for i := uint64(0); i < m.Cols(); i++ {
		for j := uint64(0); j < m.Rows(); j++ {
			fmt.Printf("%7.4f ", m.M[i*m.Cols()+j])
		}
		fmt.Printf("\n")
	}
	fmt.Printf("\n")
}

func TestMatrices(t *testing.T) {

	rand.Seed(time.Now().UnixNano())

	LogN := uint64(13)
	LogSlots := uint64(12)

	LogModuli := LogModuli{
		LogQi: []uint64{50, 35, 35, 35},
		LogPi: []uint64{55},
	}

	Scale := float64(1 << 35)

	params, err := NewParametersFromLogModuli(LogN, &LogModuli)
	if err != nil {
		panic(err)
	}
	params.SetScale(Scale)
	params.SetLogSlots(LogSlots)

	encoder := NewEncoder(params)
	kgen := NewKeyGenerator(params)
	sk := kgen.GenSecretKey()
	rlk := kgen.GenRelinKey(sk)
	rotKeys := NewRotationKeys()
	encryptor := NewEncryptorFromSk(params, sk)
	decryptor := NewDecryptor(params, sk)
	eval := NewEvaluator(params)

	// Size of the matrices (dxd)
	rows := uint64(64)
	cols := uint64(64)

	t.Run("Transpose/", func(t *testing.T) {

		m, _, ct := GenTestVectors(rows, cols, params, encoder, encryptor)

		diagMatrix, _ := GenTransposeDiagMatrix(params.MaxLevel(), float64(params.Qi()[params.MaxLevel()]), 16.0, rows, params.LogSlots(), encoder)

		kgen.GenRotKeysForDiagMatrix(diagMatrix, sk, rotKeys)

		for i := range m {
			m[i] = m[i].Transpose()
		}

		//PrintDebug(ct, rows, cols, params, encoder, decryptor)

		ct = eval.LinearTransform(ct, diagMatrix, rotKeys)[0]

		VerifyTestVectors(params, encoder, decryptor, m, ct, t)

	})

	t.Run("RotateCols", func(t *testing.T) {

		m, _, ct := GenTestVectors(rows, cols, params, encoder, encryptor)

		alpha := params.Alpha()
		beta := uint64(math.Ceil(float64(ct.Level()+1) / float64(alpha)))

		c2QiQDecompA := make([]*ring.Poly, beta)
		c2QiPDecompA := make([]*ring.Poly, beta)

		for i := uint64(0); i < beta; i++ {
			c2QiQDecompA[i] = eval.(*evaluator).ringQ.NewPolyLvl(ct.Level())
			c2QiPDecompA[i] = eval.(*evaluator).ringP.NewPoly()
		}

		eval.DecompInternal(ct.Level(), ct.value[1], c2QiQDecompA, c2QiPDecompA)

		res := NewCiphertext(params, 1, ct.Level(), ct.Scale())

		for k := uint64(1); k < rows; k++ {

			t.Run(fmt.Sprintf("k=%d/", k), func(t *testing.T) {

				level := params.MaxLevel()

				diagMatrix, _ := GenSubVectorRotationMatrix(level, float64(params.Qi()[level]), rows, k, params.LogSlots(), encoder)

				kgen.GenRotKeysForDiagMatrix(diagMatrix, sk, rotKeys)

				for j := range m {
					m[j].RotateCols(1)
				}

				eval.(*evaluator).multiplyByDiabMatrix(ct, res, diagMatrix, rotKeys, c2QiQDecompA, c2QiPDecompA)

				VerifyTestVectors(params, encoder, decryptor, m, res, t)
			})
		}
	})

	t.Run("RotateRows", func(t *testing.T) {
		m, _, ct := GenTestVectors(rows, cols, params, encoder, encryptor)

		alpha := params.Alpha()
		beta := uint64(math.Ceil(float64(ct.Level()+1) / float64(alpha)))

		c2QiQDecompA := make([]*ring.Poly, beta)
		c2QiPDecompA := make([]*ring.Poly, beta)

		for i := uint64(0); i < beta; i++ {
			c2QiQDecompA[i] = eval.(*evaluator).ringQ.NewPolyLvl(ct.Level())
			c2QiPDecompA[i] = eval.(*evaluator).ringP.NewPoly()
		}

		eval.DecompInternal(ct.Level(), ct.value[1], c2QiQDecompA, c2QiPDecompA)

		res := NewCiphertext(params, 1, ct.Level(), ct.Scale())

		for k := uint64(1); k < rows; k++ {

			t.Run(fmt.Sprintf("k=%d/", k), func(t *testing.T) {

				level := params.MaxLevel()

				diagMatrix, _ := GenSubVectorRotationMatrix(level, float64(params.Qi()[level]), rows*rows, k*rows, params.LogSlots(), encoder)

				kgen.GenRotKeysForDiagMatrix(diagMatrix, sk, rotKeys)

				for j := range m {
					m[j].RotateRows(1)
				}

				eval.(*evaluator).multiplyByDiabMatrix(ct, res, diagMatrix, rotKeys, c2QiQDecompA, c2QiPDecompA)

				VerifyTestVectors(params, encoder, decryptor, m, res, t)
			})
		}
	})

	t.Run("PermuteA/", func(t *testing.T) {
		m, _, ct := GenTestVectors(rows, cols, params, encoder, encryptor)

		level := params.MaxLevel()

		diagMatrix, _ := GenPermuteAMatrix(level, float64(params.Qi()[level]), 16.0, rows, params.LogSlots(), encoder)

		kgen.GenRotKeysForDiagMatrix(diagMatrix, sk, rotKeys)

		for j := range m {
			m[j].PermuteA()
		}

		ct = eval.LinearTransform(ct, diagMatrix, rotKeys)[0]

		//PrintDebug(ct, rows, cols, params, encoder, decryptor)

		VerifyTestVectors(params, encoder, decryptor, m, ct, t)

	})

	t.Run("PermuteB/", func(t *testing.T) {
		m, _, ct := GenTestVectors(rows, cols, params, encoder, encryptor)

		level := params.MaxLevel()

		diagMatrix, _ := GenPermuteBMatrix(level, float64(params.Qi()[level]), 16.0, rows, params.LogSlots(), encoder)

		kgen.GenRotKeysForDiagMatrix(diagMatrix, sk, rotKeys)

		for j := range m {
			m[j].PermuteB()
		}

		//PrintDebug(ct, d, params, encoder, decryptor)

		ct = eval.LinearTransform(ct, diagMatrix, rotKeys)[0]

		//PrintDebug(ct, d, params, encoder, decryptor)

		VerifyTestVectors(params, encoder, decryptor, m, ct, t)

	})

	t.Run("Multiply/Square", func(t *testing.T) {
		mA, _, ctA := GenTestVectors(rows, cols, params, encoder, encryptor)
		mB, _, ctB := GenTestVectors(rows, cols, params, encoder, encryptor)

		mmpt := GenPlaintextMatrices(params, params.MaxLevel(), rows, encoder)
		GenRotationKeys(mmpt, kgen, sk, rotKeys)

		for j := range mA {
			mA[j] = MulMat(mA[j], mB[j])
		}

		start := time.Now()
		ctAB := eval.MulMatrixAB(ctA, ctB, mmpt, rlk, rotKeys)
		fmt.Println("Done :", time.Since(start))

		//PrintDebug(ctAB, d, params, encoder, decryptor)

		VerifyTestVectors(params, encoder, decryptor, mA, ctAB, t)

	})

}

func MatricesToVector(m []*Matrix, params *Parameters) (values []complex128) {
	values = make([]complex128, params.Slots())

	d := m[0].Rows() * m[0].Cols()

	for i := uint64(0); i < params.Slots()/d; i++ {

		for j, c := range m[i].M {
			values[i*d+uint64(j)] = c
		}
	}
	return
}

func GenTestVectors(rows, cols uint64, params *Parameters, encoder Encoder, encryptor Encryptor) (m []*Matrix, pt *Plaintext, ct *Ciphertext) {

	m = GenMatrices(rows, cols, params.Slots()/(rows*cols))

	values := MatricesToVector(m, params)

	pt = encoder.EncodeNew(values, params.LogSlots())

	ct = encryptor.EncryptNew(pt)

	return m, pt, ct
}

func VerifyTestVectors(params *Parameters, encoder Encoder, decryptor Decryptor, m []*Matrix, element interface{}, t *testing.T) {

	precStats := GetPrecisionStats(params, encoder, decryptor, MatricesToVector(m, params), element, 0)

	if *printPrecisionStats {
		t.Log(precStats.String())
	}

	require.GreaterOrEqual(t, real(precStats.MeanPrecision), minPrec)
	require.GreaterOrEqual(t, imag(precStats.MeanPrecision), minPrec)
}

func PrintDebug(ct *Ciphertext, rows, cols uint64, params *Parameters, encoder Encoder, decryptor Decryptor) {

	valuesHave := encoder.Decode(decryptor.DecryptNew(ct), params.LogSlots())

	maxPrint := params.Slots() / (rows * cols)

	if maxPrint > 4 {
		maxPrint = 4
	}
	for k := uint64(0); k < maxPrint; k++ {

		index := k * rows * cols

		for i := uint64(0); i < rows; i++ {
			for j := uint64(0); j < cols; j++ {
				fmt.Printf("%7.4f ", real(valuesHave[index+i*rows+j]))
			}
			fmt.Printf("\n")
		}
		fmt.Printf("\n")
	}
}
