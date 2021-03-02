package ckks

import (
	"fmt"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/stretchr/testify/require"
	"math"
	"math/rand"
	"testing"
	"time"
)

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
	rlk := kgen.GenRelinearizationKey(sk)
	encryptor := NewEncryptorFromSk(params, sk)
	decryptor := NewDecryptor(params, sk)
	eval := NewEvaluator(params, EvaluationKey{rlk, nil})

	// Size of the matrices (dxd)
	rows := 16
	cols := 16

	t.Run("Transpose/", func(t *testing.T) {

		m, _, ct := GenTestVectors(rows, cols, params, encoder, encryptor)

		diagMatrix, _ := GenTransposeDiagMatrix(params.MaxLevel(), float64(params.Qi()[params.MaxLevel()]), 16.0, rows, params.LogSlots(), encoder)

		rotKeys := kgen.GenRotKeysForDiagMatrix(diagMatrix, false, sk)

		eval0 := eval.ShallowCopyWithKey(EvaluationKey{rlk, rotKeys})

		for i := range m {
			m[i] = m[i].Transpose()
		}

		//PrintDebug(ct, rows, cols, params, encoder, decryptor)

		start := time.Now()
		ct = eval0.LinearTransform(ct, diagMatrix)[0]
		fmt.Println("Done :", time.Since(start))

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

		for k := 1; k < rows; k++ {

			t.Run(fmt.Sprintf("k=%d/", k), func(t *testing.T) {

				level := params.MaxLevel()

				diagMatrix, _ := GenSubVectorRotationMatrix(level, float64(params.Qi()[level]), rows, k, params.LogSlots(), encoder)

				rotKeys := kgen.GenRotKeysForDiagMatrix(diagMatrix, false, sk)

				eval0 := eval.ShallowCopyWithKey(EvaluationKey{rlk, rotKeys})

				for j := range m {
					m[j].RotateCols(1)
				}

				eval0.(*evaluator).multiplyByDiabMatrix(ct, res, diagMatrix, c2QiQDecompA, c2QiPDecompA)

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

		for k := 1; k < rows; k++ {

			t.Run(fmt.Sprintf("k=%d/", k), func(t *testing.T) {

				level := params.MaxLevel()

				diagMatrix, _ := GenSubVectorRotationMatrix(level, float64(params.Qi()[level]), rows*rows, k*rows, params.LogSlots(), encoder)

				rotKeys := kgen.GenRotKeysForDiagMatrix(diagMatrix, false, sk)

				eval0 := eval.ShallowCopyWithKey(EvaluationKey{rlk, rotKeys})

				for j := range m {
					m[j].RotateRows(1)
				}

				eval0.(*evaluator).multiplyByDiabMatrix(ct, res, diagMatrix, c2QiQDecompA, c2QiPDecompA)

				VerifyTestVectors(params, encoder, decryptor, m, res, t)
			})
		}
	})

	t.Run("PermuteRows/", func(t *testing.T) {
		m, _, ct := GenTestVectors(rows, cols, params, encoder, encryptor)

		level := params.MaxLevel()

		diagMatrix, _ := genPermuteRowsMatrix(level, float64(params.Qi()[level]), 16.0, rows, params.LogSlots(), encoder)

		rotKeys := kgen.GenRotKeysForDiagMatrix(diagMatrix, false, sk)

		eval0 := eval.ShallowCopyWithKey(EvaluationKey{rlk, rotKeys})

		for j := range m {
			m[j].PermuteRows()
		}

		ct = eval0.LinearTransform(ct, diagMatrix)[0]

		//PrintDebug(ct, rows, cols, params, encoder, decryptor)

		VerifyTestVectors(params, encoder, decryptor, m, ct, t)

	})

	t.Run("PermuteCols/", func(t *testing.T) {
		m, _, ct := GenTestVectors(rows, cols, params, encoder, encryptor)

		level := params.MaxLevel()

		diagMatrix, _ := genPermuteColsMatrix(level, float64(params.Qi()[level]), 16.0, rows, params.LogSlots(), encoder)

		rotKeys := kgen.GenRotKeysForDiagMatrix(diagMatrix, false, sk)

		eval0 := eval.ShallowCopyWithKey(EvaluationKey{rlk, rotKeys})

		for j := range m {
			m[j].PermuteCols()
		}

		//PrintDebug(ct, d, params, encoder, decryptor)

		ct = eval0.LinearTransform(ct, diagMatrix)[0]

		//PrintDebug(ct, d, params, encoder, decryptor)

		VerifyTestVectors(params, encoder, decryptor, m, ct, t)

	})

	t.Run("Multiply/", func(t *testing.T) {
		mA, _, ctA := GenTestVectors(rows, cols, params, encoder, encryptor)
		mB, _, ctB := GenTestVectors(rows, cols, params, encoder, encryptor)

		mmpt := GenMatMulLinTrans(params, params.MaxLevel(), rows, encoder)
		rotKeys := kgen.GenMatMulRotKeys(mmpt, sk)

		eval0 := eval.ShallowCopyWithKey(EvaluationKey{rlk, rotKeys})

		for j := range mA {
			mA[j].MulMat(mA[j], mB[j])
		}

		start := time.Now()
		ctAB := eval0.MulMatrix(ctA, ctB, mmpt)
		fmt.Println("Done :", time.Since(start))

		//PrintDebug(ctAB, d, params, encoder, decryptor)

		VerifyTestVectors(params, encoder, decryptor, mA, ctAB, t)

	})

}

func MatricesToVector(m []*Matrix, params *Parameters) (values []complex128) {
	values = make([]complex128, params.Slots())

	d := m[0].Rows() * m[0].Cols()

	for i := 0; i < int(params.Slots())/d; i++ {

		for j, c := range m[i].M {
			values[i*d+j] = c
		}
	}
	return
}

func GenTestVectors(rows, cols int, params *Parameters, encoder Encoder, encryptor Encryptor) (m []*Matrix, pt *Plaintext, ct *Ciphertext) {

	m = GenRandomComplexMatrices(rows, cols, int(params.Slots())/(rows*cols))

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
