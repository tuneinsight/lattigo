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
	M [][]complex128
}

func (m *Matrix) N() uint64 {
	return uint64(len(m.M))
}

func MulMat(A, B *Matrix) (AB *Matrix) {
	n := A.N()
	AB = new(Matrix)
	m := make([][]complex128, n)

	for i := uint64(0); i < n; i++ {
		m[i] = make([]complex128, n)
		for j := uint64(0); j < n; j++ {
			for k := uint64(0); k < n; k++ {
				m[i][j] += A.M[i][k] * B.M[k][j]
			}
		}
	}

	AB.M = m
	return
}

func GenMatrices(d, n uint64) (Matrices []*Matrix) {

	Matrices = make([]*Matrix, n)

	for k := range Matrices {
		m := new(Matrix)
		m.M = make([][]complex128, d)
		for i := range m.M {
			m.M[i] = make([]complex128, d)
			for j := range m.M[i] {
				m.M[i][j] = complex(utils.RandFloat64(-1, 1), 0) //float64(i*int(d) + j+1  + k*int(d*d)), 0) //randomComplex(-1, 1)float64(i*int(d) + j+1), 0)) // //
			}
		}
		Matrices[k] = m
	}

	return
}

func (m *Matrix) PermuteA() {
	for i := range m.M {
		m.M[i] = append(m.M[i][i:], m.M[i][:i]...)
	}
}

func (m *Matrix) PermuteB() {
	m.Transpose()
	m.PermuteA()
	m.Transpose()
}

func (m *Matrix) RotateCols(k uint64) {
	k %= m.N()
	for i := range m.M {
		m.M[i] = append(m.M[i][k:], m.M[i][:k]...)
	}
}

func (m *Matrix) Transpose() {
	for i := uint64(0); i < m.N(); i++ {
		for j := 1 + i; j < m.N(); j++ {
			m.M[i][j], m.M[j][i] = m.M[j][i], m.M[i][j]
		}
	}
}

func (m *Matrix) RotateRows(k uint64) {
	k %= m.N()
	m.M = append(m.M[k:], m.M[:k]...)
}

func (m *Matrix) Print() {
	for i := range m.M {
		for j := range m.M[i] {
			fmt.Printf("%7.4f ", m.M[i][j])
		}
		fmt.Printf("\n")
	}
	fmt.Printf("\n")
}

func (m *Matrix) Flatten() (vec []complex128) {
	n := int(m.N())
	vec = make([]complex128, n*n)
	for i := range m.M {
		for j := range m.M[i] {
			vec[i*n+j] = m.M[i][j]
		}
	}
	return
}

func TestMatrices(t *testing.T) {

	rand.Seed(time.Now().UnixNano())

	LogN := uint64(13)
	LogSlots := uint64(12)

	LogModuli := LogModuli{
		LogQi: []uint64{60, 35, 35, 35},
		LogPi: []uint64{55, 55},
	}

	Scale := float64(1 << 35)

	params, err := NewParametersFromLogModuli(LogN, &LogModuli)
	if err != nil {
		panic(err)
	}
	params.SetScale(Scale)
	params.SetLogSlots(LogSlots)

	MM := NewMatrixMultiplier(params)
	encoder := NewEncoder(params)
	kgen := NewKeyGenerator(params)
	sk := kgen.GenSecretKey()
	rlk := kgen.GenRelinKey(sk)
	rotKeys := NewRotationKeys()
	encryptor := NewEncryptorFromSk(params, sk)
	decryptor := NewDecryptor(params, sk)
	eval := NewEvaluator(params)

	// Size of the matrices (dxd)
	d := uint64(64)

	t.Run("Transpose/", func(t *testing.T) {

		m, _, ct := GenTestVectors(d, params, encoder, encryptor)

		level := params.MaxLevel()

		diagMatrix := MM.GenTransposeDiagMatrix(d, params.LogSlots())

		scale := float64(params.Qi()[level])

		mTranspose := encoder.EncodeDiagMatrixAtLvl(level, diagMatrix, scale, 16.0, params.LogSlots())

		kgen.GenRotKeysForDiagMatrix(mTranspose, sk, rotKeys)

		for i := range m {
			m[i].Transpose()
		}

		VerifyTestVectors(d, params, encoder, decryptor, m, eval.LinearTransform(ct, mTranspose, rotKeys)[0], t)

	})

	t.Run("RotateCols", func(t *testing.T) {

		m, _, ct := GenTestVectors(d, params, encoder, encryptor)

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

		for k := uint64(1); k < d; k++ {

			t.Run(fmt.Sprintf("k=%d/", k), func(t *testing.T) {

				level := params.MaxLevel()

				diagMatrix, _ := MM.GenSubVectorRotationMatrix(level, float64(params.Qi()[level]), d, k, params.LogSlots(), encoder)

				kgen.GenRotKeysForDiagMatrix(diagMatrix, sk, rotKeys)

				for j := range m {
					m[j].RotateCols(1)
				}

				eval.(*evaluator).multiplyByDiabMatrix(ct, res, diagMatrix, rotKeys, c2QiQDecompA, c2QiPDecompA)

				VerifyTestVectors(d, params, encoder, decryptor, m, res, t)
			})
		}
	})

	t.Run("RotateRows", func(t *testing.T) {
		m, _, ct := GenTestVectors(d, params, encoder, encryptor)

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

		for k := uint64(1); k < d; k++ {

			t.Run(fmt.Sprintf("k=%d/", k), func(t *testing.T) {

				level := params.MaxLevel()

				diagMatrix, _ := MM.GenSubVectorRotationMatrix(level, float64(params.Qi()[level]), d*d, k*d, params.LogSlots(), encoder)

				kgen.GenRotKeysForDiagMatrix(diagMatrix, sk, rotKeys)

				for j := range m {
					m[j].RotateRows(1)
				}

				eval.(*evaluator).multiplyByDiabMatrix(ct, res, diagMatrix, rotKeys, c2QiQDecompA, c2QiPDecompA)

				VerifyTestVectors(d, params, encoder, decryptor, m, res, t)
			})
		}
	})

	t.Run("PermuteA/", func(t *testing.T) {
		m, _, ct := GenTestVectors(d, params, encoder, encryptor)

		level := params.MaxLevel()

		diagMatrix := MM.GenPermuteAMatrix(d, params.LogSlots())

		scale := float64(params.Qi()[level])

		mPermuteA := encoder.EncodeDiagMatrixAtLvl(level, diagMatrix, scale, 16.0, params.LogSlots())

		kgen.GenRotKeysForDiagMatrix(mPermuteA, sk, rotKeys)

		for j := range m {
			m[j].PermuteA()
		}

		ct = eval.LinearTransform(ct, mPermuteA, rotKeys)[0]

		//PrintDebug(ctAB, d, params, encoder, decryptor)

		VerifyTestVectors(d, params, encoder, decryptor, m, ct, t)

	})

	t.Run("PermuteB/", func(t *testing.T) {
		m, _, ct := GenTestVectors(d, params, encoder, encryptor)

		level := params.MaxLevel()

		diagMatrix := MM.GenPermuteBMatrix(d, params.LogSlots())

		scale := float64(params.Qi()[level])

		mPermuteB := encoder.EncodeDiagMatrixAtLvl(level, diagMatrix, scale, 16.0, params.LogSlots())

		kgen.GenRotKeysForDiagMatrix(mPermuteB, sk, rotKeys)

		for j := range m {
			m[j].PermuteB()
		}

		//PrintDebug(ct, d, params, encoder, decryptor)

		ct = eval.LinearTransform(ct, mPermuteB, rotKeys)[0]

		//PrintDebug(ct, d, params, encoder, decryptor)

		VerifyTestVectors(d, params, encoder, decryptor, m, ct, t)

	})

	t.Run("Multiply/", func(t *testing.T) {
		mA, _, ctA := GenTestVectors(d, params, encoder, encryptor)
		mB, _, ctB := GenTestVectors(d, params, encoder, encryptor)

		mmpt := MM.GenPlaintextMatrices(params, params.MaxLevel(), d, encoder)
		MM.GenRotationKeys(mmpt, kgen, sk, rotKeys)

		for j := range mA {
			mA[j] = MulMat(mA[j], mB[j])
		}

		start := time.Now()
		ctAB := eval.MulMatrixAB(ctA, ctB, mmpt, rlk, rotKeys)
		fmt.Println("Done :", time.Since(start))

		//PrintDebug(ctAB, d, params, encoder, decryptor)

		VerifyTestVectors(d, params, encoder, decryptor, mA, ctAB, t)

	})
}

func MatricesToVector(m []*Matrix, d uint64, params *Parameters) (values []complex128) {
	values = make([]complex128, params.Slots())

	for i := uint64(0); i < params.Slots()/(d*d); i++ {

		vec := m[i].Flatten()

		for j := uint64(0); j < d*d; j++ {
			values[i*d*d+j] = vec[j]
		}
	}
	return
}

func GenTestVectors(d uint64, params *Parameters, encoder Encoder, encryptor Encryptor) (m []*Matrix, pt *Plaintext, ct *Ciphertext) {

	m = GenMatrices(d, params.Slots()/(d*d))

	values := MatricesToVector(m, d, params)

	pt = encoder.EncodeNew(values, params.LogSlots())

	ct = encryptor.EncryptNew(pt)

	return m, pt, ct
}

func VerifyTestVectors(d uint64, params *Parameters, encoder Encoder, decryptor Decryptor, m []*Matrix, element interface{}, t *testing.T) {

	precStats := GetPrecisionStats(params, encoder, decryptor, MatricesToVector(m, d, params), element, 0)

	if *printPrecisionStats {
		t.Log(precStats.String())
	}

	require.GreaterOrEqual(t, real(precStats.MeanPrecision), minPrec)
	require.GreaterOrEqual(t, imag(precStats.MeanPrecision), minPrec)
}

func PrintDebug(ct *Ciphertext, d uint64, params *Parameters, encoder Encoder, decryptor Decryptor) {

	valuesHave := encoder.Decode(decryptor.DecryptNew(ct), params.LogSlots())

	maxPrint := params.Slots() / (d * d)

	if maxPrint > 4 {
		maxPrint = 4
	}
	for k := uint64(0); k < maxPrint; k++ {

		for i := uint64(0); i < d; i++ {
			for j := uint64(0); j < d; j++ {
				fmt.Printf("%4.0f ", real(valuesHave[k*d*d+i*d+j]))
			}
			fmt.Printf("\n")
		}
		fmt.Printf("\n")
	}
}
