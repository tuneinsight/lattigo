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
	GenColumnRotationMatrix(d, k, logSlots uint64) (diagMatrix map[uint64][]complex128)
	GenRowRotationMatrix(d, k, logSlots uint64) (diagMatrix map[uint64][]complex128)
}

type matrixMultiplier struct {
}

type MMPt struct {
	d         uint64
	mPermuteA *PtDiagMatrix
	mPermuteB *PtDiagMatrix
	mRotRows  []*PtDiagMatrix
	mRotCols  []*PtDiagMatrix
}

func NewMatrixMultiplier(params *Parameters) MatrixMultiplier {

	return &matrixMultiplier{}
}

func (mm *matrixMultiplier) GenPlaintextMatrices(params *Parameters, level uint64, d uint64, encoder Encoder) (mmpt *MMPt) {

	mmpt = new(MMPt)

	mmpt.d = d

	var scale float64
	scale = float64(params.Qi()[level]) * math.Sqrt(float64(params.Qi()[level-2])/params.Scale())

	mmpt.mPermuteA = encoder.EncodeDiagMatrixAtLvl(level, mm.GenPermuteAMatrix(d, params.LogSlots()), scale, 16.0, params.LogSlots())
	mmpt.mPermuteB = encoder.EncodeDiagMatrixAtLvl(level, mm.GenPermuteBMatrix(d, params.LogSlots()), scale, 16.0, params.LogSlots())

	mmpt.mRotCols = make([]*PtDiagMatrix, d-1)
	mmpt.mRotRows = make([]*PtDiagMatrix, d-1)

	scale = float64(params.Qi()[level-1])

	for i := uint64(0); i < d-1; i++ {
		mmpt.mRotCols[i] = encoder.EncodeDiagMatrixAtLvl(level-1, mm.GenColumnRotationMatrix(d, i+1, params.LogSlots()), scale, 16.0, params.LogSlots())
		mmpt.mRotRows[i] = encoder.EncodeDiagMatrixAtLvl(level-1, mm.GenRowRotationMatrix(d, i+1, params.LogSlots()), scale, 16.0, params.LogSlots())
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

	ciphertextA := eval.MultiplyByDiagMatrix(A, mmpt.mPermuteA, rotKeys)
	eval.Rescale(ciphertextA, eval.params.Scale(), ciphertextA)

	ciphertextB := eval.MultiplyByDiagMatrix(B, mmpt.mPermuteB, rotKeys)
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
	for i := uint64(0); i < mmpt.d-1; i++ {

		tmpA := NewCiphertext(eval.params, 1, ciphertextA.Level(), ciphertextA.Scale())
		tmpB := NewCiphertext(eval.params, 1, ciphertextB.Level(), ciphertextB.Scale())

		eval.multiplyByDiabMatrixNaive(ciphertextA, tmpA, mmpt.mRotCols[i], rotKeys, eval.c2QiQDecomp, eval.c2QiPDecomp)
		eval.multiplyByDiabMatrixNaive(ciphertextB, tmpB, mmpt.mRotRows[i], rotKeys, c2QiQDecompB, c2QiPDecompB)

		eval.Rescale(tmpA, eval.params.Scale(), tmpA)
		eval.Rescale(tmpB, eval.params.Scale(), tmpB)

		eval.MulRelin(tmpA, tmpB, nil, tmpC)

		eval.Add(ciphertextAB, tmpC, ciphertextAB)
	}

	eval.Relinearize(ciphertextAB, rlk, ciphertextAB)
	eval.Rescale(ciphertextAB, eval.params.Scale(), ciphertextAB)

	return
}

func (mm *matrixMultiplier) GenPermuteAMatrix(d, logSlots uint64) (diagMatrix map[uint64][]complex128) {

	slots := uint64(1 << logSlots)

	diagMatrix = make(map[uint64][]complex128)

	d2 := int(d * d)

	for i := -int(d) + 1; i < int(d); i++ {

		m := make([]complex128, slots)

		for k := 0; k < d2; k++ {

			if i < 0 {
				for j := i; j < int(d); j++ {
					x := (d2 + k - (int(d)+i)*int(d)) % d2
					if x < int(d) && x >= -i {
						m[k] = 1
					}
				}
			} else {

				for j := i; j < int(d); j++ {
					if (d2+k-int(d)*i)%d2 < int(d)-i {
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

func (mm *matrixMultiplier) GenPermuteBMatrix(d, logSlots uint64) (diagMatrix map[uint64][]complex128) {

	slots := uint64(1 << logSlots)

	diagMatrix = make(map[uint64][]complex128)

	d2 := int(d * d)

	if d*d < 1<<logSlots {

		for i := -int((d - 1) * d); i < d2; i = i + int(d) {

			m := make([]complex128, 1<<logSlots)

			if i >= 0 {
				for j := 0; j < d2-i; j = j + int(d) {
					m[i/int(d)+j] = 1
				}
			} else {
				for j := 0; j < d2+i; j = j + int(d) {
					m[-i+int(d)+(i/int(d))+j] = 1
				}
			}

			populateVector(m, d2, logSlots)

			diagMatrix[uint64((i+int(slots)))%slots] = m
		}
	} else {
		for i := 0; i < int(d); i++ {

			m := make([]complex128, 1<<logSlots)

			for j := 0; j < d2; j = j + int(d) {
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

			diagMatrix[uint64(i)*d] = m
		}
	}

	return
}

func (mm *matrixMultiplier) GenColumnRotationMatrix(d, k, logSlots uint64) (diagMatrix map[uint64][]complex128) {

	k %= d

	k = d - k

	diagMatrix = make(map[uint64][]complex128)

	d2 := int(d * d)

	m0 := make([]complex128, 1<<logSlots)
	for i := 0; i < d2; i++ {
		if i%int(d) < int(d-k) {
			m0[i] = 1
		}
	}

	populateVector(m0, d2, logSlots)

	m1 := make([]complex128, 1<<logSlots)
	for i := 0; i < d2; i++ {
		if i%int(d) >= int(d-k) {
			m1[i] = 1
		}
	}

	populateVector(m1, d2, logSlots)

	/*
		fmt.Printf("%4d", (1<<logSlots)-k)
		for i := range m0[:d2+d2]{
			fmt.Printf("%2.f ", real(m0[i]))
		}
		fmt.Println()

		fmt.Printf("%4d", d-k)
		for i := range m1[:d2+d2]{
			fmt.Printf("%2.f ", real(m1[i]))
		}
		fmt.Println()
	*/

	diagMatrix[(1<<logSlots)-k] = m0
	diagMatrix[d-k] = m1

	return
}

func (mm *matrixMultiplier) GenRowRotationMatrix(d, k, logSlots uint64) (diagMatrix map[uint64][]complex128) {

	k %= d

	diagMatrix = make(map[uint64][]complex128)

	d2 := int(d * d)

	if d*d < 1<<logSlots {

		k = d - k

		m0 := make([]complex128, 1<<logSlots)
		for i := 0; i < d2-int(k*d); i++ {
			m0[i] = 1
		}

		m1 := make([]complex128, 1<<logSlots)
		for i := (uint64(d2) - d*k); i < uint64(d2); i++ {
			m1[i] = 1
		}

		populateVector(m0, d2, logSlots)
		populateVector(m1, d2, logSlots)

		rotatem0 := (1 << logSlots) - d*k
		rotatem1 := uint64(d2) - d*k

		/*
			fmt.Printf("%4d", rotatem0)
			for i := range m0[:d2+d2]{
				fmt.Printf("%2.f ", real(m0[i]))
			}
			fmt.Println()

			fmt.Printf("%4d", rotatem1)
			for i := range m1[:d2+d2]{
				fmt.Printf("%2.f ", real(m1[i]))
			}
			fmt.Println()
		*/

		diagMatrix[rotatem0] = m0
		diagMatrix[rotatem1] = m1

	} else {

		m0 := make([]complex128, 1<<logSlots)
		for i := 0; i < d2; i++ {
			m0[i] = 1
		}

		populateVector(m0, d2, logSlots)

		diagMatrix[d*k] = m0
	}

	return
}

func (mm *matrixMultiplier) GenTransposeDiagMatrix(d, logSlots uint64) (diagMatrix map[uint64][]complex128) {

	slots := uint64(1 << logSlots)

	diagMatrix = make(map[uint64][]complex128)

	d2 := int(d * d)

	for i := -int(d) + 1; i < int(d); i++ {

		m := make([]complex128, slots)

		if i >= 0 {
			for j := 0; j < d2-i*int(d); j = j + int(d) + 1 {
				m[i+j] = 1
			}
		} else {
			for j := -i * int(d); j < d2; j = j + int(d) + 1 {
				m[j] = 1
			}
		}

		populateVector(m, d2, logSlots)

		diagMatrix[uint64(i*int(d-1)+int(slots))%slots] = m
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
