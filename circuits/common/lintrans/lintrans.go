// Package lintrans bundles generic parts of the homomorphic linear transformation circuit.
package lintrans

import (
	"fmt"
	"sort"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
	"github.com/tuneinsight/lattigo/v6/schemes"
	"github.com/tuneinsight/lattigo/v6/utils"
)

// Parameters is a struct storing the parameterization of a
// linear transformation.
//
// A homomorphic linear transformations on a ciphertext acts as evaluating:
//
// Ciphertext([1 x n] vector) <- Ciphertext([1 x n] vector) x Plaintext([n x n] matrix)
//
// where n is the number of plaintext slots.
//
// The diagonal representation of a linear transformations is defined by first expressing
// the linear transformation through its nxn matrix and then traversing the matrix diagonally.
//
// For example, the following nxn for n=4 matrix:
//
// 0 1 2 3 (diagonal index)
// | 1 2 3 0 |
// | 0 1 2 3 |
// | 3 0 1 2 |
// | 2 3 0 1 |
//
// its diagonal traversal representation is comprised of 3 non-zero diagonals at indexes [0, 1, 2]:
// 0: [1, 1, 1, 1]
// 1: [2, 2, 2, 2]
// 2: [3, 3, 3, 3]
// 3: [0, 0, 0, 0] -> this diagonal is omitted as it is composed only of zero values.
//
// Note that negative indexes can be used and will be interpreted modulo the matrix dimension.
//
// The diagonal representation is well suited for two reasons:
//  1. It is the effective format used during the homomorphic evaluation.
//  2. It enables on average a more compact and efficient representation of linear transformations
//     than their matrix representation by being able to only store the non-zero diagonals.
//
// Finally, some metrics about the time and storage complexity of homomorphic linear transformations:
//   - Storage: #diagonals polynomials mod Q * P
//   - Evaluation: #diagonals multiplications and 2sqrt(#diagonals) ciphertexts rotations.
type Parameters struct {
	// DiagonalsIndexList is the list of the non-zero diagonals of the square matrix.
	// A non zero diagonals is a diagonal with a least one non-zero element.
	DiagonalsIndexList []int

	// LevelQ is the level at which to encode the linear transformation.
	LevelQ int

	// LevelP is the level of the auxliary prime used during the automorphisms
	// User must ensure that this value is the same as the one used to generate
	// the evaluation keys used to perform the automorphisms.
	LevelP int

	// Scale is the plaintext scale at which to encode the linear transformation.
	Scale rlwe.Scale

	// LogDimensions is the log2 dimensions of the matrix that can be SIMD packed
	// in a single plaintext polynomial.
	// This method is equivalent to params.PlaintextDimensions().
	// Note that the linear transformation is evaluated independently on each rows of
	// the SIMD packed matrix.
	LogDimensions ring.Dimensions

	// LogBabyStepGianStepRatio is the log2 of the ratio n1/n2 for n = n1 * n2 and
	// n is the dimension of the linear transformation. The number of Galois keys required
	// is minimized when this value is 0 but the overall complexity of the homomorphic evaluation
	// can be reduced by increasing the ratio (at the expanse of increasing the number of keys required).
	// If the value returned is negative, then the baby-step giant-step algorithm is not used
	// and the evaluation complexity (as well as the number of keys) becomes O(n) instead of O(sqrt(n)).
	LogBabyStepGiantStepRatio int
}

type Diagonals[T any] map[int][]T

// DiagonalsIndexList returns the list of the non-zero diagonals of the square matrix.
// A non zero diagonals is a diagonal with a least one non-zero element.
func (m Diagonals[T]) DiagonalsIndexList() (indexes []int) {
	indexes = make([]int, 0, len(m))
	for k := range m {
		indexes = append(indexes, k)
	}
	return indexes
}

// At returns the i-th non-zero diagonal.
// Method accepts negative values with the equivalency -i = n - i.
func (m Diagonals[T]) At(i, slots int) ([]T, error) {

	v, ok := m[i]

	if !ok {

		var j int
		if i > 0 {
			j = i - slots
		} else if j < 0 {
			j = i + slots
		} else {
			return nil, fmt.Errorf("cannot At[0]: diagonal does not exist")
		}

		v, ok := m[j]

		if !ok {
			return nil, fmt.Errorf("cannot At[%d or %d]: diagonal does not exist", i, j)
		}

		return v, nil
	}

	return v, nil
}

// LinearTransformation is a type for linear transformations on ciphertexts.
// It stores a plaintext matrix in diagonal form and can be evaluated on a
// ciphertext using a [Evaluator].
type LinearTransformation struct {
	*rlwe.MetaData
	LogBabyStepGiantStepRatio int
	N1                        int
	LevelQ                    int
	LevelP                    int
	Vec                       map[int]ringqp.Poly
}

// GaloisElements returns the list of Galois elements needed for the evaluation of the linear transformation.
func (lt LinearTransformation) GaloisElements(params rlwe.ParameterProvider) (galEls []uint64) {
	return GaloisElements(params, utils.GetKeys(lt.Vec), 1<<lt.LogDimensions.Cols, lt.LogBabyStepGiantStepRatio)
}

// BSGSIndex returns the BSGSIndex of the target linear transformation.
func (lt LinearTransformation) BSGSIndex() (index map[int][]int, n1, n2 []int) {
	return BSGSIndex(utils.GetKeys(lt.Vec), 1<<lt.LogDimensions.Cols, lt.N1)
}

// NewLinearTransformation allocates a new LinearTransformation with zero values according to the parameters specified
// by the [Parameters].
func NewLinearTransformation(params rlwe.ParameterProvider, ltparams Parameters) LinearTransformation {

	p := params.GetRLWEParameters()

	vec := make(map[int]ringqp.Poly)
	cols := 1 << ltparams.LogDimensions.Cols
	logBabyStepGiantStepRatio := ltparams.LogBabyStepGiantStepRatio
	levelQ := ltparams.LevelQ
	levelP := ltparams.LevelP
	ringQP := p.RingQP().AtLevel(levelQ, levelP)

	diagslislt := ltparams.DiagonalsIndexList

	var N1 int
	if logBabyStepGiantStepRatio < 0 {
		N1 = 0
		for _, i := range diagslislt {
			idx := i
			if idx < 0 {
				idx += cols
			}
			vec[idx] = ringQP.NewPoly()
		}
	} else {
		N1 = FindBestBSGSRatio(diagslislt, cols, logBabyStepGiantStepRatio)
		index, _, _ := BSGSIndex(diagslislt, cols, N1)
		for j := range index {
			for _, i := range index[j] {
				vec[j+i] = ringQP.NewPoly()
			}
		}
	}

	metadata := &rlwe.MetaData{
		PlaintextMetaData: rlwe.PlaintextMetaData{
			LogDimensions: ltparams.LogDimensions,
			Scale:         ltparams.Scale,
			IsBatched:     true,
		},
		CiphertextMetaData: rlwe.CiphertextMetaData{
			IsNTT:        true,
			IsMontgomery: true,
		},
	}

	return LinearTransformation{
		MetaData:                  metadata,
		LogBabyStepGiantStepRatio: logBabyStepGiantStepRatio,
		N1:                        N1,
		LevelQ:                    levelQ,
		LevelP:                    levelP,
		Vec:                       vec,
	}
}

// Encode encodes on a pre-allocated LinearTransformation a set of non-zero diagonals of a matrix representing
// a linear transformation.
func Encode[T any](encoder schemes.Encoder, diagonals Diagonals[T], allocated LinearTransformation) (err error) {

	rows := 1 << allocated.LogDimensions.Rows
	cols := 1 << allocated.LogDimensions.Cols
	N1 := allocated.N1

	diags := diagonals.DiagonalsIndexList()

	buf := make([]T, rows*cols)

	metaData := allocated.MetaData

	metaData.Scale = allocated.Scale

	var v []T

	if N1 == 0 {
		for _, i := range diags {

			idx := i
			if idx < 0 {
				idx += cols
			}

			if vec, ok := allocated.Vec[idx]; !ok {
				return fmt.Errorf("cannot Encode: error encoding on LinearTransformation: plaintext diagonal [%d] does not exist", idx)
			} else {

				if v, err = diagonals.At(i, cols); err != nil {
					return fmt.Errorf("cannot Encode: %w", err)
				}

				if err = rotateAndEncodeDiagonal(v, encoder, 0, metaData, buf, vec); err != nil {
					return
				}
			}
		}
	} else {

		index, _, _ := allocated.BSGSIndex()

		for j := range index {

			rot := -j & (cols - 1)

			for _, i := range index[j] {

				if vec, ok := allocated.Vec[i+j]; !ok {
					return fmt.Errorf("cannot Encode: error encoding on LinearTransformation BSGS: input does not match the same non-zero diagonals")
				} else {

					if v, err = diagonals.At(i+j, cols); err != nil {
						return fmt.Errorf("cannot Encode: %w", err)
					}

					if err = rotateAndEncodeDiagonal(v, encoder, rot, metaData, buf, vec); err != nil {
						return
					}
				}
			}
		}
	}

	return
}

func rotateAndEncodeDiagonal[T any](v []T, encoder schemes.Encoder, rot int, metaData *rlwe.MetaData, buf []T, poly ringqp.Poly) (err error) {

	rows := 1 << metaData.LogDimensions.Rows
	cols := 1 << metaData.LogDimensions.Cols

	rot &= (cols - 1)

	var values []T
	if rot != 0 {

		values = buf

		for i := 0; i < rows; i++ {
			utils.RotateSliceAllocFree(v[i*cols:(i+1)*cols], rot, values[i*cols:(i+1)*cols])
		}

	} else {
		values = v
	}

	return encoder.Embed(values, metaData, poly)
}

// GaloisElements returns the list of Galois elements needed for the evaluation of a linear transformation
// given the index of its non-zero diagonals, the number of slots in the plaintext and the LogBabyStepGiantStepRatio,
// see [Parameters].
func GaloisElements(params rlwe.ParameterProvider, diags []int, slots, logBabyStepGianStepRatio int) (galEls []uint64) {

	p := params.GetRLWEParameters()

	if logBabyStepGianStepRatio < 0 {

		_, _, rotN2 := BSGSIndex(diags, slots, slots)

		galEls = make([]uint64, len(rotN2))
		for i := range rotN2 {
			galEls[i] = p.GaloisElement(rotN2[i])
		}

		return
	}

	N1 := FindBestBSGSRatio(diags, slots, logBabyStepGianStepRatio)

	_, rotN1, rotN2 := BSGSIndex(diags, slots, N1)

	return p.GaloisElements(utils.GetDistincts(append(rotN1, rotN2...)))
}

// FindBestBSGSRatio finds the best N1*N2 = N for the baby-step giant-step algorithm for matrix multiplication.
func FindBestBSGSRatio(nonZeroDiags []int, maxN int, logMaxRatio int) (minN int) {

	maxRatio := float64(int(1 << logMaxRatio))

	for N1 := 1; N1 < maxN; N1 <<= 1 {

		_, rotN1, rotN2 := BSGSIndex(nonZeroDiags, maxN, N1)

		nbN1, nbN2 := len(rotN1)-1, len(rotN2)-1

		if float64(nbN2)/float64(nbN1) == maxRatio {
			return N1
		}

		if float64(nbN2)/float64(nbN1) > maxRatio {
			return N1 / 2
		}
	}

	return 1
}

// BSGSIndex returns the index map and needed rotation for the BSGS matrix-vector multiplication algorithm.
func BSGSIndex(nonZeroDiags []int, slots, N1 int) (index map[int][]int, rotN1, rotN2 []int) {
	index = make(map[int][]int)
	rotN1Map := make(map[int]bool)
	rotN2Map := make(map[int]bool)

	for _, rot := range nonZeroDiags {
		rot &= (slots - 1)
		idxN1 := ((rot / N1) * N1) & (slots - 1)
		idxN2 := rot & (N1 - 1)
		if index[idxN1] == nil {
			index[idxN1] = []int{idxN2}
		} else {
			index[idxN1] = append(index[idxN1], idxN2)
		}
		rotN1Map[idxN1] = true
		rotN2Map[idxN2] = true
	}

	for k := range index {
		sort.Ints(index[k])
	}

	return index, utils.GetSortedKeys(rotN1Map), utils.GetSortedKeys(rotN2Map)
}
