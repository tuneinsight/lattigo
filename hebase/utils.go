package hebase

import (
	"sort"

	"github.com/tuneinsight/lattigo/v4/utils"
)

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
