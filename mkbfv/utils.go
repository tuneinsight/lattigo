package mkbfv

import (
	"sort"

	"github.com/ldsec/lattigo/v2/ring"
)

// Dot computes the dot product of two decomposed polynomials in R_Q^d
func Dot(decpoly1 *MKDecomposedPoly, decpoly2 *MKDecomposedPoly, accPoly *ring.Poly, ringQP *ring.Ring, beta uint64) {
	for l := uint64(0); l < beta; l++ {
		ringQP.MulCoeffsMontgomeryAndAdd(decpoly1.poly[l], decpoly2.poly[l], accPoly)
	}
}

// MergeSlices merges two slices of uint64 and places the result in s3
// the resulting slice is sorted in ascending order
func MergeSlices(s1, s2 []uint64) []uint64 {

	s3 := make([]uint64, len(s1))

	copy(s3, s1)

	for _, el := range s2 {

		if !Contains(s3, el) {
			s3 = append(s3, el)
		}
	}

	sort.Slice(s3, func(i, j int) bool { return s3[i] < s3[j] })

	return s3
}

// Contains return true if the element is in the slice. False otherwise
func Contains(s []uint64, e uint64) bool {

	for _, el := range s {
		if el == e {
			return true
		}
	}
	return false
}
