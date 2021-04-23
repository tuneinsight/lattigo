package mkckks

import (
	"sort"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// Dot computes the dot product of two decomposed polynomials in ring^d and store the result in res
func Dot(p1 *MKDecomposedPoly, p2 *MKDecomposedPoly, r *ring.Ring) *ring.Poly {
	if len(p1.poly) != len(p2.poly) {
		panic("Cannot compute dot product on vectors of different size !")
	}

	res := r.NewPoly()

	tmp0 := r.NewPoly()
	tmp1 := r.NewPoly()

	for i := uint64(0); i < uint64(len(p1.poly)); i++ {
		r.MForm(p1.poly[i], tmp0)
		r.MForm(p2.poly[i], tmp1)
		r.MulCoeffsMontgomeryAndAdd(tmp0, tmp1, res)
	}

	return res
}

// DotLvl computes the dot product of two decomposed polynomials in ringQ^d up to q_level and store the result in res
func DotLvl(level uint64, p1 *MKDecomposedPoly, p2 *MKDecomposedPoly, r *ring.Ring) *ring.Poly {
	if len(p1.poly) != len(p2.poly) {
		panic("Cannot compute dot product on vectors of different size !")
	}

	res := r.NewPoly()

	tmp0 := r.NewPoly()
	tmp1 := r.NewPoly()

	for i := uint64(0); i < uint64(len(p1.poly)); i++ {
		r.MFormLvl(level, p1.poly[i], tmp0)
		r.MFormLvl(level, p2.poly[i], tmp1)
		r.MulCoeffsMontgomeryAndAddLvl(level, tmp0, tmp1, res)
	}

	return res
}

// MergeSlices merges two slices of uint64 and places the result in s3
// the resulting slice is sorted in ascending order
func MergeSlices(s1, s2 []uint64) []uint64 {

	s3 := make([]uint64, len(s1))

	copy(s3, s1)

	for _, el := range s2 {

		if Contains(s3, el) < 0 {
			s3 = append(s3, el)
		}
	}

	sort.Slice(s3, func(i, j int) bool { return s3[i] < s3[j] })

	return s3
}

// Contains return the element's index if the element is in the slice. -1 otherwise
func Contains(s []uint64, e uint64) int {

	for i, el := range s {
		if el == e {
			return i
		}
	}
	return -1
}

// GetRandomPoly samples a polynomial with a uniform distribution in the given ring
func GetRandomPoly(params *ckks.Parameters, r *ring.Ring) *ring.Poly {

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	return GetUniformSampler(params, r, prng).ReadNew()
}
