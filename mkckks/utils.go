package mkckks

import (
	"sort"
	"unsafe"

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

	for i := uint64(0); i < uint64(len(p1.poly)); i++ {
		r.MulCoeffsAndAdd(p1.poly[i], p2.poly[i], res)
	}

	return res
}

// DotLvl computes the dot product of two decomposed polynomials in ringQ^d up to q_level and store the result in res
func DotLvl(level uint64, p1 *MKDecomposedPoly, p2 *MKDecomposedPoly, r *ring.Ring) *ring.Poly {
	if len(p1.poly) != len(p2.poly) {
		panic("Cannot compute dot product on vectors of different size !")
	}

	res := r.NewPoly()

	for i := uint64(0); i < uint64(len(p1.poly)); i++ {

		MulCoeffsAndAddLvl(level, p1.poly[i], p2.poly[i], res, r)
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

// returns a polynomial with all coefficient set to 1
func getOne(r *ring.Ring) *ring.Poly {

	res := r.NewPoly()

	coeffs := res.Coeffs

	for _, c := range coeffs {
		for index := range c {
			c[index] = 1
		}
	}

	res.SetCoefficients(coeffs)

	return res
}

// MulCoeffsAndAddLvl multiplies p1 by p2 coefficient-wise with
// a Barret modular reduction up to q_leveland adds the result to p3.
func MulCoeffsAndAddLvl(level uint64, p1, p2, p3 *ring.Poly, r *ring.Ring) {

	for i, qi := range r.Modulus[:level+1] {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		bredParams := r.BredParams[i]
		for j := uint64(0); j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))

			z[0] = ring.CRed(z[0]+ring.BRed(x[0], y[0], qi, bredParams), qi)
			z[1] = ring.CRed(z[1]+ring.BRed(x[1], y[1], qi, bredParams), qi)
			z[2] = ring.CRed(z[2]+ring.BRed(x[2], y[2], qi, bredParams), qi)
			z[3] = ring.CRed(z[3]+ring.BRed(x[3], y[3], qi, bredParams), qi)
			z[4] = ring.CRed(z[4]+ring.BRed(x[4], y[4], qi, bredParams), qi)
			z[5] = ring.CRed(z[5]+ring.BRed(x[5], y[5], qi, bredParams), qi)
			z[6] = ring.CRed(z[6]+ring.BRed(x[6], y[6], qi, bredParams), qi)
			z[7] = ring.CRed(z[7]+ring.BRed(x[7], y[7], qi, bredParams), qi)
		}
	}
}

// MulCoeffsAndSubLvl multiplies p1 by p2 coefficient-wise with
// a Barett modular reduction up to q_level and subtracts the result from p3.
func MulCoeffsAndSubLvl(level uint64, p1, p2, p3 *ring.Poly, r *ring.Ring) {
	for i, qi := range r.Modulus[:level+1] {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		bredParams := r.BredParams[i]
		for j := uint64(0); j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))

			z[0] = ring.CRed(z[0]+(qi-ring.BRed(x[0], y[0], qi, bredParams)), qi)
			z[1] = ring.CRed(z[1]+(qi-ring.BRed(x[1], y[1], qi, bredParams)), qi)
			z[2] = ring.CRed(z[2]+(qi-ring.BRed(x[2], y[2], qi, bredParams)), qi)
			z[3] = ring.CRed(z[3]+(qi-ring.BRed(x[3], y[3], qi, bredParams)), qi)
			z[4] = ring.CRed(z[4]+(qi-ring.BRed(x[4], y[4], qi, bredParams)), qi)
			z[5] = ring.CRed(z[5]+(qi-ring.BRed(x[5], y[5], qi, bredParams)), qi)
			z[6] = ring.CRed(z[6]+(qi-ring.BRed(x[6], y[6], qi, bredParams)), qi)
			z[7] = ring.CRed(z[7]+(qi-ring.BRed(x[7], y[7], qi, bredParams)), qi)
		}
	}
}

// MulCoeffsLvl multiplies p1 by p2 coefficient-wise at a certain level, performs a
// Barrett modular reduction and writes the result on p3.
func MulCoeffsLvl(level uint64, p1, p2, p3 *ring.Poly, r *ring.Ring) {
	for i, qi := range r.Modulus[:level+1] {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		bredParams := r.BredParams[i]
		for j := uint64(0); j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))

			z[0] = ring.BRed(x[0], y[0], qi, bredParams)
			z[1] = ring.BRed(x[1], y[1], qi, bredParams)
			z[2] = ring.BRed(x[2], y[2], qi, bredParams)
			z[3] = ring.BRed(x[3], y[3], qi, bredParams)
			z[4] = ring.BRed(x[4], y[4], qi, bredParams)
			z[5] = ring.BRed(x[5], y[5], qi, bredParams)
			z[6] = ring.BRed(x[6], y[6], qi, bredParams)
			z[7] = ring.BRed(x[7], y[7], qi, bredParams)
		}
	}
}
