package ring

import (
	"math/bits"
	"unsafe"

	"github.com/tuneinsight/lattigo/v3/utils"
)

// GenGaloisParams generates the generators for the Galois endomorphisms.
func GenGaloisParams(n, gen uint64) (galElRotCol []uint64) {

	var m, mask uint64

	m = n << 1

	mask = m - 1

	galElRotCol = make([]uint64, n>>1)

	galElRotCol[0] = 1

	for i := uint64(1); i < n>>1; i++ {
		galElRotCol[i] = (galElRotCol[i-1] * gen) & mask
	}

	return
}

// PermuteNTTIndex computes the index table for PermuteNTT.
func (r *Ring) PermuteNTTIndex(galEl uint64) (index []uint64) {

	var mask, tmp1, tmp2, logNthRoot uint64
	logNthRoot = uint64(bits.Len64(r.NthRoot) - 2)
	mask = r.NthRoot - 1
	index = make([]uint64, r.N)

	for i := uint64(0); i < uint64(r.N); i++ {
		tmp1 = 2*utils.BitReverse64(i, logNthRoot) + 1
		tmp2 = ((galEl * tmp1 & mask) - 1) >> 1
		index[i] = utils.BitReverse64(tmp2, logNthRoot)
	}

	return
}

// PermuteNTT applies the Galois transform on a polynomial in the NTT domain.
// It maps the coefficients x^i to x^(gen*i)
// It must be noted that the result cannot be in-place.
func (r *Ring) PermuteNTT(polIn *Poly, gen uint64, polOut *Poly) {
	r.PermuteNTTLvl(utils.MinInt(polIn.Level(), polOut.Level()), polIn, gen, polOut)
}

// PermuteNTTLvl applies the Galois transform on a polynomial in the NTT domain, up to a given level.
// It maps the coefficients x^i to x^(gen*i)
// It must be noted that the result cannot be in-place.
func (r *Ring) PermuteNTTLvl(level int, polIn *Poly, gen uint64, polOut *Poly) {
	r.PermuteNTTWithIndexLvl(level, polIn, r.PermuteNTTIndex(gen), polOut)
}

// PermuteNTTWithIndexLvl applies the Galois transform on a polynomial in the NTT domain, up to a given level.
// It maps the coefficients x^i to x^(gen*i) using the PermuteNTTIndex table.
// It must be noted that the result cannot be in-place.
func (r *Ring) PermuteNTTWithIndexLvl(level int, polIn *Poly, index []uint64, polOut *Poly) {

	for j := 0; j < r.N; j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&index[j]))

		for i := 0; i < level+1; i++ {

			z := (*[8]uint64)(unsafe.Pointer(&polOut.Coeffs[i][j]))
			y := polIn.Coeffs[i]

			z[0] = y[x[0]]
			z[1] = y[x[1]]
			z[2] = y[x[2]]
			z[3] = y[x[3]]
			z[4] = y[x[4]]
			z[5] = y[x[5]]
			z[6] = y[x[6]]
			z[7] = y[x[7]]
		}
	}
}

// PermuteNTTWithIndexAndAddNoModLvl applies the Galois transform on a polynomial in the NTT domain, up to a given level,
// and adds the result to the output polynomial without modular reduction.
// It maps the coefficients x^i to x^(gen*i) using the PermuteNTTIndex table.
// It must be noted that the result cannot be in-place.
func (r *Ring) PermuteNTTWithIndexAndAddNoModLvl(level int, polIn *Poly, index []uint64, polOut *Poly) {

	for j := 0; j < r.N; j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&index[j]))

		for i := 0; i < level+1; i++ {

			z := (*[8]uint64)(unsafe.Pointer(&polOut.Coeffs[i][j]))
			y := polIn.Coeffs[i]

			z[0] += y[x[0]]
			z[1] += y[x[1]]
			z[2] += y[x[2]]
			z[3] += y[x[3]]
			z[4] += y[x[4]]
			z[5] += y[x[5]]
			z[6] += y[x[6]]
			z[7] += y[x[7]]
		}
	}
}

// Permute applies the Galois transform on a polynomial outside of the NTT domain.
// It maps the coefficients x^i to x^(gen*i).
// It must be noted that the result cannot be in-place.
func (r *Ring) Permute(polIn *Poly, gen uint64, polOut *Poly) {
	r.PermuteLvl(utils.MinInt(polIn.Level(), polOut.Level()), polIn, gen, polOut)
}

// PermuteLvl applies the Galois transform on a polynomial outside of the NTT domain.
// It maps the coefficients x^i to x^(gen*i).
// It must be noted that the result cannot be in-place.
func (r *Ring) PermuteLvl(level int, polIn *Poly, gen uint64, polOut *Poly) {

	var mask, index, indexRaw, logN, tmp uint64

	mask = uint64(r.N - 1)

	logN = uint64(bits.Len64(mask))

	for i := uint64(0); i < uint64(r.N); i++ {

		indexRaw = i * gen

		index = indexRaw & mask

		tmp = (indexRaw >> logN) & 1

		for j, qi := range r.Modulus[:level+1] {

			polOut.Coeffs[j][index] = polIn.Coeffs[j][i]*(tmp^1) | (qi-polIn.Coeffs[j][i])*tmp
		}
	}
}
