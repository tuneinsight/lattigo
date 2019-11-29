package ring

import (
	"github.com/ldsec/lattigo/utils"
	"math/bits"
)

// GenGaloisParams generates the generators for the galois endomorphisms.
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
func PermuteNTTIndex(gen, power, N uint64) (index []uint64) {

	genPow := ModExp(gen, power, 2*N)

	var mask, logN, tmp1, tmp2 uint64

	logN = uint64(bits.Len64(N) - 1)

	mask = (N << 1) - 1

	index = make([]uint64, N)

	for i := uint64(0); i < N; i++ {
		tmp1 = 2*utils.BitReverse64(i, logN) + 1

		tmp2 = ((genPow * tmp1 & mask) - 1) >> 1

		index[i] = utils.BitReverse64(tmp2, logN)
	}

	return
}

// PermuteNTT applies the galois transform on a polynomial in the NTT domain.
// It maps the coefficients x^i to x^(gen*i)
// Careful, not inplace!
func PermuteNTT(polIn *Poly, gen uint64, polOut *Poly) {

	var N, tmp, mask, logN, tmp1, tmp2 uint64

	N = uint64(len(polIn.Coeffs[0]))

	logN = uint64(bits.Len64(N) - 1)

	mask = (N << 1) - 1

	index := make([]uint64, N)

	for i := uint64(0); i < N; i++ {
		tmp1 = 2*utils.BitReverse64(i, logN) + 1

		tmp2 = ((gen * tmp1 & mask) - 1) >> 1

		index[i] = utils.BitReverse64(tmp2, logN)
	}

	for j := uint64(0); j < N; j++ {

		tmp = index[j]

		for i := 0; i < len(polIn.Coeffs); i++ {

			polOut.Coeffs[i][j] = polIn.Coeffs[i][tmp]
		}
	}
}

// PermuteNTTWithIndex applies the galois transform on a polynomial in the NTT domain.
// It maps the coefficients x^i to x^(gen*i) using the PermuteNTTIndex table.
// Careful, not inplace!
func PermuteNTTWithIndex(polIn *Poly, index []uint64, polOut *Poly) {

	var tmp uint64

	for j := uint64(0); j < uint64(len(polIn.Coeffs[0])); j++ {

		tmp = index[j]

		for i := 0; i < len(polIn.Coeffs); i++ {
			polOut.Coeffs[i][j] = polIn.Coeffs[i][tmp]
		}
	}
}

// Permute applies the galois transform on a polynonial outside of the NTT domain.
// It maps the coefficients x^i to x^(gen*i)
// Careful, not inplace!
func (context *Context) Permute(polIn *Poly, gen uint64, polOut *Poly) {

	var mask, index, indexRaw, logN, tmp uint64

	mask = context.N - 1

	logN = uint64(bits.Len64(mask))

	for i := uint64(0); i < context.N; i++ {

		indexRaw = i * gen

		index = indexRaw & mask

		tmp = (indexRaw >> logN) & 1

		for j, qi := range context.Modulus {

			polOut.Coeffs[j][index] = polIn.Coeffs[j][i]*(tmp^1) | (qi-polIn.Coeffs[j][i])*tmp
		}
	}
}
