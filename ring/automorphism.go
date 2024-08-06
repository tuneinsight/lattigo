package ring

import (
	"fmt"
	"math/bits"
	"unsafe"

	"github.com/tuneinsight/lattigo/v6/utils"
)

// AutomorphismNTTIndex computes the look-up table for the automorphism X^{i} -> X^{i*k mod NthRoot}.
func AutomorphismNTTIndex(N int, NthRoot, GalEl uint64) (index []uint64, err error) {

	if N&(N-1) != 0 {
		return nil, fmt.Errorf("N must be a power of two")
	}

	if NthRoot&(NthRoot-1) != 0 {
		return nil, fmt.Errorf("NthRoot must be w power of two")
	}

	var mask, tmp1, tmp2 uint64
	logNthRoot := int(bits.Len64(NthRoot-1) - 1)
	mask = NthRoot - 1
	index = make([]uint64, N)

	for i := 0; i < N; i++ {
		tmp1 = 2*utils.BitReverse64(i, logNthRoot) + 1
		tmp2 = ((GalEl * tmp1 & mask) - 1) >> 1
		index[i] = utils.BitReverse64(tmp2, logNthRoot)
	}

	return
}

// AutomorphismNTT applies the automorphism X^{i} -> X^{i*gen} on a polynomial in the NTT domain.
// It must be noted that the result cannot be in-place.
func (r Ring) AutomorphismNTT(polIn Poly, gen uint64, polOut Poly) {
	index, err := AutomorphismNTTIndex(r.N(), r.NthRoot(), gen)
	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}
	r.AutomorphismNTTWithIndex(polIn, index, polOut)
}

// AutomorphismNTTWithIndex applies the automorphism X^{i} -> X^{i*gen} on a polynomial in the NTT domain.
// `index` is the lookup table storing the mapping of the automorphism.
// It must be noted that the result cannot be in-place.
func (r Ring) AutomorphismNTTWithIndex(polIn Poly, index []uint64, polOut Poly) {

	level := r.level

	N := r.N()

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(index)%8 != 0  */
		x := (*[8]uint64)(unsafe.Pointer(&index[j]))

		for i := 0; i < level+1; i++ {

			/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(polOut.Coeffs)%8 != 0 */
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

// AutomorphismNTTWithIndexThenAddLazy applies the automorphism X^{i} -> X^{i*gen} on a polynomial in the NTT domain .
// `index` is the lookup table storing the mapping of the automorphism.
// The result of the automorphism is added on polOut.
func (r Ring) AutomorphismNTTWithIndexThenAddLazy(polIn Poly, index []uint64, polOut Poly) {

	level := r.level

	N := r.N()

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(index)%8 != 0 */
		x := (*[8]uint64)(unsafe.Pointer(&index[j]))

		for i := 0; i < level+1; i++ {

			/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(polOut.Coeffs)%8 != 0 */
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

// Automorphism applies the automorphism X^{i} -> X^{i*gen} on a polynomial outside of the NTT domain.
// It must be noted that the result cannot be in-place.
func (r Ring) Automorphism(polIn Poly, gen uint64, polOut Poly) {

	var mask, index, indexRaw, logN, tmp uint64

	N := uint64(r.N())

	level := r.level

	if r.Type() == ConjugateInvariant {

		mask = 2*N - 1

		logN = uint64(bits.Len64(mask))

		// TODO: find a more efficient way to do
		// the automorphism on Z[X+X^-1]
		for i := uint64(0); i < 2*N; i++ {

			indexRaw = i * gen

			index = indexRaw & mask

			tmp = (indexRaw >> logN) & 1

			// Only consider i -> index if within [0, N-1]
			if index < N {

				idx := i

				// If the starting index is within [N, 2N-1]
				if idx >= N {
					idx = 2*N - idx // Wrap back between [0, N-1]
					tmp ^= 1        // Negate
				}

				for j, s := range r.SubRings[:level+1] {
					polOut.Coeffs[j][index] = polIn.Coeffs[j][idx]*(tmp^1) | (s.Modulus-polIn.Coeffs[j][idx])*tmp
				}
			}
		}

	} else {

		mask = N - 1

		logN = uint64(bits.Len64(mask))

		for i := uint64(0); i < N; i++ {

			indexRaw = i * gen

			index = indexRaw & mask

			tmp = (indexRaw >> logN) & 1

			for j, s := range r.SubRings[:level+1] {
				polOut.Coeffs[j][index] = polIn.Coeffs[j][i]*(tmp^1) | (s.Modulus-polIn.Coeffs[j][i])*tmp
			}
		}
	}
}
