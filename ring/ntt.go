package ring

import (
	"fmt"
	"math/bits"
	"unsafe"
)

const (
	// MinimumRingDegreeForLoopUnrolledNTT is the minimum ring degree
	// necessary for memory safe loop unrolling
	MinimumRingDegreeForLoopUnrolledNTT = 16
)

// NumberTheoreticTransformer is an interface to provide
// flexibility on what type of NTT is used by the struct Ring.
type NumberTheoreticTransformer interface {
	Forward(p1, p2 []uint64)
	ForwardLazy(p1, p2 []uint64)
	Backward(p1, p2 []uint64)
	BackwardLazy(p1, p2 []uint64)
}

type numberTheoreticTransformerBase struct {
	*NTTTable
	N            int
	Modulus      uint64
	MRedConstant uint64
	BRedConstant [2]uint64
}

// NumberTheoreticTransformerStandard computes the standard nega-cyclic NTT in the ring Z[X]/(X^N+1).
type NumberTheoreticTransformerStandard struct {
	numberTheoreticTransformerBase
}

// NTTTable store all the constants that are specifically tied to the NTT.
type NTTTable struct {
	NthRoot       uint64   // Nthroot used for the NTT
	PrimitiveRoot uint64   // 2N-th primitive root
	RootsForward  []uint64 //powers of the 2N-th primitive root in Montgomery form (in bit-reversed order)
	RootsBackward []uint64 //powers of the inverse of the 2N-th primitive root in Montgomery form (in bit-reversed order)
	NInv          uint64   //[N^-1] mod Modulus in Montgomery form
}

func NewNumberTheoreticTransformerStandard(r *SubRing, n int) NumberTheoreticTransformer {
	return NumberTheoreticTransformerStandard{
		numberTheoreticTransformerBase: numberTheoreticTransformerBase{
			N:            r.N,
			Modulus:      r.Modulus,
			MRedConstant: r.MRedConstant,
			BRedConstant: r.BRedConstant,
			NTTTable:     r.NTTTable,
		},
	}
}

// Forward writes the forward NTT in Z[X]/(X^N+1) of p1 on p2.
func (rntt NumberTheoreticTransformerStandard) Forward(p1, p2 []uint64) {
	NTTStandard(p1, p2, rntt.N, rntt.Modulus, rntt.MRedConstant, rntt.BRedConstant, rntt.RootsForward)
}

// ForwardLazy writes the forward NTT in Z[X]/(X^N+1) of p1 on p2.
// Returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerStandard) ForwardLazy(p1, p2 []uint64) {
	NTTStandardLazy(p1, p2, rntt.N, rntt.Modulus, rntt.MRedConstant, rntt.RootsForward)
}

// Backward writes the backward NTT in Z[X]/(X^N+1) of p1 on p2.
func (rntt NumberTheoreticTransformerStandard) Backward(p1, p2 []uint64) {
	INTTStandard(p1, p2, rntt.N, rntt.NInv, rntt.Modulus, rntt.MRedConstant, rntt.RootsBackward)
}

// BackwardLazy writes the backward NTT in Z[X]/(X^N+1) p1 on p2.
// Returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerStandard) BackwardLazy(p1, p2 []uint64) {
	INTTStandardLazy(p1, p2, rntt.N, rntt.NInv, rntt.Modulus, rntt.MRedConstant, rntt.RootsBackward)
}

// NumberTheoreticTransformerConjugateInvariant computes the NTT in the ring Z[X+X^-1]/(X^2N+1).
// Z[X+X^-1]/(X^2N+1) is a closed sub-ring of Z[X]/(X^2N+1). Note that the input polynomial only needs to be size N
// since the right half does not provide any additional information.
// See "Approximate Homomorphic Encryption over the Conjugate-invariant Ring", https://eprint.iacr.org/2018/952.
// The implemented approach is more efficient than the one proposed in the referenced work.
// It avoids the linear map Z[X + X^-1]/(X^2N + 1) <-> Z[X]/(X^N - 1) by instead directly computing the left
// half of the NTT of Z[X + X^-1]/(X^2N + 1) since the right half provides no additional information, which
// allows to (re)use nega-cyclic NTT.
type NumberTheoreticTransformerConjugateInvariant struct {
	numberTheoreticTransformerBase
}

func NewNumberTheoreticTransformerConjugateInvariant(r *SubRing, n int) NumberTheoreticTransformer {
	return NumberTheoreticTransformerConjugateInvariant{
		numberTheoreticTransformerBase: numberTheoreticTransformerBase{
			N:            r.N,
			Modulus:      r.Modulus,
			MRedConstant: r.MRedConstant,
			BRedConstant: r.BRedConstant,
			NTTTable:     r.NTTTable,
		},
	}
}

// Forward writes the forward NTT in Z[X+X^-1]/(X^2N+1) of p1 on p2.
func (rntt NumberTheoreticTransformerConjugateInvariant) Forward(p1, p2 []uint64) {
	NTTConjugateInvariant(p1, p2, rntt.N, rntt.Modulus, rntt.MRedConstant, rntt.BRedConstant, rntt.RootsForward)
}

// ForwardLazy writes the forward NTT in Z[X+X^-1]/(X^2N+1) of p1 on p2.
// Returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerConjugateInvariant) ForwardLazy(p1, p2 []uint64) {
	NTTConjugateInvariantLazy(p1, p2, rntt.N, rntt.Modulus, rntt.MRedConstant, rntt.RootsForward)
}

// Backward writes the backward NTT in Z[X+X^-1]/(X^2N+1) of p1 on p2.
func (rntt NumberTheoreticTransformerConjugateInvariant) Backward(p1, p2 []uint64) {
	INTTConjugateInvariant(p1, p2, rntt.N, rntt.NInv, rntt.Modulus, rntt.MRedConstant, rntt.RootsBackward)
}

// BackwardLazy writes the backward NTT in Z[X+X^-1]/(X^2N+1) of p1 on p2.
// Returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerConjugateInvariant) BackwardLazy(p1, p2 []uint64) {
	INTTConjugateInvariantLazy(p1, p2, rntt.N, rntt.NInv, rntt.Modulus, rntt.MRedConstant, rntt.RootsBackward)
}

// NTT evaluates p2 = NTT(P1).
func (r Ring) NTT(p1, p2 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.NTT(p1.Coeffs[i], p2.Coeffs[i])
	}
}

// NTTLazy evaluates p2 = NTT(p1) with p2 in [0, 2*modulus-1].
func (r Ring) NTTLazy(p1, p2 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.NTTLazy(p1.Coeffs[i], p2.Coeffs[i])
	}
}

// INTT evaluates p2 = INTT(p1).
func (r Ring) INTT(p1, p2 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.INTT(p1.Coeffs[i], p2.Coeffs[i])
	}
}

// INTTLazy evaluates p2 = INTT(p1) with p2 in [0, 2*modulus-1].
func (r Ring) INTTLazy(p1, p2 Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.INTTLazy(p1.Coeffs[i], p2.Coeffs[i])
	}
}

// butterfly computes X, Y = U + V*Psi, U - V*Psi mod Q.
func butterfly(U, V, Psi, twoQ, fourQ, Q, MRedConstant uint64) (uint64, uint64) {
	if U >= fourQ {
		U -= fourQ
	}
	V = MRedLazy(V, Psi, Q, MRedConstant)
	return U + V, U + twoQ - V
}

// invbutterfly computes X, Y = U + V, (U - V) * Psi mod Q.
func invbutterfly(U, V, Psi, twoQ, fourQ, Q, MRedConstant uint64) (X, Y uint64) {
	X = U + V
	if X >= twoQ {
		X -= twoQ
	}
	Y = MRedLazy(U+fourQ-V, Psi, Q, MRedConstant) // At the moment it is not possible to use MRedLazy if Q > 61 bits
	return
}

// NTTStandard computes the NTTStandard in the given SubRing.
func NTTStandard(p1, p2 []uint64, N int, Q, MRedConstant uint64, BRedConstant [2]uint64, roots []uint64) {
	nttCoreLazy(p1, p2, N, Q, MRedConstant, roots)
	reducevec(p2, p2, Q, BRedConstant)
}

// NTTStandardLazy computes the NTTStandard in the given SubRing with p2 in [0, 2*modulus-1].
func NTTStandardLazy(p1, p2 []uint64, N int, Q, MRedConstant uint64, roots []uint64) {
	nttCoreLazy(p1, p2, N, Q, MRedConstant, roots)
}

// INTTStandard evaluates p2 = INTTStandard(p1) in the given SubRing.
func INTTStandard(p1, p2 []uint64, N int, NInv, Q, MRedConstant uint64, roots []uint64) {
	inttCoreLazy(p1, p2, N, Q, MRedConstant, roots)
	if N < MinimumRingDegreeForLoopUnrolledNTT {
		for i := 0; i < N; i++ {
			p2[i] = MRed(p2[i], NInv, Q, MRedConstant)
		}
	} else {
		mulscalarmontgomeryvec(p2, NInv, p2, Q, MRedConstant)
	}
}

// INTTStandardLazy evaluates p2 = INTT(p1) in the given SubRing with p2 in [0, 2*modulus-1].
func INTTStandardLazy(p1, p2 []uint64, N int, NInv, Q, MRedConstant uint64, roots []uint64) {
	inttCoreLazy(p1, p2, N, Q, MRedConstant, roots)
	if N < MinimumRingDegreeForLoopUnrolledNTT {
		for i := 0; i < N; i++ {
			p2[i] = MRedLazy(p2[i], NInv, Q, MRedConstant)
		}
	} else {
		mulscalarmontgomeryvec(p2, NInv, p2, Q, MRedConstant)
	}
}

// nttCoreLazy computes the NTT on the input coefficients using the input parameters with output values in the range [0, 2*modulus-1].
func nttCoreLazy(p1, p2 []uint64, N int, Q, MRedConstant uint64, roots []uint64) {

	// Sanity check
	if len(p1) < N || len(p2) < N || len(roots) < N {
		panic(fmt.Sprintf("cannot nttCoreLazy: ensure that len(p1)=%d, len(p2)=%d and len(roots)=%d >= N=%d", len(p1), len(p2), len(roots), N))
	}

	if N < MinimumRingDegreeForLoopUnrolledNTT {
		nttLazy(p1, p2, N, Q, MRedConstant, roots)
	} else {
		nttUnrolled16Lazy(p1, p2, N, Q, MRedConstant, roots)
	}
}

func nttLazy(p1, p2 []uint64, N int, Q, MRedConstant uint64, roots []uint64) {

	var j1, j2, t int
	var F uint64

	fourQ := 4 * Q
	twoQ := 2 * Q

	t = N >> 1
	F = roots[1]
	j1 = 0
	j2 = j1 + t

	for jx, jy := j1, j1+t; jx < j2; jx, jy = jx+1, jy+1 {
		p2[jx], p2[jy] = butterfly(p1[jx], p1[jy], F, twoQ, fourQ, Q, MRedConstant)
	}

	for m := 2; m < N; m <<= 1 {

		t >>= 1

		for i := 0; i < m; i++ {

			j1 = (i * t) << 1

			j2 = j1 + t

			F = roots[m+i]

			for jx, jy := j1, j1+t; jx < j2; jx, jy = jx+1, jy+1 {
				p2[jx], p2[jy] = butterfly(p2[jx], p2[jy], F, twoQ, fourQ, Q, MRedConstant)
			}
		}
	}
}
func nttUnrolled16Lazy(p1, p2 []uint64, N int, Q, MRedConstant uint64, roots []uint64) {

	// Sanity check
	if len(p2) < MinimumRingDegreeForLoopUnrolledNTT {
		panic(fmt.Sprintf("unsafe call of nttUnrolled16Lazy: receiver len(p2)=%d < %d", len(p2), MinimumRingDegreeForLoopUnrolledNTT))
	}

	var j1, j2, t int
	var F, V uint64

	fourQ := 4 * Q
	twoQ := 2 * Q

	// Copy the result of the first round of butterflies on p2 with approximate reduction
	t = N >> 1
	F = roots[1]

	for jx, jy := 0, t; jx < t; jx, jy = jx+8, jy+8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 != 0 */
		xin := (*[8]uint64)(unsafe.Pointer(&p1[jx]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 != 0 */
		yin := (*[8]uint64)(unsafe.Pointer(&p1[jy]))

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 != 0 */
		xout := (*[8]uint64)(unsafe.Pointer(&p2[jx]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 != 0 */
		yout := (*[8]uint64)(unsafe.Pointer(&p2[jy]))

		V = MRedLazy(yin[0], F, Q, MRedConstant)
		xout[0], yout[0] = xin[0]+V, xin[0]+twoQ-V

		V = MRedLazy(yin[1], F, Q, MRedConstant)
		xout[1], yout[1] = xin[1]+V, xin[1]+twoQ-V

		V = MRedLazy(yin[2], F, Q, MRedConstant)
		xout[2], yout[2] = xin[2]+V, xin[2]+twoQ-V

		V = MRedLazy(yin[3], F, Q, MRedConstant)
		xout[3], yout[3] = xin[3]+V, xin[3]+twoQ-V

		V = MRedLazy(yin[4], F, Q, MRedConstant)
		xout[4], yout[4] = xin[4]+V, xin[4]+twoQ-V

		V = MRedLazy(yin[5], F, Q, MRedConstant)
		xout[5], yout[5] = xin[5]+V, xin[5]+twoQ-V

		V = MRedLazy(yin[6], F, Q, MRedConstant)
		xout[6], yout[6] = xin[6]+V, xin[6]+twoQ-V

		V = MRedLazy(yin[7], F, Q, MRedConstant)
		xout[7], yout[7] = xin[7]+V, xin[7]+twoQ-V
	}

	// Continue the rest of the second to the n-1 butterflies on p2 with approximate reduction
	var reduce bool

	for m := 2; m < N; m <<= 1 {

		reduce = (bits.Len64(uint64(m))&1 == 1)

		t >>= 1

		if t >= 8 {

			for i := 0; i < m; i++ {

				j1 = (i * t) << 1

				j2 = j1 + t

				F = roots[m+i]

				if reduce {

					for jx, jy := j1, j1+t; jx < j2; jx, jy = jx+8, jy+8 {

						/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 != 0 */
						x := (*[8]uint64)(unsafe.Pointer(&p2[jx]))
						/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 != 0 */
						y := (*[8]uint64)(unsafe.Pointer(&p2[jy]))

						x[0], y[0] = butterfly(x[0], y[0], F, twoQ, fourQ, Q, MRedConstant)
						x[1], y[1] = butterfly(x[1], y[1], F, twoQ, fourQ, Q, MRedConstant)
						x[2], y[2] = butterfly(x[2], y[2], F, twoQ, fourQ, Q, MRedConstant)
						x[3], y[3] = butterfly(x[3], y[3], F, twoQ, fourQ, Q, MRedConstant)
						x[4], y[4] = butterfly(x[4], y[4], F, twoQ, fourQ, Q, MRedConstant)
						x[5], y[5] = butterfly(x[5], y[5], F, twoQ, fourQ, Q, MRedConstant)
						x[6], y[6] = butterfly(x[6], y[6], F, twoQ, fourQ, Q, MRedConstant)
						x[7], y[7] = butterfly(x[7], y[7], F, twoQ, fourQ, Q, MRedConstant)
					}

				} else {

					for jx, jy := j1, j1+t; jx < j2; jx, jy = jx+8, jy+8 {

						/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 != 0 */
						x := (*[8]uint64)(unsafe.Pointer(&p2[jx]))
						/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 != 0 */
						y := (*[8]uint64)(unsafe.Pointer(&p2[jy]))

						V = MRedLazy(y[0], F, Q, MRedConstant)
						x[0], y[0] = x[0]+V, x[0]+twoQ-V

						V = MRedLazy(y[1], F, Q, MRedConstant)
						x[1], y[1] = x[1]+V, x[1]+twoQ-V

						V = MRedLazy(y[2], F, Q, MRedConstant)
						x[2], y[2] = x[2]+V, x[2]+twoQ-V

						V = MRedLazy(y[3], F, Q, MRedConstant)
						x[3], y[3] = x[3]+V, x[3]+twoQ-V

						V = MRedLazy(y[4], F, Q, MRedConstant)
						x[4], y[4] = x[4]+V, x[4]+twoQ-V

						V = MRedLazy(y[5], F, Q, MRedConstant)
						x[5], y[5] = x[5]+V, x[5]+twoQ-V

						V = MRedLazy(y[6], F, Q, MRedConstant)
						x[6], y[6] = x[6]+V, x[6]+twoQ-V

						V = MRedLazy(y[7], F, Q, MRedConstant)
						x[7], y[7] = x[7]+V, x[7]+twoQ-V
					}
				}
			}

		} else if t == 4 {

			if reduce {

				for i, j1 := m, 0; i < 2*m; i, j1 = i+2, j1+4*t {

					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(roots)%2 != 0 */
					psi := (*[2]uint64)(unsafe.Pointer(&roots[i]))
					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%16 != 0 */
					x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

					x[0], x[4] = butterfly(x[0], x[4], psi[0], twoQ, fourQ, Q, MRedConstant)
					x[1], x[5] = butterfly(x[1], x[5], psi[0], twoQ, fourQ, Q, MRedConstant)
					x[2], x[6] = butterfly(x[2], x[6], psi[0], twoQ, fourQ, Q, MRedConstant)
					x[3], x[7] = butterfly(x[3], x[7], psi[0], twoQ, fourQ, Q, MRedConstant)
					x[8], x[12] = butterfly(x[8], x[12], psi[1], twoQ, fourQ, Q, MRedConstant)
					x[9], x[13] = butterfly(x[9], x[13], psi[1], twoQ, fourQ, Q, MRedConstant)
					x[10], x[14] = butterfly(x[10], x[14], psi[1], twoQ, fourQ, Q, MRedConstant)
					x[11], x[15] = butterfly(x[11], x[15], psi[1], twoQ, fourQ, Q, MRedConstant)

				}
			} else {

				for i, j1 := m, 0; i < 2*m; i, j1 = i+2, j1+4*t {

					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(roots)%2 != 0 */
					psi := (*[2]uint64)(unsafe.Pointer(&roots[i]))
					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%16 != 0 */
					x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

					V = MRedLazy(x[4], psi[0], Q, MRedConstant)
					x[0], x[4] = x[0]+V, x[0]+twoQ-V

					V = MRedLazy(x[5], psi[0], Q, MRedConstant)
					x[1], x[5] = x[1]+V, x[1]+twoQ-V

					V = MRedLazy(x[6], psi[0], Q, MRedConstant)
					x[2], x[6] = x[2]+V, x[2]+twoQ-V

					V = MRedLazy(x[7], psi[0], Q, MRedConstant)
					x[3], x[7] = x[3]+V, x[3]+twoQ-V

					V = MRedLazy(x[12], psi[1], Q, MRedConstant)
					x[8], x[12] = x[8]+V, x[8]+twoQ-V

					V = MRedLazy(x[13], psi[1], Q, MRedConstant)
					x[9], x[13] = x[9]+V, x[9]+twoQ-V

					V = MRedLazy(x[14], psi[1], Q, MRedConstant)
					x[10], x[14] = x[10]+V, x[10]+twoQ-V

					V = MRedLazy(x[15], psi[1], Q, MRedConstant)
					x[11], x[15] = x[11]+V, x[11]+twoQ-V

				}

			}

		} else if t == 2 {

			if reduce {

				for i, j1 := m, 0; i < 2*m; i, j1 = i+4, j1+8*t {

					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(roots)%4 != 0 */
					psi := (*[4]uint64)(unsafe.Pointer(&roots[i]))
					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%16 != 0 */
					x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

					x[0], x[2] = butterfly(x[0], x[2], psi[0], twoQ, fourQ, Q, MRedConstant)
					x[1], x[3] = butterfly(x[1], x[3], psi[0], twoQ, fourQ, Q, MRedConstant)
					x[4], x[6] = butterfly(x[4], x[6], psi[1], twoQ, fourQ, Q, MRedConstant)
					x[5], x[7] = butterfly(x[5], x[7], psi[1], twoQ, fourQ, Q, MRedConstant)
					x[8], x[10] = butterfly(x[8], x[10], psi[2], twoQ, fourQ, Q, MRedConstant)
					x[9], x[11] = butterfly(x[9], x[11], psi[2], twoQ, fourQ, Q, MRedConstant)
					x[12], x[14] = butterfly(x[12], x[14], psi[3], twoQ, fourQ, Q, MRedConstant)
					x[13], x[15] = butterfly(x[13], x[15], psi[3], twoQ, fourQ, Q, MRedConstant)
				}
			} else {

				for i, j1 := m, 0; i < 2*m; i, j1 = i+4, j1+8*t {

					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(roots)%4 != 0 */
					psi := (*[4]uint64)(unsafe.Pointer(&roots[i]))
					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%16 != 0 */
					x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

					V = MRedLazy(x[2], psi[0], Q, MRedConstant)
					x[0], x[2] = x[0]+V, x[0]+twoQ-V

					V = MRedLazy(x[3], psi[0], Q, MRedConstant)
					x[1], x[3] = x[1]+V, x[1]+twoQ-V

					V = MRedLazy(x[6], psi[1], Q, MRedConstant)
					x[4], x[6] = x[4]+V, x[4]+twoQ-V

					V = MRedLazy(x[7], psi[1], Q, MRedConstant)
					x[5], x[7] = x[5]+V, x[5]+twoQ-V

					V = MRedLazy(x[10], psi[2], Q, MRedConstant)
					x[8], x[10] = x[8]+V, x[8]+twoQ-V

					V = MRedLazy(x[11], psi[2], Q, MRedConstant)
					x[9], x[11] = x[9]+V, x[9]+twoQ-V

					V = MRedLazy(x[14], psi[3], Q, MRedConstant)
					x[12], x[14] = x[12]+V, x[12]+twoQ-V

					V = MRedLazy(x[15], psi[3], Q, MRedConstant)
					x[13], x[15] = x[13]+V, x[13]+twoQ-V
				}
			}

		} else {

			for i, j1 := m, 0; i < 2*m; i, j1 = i+8, j1+16 {

				/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(roots)%8 != 0 */
				psi := (*[8]uint64)(unsafe.Pointer(&roots[i]))
				/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%16 != 0 */
				x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

				x[0], x[1] = butterfly(x[0], x[1], psi[0], twoQ, fourQ, Q, MRedConstant)
				x[2], x[3] = butterfly(x[2], x[3], psi[1], twoQ, fourQ, Q, MRedConstant)
				x[4], x[5] = butterfly(x[4], x[5], psi[2], twoQ, fourQ, Q, MRedConstant)
				x[6], x[7] = butterfly(x[6], x[7], psi[3], twoQ, fourQ, Q, MRedConstant)
				x[8], x[9] = butterfly(x[8], x[9], psi[4], twoQ, fourQ, Q, MRedConstant)
				x[10], x[11] = butterfly(x[10], x[11], psi[5], twoQ, fourQ, Q, MRedConstant)
				x[12], x[13] = butterfly(x[12], x[13], psi[6], twoQ, fourQ, Q, MRedConstant)
				x[14], x[15] = butterfly(x[14], x[15], psi[7], twoQ, fourQ, Q, MRedConstant)
			}

			/*
				for i := uint64(0); i < m; i = i + 8 {

					psi := (*[8]uint64)(unsafe.Pointer(&roots[m+i]))
					x := (*[16]uint64)(unsafe.Pointer(&p2[2*i]))

					V = MRedLazy(x[1], psi[0], Q, MRedConstant)
					x[0], x[1] = x[0]+V, x[0]+twoQ-V

					V = MRedLazy(x[3], psi[1], Q, MRedConstant)
					x[2], x[3] = x[2]+V, x[2]+twoQ-V

					V = MRedLazy(x[5], psi[2], Q, MRedConstant)
					x[4], x[5] = x[4]+V, x[4]+twoQ-V

					V = MRedLazy(x[7], psi[3], Q, MRedConstant)
					x[6], x[7] = x[6]+V, x[6]+twoQ-V

					V = MRedLazy(x[9], psi[4], Q, MRedConstant)
					x[8], x[9] = x[8]+V, x[8]+twoQ-V

					V = MRedLazy(x[11], psi[5], Q, MRedConstant)
					x[10], x[11] = x[10]+V, x[10]+twoQ-V

					V = MRedLazy(x[13], psi[6], Q, MRedConstant)
					x[12], x[13] = x[12]+V, x[12]+twoQ-V

					V = MRedLazy(x[15], psi[7], Q, MRedConstant)
					x[14], x[15] = x[14]+V, x[14]+twoQ-V
				}
			*/
		}
	}
}

func inttCoreLazy(p1, p2 []uint64, N int, Q, MRedConstant uint64, roots []uint64) {

	// Sanity check
	if len(p1) < N || len(p2) < N || len(roots) < N {
		panic(fmt.Sprintf("cannot inttCoreLazy: ensure that len(p1)=%d, len(p2)=%d and len(roots)=%d >= N=%d", len(p1), len(p2), len(roots), N))
	}

	if N < MinimumRingDegreeForLoopUnrolledNTT {
		inttLazy(p1, p2, N, Q, MRedConstant, roots)
	} else {
		inttLazyUnrolled16(p1, p2, N, Q, MRedConstant, roots)
	}
}

func inttLazy(p1, p2 []uint64, N int, Q, MRedConstant uint64, roots []uint64) {
	var h, t int
	var F uint64

	// Copy the result of the first round of butterflies on p2 with approximate reduction
	t = 1
	h = N >> 1
	twoQ := Q << 1
	fourQ := Q << 2

	for i, j1, j2 := 0, 0, t; i < h; i, j1, j2 = i+1, j1+2*t, j2+2*t {

		F = roots[h+i]

		for jx, jy := j1, j1+t; jx < j2; jx, jy = jx+1, jy+1 {
			p2[jx], p2[jy] = invbutterfly(p1[jx], p1[jy], F, twoQ, fourQ, Q, MRedConstant)

		}
	}

	t <<= 1

	for m := N >> 1; m > 1; m >>= 1 {

		h = m >> 1

		for i, j1, j2 := 0, 0, t; i < h; i, j1, j2 = i+1, j1+2*t, j2+2*t {

			F = roots[h+i]

			for jx, jy := j1, j1+t; jx < j2; jx, jy = jx+1, jy+1 {
				p2[jx], p2[jy] = invbutterfly(p2[jx], p2[jy], F, twoQ, fourQ, Q, MRedConstant)

			}
		}

		t <<= 1
	}
}

func inttLazyUnrolled16(p1, p2 []uint64, N int, Q, MRedConstant uint64, roots []uint64) {

	// Sanity check
	if len(p2) < MinimumRingDegreeForLoopUnrolledNTT {
		panic(fmt.Sprintf("unsafe call of inttCoreUnrolled16Lazy: receiver len(p2)=%d < %d", len(p2), MinimumRingDegreeForLoopUnrolledNTT))
	}

	var h, t int
	var F uint64

	// Copy the result of the first round of butterflies on p2 with approximate reduction
	t = 1
	h = N >> 1
	twoQ := Q << 1
	fourQ := Q << 2

	for i, j := h, 0; i < 2*h; i, j = i+8, j+16 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(roots)%8 != 0 */
		psi := (*[8]uint64)(unsafe.Pointer(&roots[i]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%16 != 0 */
		xin := (*[16]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%16 != 0 */
		xout := (*[16]uint64)(unsafe.Pointer(&p2[j]))

		xout[0], xout[1] = invbutterfly(xin[0], xin[1], psi[0], twoQ, fourQ, Q, MRedConstant)
		xout[2], xout[3] = invbutterfly(xin[2], xin[3], psi[1], twoQ, fourQ, Q, MRedConstant)
		xout[4], xout[5] = invbutterfly(xin[4], xin[5], psi[2], twoQ, fourQ, Q, MRedConstant)
		xout[6], xout[7] = invbutterfly(xin[6], xin[7], psi[3], twoQ, fourQ, Q, MRedConstant)
		xout[8], xout[9] = invbutterfly(xin[8], xin[9], psi[4], twoQ, fourQ, Q, MRedConstant)
		xout[10], xout[11] = invbutterfly(xin[10], xin[11], psi[5], twoQ, fourQ, Q, MRedConstant)
		xout[12], xout[13] = invbutterfly(xin[12], xin[13], psi[6], twoQ, fourQ, Q, MRedConstant)
		xout[14], xout[15] = invbutterfly(xin[14], xin[15], psi[7], twoQ, fourQ, Q, MRedConstant)
	}

	// Continue the rest of the second to the n-1 butterflies on p2 with approximate reduction
	t <<= 1
	for m := N >> 1; m > 1; m >>= 1 {

		h = m >> 1

		if t >= 8 {

			for i, j1, j2 := 0, 0, t; i < h; i, j1, j2 = i+1, j1+2*t, j2+2*t {

				F = roots[h+i]

				for jx, jy := j1, j1+t; jx < j2; jx, jy = jx+8, jy+8 {

					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 != 0 */
					x := (*[8]uint64)(unsafe.Pointer(&p2[jx]))
					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 != 0 */
					y := (*[8]uint64)(unsafe.Pointer(&p2[jy]))

					x[0], y[0] = invbutterfly(x[0], y[0], F, twoQ, fourQ, Q, MRedConstant)
					x[1], y[1] = invbutterfly(x[1], y[1], F, twoQ, fourQ, Q, MRedConstant)
					x[2], y[2] = invbutterfly(x[2], y[2], F, twoQ, fourQ, Q, MRedConstant)
					x[3], y[3] = invbutterfly(x[3], y[3], F, twoQ, fourQ, Q, MRedConstant)
					x[4], y[4] = invbutterfly(x[4], y[4], F, twoQ, fourQ, Q, MRedConstant)
					x[5], y[5] = invbutterfly(x[5], y[5], F, twoQ, fourQ, Q, MRedConstant)
					x[6], y[6] = invbutterfly(x[6], y[6], F, twoQ, fourQ, Q, MRedConstant)
					x[7], y[7] = invbutterfly(x[7], y[7], F, twoQ, fourQ, Q, MRedConstant)
				}
			}

		} else if t == 4 {

			for i, j1 := h, 0; i < 2*h; i, j1 = i+2, j1+4*t {

				/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(roots)%2 != 0 */
				psi := (*[2]uint64)(unsafe.Pointer(&roots[i]))
				/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%16 != 0 */
				x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

				x[0], x[4] = invbutterfly(x[0], x[4], psi[0], twoQ, fourQ, Q, MRedConstant)
				x[1], x[5] = invbutterfly(x[1], x[5], psi[0], twoQ, fourQ, Q, MRedConstant)
				x[2], x[6] = invbutterfly(x[2], x[6], psi[0], twoQ, fourQ, Q, MRedConstant)
				x[3], x[7] = invbutterfly(x[3], x[7], psi[0], twoQ, fourQ, Q, MRedConstant)
				x[8], x[12] = invbutterfly(x[8], x[12], psi[1], twoQ, fourQ, Q, MRedConstant)
				x[9], x[13] = invbutterfly(x[9], x[13], psi[1], twoQ, fourQ, Q, MRedConstant)
				x[10], x[14] = invbutterfly(x[10], x[14], psi[1], twoQ, fourQ, Q, MRedConstant)
				x[11], x[15] = invbutterfly(x[11], x[15], psi[1], twoQ, fourQ, Q, MRedConstant)
			}

		} else {

			for i, j1 := h, 0; i < 2*h; i, j1 = i+4, j1+8*t {

				/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(roots)%4 != 0 */
				psi := (*[4]uint64)(unsafe.Pointer(&roots[i]))
				/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%16 != 0 */
				x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

				x[0], x[2] = invbutterfly(x[0], x[2], psi[0], twoQ, fourQ, Q, MRedConstant)
				x[1], x[3] = invbutterfly(x[1], x[3], psi[0], twoQ, fourQ, Q, MRedConstant)
				x[4], x[6] = invbutterfly(x[4], x[6], psi[1], twoQ, fourQ, Q, MRedConstant)
				x[5], x[7] = invbutterfly(x[5], x[7], psi[1], twoQ, fourQ, Q, MRedConstant)
				x[8], x[10] = invbutterfly(x[8], x[10], psi[2], twoQ, fourQ, Q, MRedConstant)
				x[9], x[11] = invbutterfly(x[9], x[11], psi[2], twoQ, fourQ, Q, MRedConstant)
				x[12], x[14] = invbutterfly(x[12], x[14], psi[3], twoQ, fourQ, Q, MRedConstant)
				x[13], x[15] = invbutterfly(x[13], x[15], psi[3], twoQ, fourQ, Q, MRedConstant)
			}
		}

		t <<= 1
	}
}

// NTTConjugateInvariant evaluates p2 = NTT(p1) in the sub-ring Z[X + X^-1]/(X^2N +1) of Z[X]/(X^2N+1).
func NTTConjugateInvariant(p1, p2 []uint64, N int, Q, MRedConstant uint64, BRedConstant [2]uint64, roots []uint64) {
	nttCoreConjugateInvariantLazy(p1, p2, N, Q, MRedConstant, roots)
	reducevec(p2, p2, Q, BRedConstant)
}

// NTTConjugateInvariantLazy evaluates p2 = NTT(p1) in the sub-ring Z[X + X^-1]/(X^2N +1) of Z[X]/(X^2N+1) with p2 in the range [0, 2*modulus-1].
func NTTConjugateInvariantLazy(p1, p2 []uint64, N int, Q, MRedConstant uint64, roots []uint64) {
	nttCoreConjugateInvariantLazy(p1, p2, N, Q, MRedConstant, roots)
}

// INTTConjugateInvariant evaluates p2 = INTT(p1) in the closed sub-ring Z[X + X^-1]/(X^2N +1) of Z[X]/(X^2N+1).
func INTTConjugateInvariant(p1, p2 []uint64, N int, NInv, Q, MRedConstant uint64, roots []uint64) {
	inttCoreConjugateInvariantLazy(p1, p2, N, Q, MRedConstant, roots)
	mulscalarmontgomeryvec(p2, NInv, p2, Q, MRedConstant)
}

// INTTConjugateInvariantLazy evaluates p2 = INTT(p1) in the closed sub-ring Z[X + X^-1]/(X^2N +1) of Z[X]/(X^2N+1) with p2 in the range [0, 2*modulus-1].
func INTTConjugateInvariantLazy(p1, p2 []uint64, N int, NInv, Q, MRedConstant uint64, roots []uint64) {
	inttCoreConjugateInvariantLazy(p1, p2, N, Q, MRedConstant, roots)
	mulscalarmontgomerylazyvec(p2, NInv, p2, Q, MRedConstant)
}

// nttCoreConjugateInvariantLazy evaluates p2 = NTT(p1) in the sub-ring Z[X + X^-1]/(X^2N +1) of Z[X]/(X^2N+1) with p2 [0, 2*modulus-1].
func nttCoreConjugateInvariantLazy(p1, p2 []uint64, N int, Q, MRedConstant uint64, roots []uint64) {

	// Sanity check
	if len(p1) < N || len(p2) < N || len(roots) < N {
		panic(fmt.Sprintf("cannot nttCoreConjugateInvariantLazy: ensure that len(p1)=%d, len(p2)=%d and len(roots)=%d >= N=%d", len(p1), len(p2), len(roots), N))
	}

	if N < MinimumRingDegreeForLoopUnrolledNTT {
		nttConjugateInvariantLazy(p1, p2, N, Q, MRedConstant, roots)
	} else {
		nttConjugateInvariantLazyUnrolled16(p1, p2, N, Q, MRedConstant, roots)
	}
}

func nttConjugateInvariantLazy(p1, p2 []uint64, N int, Q, MRedConstant uint64, roots []uint64) {
	var t, h int
	var F uint64

	fourQ := 4 * Q
	twoQ := 2 * Q

	t = N
	F = roots[1]

	for jx, jy := 1, N-1; jx < (N >> 1); jx, jy = jx+1, jy-1 {
		p2[jx], p2[jy] = p1[jx]+twoQ-MRedLazy(p1[jy], F, Q, MRedConstant), p1[jy]+twoQ-MRedLazy(p1[jx], F, Q, MRedConstant)
	}

	p2[N>>1] = p1[N>>1] + twoQ - MRedLazy(p1[N>>1], F, Q, MRedConstant)
	p2[0] = p1[0]

	// Continue the rest of the second to the n-1 butterflies on p2 with approximate reduction
	for m := 2; m < 2*N; m <<= 1 {

		t >>= 1
		h = m >> 1

		for i, j1, j2 := 0, 0, t; i < h; i, j1, j2 = i+1, j1+2*t, j2+2*t {

			F = roots[m+i]

			for jx, jy := j1, j1+t; jx < j2; jx, jy = jx+1, jy+1 {
				p2[jx], p2[jy] = butterfly(p2[jx], p2[jy], F, twoQ, fourQ, Q, MRedConstant)
			}
		}
	}
}

func nttConjugateInvariantLazyUnrolled16(p1, p2 []uint64, N int, Q, MRedConstant uint64, roots []uint64) {

	// Sanity check
	if len(p2) < MinimumRingDegreeForLoopUnrolledNTT {
		panic(fmt.Sprintf("unsafe call of nttCoreConjugateInvariantLazyUnrolled16: receiver len(p2)=%d < %d", len(p2), MinimumRingDegreeForLoopUnrolledNTT))
	}

	var t, h int
	var F, V uint64
	var reduce bool

	fourQ := 4 * Q
	twoQ := 2 * Q

	// Copy the result of the first round of butterflies on p2 with approximate reduction
	t = N
	F = roots[1]

	for jx, jy := 1, N-8; jx < (N>>1)-7; jx, jy = jx+8, jy-8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 != 0 */
		xin := (*[8]uint64)(unsafe.Pointer(&p1[jx]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 != 0 */
		yin := (*[8]uint64)(unsafe.Pointer(&p1[jy]))

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 != 0 */
		xout := (*[8]uint64)(unsafe.Pointer(&p2[jx]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 != 0 */
		yout := (*[8]uint64)(unsafe.Pointer(&p2[jy]))

		xout[0], yout[7] = xin[0]+twoQ-MRedLazy(yin[7], F, Q, MRedConstant), yin[7]+twoQ-MRedLazy(xin[0], F, Q, MRedConstant)
		xout[1], yout[6] = xin[1]+twoQ-MRedLazy(yin[6], F, Q, MRedConstant), yin[6]+twoQ-MRedLazy(xin[1], F, Q, MRedConstant)
		xout[2], yout[5] = xin[2]+twoQ-MRedLazy(yin[5], F, Q, MRedConstant), yin[5]+twoQ-MRedLazy(xin[2], F, Q, MRedConstant)
		xout[3], yout[4] = xin[3]+twoQ-MRedLazy(yin[4], F, Q, MRedConstant), yin[4]+twoQ-MRedLazy(xin[3], F, Q, MRedConstant)
		xout[4], yout[3] = xin[4]+twoQ-MRedLazy(yin[3], F, Q, MRedConstant), yin[3]+twoQ-MRedLazy(xin[4], F, Q, MRedConstant)
		xout[5], yout[2] = xin[5]+twoQ-MRedLazy(yin[2], F, Q, MRedConstant), yin[2]+twoQ-MRedLazy(xin[5], F, Q, MRedConstant)
		xout[6], yout[1] = xin[6]+twoQ-MRedLazy(yin[1], F, Q, MRedConstant), yin[1]+twoQ-MRedLazy(xin[6], F, Q, MRedConstant)
		xout[7], yout[0] = xin[7]+twoQ-MRedLazy(yin[0], F, Q, MRedConstant), yin[0]+twoQ-MRedLazy(xin[7], F, Q, MRedConstant)
	}

	j := (N >> 1) - 7
	/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 != 0 */
	xin := (*[7]uint64)(unsafe.Pointer(&p1[j]))
	/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 != 0 */
	yin := (*[7]uint64)(unsafe.Pointer(&p1[N-j-6]))
	/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 != 0 */
	xout := (*[7]uint64)(unsafe.Pointer(&p2[j]))
	/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 != 0 */
	yout := (*[7]uint64)(unsafe.Pointer(&p2[N-j-6]))

	xout[0], yout[6] = xin[0]+twoQ-MRedLazy(yin[6], F, Q, MRedConstant), yin[6]+twoQ-MRedLazy(xin[0], F, Q, MRedConstant)
	xout[1], yout[5] = xin[1]+twoQ-MRedLazy(yin[5], F, Q, MRedConstant), yin[5]+twoQ-MRedLazy(xin[1], F, Q, MRedConstant)
	xout[2], yout[4] = xin[2]+twoQ-MRedLazy(yin[4], F, Q, MRedConstant), yin[4]+twoQ-MRedLazy(xin[2], F, Q, MRedConstant)
	xout[3], yout[3] = xin[3]+twoQ-MRedLazy(yin[3], F, Q, MRedConstant), yin[3]+twoQ-MRedLazy(xin[3], F, Q, MRedConstant)
	xout[4], yout[2] = xin[4]+twoQ-MRedLazy(yin[2], F, Q, MRedConstant), yin[2]+twoQ-MRedLazy(xin[4], F, Q, MRedConstant)
	xout[5], yout[1] = xin[5]+twoQ-MRedLazy(yin[1], F, Q, MRedConstant), yin[1]+twoQ-MRedLazy(xin[5], F, Q, MRedConstant)
	xout[6], yout[0] = xin[6]+twoQ-MRedLazy(yin[0], F, Q, MRedConstant), yin[0]+twoQ-MRedLazy(xin[6], F, Q, MRedConstant)

	p2[N>>1] = p1[N>>1] + twoQ - MRedLazy(p1[N>>1], F, Q, MRedConstant)
	p2[0] = p1[0]

	// Continue the rest of the second to the n-1 butterflies on p2 with approximate reduction
	for m := 2; m < 2*N; m <<= 1 {

		reduce = (bits.Len64(uint64(m))&1 == 1)

		t >>= 1
		h = m >> 1

		if t >= 8 {

			for i, j1, j2 := 0, 0, t; i < h; i, j1, j2 = i+1, j1+2*t, j2+2*t {

				F = roots[m+i]

				if reduce {

					for jx, jy := j1, j1+t; jx < j2; jx, jy = jx+8, jy+8 {

						/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 != 0 */
						x := (*[8]uint64)(unsafe.Pointer(&p2[jx]))
						/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 != 0 */
						y := (*[8]uint64)(unsafe.Pointer(&p2[jy]))

						x[0], y[0] = butterfly(x[0], y[0], F, twoQ, fourQ, Q, MRedConstant)
						x[1], y[1] = butterfly(x[1], y[1], F, twoQ, fourQ, Q, MRedConstant)
						x[2], y[2] = butterfly(x[2], y[2], F, twoQ, fourQ, Q, MRedConstant)
						x[3], y[3] = butterfly(x[3], y[3], F, twoQ, fourQ, Q, MRedConstant)
						x[4], y[4] = butterfly(x[4], y[4], F, twoQ, fourQ, Q, MRedConstant)
						x[5], y[5] = butterfly(x[5], y[5], F, twoQ, fourQ, Q, MRedConstant)
						x[6], y[6] = butterfly(x[6], y[6], F, twoQ, fourQ, Q, MRedConstant)
						x[7], y[7] = butterfly(x[7], y[7], F, twoQ, fourQ, Q, MRedConstant)
					}

				} else {

					for jx, jy := j1, j1+t; jx < j2; jx, jy = jx+8, jy+8 {

						/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 != 0 */
						x := (*[8]uint64)(unsafe.Pointer(&p2[jx]))
						/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 != 0 */
						y := (*[8]uint64)(unsafe.Pointer(&p2[jy]))

						V = MRedLazy(y[0], F, Q, MRedConstant)
						x[0], y[0] = x[0]+V, x[0]+twoQ-V

						V = MRedLazy(y[1], F, Q, MRedConstant)
						x[1], y[1] = x[1]+V, x[1]+twoQ-V

						V = MRedLazy(y[2], F, Q, MRedConstant)
						x[2], y[2] = x[2]+V, x[2]+twoQ-V

						V = MRedLazy(y[3], F, Q, MRedConstant)
						x[3], y[3] = x[3]+V, x[3]+twoQ-V

						V = MRedLazy(y[4], F, Q, MRedConstant)
						x[4], y[4] = x[4]+V, x[4]+twoQ-V

						V = MRedLazy(y[5], F, Q, MRedConstant)
						x[5], y[5] = x[5]+V, x[5]+twoQ-V

						V = MRedLazy(y[6], F, Q, MRedConstant)
						x[6], y[6] = x[6]+V, x[6]+twoQ-V

						V = MRedLazy(y[7], F, Q, MRedConstant)
						x[7], y[7] = x[7]+V, x[7]+twoQ-V
					}
				}
			}

		} else if t == 4 {

			if reduce {

				for i, j1 := m, 0; i < h+m; i, j1 = i+2, j1+4*t {

					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(roots)%2 != 0 */
					psi := (*[2]uint64)(unsafe.Pointer(&roots[i]))
					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%16 != 0 */
					x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

					x[0], x[4] = butterfly(x[0], x[4], psi[0], twoQ, fourQ, Q, MRedConstant)
					x[1], x[5] = butterfly(x[1], x[5], psi[0], twoQ, fourQ, Q, MRedConstant)
					x[2], x[6] = butterfly(x[2], x[6], psi[0], twoQ, fourQ, Q, MRedConstant)
					x[3], x[7] = butterfly(x[3], x[7], psi[0], twoQ, fourQ, Q, MRedConstant)
					x[8], x[12] = butterfly(x[8], x[12], psi[1], twoQ, fourQ, Q, MRedConstant)
					x[9], x[13] = butterfly(x[9], x[13], psi[1], twoQ, fourQ, Q, MRedConstant)
					x[10], x[14] = butterfly(x[10], x[14], psi[1], twoQ, fourQ, Q, MRedConstant)
					x[11], x[15] = butterfly(x[11], x[15], psi[1], twoQ, fourQ, Q, MRedConstant)

				}
			} else {

				for i, j1 := m, 0; i < h+m; i, j1 = i+2, j1+4*t {

					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(roots)%2 != 0 */
					psi := (*[2]uint64)(unsafe.Pointer(&roots[i]))
					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%16 != 0 */
					x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

					V = MRedLazy(x[4], psi[0], Q, MRedConstant)
					x[0], x[4] = x[0]+V, x[0]+twoQ-V

					V = MRedLazy(x[5], psi[0], Q, MRedConstant)
					x[1], x[5] = x[1]+V, x[1]+twoQ-V

					V = MRedLazy(x[6], psi[0], Q, MRedConstant)
					x[2], x[6] = x[2]+V, x[2]+twoQ-V

					V = MRedLazy(x[7], psi[0], Q, MRedConstant)
					x[3], x[7] = x[3]+V, x[3]+twoQ-V

					V = MRedLazy(x[12], psi[1], Q, MRedConstant)
					x[8], x[12] = x[8]+V, x[8]+twoQ-V

					V = MRedLazy(x[13], psi[1], Q, MRedConstant)
					x[9], x[13] = x[9]+V, x[9]+twoQ-V

					V = MRedLazy(x[14], psi[1], Q, MRedConstant)
					x[10], x[14] = x[10]+V, x[10]+twoQ-V

					V = MRedLazy(x[15], psi[1], Q, MRedConstant)
					x[11], x[15] = x[11]+V, x[11]+twoQ-V

				}
			}

		} else if t == 2 {

			if reduce {

				for i, j1 := m, 0; i < h+m; i, j1 = i+4, j1+8*t {

					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(roots)%4 != 0 */
					psi := (*[4]uint64)(unsafe.Pointer(&roots[i]))
					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%16 != 0 */
					x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

					x[0], x[2] = butterfly(x[0], x[2], psi[0], twoQ, fourQ, Q, MRedConstant)
					x[1], x[3] = butterfly(x[1], x[3], psi[0], twoQ, fourQ, Q, MRedConstant)
					x[4], x[6] = butterfly(x[4], x[6], psi[1], twoQ, fourQ, Q, MRedConstant)
					x[5], x[7] = butterfly(x[5], x[7], psi[1], twoQ, fourQ, Q, MRedConstant)
					x[8], x[10] = butterfly(x[8], x[10], psi[2], twoQ, fourQ, Q, MRedConstant)
					x[9], x[11] = butterfly(x[9], x[11], psi[2], twoQ, fourQ, Q, MRedConstant)
					x[12], x[14] = butterfly(x[12], x[14], psi[3], twoQ, fourQ, Q, MRedConstant)
					x[13], x[15] = butterfly(x[13], x[15], psi[3], twoQ, fourQ, Q, MRedConstant)
				}
			} else {

				for i, j1 := m, 0; i < h+m; i, j1 = i+4, j1+8*t {

					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(roots)%4 != 0 */
					psi := (*[4]uint64)(unsafe.Pointer(&roots[i]))
					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%16 != 0 */
					x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

					V = MRedLazy(x[2], psi[0], Q, MRedConstant)
					x[0], x[2] = x[0]+V, x[0]+twoQ-V

					V = MRedLazy(x[3], psi[0], Q, MRedConstant)
					x[1], x[3] = x[1]+V, x[1]+twoQ-V

					V = MRedLazy(x[6], psi[1], Q, MRedConstant)
					x[4], x[6] = x[4]+V, x[4]+twoQ-V

					V = MRedLazy(x[7], psi[1], Q, MRedConstant)
					x[5], x[7] = x[5]+V, x[5]+twoQ-V

					V = MRedLazy(x[10], psi[2], Q, MRedConstant)
					x[8], x[10] = x[8]+V, x[8]+twoQ-V

					V = MRedLazy(x[11], psi[2], Q, MRedConstant)
					x[9], x[11] = x[9]+V, x[9]+twoQ-V

					V = MRedLazy(x[14], psi[3], Q, MRedConstant)
					x[12], x[14] = x[12]+V, x[12]+twoQ-V

					V = MRedLazy(x[15], psi[3], Q, MRedConstant)
					x[13], x[15] = x[13]+V, x[13]+twoQ-V
				}
			}

		} else {

			if reduce {

				for i, j1 := m, 0; i < h+m; i, j1 = i+8, j1+16 {

					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(roots)%8 != 0 */
					psi := (*[8]uint64)(unsafe.Pointer(&roots[i]))
					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%16 != 0 */
					x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

					x[0], x[1] = butterfly(x[0], x[1], psi[0], twoQ, fourQ, Q, MRedConstant)
					x[2], x[3] = butterfly(x[2], x[3], psi[1], twoQ, fourQ, Q, MRedConstant)
					x[4], x[5] = butterfly(x[4], x[5], psi[2], twoQ, fourQ, Q, MRedConstant)
					x[6], x[7] = butterfly(x[6], x[7], psi[3], twoQ, fourQ, Q, MRedConstant)
					x[8], x[9] = butterfly(x[8], x[9], psi[4], twoQ, fourQ, Q, MRedConstant)
					x[10], x[11] = butterfly(x[10], x[11], psi[5], twoQ, fourQ, Q, MRedConstant)
					x[12], x[13] = butterfly(x[12], x[13], psi[6], twoQ, fourQ, Q, MRedConstant)
					x[14], x[15] = butterfly(x[14], x[15], psi[7], twoQ, fourQ, Q, MRedConstant)
				}
			} else {

				for i, j1 := m, 0; i < h+m; i, j1 = i+8, j1+16 {

					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(roots)%16 != 0 */
					psi := (*[8]uint64)(unsafe.Pointer(&roots[i]))
					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%16 != 0 */
					x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

					V = MRedLazy(x[1], psi[0], Q, MRedConstant)
					x[0], x[1] = x[0]+V, x[0]+twoQ-V

					V = MRedLazy(x[3], psi[1], Q, MRedConstant)
					x[2], x[3] = x[2]+V, x[2]+twoQ-V

					V = MRedLazy(x[5], psi[2], Q, MRedConstant)
					x[4], x[5] = x[4]+V, x[4]+twoQ-V

					V = MRedLazy(x[7], psi[3], Q, MRedConstant)
					x[6], x[7] = x[6]+V, x[6]+twoQ-V

					V = MRedLazy(x[9], psi[4], Q, MRedConstant)
					x[8], x[9] = x[8]+V, x[8]+twoQ-V

					V = MRedLazy(x[11], psi[5], Q, MRedConstant)
					x[10], x[11] = x[10]+V, x[10]+twoQ-V

					V = MRedLazy(x[13], psi[6], Q, MRedConstant)
					x[12], x[13] = x[12]+V, x[12]+twoQ-V

					V = MRedLazy(x[15], psi[7], Q, MRedConstant)
					x[14], x[15] = x[14]+V, x[14]+twoQ-V
				}
			}
		}
	}
}

// inttCoreConjugateInvariantLazy evaluates p2 = INTT(p1) in the sub-ring Z[X + X^-1]/(X^2N +1) of Z[X]/(X^2N+1) with p2 [0, 2*modulus-1].
func inttCoreConjugateInvariantLazy(p1, p2 []uint64, N int, Q, MRedConstant uint64, roots []uint64) {

	// Sanity check
	if len(p1) < N || len(p2) < N || len(roots) < N {
		panic(fmt.Sprintf("cannot inttCoreConjugateInvariantLazy: ensure that len(p1)=%d, len(p2)=%d and len(roots)=%d >= N=%d", len(p1), len(p2), len(roots), N))
	}

	if N < MinimumRingDegreeForLoopUnrolledNTT {
		inttConjugateInvariantLazy(p1, p2, N, Q, MRedConstant, roots)
	} else {
		inttConjugateInvariantLazyUnrolled16(p1, p2, N, Q, MRedConstant, roots)
	}
}

func inttConjugateInvariantLazy(p1, p2 []uint64, N int, Q, MRedConstant uint64, roots []uint64) {
	var j1, j2, h, t int
	var F uint64

	twoQ := Q << 1
	fourQ := Q << 2

	t = 1
	h = N >> 1
	j1 = 0
	for i := 0; i < h; i++ {

		j2 = j1 + t

		F = roots[N+i]

		for jx, jy := j1, j1+t; jx < j2; jx, jy = jx+1, jy+1 {
			p2[jx], p2[jy] = invbutterfly(p1[jx], p1[jy], F, twoQ, fourQ, Q, MRedConstant)
		}

		j1 = j1 + (t << 1)
	}

	t <<= 1

	for m := N >> 1; m > 1; m >>= 1 {

		j1 = 0
		h = m >> 1

		for i := 0; i < h; i++ {

			j2 = j1 + t

			F = roots[m+i]

			for jx, jy := j1, j1+t; jx < j2; jx, jy = jx+1, jy+1 {
				p2[jx], p2[jy] = invbutterfly(p2[jx], p2[jy], F, twoQ, fourQ, Q, MRedConstant)
			}

			j1 = j1 + (t << 1)
		}

		t <<= 1
	}

	F = roots[1]

	for jx, jy := 1, N-1; jx < (N >> 1); jx, jy = jx+1, jy-1 {
		p2[jx], p2[jy] = p2[jx]+twoQ-MRedLazy(p2[jy], F, Q, MRedConstant), p2[jy]+twoQ-MRedLazy(p2[jx], F, Q, MRedConstant)
	}

	p2[N>>1] = p2[N>>1] + twoQ - MRedLazy(p2[N>>1], F, Q, MRedConstant)
	p2[0] = CRed(p2[0]<<1, Q)
}

func inttConjugateInvariantLazyUnrolled16(p1, p2 []uint64, N int, Q, MRedConstant uint64, roots []uint64) {

	// Sanity check
	if len(p2) < MinimumRingDegreeForLoopUnrolledNTT {
		panic(fmt.Sprintf("unsafe call of inttConjugateInvariantLazyUnrolled16: receiver len(p2)=%d < %d", len(p2), MinimumRingDegreeForLoopUnrolledNTT))
	}

	var j1, j2, h, t int
	var F uint64

	// Copy the result of the first round of butterflies on p2 with approximate reduction
	t = 1
	h = N >> 1
	twoQ := Q << 1
	fourQ := Q << 2

	for i, j := N, 0; i < h+N; i, j = i+8, j+16 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(roots)%8 != 0 */
		psi := (*[8]uint64)(unsafe.Pointer(&roots[i]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%16 != 0 */
		xin := (*[16]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%16 != 0 */
		xout := (*[16]uint64)(unsafe.Pointer(&p2[j]))

		xout[0], xout[1] = invbutterfly(xin[0], xin[1], psi[0], twoQ, fourQ, Q, MRedConstant)
		xout[2], xout[3] = invbutterfly(xin[2], xin[3], psi[1], twoQ, fourQ, Q, MRedConstant)
		xout[4], xout[5] = invbutterfly(xin[4], xin[5], psi[2], twoQ, fourQ, Q, MRedConstant)
		xout[6], xout[7] = invbutterfly(xin[6], xin[7], psi[3], twoQ, fourQ, Q, MRedConstant)
		xout[8], xout[9] = invbutterfly(xin[8], xin[9], psi[4], twoQ, fourQ, Q, MRedConstant)
		xout[10], xout[11] = invbutterfly(xin[10], xin[11], psi[5], twoQ, fourQ, Q, MRedConstant)
		xout[12], xout[13] = invbutterfly(xin[12], xin[13], psi[6], twoQ, fourQ, Q, MRedConstant)
		xout[14], xout[15] = invbutterfly(xin[14], xin[15], psi[7], twoQ, fourQ, Q, MRedConstant)
	}

	// Continue the rest of the second to the n-1 butterflies on p2 with approximate reduction
	t <<= 1
	for m := N >> 1; m > 1; m >>= 1 {

		j1 = 0
		h = m >> 1

		if t >= 8 {

			for i := 0; i < h; i++ {

				j2 = j1 + t

				F = roots[m+i]

				for jx, jy := j1, j1+t; jx < j2; jx, jy = jx+8, jy+8 {

					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
					x := (*[8]uint64)(unsafe.Pointer(&p2[jx]))
					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
					y := (*[8]uint64)(unsafe.Pointer(&p2[jy]))

					x[0], y[0] = invbutterfly(x[0], y[0], F, twoQ, fourQ, Q, MRedConstant)
					x[1], y[1] = invbutterfly(x[1], y[1], F, twoQ, fourQ, Q, MRedConstant)
					x[2], y[2] = invbutterfly(x[2], y[2], F, twoQ, fourQ, Q, MRedConstant)
					x[3], y[3] = invbutterfly(x[3], y[3], F, twoQ, fourQ, Q, MRedConstant)
					x[4], y[4] = invbutterfly(x[4], y[4], F, twoQ, fourQ, Q, MRedConstant)
					x[5], y[5] = invbutterfly(x[5], y[5], F, twoQ, fourQ, Q, MRedConstant)
					x[6], y[6] = invbutterfly(x[6], y[6], F, twoQ, fourQ, Q, MRedConstant)
					x[7], y[7] = invbutterfly(x[7], y[7], F, twoQ, fourQ, Q, MRedConstant)
				}

				j1 = j1 + (t << 1)
			}

		} else if t == 4 {

			for i := m; i < h+m; i = i + 2 {

				/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(roots)%2 */
				psi := (*[2]uint64)(unsafe.Pointer(&roots[i]))
				/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%16 */
				x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

				x[0], x[4] = invbutterfly(x[0], x[4], psi[0], twoQ, fourQ, Q, MRedConstant)
				x[1], x[5] = invbutterfly(x[1], x[5], psi[0], twoQ, fourQ, Q, MRedConstant)
				x[2], x[6] = invbutterfly(x[2], x[6], psi[0], twoQ, fourQ, Q, MRedConstant)
				x[3], x[7] = invbutterfly(x[3], x[7], psi[0], twoQ, fourQ, Q, MRedConstant)
				x[8], x[12] = invbutterfly(x[8], x[12], psi[1], twoQ, fourQ, Q, MRedConstant)
				x[9], x[13] = invbutterfly(x[9], x[13], psi[1], twoQ, fourQ, Q, MRedConstant)
				x[10], x[14] = invbutterfly(x[10], x[14], psi[1], twoQ, fourQ, Q, MRedConstant)
				x[11], x[15] = invbutterfly(x[11], x[15], psi[1], twoQ, fourQ, Q, MRedConstant)

				j1 = j1 + (t << 2)
			}

		} else {

			for i := m; i < h+m; i = i + 4 {

				/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(roots)%4 */
				psi := (*[4]uint64)(unsafe.Pointer(&roots[i]))
				/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%16 */
				x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

				x[0], x[2] = invbutterfly(x[0], x[2], psi[0], twoQ, fourQ, Q, MRedConstant)
				x[1], x[3] = invbutterfly(x[1], x[3], psi[0], twoQ, fourQ, Q, MRedConstant)
				x[4], x[6] = invbutterfly(x[4], x[6], psi[1], twoQ, fourQ, Q, MRedConstant)
				x[5], x[7] = invbutterfly(x[5], x[7], psi[1], twoQ, fourQ, Q, MRedConstant)
				x[8], x[10] = invbutterfly(x[8], x[10], psi[2], twoQ, fourQ, Q, MRedConstant)
				x[9], x[11] = invbutterfly(x[9], x[11], psi[2], twoQ, fourQ, Q, MRedConstant)
				x[12], x[14] = invbutterfly(x[12], x[14], psi[3], twoQ, fourQ, Q, MRedConstant)
				x[13], x[15] = invbutterfly(x[13], x[15], psi[3], twoQ, fourQ, Q, MRedConstant)

				j1 = j1 + (t << 3)
			}
		}

		t <<= 1
	}

	F = roots[1]

	for jx, jy := 1, N-8; jx < (N>>1)-7; jx, jy = jx+8, jy-8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		xout := (*[8]uint64)(unsafe.Pointer(&p2[jx]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		yout := (*[8]uint64)(unsafe.Pointer(&p2[jy]))

		xout[0], yout[7] = xout[0]+twoQ-MRedLazy(yout[7], F, Q, MRedConstant), yout[7]+twoQ-MRedLazy(xout[0], F, Q, MRedConstant)
		xout[1], yout[6] = xout[1]+twoQ-MRedLazy(yout[6], F, Q, MRedConstant), yout[6]+twoQ-MRedLazy(xout[1], F, Q, MRedConstant)
		xout[2], yout[5] = xout[2]+twoQ-MRedLazy(yout[5], F, Q, MRedConstant), yout[5]+twoQ-MRedLazy(xout[2], F, Q, MRedConstant)
		xout[3], yout[4] = xout[3]+twoQ-MRedLazy(yout[4], F, Q, MRedConstant), yout[4]+twoQ-MRedLazy(xout[3], F, Q, MRedConstant)
		xout[4], yout[3] = xout[4]+twoQ-MRedLazy(yout[3], F, Q, MRedConstant), yout[3]+twoQ-MRedLazy(xout[4], F, Q, MRedConstant)
		xout[5], yout[2] = xout[5]+twoQ-MRedLazy(yout[2], F, Q, MRedConstant), yout[2]+twoQ-MRedLazy(xout[5], F, Q, MRedConstant)
		xout[6], yout[1] = xout[6]+twoQ-MRedLazy(yout[1], F, Q, MRedConstant), yout[1]+twoQ-MRedLazy(xout[6], F, Q, MRedConstant)
		xout[7], yout[0] = xout[7]+twoQ-MRedLazy(yout[0], F, Q, MRedConstant), yout[0]+twoQ-MRedLazy(xout[7], F, Q, MRedConstant)
	}

	j := (N >> 1) - 7
	/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
	xout := (*[7]uint64)(unsafe.Pointer(&p2[j]))
	/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
	yout := (*[7]uint64)(unsafe.Pointer(&p2[N-j-6]))

	xout[0], yout[6] = xout[0]+twoQ-MRedLazy(yout[6], F, Q, MRedConstant), yout[6]+twoQ-MRedLazy(xout[0], F, Q, MRedConstant)
	xout[1], yout[5] = xout[1]+twoQ-MRedLazy(yout[5], F, Q, MRedConstant), yout[5]+twoQ-MRedLazy(xout[1], F, Q, MRedConstant)
	xout[2], yout[4] = xout[2]+twoQ-MRedLazy(yout[4], F, Q, MRedConstant), yout[4]+twoQ-MRedLazy(xout[2], F, Q, MRedConstant)
	xout[3], yout[3] = xout[3]+twoQ-MRedLazy(yout[3], F, Q, MRedConstant), yout[3]+twoQ-MRedLazy(xout[3], F, Q, MRedConstant)
	xout[4], yout[2] = xout[4]+twoQ-MRedLazy(yout[2], F, Q, MRedConstant), yout[2]+twoQ-MRedLazy(xout[4], F, Q, MRedConstant)
	xout[5], yout[1] = xout[5]+twoQ-MRedLazy(yout[1], F, Q, MRedConstant), yout[1]+twoQ-MRedLazy(xout[5], F, Q, MRedConstant)
	xout[6], yout[0] = xout[6]+twoQ-MRedLazy(yout[0], F, Q, MRedConstant), yout[0]+twoQ-MRedLazy(xout[6], F, Q, MRedConstant)

	p2[N>>1] = p2[N>>1] + twoQ - MRedLazy(p2[N>>1], F, Q, MRedConstant)
	p2[0] = CRed(p2[0]<<1, Q)
}
