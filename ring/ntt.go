package ring

import (
	"math/bits"
	"unsafe"
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
	n                                          int
	nInv, modulus, mredConstant                uint64
	bredConstants, rootsForward, rootsBackward []uint64
}

// NumberTheoreticTransformerStandard computes the standard nega-cyclic NTT in the ring Z[X]/(X^N+1).
type NumberTheoreticTransformerStandard struct {
	numberTheoreticTransformerBase
}

func NewNumberTheoreticTransformerStandard(r *SubRing, n int) NumberTheoreticTransformer {
	return NumberTheoreticTransformerStandard{
		numberTheoreticTransformerBase: numberTheoreticTransformerBase{
			n: n, nInv: r.NInv, modulus: r.Modulus, mredConstant: r.MRedConstant, bredConstants: r.BRedConstant, rootsForward: r.RootsForward, rootsBackward: r.RootsBackward,
		},
	}
}

// Forward writes the forward NTT in Z[X]/(X^N+1) of p1 on p2.
func (rntt NumberTheoreticTransformerStandard) Forward(p1, p2 []uint64) {
	nttStandard(p1, p2, rntt.n, rntt.modulus, rntt.mredConstant, rntt.bredConstants, rntt.rootsForward)
}

// ForwardLazy writes the forward NTT in Z[X]/(X^N+1) of p1 on p2.
// Returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerStandard) ForwardLazy(p1, p2 []uint64) {
	nttStandardLazy(p1, p2, rntt.n, rntt.modulus, rntt.mredConstant, rntt.rootsBackward)
}

// Backward writes the backward NTT in Z[X]/(X^N+1) of p1 on p2.
func (rntt NumberTheoreticTransformerStandard) Backward(p1, p2 []uint64) {
	inttStandard(p1, p2, rntt.n, rntt.modulus, rntt.nInv, rntt.mredConstant, rntt.rootsBackward)
}

// BackwardLazy writes the backward NTT in Z[X]/(X^N+1) p1 on p2.
// Returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerStandard) BackwardLazy(p1, p2 []uint64) {
	inttStandardLazy(p1, p2, rntt.n, rntt.modulus, rntt.nInv, rntt.mredConstant, rntt.rootsBackward)
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
			n: n, nInv: r.NInv, modulus: r.Modulus, mredConstant: r.MRedConstant, bredConstants: r.BRedConstant, rootsForward: r.RootsForward, rootsBackward: r.RootsBackward,
		},
	}
}

// Forward writes the forward NTT in Z[X+X^-1]/(X^2N+1) of p1 on p2.
func (rntt NumberTheoreticTransformerConjugateInvariant) Forward(p1, p2 []uint64) {
	nttConjugateInvariant(p1, p2, rntt.n, rntt.modulus, rntt.mredConstant, rntt.bredConstants, rntt.rootsForward)
}

// ForwardLazy writes the forward NTT in Z[X+X^-1]/(X^2N+1) of p1 on p2.
// Returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerConjugateInvariant) ForwardLazy(p1, p2 []uint64) {
	nttConjugateInvariantLazy(p1, p2, rntt.n, rntt.modulus, rntt.mredConstant, rntt.rootsForward)
}

// Backward writes the backward NTT in Z[X+X^-1]/(X^2N+1) of p1 on p2.
func (rntt NumberTheoreticTransformerConjugateInvariant) Backward(p1, p2 []uint64) {
	inttConjugateInvariant(p1, p2, rntt.n, rntt.modulus, rntt.nInv, rntt.mredConstant, rntt.rootsBackward)
}

// BackwardLazy writes the backward NTT in Z[X+X^-1]/(X^2N+1) of p1 on p2.
// Returns values in the range [0, 2q-1].
func (rntt NumberTheoreticTransformerConjugateInvariant) BackwardLazy(p1, p2 []uint64) {
	inttConjugateInvariantLazy(p1, p2, rntt.n, rntt.modulus, rntt.nInv, rntt.mredConstant, rntt.mredConstant, rntt.rootsBackward)
}

// NTT evaluates p2 = NTT(P1).
func (r *Ring) NTT(p1, p2 *Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.NTT(p1.Coeffs[i], p2.Coeffs[i])
	}
}

// NTTLazy evaluates p2 = NTT(p1) with p2 in [0, 2*modulus-1].
func (r *Ring) NTTLazy(p1, p2 *Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.NTTLazy(p1.Coeffs[i], p2.Coeffs[i])
	}
}

// INTT evaluates p2 = INTT(p1).
func (r *Ring) INTT(p1, p2 *Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.INTT(p1.Coeffs[i], p2.Coeffs[i])
	}
}

// INTTLazy evaluates p2 = INTT(p1) with p2 in [0, 2*modulus-1].
func (r *Ring) INTTLazy(p1, p2 *Poly) {
	for i, s := range r.SubRings[:r.level+1] {
		s.INTTLazy(p1.Coeffs[i], p2.Coeffs[i])
	}
}

// NttSparseAndMontgomery takes the polynomial polIn Z[Y] outside of the NTT domain to the polynomial Z[X] in the NTT domain where Y = X^(gap).
// This method is used to accelerate the NTT of polynomials that encode sparse plaintexts.
func (r *Ring) NttSparseAndMontgomery(logSlots int, montgomery bool, pol *Poly) {

	if 1<<logSlots == r.NthRoot()>>2 {
		r.NTT(pol, pol)
		if montgomery {
			r.MForm(pol, pol)
		}
	} else {

		var n int
		var ntt func(p1, p2 []uint64, N int, Q, QInv uint64, BRedConstant, nttPsi []uint64)
		switch r.Type() {
		case Standard:
			n = 2 << logSlots
			ntt = nttStandard
		case ConjugateInvariant:
			n = 1 << logSlots
			ntt = nttConjugateInvariant
		}

		N := r.N()
		gap := N / n
		for i, s := range r.SubRings[:r.Level()+1] {

			coeffs := pol.Coeffs[i]

			// Hack!
			// NTT in dimension n but with roots of N
			ntt(coeffs[:n], coeffs[:n], n, s.Modulus, s.MRedConstant, s.BRedConstant, s.RootsForward)

			if montgomery {
				s.MForm(coeffs[:n], coeffs[:n])
			}

			// Maps NTT in dimension n to NTT in dimension N
			for j := n - 1; j >= 0; j-- {
				c := coeffs[j]
				for w := 0; w < gap; w++ {
					coeffs[j*gap+w] = c
				}
			}
		}
	}
}

// butterfly computes X, Y = U + V*Psi, U - V*Psi mod Q.
func butterfly(U, V, Psi, twoQ, fourQ, Q, Qinv uint64) (uint64, uint64) {
	if U >= fourQ {
		U -= fourQ
	}
	V = MRedLazy(V, Psi, Q, Qinv)
	return U + V, U + twoQ - V
}

// invbutterfly computes X, Y = U + V, (U - V) * Psi mod Q.
func invbutterfly(U, V, Psi, twoQ, fourQ, Q, Qinv uint64) (X, Y uint64) {
	X = U + V
	if X >= twoQ {
		X -= twoQ
	}
	Y = MRedLazy(U+fourQ-V, Psi, Q, Qinv) // At the moment it is not possible to use MRedLazy if Q > 61 bits
	return
}

// nttStandard computes the nttStandard on the input coefficients using the input parameters.
func nttStandard(p1, p2 []uint64, N int, Q, QInv uint64, BRedConstant, nttPsi []uint64) {
	nttStandardLazy(p1, p2, N, Q, QInv, nttPsi)
	reducevec(p2, p2, Q, BRedConstant)
}

// inttStandard evalues p2 = inttStandard(p1) in the given SubRing.
func inttStandard(p1, p2 []uint64, N int, NInv, Q, MRedConstant uint64, nttPsiInv []uint64) {
	iNTTCore(p1, p2, N, Q, MRedConstant, nttPsiInv)
	mulscalarmontgomeryvec(p2, NInv, p2, Q, MRedConstant)
}

// inttStandardLazy evalues p2 = INTT(p1) in the given SubRing with p2 in [0, 2*modulus-1].
func inttStandardLazy(p1, p2 []uint64, N int, NInv, Q, MRedConstant uint64, nttPsiInv []uint64) {
	iNTTCore(p1, p2, N, Q, MRedConstant, nttPsiInv)
	mulscalarmontgomerylazyvec(p2, NInv, p2, Q, MRedConstant)
}

// nttConjugateInvariant evaluates p2 = NTT(p1) in the sub-ring Z[X + X^-1]/(X^2N +1) of Z[X]/(X^2N+1).
func nttConjugateInvariant(p1, p2 []uint64, N int, Q, MRedConstant uint64, BRedConstant, nttPsi []uint64) {
	nttConjugateInvariantLazy(p1, p2, N, Q, MRedConstant, nttPsi)
	reducevec(p2, p2, Q, BRedConstant)
}

// inttConjugateInvariant evaluates p2 = INTT(p1) in the closed sub-ring Z[X + X^-1]/(X^2N +1) of Z[X]/(X^2N+1).
func inttConjugateInvariant(p1, p2 []uint64, N int, NInv, Q, MRedConstant uint64, nttPsiInv []uint64) {
	iNTTConjugateInvariantCore(p1, p2, N, Q, MRedConstant, nttPsiInv)
	mulscalarmontgomeryvec(p2, NInv, p2, Q, MRedConstant)
}

// inttConjugateInvariantLazy evaluates p2 = INTT(p1) in the closed sub-ring Z[X + X^-1]/(X^2N +1) of Z[X]/(X^2N+1) with p2 in the range [0, 2*modulus-1].
func inttConjugateInvariantLazy(p1, p2 []uint64, N int, NInv, Q, QInv, MRedConstant uint64, nttPsiInv []uint64) {
	iNTTConjugateInvariantCore(p1, p2, N, Q, QInv, nttPsiInv)
	mulscalarmontgomerylazyvec(p2, NInv, p2, Q, MRedConstant)
}

// nttStandardLazy computes the NTT on the input coefficients using the input parameters with output values in the range [0, 2*modulus-1].
func nttStandardLazy(p1, p2 []uint64, N int, Q, QInv uint64, nttPsi []uint64) {

	var j1, j2, t int
	var F, V uint64

	fourQ := 4 * Q
	twoQ := 2 * Q

	// Copy the result of the first round of butterflies on p2 with approximate reduction
	t = N >> 1
	F = nttPsi[1]

	for jx, jy := 0, t; jx <= t-1; jx, jy = jx+8, jy+8 {

		xin := (*[8]uint64)(unsafe.Pointer(&p1[jx]))
		yin := (*[8]uint64)(unsafe.Pointer(&p1[jy]))

		xout := (*[8]uint64)(unsafe.Pointer(&p2[jx]))
		yout := (*[8]uint64)(unsafe.Pointer(&p2[jy]))

		V = MRedLazy(yin[0], F, Q, QInv)
		xout[0], yout[0] = xin[0]+V, xin[0]+twoQ-V

		V = MRedLazy(yin[1], F, Q, QInv)
		xout[1], yout[1] = xin[1]+V, xin[1]+twoQ-V

		V = MRedLazy(yin[2], F, Q, QInv)
		xout[2], yout[2] = xin[2]+V, xin[2]+twoQ-V

		V = MRedLazy(yin[3], F, Q, QInv)
		xout[3], yout[3] = xin[3]+V, xin[3]+twoQ-V

		V = MRedLazy(yin[4], F, Q, QInv)
		xout[4], yout[4] = xin[4]+V, xin[4]+twoQ-V

		V = MRedLazy(yin[5], F, Q, QInv)
		xout[5], yout[5] = xin[5]+V, xin[5]+twoQ-V

		V = MRedLazy(yin[6], F, Q, QInv)
		xout[6], yout[6] = xin[6]+V, xin[6]+twoQ-V

		V = MRedLazy(yin[7], F, Q, QInv)
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

				j2 = j1 + t - 1

				F = nttPsi[m+i]

				if reduce {

					for jx, jy := j1, j1+t; jx <= j2; jx, jy = jx+8, jy+8 {

						x := (*[8]uint64)(unsafe.Pointer(&p2[jx]))
						y := (*[8]uint64)(unsafe.Pointer(&p2[jy]))

						x[0], y[0] = butterfly(x[0], y[0], F, twoQ, fourQ, Q, QInv)
						x[1], y[1] = butterfly(x[1], y[1], F, twoQ, fourQ, Q, QInv)
						x[2], y[2] = butterfly(x[2], y[2], F, twoQ, fourQ, Q, QInv)
						x[3], y[3] = butterfly(x[3], y[3], F, twoQ, fourQ, Q, QInv)
						x[4], y[4] = butterfly(x[4], y[4], F, twoQ, fourQ, Q, QInv)
						x[5], y[5] = butterfly(x[5], y[5], F, twoQ, fourQ, Q, QInv)
						x[6], y[6] = butterfly(x[6], y[6], F, twoQ, fourQ, Q, QInv)
						x[7], y[7] = butterfly(x[7], y[7], F, twoQ, fourQ, Q, QInv)
					}

				} else {

					for jx, jy := j1, j1+t; jx <= j2; jx, jy = jx+8, jy+8 {

						x := (*[8]uint64)(unsafe.Pointer(&p2[jx]))
						y := (*[8]uint64)(unsafe.Pointer(&p2[jy]))

						V = MRedLazy(y[0], F, Q, QInv)
						x[0], y[0] = x[0]+V, x[0]+twoQ-V

						V = MRedLazy(y[1], F, Q, QInv)
						x[1], y[1] = x[1]+V, x[1]+twoQ-V

						V = MRedLazy(y[2], F, Q, QInv)
						x[2], y[2] = x[2]+V, x[2]+twoQ-V

						V = MRedLazy(y[3], F, Q, QInv)
						x[3], y[3] = x[3]+V, x[3]+twoQ-V

						V = MRedLazy(y[4], F, Q, QInv)
						x[4], y[4] = x[4]+V, x[4]+twoQ-V

						V = MRedLazy(y[5], F, Q, QInv)
						x[5], y[5] = x[5]+V, x[5]+twoQ-V

						V = MRedLazy(y[6], F, Q, QInv)
						x[6], y[6] = x[6]+V, x[6]+twoQ-V

						V = MRedLazy(y[7], F, Q, QInv)
						x[7], y[7] = x[7]+V, x[7]+twoQ-V
					}
				}
			}

		} else if t == 4 {

			if reduce {

				for i, j1 := m, 0; i < 2*m; i, j1 = i+2, j1+4*t {

					psi := (*[2]uint64)(unsafe.Pointer(&nttPsi[i]))
					x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

					x[0], x[4] = butterfly(x[0], x[4], psi[0], twoQ, fourQ, Q, QInv)
					x[1], x[5] = butterfly(x[1], x[5], psi[0], twoQ, fourQ, Q, QInv)
					x[2], x[6] = butterfly(x[2], x[6], psi[0], twoQ, fourQ, Q, QInv)
					x[3], x[7] = butterfly(x[3], x[7], psi[0], twoQ, fourQ, Q, QInv)
					x[8], x[12] = butterfly(x[8], x[12], psi[1], twoQ, fourQ, Q, QInv)
					x[9], x[13] = butterfly(x[9], x[13], psi[1], twoQ, fourQ, Q, QInv)
					x[10], x[14] = butterfly(x[10], x[14], psi[1], twoQ, fourQ, Q, QInv)
					x[11], x[15] = butterfly(x[11], x[15], psi[1], twoQ, fourQ, Q, QInv)

				}
			} else {

				for i, j1 := m, 0; i < 2*m; i, j1 = i+2, j1+4*t {

					psi := (*[2]uint64)(unsafe.Pointer(&nttPsi[i]))
					x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

					V = MRedLazy(x[4], psi[0], Q, QInv)
					x[0], x[4] = x[0]+V, x[0]+twoQ-V

					V = MRedLazy(x[5], psi[0], Q, QInv)
					x[1], x[5] = x[1]+V, x[1]+twoQ-V

					V = MRedLazy(x[6], psi[0], Q, QInv)
					x[2], x[6] = x[2]+V, x[2]+twoQ-V

					V = MRedLazy(x[7], psi[0], Q, QInv)
					x[3], x[7] = x[3]+V, x[3]+twoQ-V

					V = MRedLazy(x[12], psi[1], Q, QInv)
					x[8], x[12] = x[8]+V, x[8]+twoQ-V

					V = MRedLazy(x[13], psi[1], Q, QInv)
					x[9], x[13] = x[9]+V, x[9]+twoQ-V

					V = MRedLazy(x[14], psi[1], Q, QInv)
					x[10], x[14] = x[10]+V, x[10]+twoQ-V

					V = MRedLazy(x[15], psi[1], Q, QInv)
					x[11], x[15] = x[11]+V, x[11]+twoQ-V

				}

			}

		} else if t == 2 {

			if reduce {

				for i, j1 := m, 0; i < 2*m; i, j1 = i+4, j1+8*t {

					psi := (*[4]uint64)(unsafe.Pointer(&nttPsi[i]))
					x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

					x[0], x[2] = butterfly(x[0], x[2], psi[0], twoQ, fourQ, Q, QInv)
					x[1], x[3] = butterfly(x[1], x[3], psi[0], twoQ, fourQ, Q, QInv)
					x[4], x[6] = butterfly(x[4], x[6], psi[1], twoQ, fourQ, Q, QInv)
					x[5], x[7] = butterfly(x[5], x[7], psi[1], twoQ, fourQ, Q, QInv)
					x[8], x[10] = butterfly(x[8], x[10], psi[2], twoQ, fourQ, Q, QInv)
					x[9], x[11] = butterfly(x[9], x[11], psi[2], twoQ, fourQ, Q, QInv)
					x[12], x[14] = butterfly(x[12], x[14], psi[3], twoQ, fourQ, Q, QInv)
					x[13], x[15] = butterfly(x[13], x[15], psi[3], twoQ, fourQ, Q, QInv)
				}
			} else {

				for i, j1 := m, 0; i < 2*m; i, j1 = i+4, j1+8*t {

					psi := (*[4]uint64)(unsafe.Pointer(&nttPsi[i]))
					x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

					V = MRedLazy(x[2], psi[0], Q, QInv)
					x[0], x[2] = x[0]+V, x[0]+twoQ-V

					V = MRedLazy(x[3], psi[0], Q, QInv)
					x[1], x[3] = x[1]+V, x[1]+twoQ-V

					V = MRedLazy(x[6], psi[1], Q, QInv)
					x[4], x[6] = x[4]+V, x[4]+twoQ-V

					V = MRedLazy(x[7], psi[1], Q, QInv)
					x[5], x[7] = x[5]+V, x[5]+twoQ-V

					V = MRedLazy(x[10], psi[2], Q, QInv)
					x[8], x[10] = x[8]+V, x[8]+twoQ-V

					V = MRedLazy(x[11], psi[2], Q, QInv)
					x[9], x[11] = x[9]+V, x[9]+twoQ-V

					V = MRedLazy(x[14], psi[3], Q, QInv)
					x[12], x[14] = x[12]+V, x[12]+twoQ-V

					V = MRedLazy(x[15], psi[3], Q, QInv)
					x[13], x[15] = x[13]+V, x[13]+twoQ-V
				}
			}

		} else {

			for i, j1 := m, 0; i < 2*m; i, j1 = i+8, j1+16 {

				psi := (*[8]uint64)(unsafe.Pointer(&nttPsi[i]))
				x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

				x[0], x[1] = butterfly(x[0], x[1], psi[0], twoQ, fourQ, Q, QInv)
				x[2], x[3] = butterfly(x[2], x[3], psi[1], twoQ, fourQ, Q, QInv)
				x[4], x[5] = butterfly(x[4], x[5], psi[2], twoQ, fourQ, Q, QInv)
				x[6], x[7] = butterfly(x[6], x[7], psi[3], twoQ, fourQ, Q, QInv)
				x[8], x[9] = butterfly(x[8], x[9], psi[4], twoQ, fourQ, Q, QInv)
				x[10], x[11] = butterfly(x[10], x[11], psi[5], twoQ, fourQ, Q, QInv)
				x[12], x[13] = butterfly(x[12], x[13], psi[6], twoQ, fourQ, Q, QInv)
				x[14], x[15] = butterfly(x[14], x[15], psi[7], twoQ, fourQ, Q, QInv)
			}

			/*
				for i := uint64(0); i < m; i = i + 8 {

					psi := (*[8]uint64)(unsafe.Pointer(&nttPsi[m+i]))
					x := (*[16]uint64)(unsafe.Pointer(&p2[2*i]))

					V = MRedLazy(x[1], psi[0], Q, QInv)
					x[0], x[1] = x[0]+V, x[0]+twoQ-V

					V = MRedLazy(x[3], psi[1], Q, QInv)
					x[2], x[3] = x[2]+V, x[2]+twoQ-V

					V = MRedLazy(x[5], psi[2], Q, QInv)
					x[4], x[5] = x[4]+V, x[4]+twoQ-V

					V = MRedLazy(x[7], psi[3], Q, QInv)
					x[6], x[7] = x[6]+V, x[6]+twoQ-V

					V = MRedLazy(x[9], psi[4], Q, QInv)
					x[8], x[9] = x[8]+V, x[8]+twoQ-V

					V = MRedLazy(x[11], psi[5], Q, QInv)
					x[10], x[11] = x[10]+V, x[10]+twoQ-V

					V = MRedLazy(x[13], psi[6], Q, QInv)
					x[12], x[13] = x[12]+V, x[12]+twoQ-V

					V = MRedLazy(x[15], psi[7], Q, QInv)
					x[14], x[15] = x[14]+V, x[14]+twoQ-V
				}
			*/
		}
	}
}

func iNTTCore(p1, p2 []uint64, N int, Q, QInv uint64, nttPsiInv []uint64) {

	var h, t int
	var F uint64

	// Copy the result of the first round of butterflies on p2 with approximate reduction
	t = 1
	h = N >> 1
	twoQ := Q << 1
	fourQ := Q << 2

	for i, j := h, 0; i < 2*h; i, j = i+8, j+16 {

		psi := (*[8]uint64)(unsafe.Pointer(&nttPsiInv[i]))
		xin := (*[16]uint64)(unsafe.Pointer(&p1[j]))
		xout := (*[16]uint64)(unsafe.Pointer(&p2[j]))

		xout[0], xout[1] = invbutterfly(xin[0], xin[1], psi[0], twoQ, fourQ, Q, QInv)
		xout[2], xout[3] = invbutterfly(xin[2], xin[3], psi[1], twoQ, fourQ, Q, QInv)
		xout[4], xout[5] = invbutterfly(xin[4], xin[5], psi[2], twoQ, fourQ, Q, QInv)
		xout[6], xout[7] = invbutterfly(xin[6], xin[7], psi[3], twoQ, fourQ, Q, QInv)
		xout[8], xout[9] = invbutterfly(xin[8], xin[9], psi[4], twoQ, fourQ, Q, QInv)
		xout[10], xout[11] = invbutterfly(xin[10], xin[11], psi[5], twoQ, fourQ, Q, QInv)
		xout[12], xout[13] = invbutterfly(xin[12], xin[13], psi[6], twoQ, fourQ, Q, QInv)
		xout[14], xout[15] = invbutterfly(xin[14], xin[15], psi[7], twoQ, fourQ, Q, QInv)
	}

	// Continue the rest of the second to the n-1 butterflies on p2 with approximate reduction
	t <<= 1
	for m := N >> 1; m > 1; m >>= 1 {

		h = m >> 1

		if t >= 8 {

			for i, j1, j2 := 0, 0, t-1; i < h; i, j1, j2 = i+1, j1+2*t, j2+2*t {

				F = nttPsiInv[h+i]

				for jx, jy := j1, j1+t; jx <= j2; jx, jy = jx+8, jy+8 {

					x := (*[8]uint64)(unsafe.Pointer(&p2[jx]))
					y := (*[8]uint64)(unsafe.Pointer(&p2[jy]))

					x[0], y[0] = invbutterfly(x[0], y[0], F, twoQ, fourQ, Q, QInv)
					x[1], y[1] = invbutterfly(x[1], y[1], F, twoQ, fourQ, Q, QInv)
					x[2], y[2] = invbutterfly(x[2], y[2], F, twoQ, fourQ, Q, QInv)
					x[3], y[3] = invbutterfly(x[3], y[3], F, twoQ, fourQ, Q, QInv)
					x[4], y[4] = invbutterfly(x[4], y[4], F, twoQ, fourQ, Q, QInv)
					x[5], y[5] = invbutterfly(x[5], y[5], F, twoQ, fourQ, Q, QInv)
					x[6], y[6] = invbutterfly(x[6], y[6], F, twoQ, fourQ, Q, QInv)
					x[7], y[7] = invbutterfly(x[7], y[7], F, twoQ, fourQ, Q, QInv)
				}
			}

		} else if t == 4 {

			for i, j1 := h, 0; i < 2*h; i, j1 = i+2, j1+4*t {

				psi := (*[2]uint64)(unsafe.Pointer(&nttPsiInv[i]))
				x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

				x[0], x[4] = invbutterfly(x[0], x[4], psi[0], twoQ, fourQ, Q, QInv)
				x[1], x[5] = invbutterfly(x[1], x[5], psi[0], twoQ, fourQ, Q, QInv)
				x[2], x[6] = invbutterfly(x[2], x[6], psi[0], twoQ, fourQ, Q, QInv)
				x[3], x[7] = invbutterfly(x[3], x[7], psi[0], twoQ, fourQ, Q, QInv)
				x[8], x[12] = invbutterfly(x[8], x[12], psi[1], twoQ, fourQ, Q, QInv)
				x[9], x[13] = invbutterfly(x[9], x[13], psi[1], twoQ, fourQ, Q, QInv)
				x[10], x[14] = invbutterfly(x[10], x[14], psi[1], twoQ, fourQ, Q, QInv)
				x[11], x[15] = invbutterfly(x[11], x[15], psi[1], twoQ, fourQ, Q, QInv)
			}

		} else {

			for i, j1 := h, 0; i < 2*h; i, j1 = i+4, j1+8*t {

				psi := (*[4]uint64)(unsafe.Pointer(&nttPsiInv[i]))
				x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

				x[0], x[2] = invbutterfly(x[0], x[2], psi[0], twoQ, fourQ, Q, QInv)
				x[1], x[3] = invbutterfly(x[1], x[3], psi[0], twoQ, fourQ, Q, QInv)
				x[4], x[6] = invbutterfly(x[4], x[6], psi[1], twoQ, fourQ, Q, QInv)
				x[5], x[7] = invbutterfly(x[5], x[7], psi[1], twoQ, fourQ, Q, QInv)
				x[8], x[10] = invbutterfly(x[8], x[10], psi[2], twoQ, fourQ, Q, QInv)
				x[9], x[11] = invbutterfly(x[9], x[11], psi[2], twoQ, fourQ, Q, QInv)
				x[12], x[14] = invbutterfly(x[12], x[14], psi[3], twoQ, fourQ, Q, QInv)
				x[13], x[15] = invbutterfly(x[13], x[15], psi[3], twoQ, fourQ, Q, QInv)
			}
		}

		t <<= 1
	}
}

// nttConjugateInvariantLazy evaluates p2 = NTT(p1) in the sub-ring Z[X + X^-1]/(X^2N +1) of Z[X]/(X^2N+1) with p2 [0, 2*modulus-1].
func nttConjugateInvariantLazy(p1, p2 []uint64, N int, Q, QInv uint64, nttPsi []uint64) {

	var t, h int
	var F, V uint64
	var reduce bool

	fourQ := 4 * Q
	twoQ := 2 * Q

	// Copy the result of the first round of butterflies on p2 with approximate reduction
	t = N
	F = nttPsi[1]

	for jx, jy := 1, N-8; jx < (N>>1)-7; jx, jy = jx+8, jy-8 {

		xin := (*[8]uint64)(unsafe.Pointer(&p1[jx]))
		yin := (*[8]uint64)(unsafe.Pointer(&p1[jy]))

		xout := (*[8]uint64)(unsafe.Pointer(&p2[jx]))
		yout := (*[8]uint64)(unsafe.Pointer(&p2[jy]))

		xout[0], yout[7] = xin[0]+twoQ-MRedLazy(yin[7], F, Q, QInv), yin[7]+twoQ-MRedLazy(xin[0], F, Q, QInv)
		xout[1], yout[6] = xin[1]+twoQ-MRedLazy(yin[6], F, Q, QInv), yin[6]+twoQ-MRedLazy(xin[1], F, Q, QInv)
		xout[2], yout[5] = xin[2]+twoQ-MRedLazy(yin[5], F, Q, QInv), yin[5]+twoQ-MRedLazy(xin[2], F, Q, QInv)
		xout[3], yout[4] = xin[3]+twoQ-MRedLazy(yin[4], F, Q, QInv), yin[4]+twoQ-MRedLazy(xin[3], F, Q, QInv)
		xout[4], yout[3] = xin[4]+twoQ-MRedLazy(yin[3], F, Q, QInv), yin[3]+twoQ-MRedLazy(xin[4], F, Q, QInv)
		xout[5], yout[2] = xin[5]+twoQ-MRedLazy(yin[2], F, Q, QInv), yin[2]+twoQ-MRedLazy(xin[5], F, Q, QInv)
		xout[6], yout[1] = xin[6]+twoQ-MRedLazy(yin[1], F, Q, QInv), yin[1]+twoQ-MRedLazy(xin[6], F, Q, QInv)
		xout[7], yout[0] = xin[7]+twoQ-MRedLazy(yin[0], F, Q, QInv), yin[0]+twoQ-MRedLazy(xin[7], F, Q, QInv)
	}

	j := (N >> 1) - 7
	xin := (*[7]uint64)(unsafe.Pointer(&p1[j]))
	yin := (*[7]uint64)(unsafe.Pointer(&p1[N-j-6]))
	xout := (*[7]uint64)(unsafe.Pointer(&p2[j]))
	yout := (*[7]uint64)(unsafe.Pointer(&p2[N-j-6]))

	xout[0], yout[6] = xin[0]+twoQ-MRedLazy(yin[6], F, Q, QInv), yin[6]+twoQ-MRedLazy(xin[0], F, Q, QInv)
	xout[1], yout[5] = xin[1]+twoQ-MRedLazy(yin[5], F, Q, QInv), yin[5]+twoQ-MRedLazy(xin[1], F, Q, QInv)
	xout[2], yout[4] = xin[2]+twoQ-MRedLazy(yin[4], F, Q, QInv), yin[4]+twoQ-MRedLazy(xin[2], F, Q, QInv)
	xout[3], yout[3] = xin[3]+twoQ-MRedLazy(yin[3], F, Q, QInv), yin[3]+twoQ-MRedLazy(xin[3], F, Q, QInv)
	xout[4], yout[2] = xin[4]+twoQ-MRedLazy(yin[2], F, Q, QInv), yin[2]+twoQ-MRedLazy(xin[4], F, Q, QInv)
	xout[5], yout[1] = xin[5]+twoQ-MRedLazy(yin[1], F, Q, QInv), yin[1]+twoQ-MRedLazy(xin[5], F, Q, QInv)
	xout[6], yout[0] = xin[6]+twoQ-MRedLazy(yin[0], F, Q, QInv), yin[0]+twoQ-MRedLazy(xin[6], F, Q, QInv)

	p2[N>>1] = p1[N>>1] + twoQ - MRedLazy(p1[N>>1], F, Q, QInv)
	p2[0] = p1[0]

	// Continue the rest of the second to the n-1 butterflies on p2 with approximate reduction
	for m := 2; m < 2*N; m <<= 1 {

		reduce = (bits.Len64(uint64(m))&1 == 1)

		t >>= 1
		h = m >> 1

		if t >= 8 {

			for i, j1, j2 := 0, 0, t-1; i < h; i, j1, j2 = i+1, j1+2*t, j2+2*t {

				F = nttPsi[m+i]

				if reduce {

					for jx, jy := j1, j1+t; jx <= j2; jx, jy = jx+8, jy+8 {

						x := (*[8]uint64)(unsafe.Pointer(&p2[jx]))
						y := (*[8]uint64)(unsafe.Pointer(&p2[jy]))

						x[0], y[0] = butterfly(x[0], y[0], F, twoQ, fourQ, Q, QInv)
						x[1], y[1] = butterfly(x[1], y[1], F, twoQ, fourQ, Q, QInv)
						x[2], y[2] = butterfly(x[2], y[2], F, twoQ, fourQ, Q, QInv)
						x[3], y[3] = butterfly(x[3], y[3], F, twoQ, fourQ, Q, QInv)
						x[4], y[4] = butterfly(x[4], y[4], F, twoQ, fourQ, Q, QInv)
						x[5], y[5] = butterfly(x[5], y[5], F, twoQ, fourQ, Q, QInv)
						x[6], y[6] = butterfly(x[6], y[6], F, twoQ, fourQ, Q, QInv)
						x[7], y[7] = butterfly(x[7], y[7], F, twoQ, fourQ, Q, QInv)
					}

				} else {

					for jx, jy := j1, j1+t; jx <= j2; jx, jy = jx+8, jy+8 {

						x := (*[8]uint64)(unsafe.Pointer(&p2[jx]))
						y := (*[8]uint64)(unsafe.Pointer(&p2[jy]))

						V = MRedLazy(y[0], F, Q, QInv)
						x[0], y[0] = x[0]+V, x[0]+twoQ-V

						V = MRedLazy(y[1], F, Q, QInv)
						x[1], y[1] = x[1]+V, x[1]+twoQ-V

						V = MRedLazy(y[2], F, Q, QInv)
						x[2], y[2] = x[2]+V, x[2]+twoQ-V

						V = MRedLazy(y[3], F, Q, QInv)
						x[3], y[3] = x[3]+V, x[3]+twoQ-V

						V = MRedLazy(y[4], F, Q, QInv)
						x[4], y[4] = x[4]+V, x[4]+twoQ-V

						V = MRedLazy(y[5], F, Q, QInv)
						x[5], y[5] = x[5]+V, x[5]+twoQ-V

						V = MRedLazy(y[6], F, Q, QInv)
						x[6], y[6] = x[6]+V, x[6]+twoQ-V

						V = MRedLazy(y[7], F, Q, QInv)
						x[7], y[7] = x[7]+V, x[7]+twoQ-V
					}
				}
			}

		} else if t == 4 {

			if reduce {

				for i, j1 := m, 0; i < h+m; i, j1 = i+2, j1+4*t {

					psi := (*[2]uint64)(unsafe.Pointer(&nttPsi[i]))
					x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

					x[0], x[4] = butterfly(x[0], x[4], psi[0], twoQ, fourQ, Q, QInv)
					x[1], x[5] = butterfly(x[1], x[5], psi[0], twoQ, fourQ, Q, QInv)
					x[2], x[6] = butterfly(x[2], x[6], psi[0], twoQ, fourQ, Q, QInv)
					x[3], x[7] = butterfly(x[3], x[7], psi[0], twoQ, fourQ, Q, QInv)
					x[8], x[12] = butterfly(x[8], x[12], psi[1], twoQ, fourQ, Q, QInv)
					x[9], x[13] = butterfly(x[9], x[13], psi[1], twoQ, fourQ, Q, QInv)
					x[10], x[14] = butterfly(x[10], x[14], psi[1], twoQ, fourQ, Q, QInv)
					x[11], x[15] = butterfly(x[11], x[15], psi[1], twoQ, fourQ, Q, QInv)

				}
			} else {

				for i, j1 := m, 0; i < h+m; i, j1 = i+2, j1+4*t {

					psi := (*[2]uint64)(unsafe.Pointer(&nttPsi[i]))
					x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

					V = MRedLazy(x[4], psi[0], Q, QInv)
					x[0], x[4] = x[0]+V, x[0]+twoQ-V

					V = MRedLazy(x[5], psi[0], Q, QInv)
					x[1], x[5] = x[1]+V, x[1]+twoQ-V

					V = MRedLazy(x[6], psi[0], Q, QInv)
					x[2], x[6] = x[2]+V, x[2]+twoQ-V

					V = MRedLazy(x[7], psi[0], Q, QInv)
					x[3], x[7] = x[3]+V, x[3]+twoQ-V

					V = MRedLazy(x[12], psi[1], Q, QInv)
					x[8], x[12] = x[8]+V, x[8]+twoQ-V

					V = MRedLazy(x[13], psi[1], Q, QInv)
					x[9], x[13] = x[9]+V, x[9]+twoQ-V

					V = MRedLazy(x[14], psi[1], Q, QInv)
					x[10], x[14] = x[10]+V, x[10]+twoQ-V

					V = MRedLazy(x[15], psi[1], Q, QInv)
					x[11], x[15] = x[11]+V, x[11]+twoQ-V

				}
			}

		} else if t == 2 {

			if reduce {

				for i, j1 := m, 0; i < h+m; i, j1 = i+4, j1+8*t {

					psi := (*[4]uint64)(unsafe.Pointer(&nttPsi[i]))
					x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

					x[0], x[2] = butterfly(x[0], x[2], psi[0], twoQ, fourQ, Q, QInv)
					x[1], x[3] = butterfly(x[1], x[3], psi[0], twoQ, fourQ, Q, QInv)
					x[4], x[6] = butterfly(x[4], x[6], psi[1], twoQ, fourQ, Q, QInv)
					x[5], x[7] = butterfly(x[5], x[7], psi[1], twoQ, fourQ, Q, QInv)
					x[8], x[10] = butterfly(x[8], x[10], psi[2], twoQ, fourQ, Q, QInv)
					x[9], x[11] = butterfly(x[9], x[11], psi[2], twoQ, fourQ, Q, QInv)
					x[12], x[14] = butterfly(x[12], x[14], psi[3], twoQ, fourQ, Q, QInv)
					x[13], x[15] = butterfly(x[13], x[15], psi[3], twoQ, fourQ, Q, QInv)
				}
			} else {

				for i, j1 := m, 0; i < h+m; i, j1 = i+4, j1+8*t {

					psi := (*[4]uint64)(unsafe.Pointer(&nttPsi[i]))
					x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

					V = MRedLazy(x[2], psi[0], Q, QInv)
					x[0], x[2] = x[0]+V, x[0]+twoQ-V

					V = MRedLazy(x[3], psi[0], Q, QInv)
					x[1], x[3] = x[1]+V, x[1]+twoQ-V

					V = MRedLazy(x[6], psi[1], Q, QInv)
					x[4], x[6] = x[4]+V, x[4]+twoQ-V

					V = MRedLazy(x[7], psi[1], Q, QInv)
					x[5], x[7] = x[5]+V, x[5]+twoQ-V

					V = MRedLazy(x[10], psi[2], Q, QInv)
					x[8], x[10] = x[8]+V, x[8]+twoQ-V

					V = MRedLazy(x[11], psi[2], Q, QInv)
					x[9], x[11] = x[9]+V, x[9]+twoQ-V

					V = MRedLazy(x[14], psi[3], Q, QInv)
					x[12], x[14] = x[12]+V, x[12]+twoQ-V

					V = MRedLazy(x[15], psi[3], Q, QInv)
					x[13], x[15] = x[13]+V, x[13]+twoQ-V
				}
			}

		} else {

			if reduce {

				for i, j1 := m, 0; i < h+m; i, j1 = i+8, j1+16 {

					psi := (*[8]uint64)(unsafe.Pointer(&nttPsi[i]))
					x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

					x[0], x[1] = butterfly(x[0], x[1], psi[0], twoQ, fourQ, Q, QInv)
					x[2], x[3] = butterfly(x[2], x[3], psi[1], twoQ, fourQ, Q, QInv)
					x[4], x[5] = butterfly(x[4], x[5], psi[2], twoQ, fourQ, Q, QInv)
					x[6], x[7] = butterfly(x[6], x[7], psi[3], twoQ, fourQ, Q, QInv)
					x[8], x[9] = butterfly(x[8], x[9], psi[4], twoQ, fourQ, Q, QInv)
					x[10], x[11] = butterfly(x[10], x[11], psi[5], twoQ, fourQ, Q, QInv)
					x[12], x[13] = butterfly(x[12], x[13], psi[6], twoQ, fourQ, Q, QInv)
					x[14], x[15] = butterfly(x[14], x[15], psi[7], twoQ, fourQ, Q, QInv)
				}
			} else {

				for i, j1 := m, 0; i < h+m; i, j1 = i+8, j1+16 {

					psi := (*[8]uint64)(unsafe.Pointer(&nttPsi[i]))
					x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

					V = MRedLazy(x[1], psi[0], Q, QInv)
					x[0], x[1] = x[0]+V, x[0]+twoQ-V

					V = MRedLazy(x[3], psi[1], Q, QInv)
					x[2], x[3] = x[2]+V, x[2]+twoQ-V

					V = MRedLazy(x[5], psi[2], Q, QInv)
					x[4], x[5] = x[4]+V, x[4]+twoQ-V

					V = MRedLazy(x[7], psi[3], Q, QInv)
					x[6], x[7] = x[6]+V, x[6]+twoQ-V

					V = MRedLazy(x[9], psi[4], Q, QInv)
					x[8], x[9] = x[8]+V, x[8]+twoQ-V

					V = MRedLazy(x[11], psi[5], Q, QInv)
					x[10], x[11] = x[10]+V, x[10]+twoQ-V

					V = MRedLazy(x[13], psi[6], Q, QInv)
					x[12], x[13] = x[12]+V, x[12]+twoQ-V

					V = MRedLazy(x[15], psi[7], Q, QInv)
					x[14], x[15] = x[14]+V, x[14]+twoQ-V
				}
			}
		}
	}
}

func iNTTConjugateInvariantCore(p1, p2 []uint64, N int, Q, QInv uint64, nttPsiInv []uint64) {

	var j1, j2, h, t int
	var F uint64

	// Copy the result of the first round of butterflies on p2 with approximate reduction
	t = 1
	h = N >> 1
	twoQ := Q << 1
	fourQ := Q << 2

	for i, j := N, 0; i < h+N; i, j = i+8, j+16 {

		psi := (*[8]uint64)(unsafe.Pointer(&nttPsiInv[i]))
		xin := (*[16]uint64)(unsafe.Pointer(&p1[j]))
		xout := (*[16]uint64)(unsafe.Pointer(&p2[j]))

		xout[0], xout[1] = invbutterfly(xin[0], xin[1], psi[0], twoQ, fourQ, Q, QInv)
		xout[2], xout[3] = invbutterfly(xin[2], xin[3], psi[1], twoQ, fourQ, Q, QInv)
		xout[4], xout[5] = invbutterfly(xin[4], xin[5], psi[2], twoQ, fourQ, Q, QInv)
		xout[6], xout[7] = invbutterfly(xin[6], xin[7], psi[3], twoQ, fourQ, Q, QInv)
		xout[8], xout[9] = invbutterfly(xin[8], xin[9], psi[4], twoQ, fourQ, Q, QInv)
		xout[10], xout[11] = invbutterfly(xin[10], xin[11], psi[5], twoQ, fourQ, Q, QInv)
		xout[12], xout[13] = invbutterfly(xin[12], xin[13], psi[6], twoQ, fourQ, Q, QInv)
		xout[14], xout[15] = invbutterfly(xin[14], xin[15], psi[7], twoQ, fourQ, Q, QInv)
	}

	// Continue the rest of the second to the n-1 butterflies on p2 with approximate reduction
	t <<= 1
	for m := N >> 1; m > 1; m >>= 1 {

		j1 = 0
		h = m >> 1

		if t >= 8 {

			for i := 0; i < h; i++ {

				j2 = j1 + t - 1

				F = nttPsiInv[m+i]

				for jx, jy := j1, j1+t; jx <= j2; jx, jy = jx+8, jy+8 {

					x := (*[8]uint64)(unsafe.Pointer(&p2[jx]))
					y := (*[8]uint64)(unsafe.Pointer(&p2[jy]))

					x[0], y[0] = invbutterfly(x[0], y[0], F, twoQ, fourQ, Q, QInv)
					x[1], y[1] = invbutterfly(x[1], y[1], F, twoQ, fourQ, Q, QInv)
					x[2], y[2] = invbutterfly(x[2], y[2], F, twoQ, fourQ, Q, QInv)
					x[3], y[3] = invbutterfly(x[3], y[3], F, twoQ, fourQ, Q, QInv)
					x[4], y[4] = invbutterfly(x[4], y[4], F, twoQ, fourQ, Q, QInv)
					x[5], y[5] = invbutterfly(x[5], y[5], F, twoQ, fourQ, Q, QInv)
					x[6], y[6] = invbutterfly(x[6], y[6], F, twoQ, fourQ, Q, QInv)
					x[7], y[7] = invbutterfly(x[7], y[7], F, twoQ, fourQ, Q, QInv)
				}

				j1 = j1 + (t << 1)
			}

		} else if t == 4 {

			for i := m; i < h+m; i = i + 2 {

				psi := (*[2]uint64)(unsafe.Pointer(&nttPsiInv[i]))
				x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

				x[0], x[4] = invbutterfly(x[0], x[4], psi[0], twoQ, fourQ, Q, QInv)
				x[1], x[5] = invbutterfly(x[1], x[5], psi[0], twoQ, fourQ, Q, QInv)
				x[2], x[6] = invbutterfly(x[2], x[6], psi[0], twoQ, fourQ, Q, QInv)
				x[3], x[7] = invbutterfly(x[3], x[7], psi[0], twoQ, fourQ, Q, QInv)
				x[8], x[12] = invbutterfly(x[8], x[12], psi[1], twoQ, fourQ, Q, QInv)
				x[9], x[13] = invbutterfly(x[9], x[13], psi[1], twoQ, fourQ, Q, QInv)
				x[10], x[14] = invbutterfly(x[10], x[14], psi[1], twoQ, fourQ, Q, QInv)
				x[11], x[15] = invbutterfly(x[11], x[15], psi[1], twoQ, fourQ, Q, QInv)

				j1 = j1 + (t << 2)
			}

		} else {

			for i := m; i < h+m; i = i + 4 {

				psi := (*[4]uint64)(unsafe.Pointer(&nttPsiInv[i]))
				x := (*[16]uint64)(unsafe.Pointer(&p2[j1]))

				x[0], x[2] = invbutterfly(x[0], x[2], psi[0], twoQ, fourQ, Q, QInv)
				x[1], x[3] = invbutterfly(x[1], x[3], psi[0], twoQ, fourQ, Q, QInv)
				x[4], x[6] = invbutterfly(x[4], x[6], psi[1], twoQ, fourQ, Q, QInv)
				x[5], x[7] = invbutterfly(x[5], x[7], psi[1], twoQ, fourQ, Q, QInv)
				x[8], x[10] = invbutterfly(x[8], x[10], psi[2], twoQ, fourQ, Q, QInv)
				x[9], x[11] = invbutterfly(x[9], x[11], psi[2], twoQ, fourQ, Q, QInv)
				x[12], x[14] = invbutterfly(x[12], x[14], psi[3], twoQ, fourQ, Q, QInv)
				x[13], x[15] = invbutterfly(x[13], x[15], psi[3], twoQ, fourQ, Q, QInv)

				j1 = j1 + (t << 3)
			}
		}

		t <<= 1
	}

	F = nttPsiInv[1]

	for jx, jy := 1, N-8; jx < (N>>1)-7; jx, jy = jx+8, jy-8 {

		xout := (*[8]uint64)(unsafe.Pointer(&p2[jx]))
		yout := (*[8]uint64)(unsafe.Pointer(&p2[jy]))

		xout[0], yout[7] = xout[0]+twoQ-MRedLazy(yout[7], F, Q, QInv), yout[7]+twoQ-MRedLazy(xout[0], F, Q, QInv)
		xout[1], yout[6] = xout[1]+twoQ-MRedLazy(yout[6], F, Q, QInv), yout[6]+twoQ-MRedLazy(xout[1], F, Q, QInv)
		xout[2], yout[5] = xout[2]+twoQ-MRedLazy(yout[5], F, Q, QInv), yout[5]+twoQ-MRedLazy(xout[2], F, Q, QInv)
		xout[3], yout[4] = xout[3]+twoQ-MRedLazy(yout[4], F, Q, QInv), yout[4]+twoQ-MRedLazy(xout[3], F, Q, QInv)
		xout[4], yout[3] = xout[4]+twoQ-MRedLazy(yout[3], F, Q, QInv), yout[3]+twoQ-MRedLazy(xout[4], F, Q, QInv)
		xout[5], yout[2] = xout[5]+twoQ-MRedLazy(yout[2], F, Q, QInv), yout[2]+twoQ-MRedLazy(xout[5], F, Q, QInv)
		xout[6], yout[1] = xout[6]+twoQ-MRedLazy(yout[1], F, Q, QInv), yout[1]+twoQ-MRedLazy(xout[6], F, Q, QInv)
		xout[7], yout[0] = xout[7]+twoQ-MRedLazy(yout[0], F, Q, QInv), yout[0]+twoQ-MRedLazy(xout[7], F, Q, QInv)
	}

	j := (N >> 1) - 7
	xout := (*[7]uint64)(unsafe.Pointer(&p2[j]))
	yout := (*[7]uint64)(unsafe.Pointer(&p2[N-j-6]))

	xout[0], yout[6] = xout[0]+twoQ-MRedLazy(yout[6], F, Q, QInv), yout[6]+twoQ-MRedLazy(xout[0], F, Q, QInv)
	xout[1], yout[5] = xout[1]+twoQ-MRedLazy(yout[5], F, Q, QInv), yout[5]+twoQ-MRedLazy(xout[1], F, Q, QInv)
	xout[2], yout[4] = xout[2]+twoQ-MRedLazy(yout[4], F, Q, QInv), yout[4]+twoQ-MRedLazy(xout[2], F, Q, QInv)
	xout[3], yout[3] = xout[3]+twoQ-MRedLazy(yout[3], F, Q, QInv), yout[3]+twoQ-MRedLazy(xout[3], F, Q, QInv)
	xout[4], yout[2] = xout[4]+twoQ-MRedLazy(yout[2], F, Q, QInv), yout[2]+twoQ-MRedLazy(xout[4], F, Q, QInv)
	xout[5], yout[1] = xout[5]+twoQ-MRedLazy(yout[1], F, Q, QInv), yout[1]+twoQ-MRedLazy(xout[5], F, Q, QInv)
	xout[6], yout[0] = xout[6]+twoQ-MRedLazy(yout[0], F, Q, QInv), yout[0]+twoQ-MRedLazy(xout[6], F, Q, QInv)

	p2[N>>1] = p2[N>>1] + twoQ - MRedLazy(p2[N>>1], F, Q, QInv)
	p2[0] = CRed(p2[0]<<1, Q)
}
