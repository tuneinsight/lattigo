package ring

import (
	"math/bits"
	"unsafe"
)

// NTT computes the NTT of p1 and returns the result on p2.
func (r *Ring) NTT(p1, p2 *Poly) {
	r.NumberTheoreticTransformer.Forward(r, p1, p2)
}

// NTTLvl computes the NTT of p1 and returns the result on p2.
// The value level defines the number of moduli of the input polynomials.
func (r *Ring) NTTLvl(level int, p1, p2 *Poly) {
	r.NumberTheoreticTransformer.ForwardLvl(r, level, p1, p2)
}

// NTTLazy computes the NTT of p1 and returns the result on p2.
// Output values are in the range [0, 2q-1]
func (r *Ring) NTTLazy(p1, p2 *Poly) {
	r.NumberTheoreticTransformer.ForwardLazy(r, p1, p2)
}

// NTTLazyLvl computes the NTT of p1 and returns the result on p2.
// The value level defines the number of moduli of the input polynomials.
// Output values are in the range [0, 2q-1]
func (r *Ring) NTTLazyLvl(level int, p1, p2 *Poly) {
	r.NumberTheoreticTransformer.ForwardLazyLvl(r, level, p1, p2)
}

// NTTSingle computes the NTT of p1 and returns the result on p2.
// The level-th moduli of the ring NTT params are used.
func (r *Ring) NTTSingle(level int, p1, p2 []uint64) {
	r.NumberTheoreticTransformer.ForwardVec(r, level, p1, p2)
}

// NTTSingleLazy computes the NTT of p1 and returns the result on p2.
// The level-th moduli of the ring NTT params are used.
// Output values are in the range [0, 2q-1]
func (r *Ring) NTTSingleLazy(level int, p1, p2 []uint64) {
	r.NumberTheoreticTransformer.ForwardLazyVec(r, level, p1, p2)
}

// InvNTT computes the inverse-NTT of p1 and returns the result on p2.
func (r *Ring) InvNTT(p1, p2 *Poly) {
	r.NumberTheoreticTransformer.Backward(r, p1, p2)
}

// InvNTTLvl computes the inverse-NTT of p1 and returns the result on p2.
// The value level defines the number of moduli of the input polynomials.
func (r *Ring) InvNTTLvl(level int, p1, p2 *Poly) {
	r.NumberTheoreticTransformer.BackwardLvl(r, level, p1, p2)
}

// InvNTTLazy computes the inverse-NTT of p1 and returns the result on p2.
// Output values are in the range [0, 2q-1]
func (r *Ring) InvNTTLazy(p1, p2 *Poly) {
	r.NumberTheoreticTransformer.BackwardLazy(r, p1, p2)
}

// InvNTTLazyLvl computes the inverse-NTT of p1 and returns the result on p2.
// The value level defines the number of moduli of the input polynomials.
// Output values are in the range [0, 2q-1]
func (r *Ring) InvNTTLazyLvl(level int, p1, p2 *Poly) {
	r.NumberTheoreticTransformer.BackwardLazyLvl(r, level, p1, p2)
}

// InvNTTSingle computes the InvNTT of p1 and returns the result on p2.
// The level-th moduli of the ring InvNTT params are used.
// Only computes the InvNTT for the i-th level.
func (r *Ring) InvNTTSingle(level int, p1, p2 []uint64) {
	r.NumberTheoreticTransformer.BackwardVec(r, level, p1, p2)
}

// InvNTTSingleLazy computes the InvNTT of p1 and returns the result on p2.
// The level-th moduli of the ring InvNTT params are used.
// Output values are in the range [0, 2q-1]
func (r *Ring) InvNTTSingleLazy(level int, p1, p2 []uint64) {
	r.NumberTheoreticTransformer.BackwardLazyVec(r, level, p1, p2)
}

// butterfly computes X, Y = U + V*Psi, U - V*Psi mod Q.
func butterfly(U, V, Psi, twoQ, fourQ, Q, Qinv uint64) (uint64, uint64) {
	if U >= fourQ {
		U -= fourQ
	}
	V = MRedConstant(V, Psi, Q, Qinv)
	return U + V, U + twoQ - V
}

// invbutterfly computes X, Y = U + V, (U - V) * Psi mod Q.
func invbutterfly(U, V, Psi, twoQ, fourQ, Q, Qinv uint64) (X, Y uint64) {
	X = U + V
	if X >= twoQ {
		X -= twoQ
	}
	Y = MRedConstant(U+fourQ-V, Psi, Q, Qinv) // At the moment it is not possible to use MRedConstant if Q > 61 bits
	return
}

// NTT computes the NTT on the input coefficients using the input parameters.
func NTT(coeffsIn, coeffsOut []uint64, N int, nttPsi []uint64, Q, mredParams uint64, bredParams []uint64) {
	NTTLazy(coeffsIn, coeffsOut, N, nttPsi, Q, mredParams, bredParams)
	ReduceVec(coeffsOut, coeffsOut, Q, bredParams)
}

// NTTLazy computes the NTT on the input coefficients using the input parameters with output values in the range [0, 2q-1].
func NTTLazy(coeffsIn, coeffsOut []uint64, N int, nttPsi []uint64, Q, QInv uint64, bredParams []uint64) {
	var j1, j2, t int
	var F, V uint64

	fourQ := 4 * Q
	twoQ := 2 * Q

	// Copy the result of the first round of butterflies on p2 with approximate reduction
	t = N >> 1
	F = nttPsi[1]

	for jx, jy := 0, t; jx <= t-1; jx, jy = jx+8, jy+8 {

		xin := (*[8]uint64)(unsafe.Pointer(&coeffsIn[jx]))
		yin := (*[8]uint64)(unsafe.Pointer(&coeffsIn[jy]))

		xout := (*[8]uint64)(unsafe.Pointer(&coeffsOut[jx]))
		yout := (*[8]uint64)(unsafe.Pointer(&coeffsOut[jy]))

		V = MRedConstant(yin[0], F, Q, QInv)
		xout[0], yout[0] = xin[0]+V, xin[0]+twoQ-V

		V = MRedConstant(yin[1], F, Q, QInv)
		xout[1], yout[1] = xin[1]+V, xin[1]+twoQ-V

		V = MRedConstant(yin[2], F, Q, QInv)
		xout[2], yout[2] = xin[2]+V, xin[2]+twoQ-V

		V = MRedConstant(yin[3], F, Q, QInv)
		xout[3], yout[3] = xin[3]+V, xin[3]+twoQ-V

		V = MRedConstant(yin[4], F, Q, QInv)
		xout[4], yout[4] = xin[4]+V, xin[4]+twoQ-V

		V = MRedConstant(yin[5], F, Q, QInv)
		xout[5], yout[5] = xin[5]+V, xin[5]+twoQ-V

		V = MRedConstant(yin[6], F, Q, QInv)
		xout[6], yout[6] = xin[6]+V, xin[6]+twoQ-V

		V = MRedConstant(yin[7], F, Q, QInv)
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

						x := (*[8]uint64)(unsafe.Pointer(&coeffsOut[jx]))
						y := (*[8]uint64)(unsafe.Pointer(&coeffsOut[jy]))

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

						x := (*[8]uint64)(unsafe.Pointer(&coeffsOut[jx]))
						y := (*[8]uint64)(unsafe.Pointer(&coeffsOut[jy]))

						V = MRedConstant(y[0], F, Q, QInv)
						x[0], y[0] = x[0]+V, x[0]+twoQ-V

						V = MRedConstant(y[1], F, Q, QInv)
						x[1], y[1] = x[1]+V, x[1]+twoQ-V

						V = MRedConstant(y[2], F, Q, QInv)
						x[2], y[2] = x[2]+V, x[2]+twoQ-V

						V = MRedConstant(y[3], F, Q, QInv)
						x[3], y[3] = x[3]+V, x[3]+twoQ-V

						V = MRedConstant(y[4], F, Q, QInv)
						x[4], y[4] = x[4]+V, x[4]+twoQ-V

						V = MRedConstant(y[5], F, Q, QInv)
						x[5], y[5] = x[5]+V, x[5]+twoQ-V

						V = MRedConstant(y[6], F, Q, QInv)
						x[6], y[6] = x[6]+V, x[6]+twoQ-V

						V = MRedConstant(y[7], F, Q, QInv)
						x[7], y[7] = x[7]+V, x[7]+twoQ-V
					}
				}
			}

		} else if t == 4 {

			if reduce {

				for i, j1 := m, 0; i < 2*m; i, j1 = i+2, j1+4*t {

					psi := (*[2]uint64)(unsafe.Pointer(&nttPsi[i]))
					x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

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
					x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

					V = MRedConstant(x[4], psi[0], Q, QInv)
					x[0], x[4] = x[0]+V, x[0]+twoQ-V

					V = MRedConstant(x[5], psi[0], Q, QInv)
					x[1], x[5] = x[1]+V, x[1]+twoQ-V

					V = MRedConstant(x[6], psi[0], Q, QInv)
					x[2], x[6] = x[2]+V, x[2]+twoQ-V

					V = MRedConstant(x[7], psi[0], Q, QInv)
					x[3], x[7] = x[3]+V, x[3]+twoQ-V

					V = MRedConstant(x[12], psi[1], Q, QInv)
					x[8], x[12] = x[8]+V, x[8]+twoQ-V

					V = MRedConstant(x[13], psi[1], Q, QInv)
					x[9], x[13] = x[9]+V, x[9]+twoQ-V

					V = MRedConstant(x[14], psi[1], Q, QInv)
					x[10], x[14] = x[10]+V, x[10]+twoQ-V

					V = MRedConstant(x[15], psi[1], Q, QInv)
					x[11], x[15] = x[11]+V, x[11]+twoQ-V

				}

			}

		} else if t == 2 {

			if reduce {

				for i, j1 := m, 0; i < 2*m; i, j1 = i+4, j1+8*t {

					psi := (*[4]uint64)(unsafe.Pointer(&nttPsi[i]))
					x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

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
					x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

					V = MRedConstant(x[2], psi[0], Q, QInv)
					x[0], x[2] = x[0]+V, x[0]+twoQ-V

					V = MRedConstant(x[3], psi[0], Q, QInv)
					x[1], x[3] = x[1]+V, x[1]+twoQ-V

					V = MRedConstant(x[6], psi[1], Q, QInv)
					x[4], x[6] = x[4]+V, x[4]+twoQ-V

					V = MRedConstant(x[7], psi[1], Q, QInv)
					x[5], x[7] = x[5]+V, x[5]+twoQ-V

					V = MRedConstant(x[10], psi[2], Q, QInv)
					x[8], x[10] = x[8]+V, x[8]+twoQ-V

					V = MRedConstant(x[11], psi[2], Q, QInv)
					x[9], x[11] = x[9]+V, x[9]+twoQ-V

					V = MRedConstant(x[14], psi[3], Q, QInv)
					x[12], x[14] = x[12]+V, x[12]+twoQ-V

					V = MRedConstant(x[15], psi[3], Q, QInv)
					x[13], x[15] = x[13]+V, x[13]+twoQ-V
				}
			}

		} else {

			for i, j1 := m, 0; i < 2*m; i, j1 = i+8, j1+16 {

				psi := (*[8]uint64)(unsafe.Pointer(&nttPsi[i]))
				x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

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
					x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[2*i]))

					V = MRedConstant(x[1], psi[0], Q, QInv)
					x[0], x[1] = x[0]+V, x[0]+twoQ-V

					V = MRedConstant(x[3], psi[1], Q, QInv)
					x[2], x[3] = x[2]+V, x[2]+twoQ-V

					V = MRedConstant(x[5], psi[2], Q, QInv)
					x[4], x[5] = x[4]+V, x[4]+twoQ-V

					V = MRedConstant(x[7], psi[3], Q, QInv)
					x[6], x[7] = x[6]+V, x[6]+twoQ-V

					V = MRedConstant(x[9], psi[4], Q, QInv)
					x[8], x[9] = x[8]+V, x[8]+twoQ-V

					V = MRedConstant(x[11], psi[5], Q, QInv)
					x[10], x[11] = x[10]+V, x[10]+twoQ-V

					V = MRedConstant(x[13], psi[6], Q, QInv)
					x[12], x[13] = x[12]+V, x[12]+twoQ-V

					V = MRedConstant(x[15], psi[7], Q, QInv)
					x[14], x[15] = x[14]+V, x[14]+twoQ-V
				}
			*/
		}
	}
}

func invNTTCore(coeffsIn, coeffsOut []uint64, N int, nttPsiInv []uint64, Q, QInv uint64) {
	var h, t int
	var F uint64

	// Copy the result of the first round of butterflies on p2 with approximate reduction
	t = 1
	h = N >> 1
	twoQ := Q << 1
	fourQ := Q << 2

	for i, j := h, 0; i < 2*h; i, j = i+8, j+16 {

		psi := (*[8]uint64)(unsafe.Pointer(&nttPsiInv[i]))
		xin := (*[16]uint64)(unsafe.Pointer(&coeffsIn[j]))
		xout := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j]))

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

					x := (*[8]uint64)(unsafe.Pointer(&coeffsOut[jx]))
					y := (*[8]uint64)(unsafe.Pointer(&coeffsOut[jy]))

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
				x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

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
				x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

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

// InvNTT computes the InvNTT transformation on the input coefficients using the input parameters.
func InvNTT(coeffsIn, coeffsOut []uint64, N int, nttPsiInv []uint64, nttNInv, Q, QInv uint64) {
	invNTTCore(coeffsIn, coeffsOut, N, nttPsiInv, Q, QInv)
	MulScalarMontgomeryVec(coeffsOut, coeffsOut, nttNInv, Q, QInv)
}

// InvNTTLazy computes the InvNTT transformation on the input coefficients using the input parameters with output values in the range [0, 2q-1].
func InvNTTLazy(coeffsIn, coeffsOut []uint64, N int, nttPsiInv []uint64, nttNInv, Q, QInv uint64) {
	invNTTCore(coeffsIn, coeffsOut, N, nttPsiInv, Q, QInv)
	MulScalarMontgomeryConstantVec(coeffsOut, coeffsOut, nttNInv, Q, QInv)
}

// NTTConjugateInvariant computes the NTT in the closed sub-ring Z[X + X^-1]/(X^2N +1) of Z[X]/(X^2N+1).
func NTTConjugateInvariant(coeffsIn, coeffsOut []uint64, N int, nttPsi []uint64, Q, QInv uint64, bredParams []uint64) {
	NTTConjugateInvariantLazy(coeffsIn, coeffsOut, N, nttPsi, Q, QInv, bredParams)
	ReduceVec(coeffsOut, coeffsOut, Q, bredParams)
}

// NTTConjugateInvariantLazy computes the NTT in the closed sub-ring Z[X + X^-1]/(X^2N +1) of Z[X]/(X^2N+1) with output values in the range [0, 2q-1].
func NTTConjugateInvariantLazy(coeffsIn, coeffsOut []uint64, N int, nttPsi []uint64, Q, QInv uint64, bredParams []uint64) {
	var t, h int
	var F, V uint64
	var reduce bool

	fourQ := 4 * Q
	twoQ := 2 * Q

	// Copy the result of the first round of butterflies on p2 with approximate reduction
	t = N
	F = nttPsi[1]

	for jx, jy := 1, N-8; jx < (N>>1)-7; jx, jy = jx+8, jy-8 {

		xin := (*[8]uint64)(unsafe.Pointer(&coeffsIn[jx]))
		yin := (*[8]uint64)(unsafe.Pointer(&coeffsIn[jy]))

		xout := (*[8]uint64)(unsafe.Pointer(&coeffsOut[jx]))
		yout := (*[8]uint64)(unsafe.Pointer(&coeffsOut[jy]))

		xout[0], yout[7] = xin[0]+twoQ-MRedConstant(yin[7], F, Q, QInv), yin[7]+twoQ-MRedConstant(xin[0], F, Q, QInv)
		xout[1], yout[6] = xin[1]+twoQ-MRedConstant(yin[6], F, Q, QInv), yin[6]+twoQ-MRedConstant(xin[1], F, Q, QInv)
		xout[2], yout[5] = xin[2]+twoQ-MRedConstant(yin[5], F, Q, QInv), yin[5]+twoQ-MRedConstant(xin[2], F, Q, QInv)
		xout[3], yout[4] = xin[3]+twoQ-MRedConstant(yin[4], F, Q, QInv), yin[4]+twoQ-MRedConstant(xin[3], F, Q, QInv)
		xout[4], yout[3] = xin[4]+twoQ-MRedConstant(yin[3], F, Q, QInv), yin[3]+twoQ-MRedConstant(xin[4], F, Q, QInv)
		xout[5], yout[2] = xin[5]+twoQ-MRedConstant(yin[2], F, Q, QInv), yin[2]+twoQ-MRedConstant(xin[5], F, Q, QInv)
		xout[6], yout[1] = xin[6]+twoQ-MRedConstant(yin[1], F, Q, QInv), yin[1]+twoQ-MRedConstant(xin[6], F, Q, QInv)
		xout[7], yout[0] = xin[7]+twoQ-MRedConstant(yin[0], F, Q, QInv), yin[0]+twoQ-MRedConstant(xin[7], F, Q, QInv)
	}

	j := (N >> 1) - 7
	xin := (*[7]uint64)(unsafe.Pointer(&coeffsIn[j]))
	yin := (*[7]uint64)(unsafe.Pointer(&coeffsIn[N-j-6]))
	xout := (*[7]uint64)(unsafe.Pointer(&coeffsOut[j]))
	yout := (*[7]uint64)(unsafe.Pointer(&coeffsOut[N-j-6]))

	xout[0], yout[6] = xin[0]+twoQ-MRedConstant(yin[6], F, Q, QInv), yin[6]+twoQ-MRedConstant(xin[0], F, Q, QInv)
	xout[1], yout[5] = xin[1]+twoQ-MRedConstant(yin[5], F, Q, QInv), yin[5]+twoQ-MRedConstant(xin[1], F, Q, QInv)
	xout[2], yout[4] = xin[2]+twoQ-MRedConstant(yin[4], F, Q, QInv), yin[4]+twoQ-MRedConstant(xin[2], F, Q, QInv)
	xout[3], yout[3] = xin[3]+twoQ-MRedConstant(yin[3], F, Q, QInv), yin[3]+twoQ-MRedConstant(xin[3], F, Q, QInv)
	xout[4], yout[2] = xin[4]+twoQ-MRedConstant(yin[2], F, Q, QInv), yin[2]+twoQ-MRedConstant(xin[4], F, Q, QInv)
	xout[5], yout[1] = xin[5]+twoQ-MRedConstant(yin[1], F, Q, QInv), yin[1]+twoQ-MRedConstant(xin[5], F, Q, QInv)
	xout[6], yout[0] = xin[6]+twoQ-MRedConstant(yin[0], F, Q, QInv), yin[0]+twoQ-MRedConstant(xin[6], F, Q, QInv)

	coeffsOut[N>>1] = coeffsIn[N>>1] + twoQ - MRedConstant(coeffsIn[N>>1], F, Q, QInv)
	coeffsOut[0] = coeffsIn[0]

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

						x := (*[8]uint64)(unsafe.Pointer(&coeffsOut[jx]))
						y := (*[8]uint64)(unsafe.Pointer(&coeffsOut[jy]))

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

						x := (*[8]uint64)(unsafe.Pointer(&coeffsOut[jx]))
						y := (*[8]uint64)(unsafe.Pointer(&coeffsOut[jy]))

						V = MRedConstant(y[0], F, Q, QInv)
						x[0], y[0] = x[0]+V, x[0]+twoQ-V

						V = MRedConstant(y[1], F, Q, QInv)
						x[1], y[1] = x[1]+V, x[1]+twoQ-V

						V = MRedConstant(y[2], F, Q, QInv)
						x[2], y[2] = x[2]+V, x[2]+twoQ-V

						V = MRedConstant(y[3], F, Q, QInv)
						x[3], y[3] = x[3]+V, x[3]+twoQ-V

						V = MRedConstant(y[4], F, Q, QInv)
						x[4], y[4] = x[4]+V, x[4]+twoQ-V

						V = MRedConstant(y[5], F, Q, QInv)
						x[5], y[5] = x[5]+V, x[5]+twoQ-V

						V = MRedConstant(y[6], F, Q, QInv)
						x[6], y[6] = x[6]+V, x[6]+twoQ-V

						V = MRedConstant(y[7], F, Q, QInv)
						x[7], y[7] = x[7]+V, x[7]+twoQ-V
					}
				}
			}

		} else if t == 4 {

			if reduce {

				for i, j1 := m, 0; i < h+m; i, j1 = i+2, j1+4*t {

					psi := (*[2]uint64)(unsafe.Pointer(&nttPsi[i]))
					x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

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
					x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

					V = MRedConstant(x[4], psi[0], Q, QInv)
					x[0], x[4] = x[0]+V, x[0]+twoQ-V

					V = MRedConstant(x[5], psi[0], Q, QInv)
					x[1], x[5] = x[1]+V, x[1]+twoQ-V

					V = MRedConstant(x[6], psi[0], Q, QInv)
					x[2], x[6] = x[2]+V, x[2]+twoQ-V

					V = MRedConstant(x[7], psi[0], Q, QInv)
					x[3], x[7] = x[3]+V, x[3]+twoQ-V

					V = MRedConstant(x[12], psi[1], Q, QInv)
					x[8], x[12] = x[8]+V, x[8]+twoQ-V

					V = MRedConstant(x[13], psi[1], Q, QInv)
					x[9], x[13] = x[9]+V, x[9]+twoQ-V

					V = MRedConstant(x[14], psi[1], Q, QInv)
					x[10], x[14] = x[10]+V, x[10]+twoQ-V

					V = MRedConstant(x[15], psi[1], Q, QInv)
					x[11], x[15] = x[11]+V, x[11]+twoQ-V

				}
			}

		} else if t == 2 {

			if reduce {

				for i, j1 := m, 0; i < h+m; i, j1 = i+4, j1+8*t {

					psi := (*[4]uint64)(unsafe.Pointer(&nttPsi[i]))
					x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

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
					x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

					V = MRedConstant(x[2], psi[0], Q, QInv)
					x[0], x[2] = x[0]+V, x[0]+twoQ-V

					V = MRedConstant(x[3], psi[0], Q, QInv)
					x[1], x[3] = x[1]+V, x[1]+twoQ-V

					V = MRedConstant(x[6], psi[1], Q, QInv)
					x[4], x[6] = x[4]+V, x[4]+twoQ-V

					V = MRedConstant(x[7], psi[1], Q, QInv)
					x[5], x[7] = x[5]+V, x[5]+twoQ-V

					V = MRedConstant(x[10], psi[2], Q, QInv)
					x[8], x[10] = x[8]+V, x[8]+twoQ-V

					V = MRedConstant(x[11], psi[2], Q, QInv)
					x[9], x[11] = x[9]+V, x[9]+twoQ-V

					V = MRedConstant(x[14], psi[3], Q, QInv)
					x[12], x[14] = x[12]+V, x[12]+twoQ-V

					V = MRedConstant(x[15], psi[3], Q, QInv)
					x[13], x[15] = x[13]+V, x[13]+twoQ-V
				}
			}

		} else {

			if reduce {

				for i, j1 := m, 0; i < h+m; i, j1 = i+8, j1+16 {

					psi := (*[8]uint64)(unsafe.Pointer(&nttPsi[i]))
					x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

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
					x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

					V = MRedConstant(x[1], psi[0], Q, QInv)
					x[0], x[1] = x[0]+V, x[0]+twoQ-V

					V = MRedConstant(x[3], psi[1], Q, QInv)
					x[2], x[3] = x[2]+V, x[2]+twoQ-V

					V = MRedConstant(x[5], psi[2], Q, QInv)
					x[4], x[5] = x[4]+V, x[4]+twoQ-V

					V = MRedConstant(x[7], psi[3], Q, QInv)
					x[6], x[7] = x[6]+V, x[6]+twoQ-V

					V = MRedConstant(x[9], psi[4], Q, QInv)
					x[8], x[9] = x[8]+V, x[8]+twoQ-V

					V = MRedConstant(x[11], psi[5], Q, QInv)
					x[10], x[11] = x[10]+V, x[10]+twoQ-V

					V = MRedConstant(x[13], psi[6], Q, QInv)
					x[12], x[13] = x[12]+V, x[12]+twoQ-V

					V = MRedConstant(x[15], psi[7], Q, QInv)
					x[14], x[15] = x[14]+V, x[14]+twoQ-V
				}
			}
		}
	}
}

// InvNTTConjugateInvariant computes the InvNTT in the closed sub-ring Z[X + X^-1]/(X^2N +1) of Z[X]/(X^2N+1).
func InvNTTConjugateInvariant(coeffsIn, coeffsOut []uint64, N int, nttPsiInv []uint64, nttNInv, Q, QInv uint64) {
	invNTTConjugateInvariantCore(coeffsIn, coeffsOut, N, nttPsiInv, Q, QInv)
	MulScalarMontgomeryVec(coeffsOut, coeffsOut, nttNInv, Q, QInv)
}

// InvNTTConjugateInvariantLazy computes the InvNTT in the closed sub-ring Z[X + X^-1]/(X^2N +1) of Z[X]/(X^2N+1) with output values in the range [0, 2q-1].
func InvNTTConjugateInvariantLazy(coeffsIn, coeffsOut []uint64, N int, nttPsiInv []uint64, nttNInv, Q, QInv uint64) {
	invNTTConjugateInvariantCore(coeffsIn, coeffsOut, N, nttPsiInv, Q, QInv)
	MulScalarMontgomeryConstantVec(coeffsOut, coeffsOut, nttNInv, Q, QInv)
}

func invNTTConjugateInvariantCore(coeffsIn, coeffsOut []uint64, N int, nttPsiInv []uint64, Q, QInv uint64) {
	var j1, j2, h, t int
	var F uint64

	// Copy the result of the first round of butterflies on p2 with approximate reduction
	t = 1
	h = N >> 1
	twoQ := Q << 1
	fourQ := Q << 2

	for i, j := N, 0; i < h+N; i, j = i+8, j+16 {

		psi := (*[8]uint64)(unsafe.Pointer(&nttPsiInv[i]))
		xin := (*[16]uint64)(unsafe.Pointer(&coeffsIn[j]))
		xout := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j]))

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

					x := (*[8]uint64)(unsafe.Pointer(&coeffsOut[jx]))
					y := (*[8]uint64)(unsafe.Pointer(&coeffsOut[jy]))

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
				x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

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
				x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

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

		xout := (*[8]uint64)(unsafe.Pointer(&coeffsOut[jx]))
		yout := (*[8]uint64)(unsafe.Pointer(&coeffsOut[jy]))

		xout[0], yout[7] = xout[0]+twoQ-MRedConstant(yout[7], F, Q, QInv), yout[7]+twoQ-MRedConstant(xout[0], F, Q, QInv)
		xout[1], yout[6] = xout[1]+twoQ-MRedConstant(yout[6], F, Q, QInv), yout[6]+twoQ-MRedConstant(xout[1], F, Q, QInv)
		xout[2], yout[5] = xout[2]+twoQ-MRedConstant(yout[5], F, Q, QInv), yout[5]+twoQ-MRedConstant(xout[2], F, Q, QInv)
		xout[3], yout[4] = xout[3]+twoQ-MRedConstant(yout[4], F, Q, QInv), yout[4]+twoQ-MRedConstant(xout[3], F, Q, QInv)
		xout[4], yout[3] = xout[4]+twoQ-MRedConstant(yout[3], F, Q, QInv), yout[3]+twoQ-MRedConstant(xout[4], F, Q, QInv)
		xout[5], yout[2] = xout[5]+twoQ-MRedConstant(yout[2], F, Q, QInv), yout[2]+twoQ-MRedConstant(xout[5], F, Q, QInv)
		xout[6], yout[1] = xout[6]+twoQ-MRedConstant(yout[1], F, Q, QInv), yout[1]+twoQ-MRedConstant(xout[6], F, Q, QInv)
		xout[7], yout[0] = xout[7]+twoQ-MRedConstant(yout[0], F, Q, QInv), yout[0]+twoQ-MRedConstant(xout[7], F, Q, QInv)
	}

	j := (N >> 1) - 7
	xout := (*[7]uint64)(unsafe.Pointer(&coeffsOut[j]))
	yout := (*[7]uint64)(unsafe.Pointer(&coeffsOut[N-j-6]))

	xout[0], yout[6] = xout[0]+twoQ-MRedConstant(yout[6], F, Q, QInv), yout[6]+twoQ-MRedConstant(xout[0], F, Q, QInv)
	xout[1], yout[5] = xout[1]+twoQ-MRedConstant(yout[5], F, Q, QInv), yout[5]+twoQ-MRedConstant(xout[1], F, Q, QInv)
	xout[2], yout[4] = xout[2]+twoQ-MRedConstant(yout[4], F, Q, QInv), yout[4]+twoQ-MRedConstant(xout[2], F, Q, QInv)
	xout[3], yout[3] = xout[3]+twoQ-MRedConstant(yout[3], F, Q, QInv), yout[3]+twoQ-MRedConstant(xout[3], F, Q, QInv)
	xout[4], yout[2] = xout[4]+twoQ-MRedConstant(yout[2], F, Q, QInv), yout[2]+twoQ-MRedConstant(xout[4], F, Q, QInv)
	xout[5], yout[1] = xout[5]+twoQ-MRedConstant(yout[1], F, Q, QInv), yout[1]+twoQ-MRedConstant(xout[5], F, Q, QInv)
	xout[6], yout[0] = xout[6]+twoQ-MRedConstant(yout[0], F, Q, QInv), yout[0]+twoQ-MRedConstant(xout[6], F, Q, QInv)

	coeffsOut[N>>1] = coeffsOut[N>>1] + twoQ - MRedConstant(coeffsOut[N>>1], F, Q, QInv)
	coeffsOut[0] = CRed(coeffsOut[0]<<1, Q)
}
