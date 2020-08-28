package ring

import(
	"unsafe"
)

// NTT computes the NTT transformation of p1 and returns the result on p2.
func (context *Context) NTT(p1, p2 *Poly) {
	for x := range context.Modulus {
		NTT(p1.Coeffs[x], p2.Coeffs[x], context.N, context.nttPsi[x], context.Modulus[x], context.mredParams[x], context.bredParams[x])
	}
}

// NTTLvl computes the NTT transformation of p1 and returns the result on p2.
// The value level defines the number of moduli of the input polynomials.
func (context *Context) NTTLvl(level uint64, p1, p2 *Poly) {
	for x := uint64(0); x < level+1; x++ {
		NTT(p1.Coeffs[x], p2.Coeffs[x], context.N, context.nttPsi[x], context.Modulus[x], context.mredParams[x], context.bredParams[x])
	}
}

// InvNTT computes the inverse NTT transformation of p1 and returns the result on p2.
func (context *Context) InvNTT(p1, p2 *Poly) {
	for x := range context.Modulus {
		InvNTT(p1.Coeffs[x], p2.Coeffs[x], context.N, context.nttPsiInv[x], context.nttNInv[x], context.Modulus[x], context.mredParams[x])
	}
}

// InvNTTLvl computes the inverse NTT transformation of p1 and returns the result on p2.
// The value level defines the number of moduli of the input polynomials.
func (context *Context) InvNTTLvl(level uint64, p1, p2 *Poly) {
	for x := uint64(0); x < level+1; x++ {
		InvNTT(p1.Coeffs[x], p2.Coeffs[x], context.N, context.nttPsiInv[x], context.nttNInv[x], context.Modulus[x], context.mredParams[x])
	}
}

// butterfly computes X, Y = U + V*Psi, U - V*Psi mod Q.
func butterfly(U, V, Psi, Q, Qinv uint64) (X, Y uint64) {
	V = MRed(V, Psi, Q, Qinv)
	X = U + V
	Y = U + Q - V
	return
}

// invbutterfly computes X, Y = U + V, (U - V) * Psi mod Q.
func invbutterfly(U, V, Psi, Q, Qinv uint64) (X, Y uint64) {
	X = U + V
	if X > 2*Q {
		X -= 2 * Q
	}
	Y = MRedConstant(U+2*Q-V, Psi, Q, Qinv) // At the moment it is not possible to use MRedConstant if Q > 61 bits
	return
}

// NTT computes the NTT transformation on the input coefficients using the input parameters.
func NTT(coeffsIn, coeffsOut []uint64, N uint64, nttPsi []uint64, Q, mredParams uint64, bredParams []uint64) {
	var j1, j2, t uint64
	var F uint64

	// Copies the result of the first round of butterflies on p2 with approximate reduction
	t = N >> 1
	j2 = t - 1
	F = nttPsi[1]
	for j := uint64(0); j <= j2; j = j+8 {

		x_in := (*[8]uint64)(unsafe.Pointer(&coeffsIn[j]))
		y_in := (*[8]uint64)(unsafe.Pointer(&coeffsIn[j+t]))

		x_out := (*[8]uint64)(unsafe.Pointer(&coeffsOut[j]))
		y_out := (*[8]uint64)(unsafe.Pointer(&coeffsOut[j+t]))

		x_out[0], y_out[0] = butterfly(x_in[0], y_in[0], F, Q, mredParams)
		x_out[1], y_out[1] = butterfly(x_in[1], y_in[1], F, Q, mredParams)
		x_out[2], y_out[2] = butterfly(x_in[2], y_in[2], F, Q, mredParams)
		x_out[3], y_out[3] = butterfly(x_in[3], y_in[3], F, Q, mredParams)
		x_out[4], y_out[4] = butterfly(x_in[4], y_in[4], F, Q, mredParams)
		x_out[5], y_out[5] = butterfly(x_in[5], y_in[5], F, Q, mredParams)
		x_out[6], y_out[6] = butterfly(x_in[6], y_in[6], F, Q, mredParams)
		x_out[7], y_out[7] = butterfly(x_in[7], y_in[7], F, Q, mredParams)
	}

	// Continues the rest of the second to the n-1 butterflies on p2 with approximate reduction
	for m := uint64(2); m < N; m <<= 1 {

		t >>= 1

		if t >= 8 {

			for i := uint64(0); i < m; i++ {

				j1 = (i * t) << 1

				j2 = j1 + t - 1

				F = nttPsi[m+i]

				for j := j1; j <= j2; j = j+8 {

					x := (*[8]uint64)(unsafe.Pointer(&coeffsOut[j]))
					y := (*[8]uint64)(unsafe.Pointer(&coeffsOut[j+t]))

					x[0], y[0] = butterfly(x[0], y[0], F, Q, mredParams)
					x[1], y[1] = butterfly(x[1], y[1], F, Q, mredParams)
					x[2], y[2] = butterfly(x[2], y[2], F, Q, mredParams)
					x[3], y[3] = butterfly(x[3], y[3], F, Q, mredParams)
					x[4], y[4] = butterfly(x[4], y[4], F, Q, mredParams)
					x[5], y[5] = butterfly(x[5], y[5], F, Q, mredParams)
					x[6], y[6] = butterfly(x[6], y[6], F, Q, mredParams)
					x[7], y[7] = butterfly(x[7], y[7], F, Q, mredParams)
				}
			}

		}else if t == 4 {

			for i := uint64(0); i < m; i++ {

				j1 = (i * t) << 1

				j2 = j1 + t - 1

				F = nttPsi[m+i]

				for j := j1; j <= j2; j = j+4 {

					x := (*[4]uint64)(unsafe.Pointer(&coeffsOut[j]))
					y := (*[4]uint64)(unsafe.Pointer(&coeffsOut[j+t]))

					x[0], y[0] = butterfly(x[0], y[0], F, Q, mredParams)
					x[1], y[1] = butterfly(x[1], y[1], F, Q, mredParams)
					x[2], y[2] = butterfly(x[2], y[2], F, Q, mredParams)
					x[3], y[3] = butterfly(x[3], y[3], F, Q, mredParams)
				}
			}

		} else {

			for i := uint64(0); i < m; i++ {

				j1 = (i * t) << 1

				j2 = j1 + t - 1

				F = nttPsi[m+i]

				for j := j1; j <= j2; j++ {
					coeffsOut[j], coeffsOut[j+t] = butterfly(coeffsOut[j], coeffsOut[j+t], F, Q, mredParams)
				}
			}
		}
	}


	// Finishes with an exact reduction
	for i := uint64(0); i < N; i = i + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&coeffsOut[i]))

		x[0] = BRedAdd(x[0], Q, bredParams)
		x[1] = BRedAdd(x[1], Q, bredParams)
		x[2] = BRedAdd(x[2], Q, bredParams)
		x[3] = BRedAdd(x[3], Q, bredParams)
		x[4] = BRedAdd(x[4], Q, bredParams)
		x[5] = BRedAdd(x[5], Q, bredParams)
		x[6] = BRedAdd(x[6], Q, bredParams)
		x[7] = BRedAdd(x[7], Q, bredParams)
	}

}

// InvNTT computes the InvNTT transformation on the input coefficients using the input parameters.
func InvNTT(coeffsIn, coeffsOut []uint64, N uint64, nttPsiInv []uint64, nttNInv, Q, mredParams uint64) {

	var j1, j2, h, t uint64
	var F uint64

	// Copies the result of the first round of butterflies on p2 with approximate reduction
	t = 1
	h = N >> 1

	for i := uint64(0); i < h; i = i + 8 {

		psi   := (*[8]uint64)(unsafe.Pointer(&nttPsiInv[h + i]))
		x_in  := (*[16]uint64)(unsafe.Pointer(&coeffsIn[2*i]))
		x_out := (*[16]uint64)(unsafe.Pointer(&coeffsOut[2*i]))

		x_out[ 0], x_out[ 1] = invbutterfly(x_in[ 0], x_in[ 1], psi[0], Q, mredParams)
		x_out[ 2], x_out[ 3] = invbutterfly(x_in[ 2], x_in[ 3], psi[1], Q, mredParams)
		x_out[ 4], x_out[ 5] = invbutterfly(x_in[ 4], x_in[ 5], psi[2], Q, mredParams)
		x_out[ 6], x_out[ 7] = invbutterfly(x_in[ 6], x_in[ 7], psi[3], Q, mredParams)
		x_out[ 8], x_out[ 9] = invbutterfly(x_in[ 8], x_in[ 9], psi[4], Q, mredParams)
		x_out[10], x_out[11] = invbutterfly(x_in[10], x_in[11], psi[5], Q, mredParams)
		x_out[12], x_out[13] = invbutterfly(x_in[12], x_in[13], psi[6], Q, mredParams)
		x_out[14], x_out[15] = invbutterfly(x_in[14], x_in[15], psi[7], Q, mredParams)
	}

	// Continues the rest of the second to the n-1 butterflies on p2 with approximate reduction
	t <<= 1
	for m := N >> 1; m > 1; m >>= 1 {

		j1 = 0
		h = m >> 1

		if t >= 8 {

			for i := uint64(0); i < h; i++ {

				j2 = j1 + t - 1

				F = nttPsiInv[h+i]

				for j := j1; j <= j2; j =  j + 8 {

					x := (*[8]uint64)(unsafe.Pointer(&coeffsOut[j]))
					y := (*[8]uint64)(unsafe.Pointer(&coeffsOut[j+t]))

					x[0], y[0] = invbutterfly(x[0], y[0], F, Q, mredParams)
					x[1], y[1] = invbutterfly(x[1], y[1], F, Q, mredParams)
					x[2], y[2] = invbutterfly(x[2], y[2], F, Q, mredParams)
					x[3], y[3] = invbutterfly(x[3], y[3], F, Q, mredParams)
					x[4], y[4] = invbutterfly(x[4], y[4], F, Q, mredParams)
					x[5], y[5] = invbutterfly(x[5], y[5], F, Q, mredParams)
					x[6], y[6] = invbutterfly(x[6], y[6], F, Q, mredParams)
					x[7], y[7] = invbutterfly(x[7], y[7], F, Q, mredParams)
				}

				j1 = j1 + (t << 1)
			}

		}else if t == 4 {

			for i := uint64(0); i < h; i++ {

				j2 = j1 + t - 1

				F = nttPsiInv[h+i]

				for j := j1; j <= j2; j =  j + 4 {

					x := (*[4]uint64)(unsafe.Pointer(&coeffsOut[j]))
					y := (*[4]uint64)(unsafe.Pointer(&coeffsOut[j+4]))

					x[0], y[0] = invbutterfly(x[0], y[0], F, Q, mredParams)
					x[1], y[1] = invbutterfly(x[1], y[1], F, Q, mredParams)
					x[2], y[2] = invbutterfly(x[2], y[2], F, Q, mredParams)
					x[3], y[3] = invbutterfly(x[3], y[3], F, Q, mredParams)
				}

				j1 = j1 + (t << 1)
			}

		}else{

			for i := uint64(0); i < h; i++ {

				j2 = j1 + t - 1

				F = nttPsiInv[h+i]

				for j := j1; j <= j2; j++ {
					coeffsOut[j], coeffsOut[j+t] = invbutterfly(coeffsOut[j], coeffsOut[j+t], F, Q, mredParams)
				}

				j1 = j1 + (t << 1)
			}
		}

		t <<= 1
	}

	// Finishes with an exact reduction
	for i := uint64(0); i < N; i = i + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&coeffsOut[i]))

		x[0] = MRed(x[0], nttNInv, Q, mredParams)
		x[1] = MRed(x[1], nttNInv, Q, mredParams)
		x[2] = MRed(x[2], nttNInv, Q, mredParams)
		x[3] = MRed(x[3], nttNInv, Q, mredParams)
		x[4] = MRed(x[4], nttNInv, Q, mredParams)
		x[5] = MRed(x[5], nttNInv, Q, mredParams)
		x[6] = MRed(x[6], nttNInv, Q, mredParams)
		x[7] = MRed(x[7], nttNInv, Q, mredParams)
	}
}

///////////////////////////////////
/// For benchmark purposes only ///
///////////////////////////////////

// NTTBarrett performs the NTT operation using Barrett reduction.
// For benchmark purposes only.
func (context *Context) NTTBarrett(p1, p2 *Poly) {
	for x := range context.Modulus {
		NTTBarrett(p1.Coeffs[x], p2.Coeffs[x], context.N, context.nttPsi[x], context.Modulus[x], context.bredParams[x])
	}
}

// InvNTTBarrett performs the inverse NTT operation using Barrett reduction.
// For benchmark purposes only.
func (context *Context) InvNTTBarrett(p1, p2 *Poly) {
	for x := range context.Modulus {
		InvNTTBarrett(p1.Coeffs[x], p2.Coeffs[x], context.N, context.nttPsiInv[x], context.nttNInv[x], context.Modulus[x], context.bredParams[x])
	}
}

// butterflyBarrett computes X, Y = U + V*Psi, U - V*Psi mod Q using Barrett reduction.
// For benchmark purposes only.
func butterflyBarrett(U, V, Psi, Q uint64, bredParams []uint64) (X, Y uint64) {
	if U > 2*Q {
		U -= 2 * Q
	}
	V = BRedConstant(V, Psi, Q, bredParams)
	X = U + V
	Y = U + 2*Q - V
	return
}

// invbutterflyBarrett computes X, Y = U + V, (U - V) * Psi mod Q using Barrett reduction.
// For benchmark purposes only.
func invbutterflyBarrett(U, V, Psi, Q uint64, bredParams []uint64) (X, Y uint64) {
	X = U + V
	if X > 2*Q {
		X -= 2 * Q
	}
	Y = BRedConstant(U+2*Q-V, Psi, Q, bredParams) // At the moment it is not possible to use MRedConstant if Q > 61 bits
	return
}

// NTTBarrett computes the NTT transformation using Barrett reduction.
// For benchmark purposes only.
func NTTBarrett(coeffsIn, coeffsOut []uint64, N uint64, nttPsi []uint64, Q uint64, bredParams []uint64) {
	var j1, j2, t uint64
	var F uint64

	t = N >> 1
	j2 = t - 1
	F = nttPsi[1]
	for j := uint64(0); j <= j2; j++ {
		coeffsOut[j], coeffsOut[j+t] = butterflyBarrett(coeffsIn[j], coeffsIn[j+t], F, Q, bredParams)
	}

	for m := uint64(2); m < N; m <<= 1 {
		t >>= 1
		for i := uint64(0); i < m; i++ {

			j1 = (i * t) << 1

			j2 = j1 + t - 1

			F = nttPsi[m+i]

			for j := j1; j <= j2; j++ {
				coeffsOut[j], coeffsOut[j+t] = butterflyBarrett(coeffsOut[j], coeffsOut[j+t], F, Q, bredParams)
			}
		}
	}

	for i := uint64(0); i < N; i++ {
		coeffsOut[i] = BRedAdd(coeffsOut[i], Q, bredParams)
	}
}

// InvNTTBarrett computes the InvNTT transformation using Barrett reduction.
// For benchmark purposes only.
func InvNTTBarrett(coeffsIn, coeffsOut []uint64, N uint64, nttPsiInv []uint64, nttNInv, Q uint64, bredParams []uint64) {

	var j1, j2, h, t uint64
	var F uint64

	t = 1
	j1 = 0
	h = N >> 1

	for i := uint64(0); i < h; i++ {

		j2 = j1

		F = nttPsiInv[h+i]

		for j := j1; j <= j2; j++ {
			coeffsOut[j], coeffsOut[j+t] = invbutterflyBarrett(coeffsIn[j], coeffsIn[j+t], F, Q, bredParams)
		}

		j1 = j1 + (t << 1)
	}

	t <<= 1
	for m := N >> 1; m > 1; m >>= 1 {

		j1 = 0
		h = m >> 1

		for i := uint64(0); i < h; i++ {

			j2 = j1 + t - 1

			F = nttPsiInv[h+i]

			for j := j1; j <= j2; j++ {
				coeffsOut[j], coeffsOut[j+t] = invbutterflyBarrett(coeffsOut[j], coeffsOut[j+t], F, Q, bredParams)
			}

			j1 = j1 + (t << 1)
		}

		t <<= 1
	}

	for j := uint64(0); j < N; j++ {
		coeffsOut[j] = BRed(coeffsOut[j], nttNInv, Q, bredParams)
	}
}
