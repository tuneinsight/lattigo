package ring

// NTT computes the NTT transformation of p1 and returns the result on p2.
func (context *Ring) NTT(p1, p2 *Poly) {
	for x := range context.Modulus {
		NTT(p1.Coeffs[x], p2.Coeffs[x], context.N, context.nttPsi[x], context.Modulus[x], context.mredParams[x], context.bredParams[x])
	}
}

// NTTLvl computes the NTT transformation of p1 and returns the result on p2.
// The value level defines the number of moduli of the input polynomials.
func (context *Ring) NTTLvl(level uint64, p1, p2 *Poly) {
	for x := uint64(0); x < level+1; x++ {
		NTT(p1.Coeffs[x], p2.Coeffs[x], context.N, context.nttPsi[x], context.Modulus[x], context.mredParams[x], context.bredParams[x])
	}
}

// InvNTT computes the inverse NTT transformation of p1 and returns the result on p2.
func (context *Ring) InvNTT(p1, p2 *Poly) {
	for x := range context.Modulus {
		InvNTT(p1.Coeffs[x], p2.Coeffs[x], context.N, context.nttPsiInv[x], context.nttNInv[x], context.Modulus[x], context.mredParams[x])
	}
}

// InvNTTLvl computes the inverse NTT transformation of p1 and returns the result on p2.
// The value level defines the number of moduli of the input polynomials.
func (context *Ring) InvNTTLvl(level uint64, p1, p2 *Poly) {
	for x := uint64(0); x < level+1; x++ {
		InvNTT(p1.Coeffs[x], p2.Coeffs[x], context.N, context.nttPsiInv[x], context.nttNInv[x], context.Modulus[x], context.mredParams[x])
	}
}

// butterfly computes X, Y = U + V*Psi, U - V*Psi mod Q.
func butterfly(U, V, Psi, Q, Qinv uint64) (X, Y uint64) {
	if U > 2*Q {
		U -= 2 * Q
	}
	V = MRedConstant(V, Psi, Q, Qinv)
	X = U + V
	Y = U + 2*Q - V
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
	for j := uint64(0); j <= j2; j++ {
		coeffsOut[j], coeffsOut[j+t] = butterfly(coeffsIn[j], coeffsIn[j+t], F, Q, mredParams)
	}

	// Continues the rest of the second to the n-1 butterflies on p2 with approximate reduction
	for m := uint64(2); m < N; m <<= 1 {
		t >>= 1
		for i := uint64(0); i < m; i++ {

			j1 = (i * t) << 1

			j2 = j1 + t - 1

			F = nttPsi[m+i]

			for j := j1; j <= j2; j++ {
				coeffsOut[j], coeffsOut[j+t] = butterfly(coeffsOut[j], coeffsOut[j+t], F, Q, mredParams)
			}
		}
	}

	// Finishes with an exact reduction
	for i := uint64(0); i < N; i++ {
		coeffsOut[i] = BRedAdd(coeffsOut[i], Q, bredParams)
	}
}

// InvNTT computes the InvNTT transformation on the input coefficients using the input parameters.
func InvNTT(coeffsIn, coeffsOut []uint64, N uint64, nttPsiInv []uint64, nttNInv, Q, mredParams uint64) {

	var j1, j2, h, t uint64
	var F uint64

	// Copies the result of the first round of butterflies on p2 with approximate reduction
	t = 1
	j1 = 0
	h = N >> 1

	for i := uint64(0); i < h; i++ {

		j2 = j1

		F = nttPsiInv[h+i]

		for j := j1; j <= j2; j++ {
			coeffsOut[j], coeffsOut[j+t] = invbutterfly(coeffsIn[j], coeffsIn[j+t], F, Q, mredParams)
		}

		j1 = j1 + (t << 1)
	}

	// Continues the rest of the second to the n-1 butterflies on p2 with approximate reduction
	t <<= 1
	for m := N >> 1; m > 1; m >>= 1 {

		j1 = 0
		h = m >> 1

		for i := uint64(0); i < h; i++ {

			j2 = j1 + t - 1

			F = nttPsiInv[h+i]

			for j := j1; j <= j2; j++ {
				coeffsOut[j], coeffsOut[j+t] = invbutterfly(coeffsOut[j], coeffsOut[j+t], F, Q, mredParams)
			}

			j1 = j1 + (t << 1)
		}

		t <<= 1
	}

	// Finishes with an exact reduction
	for j := uint64(0); j < N; j++ {
		coeffsOut[j] = MRed(coeffsOut[j], nttNInv, Q, mredParams)
	}
}

///////////////////////////////////
/// For benchmark purposes only ///
///////////////////////////////////

// NTTBarrett performs the NTT operation using Barrett reduction.
// For benchmark purposes only.
func (context *Ring) NTTBarrett(p1, p2 *Poly) {
	for x := range context.Modulus {
		NTTBarrett(p1.Coeffs[x], p2.Coeffs[x], context.N, context.nttPsi[x], context.Modulus[x], context.bredParams[x])
	}
}

// InvNTTBarrett performs the inverse NTT operation using Barrett reduction.
// For benchmark purposes only.
func (context *Ring) InvNTTBarrett(p1, p2 *Poly) {
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
