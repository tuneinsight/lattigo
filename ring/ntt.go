package ring

// NTT performes the NTT transformation on the coefficients a Polynomial,
// based on the Context of the Polynomial.
// https://arxiv.org/abs/1205.2926v2

// NTT performes the NTT transformation on the coefficients of a polynomial.
func (context *Context) NTT(p1, p2 *Poly) {
	for x := range context.Modulus {
		NTT(p1.Coeffs[x], p2.Coeffs[x], context.N, context.nttPsi[x], context.Modulus[x], context.mredParams[x], context.bredParams[x])
	}
}

// InvNTT performes the inverse NTT transformation on the coefficients of a polynomial,
func (context *Context) InvNTT(p1, p2 *Poly) {
	for x := range context.Modulus {
		InvNTT(p1.Coeffs[x], p2.Coeffs[x], context.N, context.nttPsiInv[x], context.nttNInv[x], context.Modulus[x], context.mredParams[x])
	}
}

func Butterfly(U, V, Psi, Q, Qinv uint64) (X, Y uint64) {
	if U > 2*Q {
		U -= 2 * Q
	}
	V = MRedConstant(V, Psi, Q, Qinv)
	X = U + V
	Y = U + 2*Q - V
	return
}

func InvButterfly(U, V, Psi, Q, Qinv uint64) (X, Y uint64) {
	X = U + V
	if X > 2*Q {
		X -= 2 * Q
	}
	Y = MRedConstant(U+2*Q-V, Psi, Q, Qinv) // At the moment it is not possible to use MRedConstant if Q > 61 bits
	return
}

func NTT(coeffs_in, coeffs_out []uint64, N uint64, nttPsi []uint64, Q, mredParams uint64, bredParams []uint64) {
	var j1, j2, t uint64
	var F uint64

	// Copies the result of the first round of butterflies on p2 with approximate reduction
	t = N >> 1
	j2 = t - 1
	F = nttPsi[1]
	for j := uint64(0); j <= j2; j++ {
		coeffs_out[j], coeffs_out[j+t] = Butterfly(coeffs_in[j], coeffs_in[j+t], F, Q, mredParams)
	}

	// Continues the rest of the second to the n-1 butterflies on p2 with approximate reduction
	for m := uint64(2); m < N; m <<= 1 {
		t >>= 1
		for i := uint64(0); i < m; i++ {

			j1 = (i * t) << 1

			j2 = j1 + t - 1

			F = nttPsi[m+i]

			for j := j1; j <= j2; j++ {
				coeffs_out[j], coeffs_out[j+t] = Butterfly(coeffs_out[j], coeffs_out[j+t], F, Q, mredParams)
			}
		}
	}

	// Finishes with an exact reduction
	for i := uint64(0); i < N; i++ {
		coeffs_out[i] = BRedAdd(coeffs_out[i], Q, bredParams)
	}
}

func InvNTT(coeffs_in, coeffs_out []uint64, N uint64, nttPsiInv []uint64, nttNInv, Q, mredParams uint64) {

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
			coeffs_out[j], coeffs_out[j+t] = InvButterfly(coeffs_in[j], coeffs_in[j+t], F, Q, mredParams)
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
				coeffs_out[j], coeffs_out[j+t] = InvButterfly(coeffs_out[j], coeffs_out[j+t], F, Q, mredParams)
			}

			j1 = j1 + (t << 1)
		}

		t <<= 1
	}

	// Finishes with an exact reduction given
	for j := uint64(0); j < N; j++ {
		coeffs_out[j] = MRed(coeffs_out[j], nttNInv, Q, mredParams)
	}
}
