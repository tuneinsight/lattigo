package ring

import (
	"math/bits"
	"unsafe"
)

// Interpolator is a struct storing the necessary
// buffer and pre-computation for polynomial interpolation
// with coefficient in finite fields.
type Interpolator struct {
	r *Ring
	x Poly
}

// NewInterpolator creates a new Interpolator. Returns an error if T is not
// prime or not congruent to 1 mod 2N, where N is the next power of two greater
// than degree+1.
func NewInterpolator(degree int, T uint64) (itp *Interpolator, err error) {
	itp = new(Interpolator)

	if itp.r, err = NewRing(1<<bits.Len64(uint64(degree)), []uint64{T}); err != nil {
		return nil, err
	}

	// NTT(x)
	itp.x = itp.r.NewPoly()
	itp.x.Coeffs[0][1] = 1
	itp.r.NTT(itp.x, itp.x)
	itp.r.MForm(itp.x, itp.x)

	return
}

// Interpolate takes a list of roots the coefficients of P(roots) = 0 mod T.
func (itp *Interpolator) Interpolate(roots []uint64) (coeffs []uint64) {

	r := itp.r
	s := r.SubRings[0]
	x := itp.x
	T := s.Modulus
	mredParams := s.MRedConstant
	bredParams := s.BRedConstant

	// res = NTT(x-root[0])
	res := *itp.x.CopyNew()
	r.SubScalar(res, MForm(roots[0], T, bredParams), res)

	// res = res * (x-root[i])
	for i := 1; i < len(roots); i++ {
		subScalarMontgomeryAndMulCoeffsMontgomery(x.Coeffs[0], MForm(roots[i], T, bredParams), res.Coeffs[0], res.Coeffs[0], T, mredParams)
	}

	r.INTT(res, res)

	return res.Coeffs[0][:len(roots)+1]
}

// Lagrange takes as input (x, y) and returns P(xi) = yi mod T.
func (itp *Interpolator) Lagrange(x, y []uint64) (coeffs []uint64, err error) {
	r := itp.r
	s := r.SubRings[0]
	X := itp.x
	T := s.Modulus
	N := r.N()
	mredParams := s.MRedConstant
	bredParams := s.BRedConstant

	// Powers of w are stored in bit-reversed order -> even powers of w are on the right n half
	roots := make(map[uint64]bool) // -> map that stores all the roots of X^{N} + 1 mod T
	for i := 0; i < N>>1; i++ {
		roots[s.RootsForward[N>>1+i]] = true
		roots[s.RootsBackward[N>>1+i]] = true
	}

	basis := r.NewPoly()
	for i := 0; i < N; i++ {
		basis.Coeffs[0][i] = 1
	}

	// Computes the Lagrange basis (X-x[0]) * (X-x[1]) * ... * (X-x[i])
	// but omits x[i] which are roots of X^{N} + 1 mod T.
	// The roots of X^{N} + 1 mod T are the even powers of w, where w is
	// is a primitive 2N-th roots of unity mod T.
	missing := make(map[uint64]bool)
	for i := 0; i < len(x); i++ {
		if _, ok := roots[x[i]]; ok {
			missing[x[i]] = true
		} else {
			subScalarMontgomeryAndMulCoeffsMontgomery(X.Coeffs[0], MForm(x[i], T, bredParams), basis.Coeffs[0], basis.Coeffs[0], T, mredParams)
		}
	}

	poly := r.NewPoly()
	tmp := r.NewPoly()
	tmp1 := r.NewPoly()

	for i := 0; i < len(x); i++ {

		tmp.Copy(basis)

		// If x[i] is a root of X^{N} + 1 mod T then it is not part
		// of the Lagrange basis pre-computation, so all we need is
		// to add the missing roots (if any), skipping x[i].
		if _, ok := missing[x[i]]; ok {

			// with the missing roots, except x[i]
			for root := range missing {
				if root != x[i] {
					subScalarMontgomeryAndMulCoeffsMontgomery(X.Coeffs[0], MForm(root, T, bredParams), tmp.Coeffs[0], tmp.Coeffs[0], T, mredParams)
				}
			}

			// If x[i] is not a root of X^{N} + 1 mod T, then we need
			// to remove it from the Lagrange basis pre-computation.
			// But first we add the missing x[i], which are the
			// roots of X^{N} + 1 mod T (if any).
		} else {
			// Continue with all the missing roots
			for root := range missing {
				subScalarMontgomeryAndMulCoeffsMontgomery(X.Coeffs[0], MForm(root, T, bredParams), tmp.Coeffs[0], tmp.Coeffs[0], T, mredParams)
			}

			// And then removes (X - x[i])
			s.SubScalar(X.Coeffs[0], x[i], tmp1.Coeffs[0])

			// TODO: unrol loop and use unsafe
			coeffs := tmp1.Coeffs[0]
			for j := 0; j < N; j++ {
				coeffs[j] = ModexpMontgomery(coeffs[j], int(T-2), T, mredParams, bredParams)
			}

			s.MulCoeffsMontgomery(tmp.Coeffs[0], tmp1.Coeffs[0], tmp.Coeffs[0])
		}

		// prod(x[i] - x[j]) i != j
		// TODO: make 2 iterations to avoid the if condition
		var den uint64 = 1
		for j := 0; j < len(x); j++ {
			if j != i {
				den = BRed(den, x[i]+T-x[j], T, bredParams)
			}
		}

		// 1 / prod(x[i] - x[j])
		den = ModExp(den, T-2, T)

		// y[i] / prod(x[i] - x[j])
		den = BRed(y[i], den, T, bredParams)

		// P(X) += (y[i] / prod(x[i] - x[j])) * prod(X-x[j])
		s.MulScalarMontgomeryThenAdd(tmp.Coeffs[0], MForm(den, T, bredParams), poly.Coeffs[0])
	}

	r.INTT(poly, poly)

	return poly.Coeffs[0][:len(x)], nil

}

// computes p3 = (p1 - a) * p2
func subScalarMontgomeryAndMulCoeffsMontgomery(p1 []uint64, a uint64, p2, p3 []uint64, t, mredParams uint64) {
	for j := 0; j < len(p1); j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 != 0 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 != 0 */
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p3)%8 != 0 */
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = MRedLazy(x[0]+t-a, y[0], t, mredParams)
		z[1] = MRedLazy(x[1]+t-a, y[1], t, mredParams)
		z[2] = MRedLazy(x[2]+t-a, y[2], t, mredParams)
		z[3] = MRedLazy(x[3]+t-a, y[3], t, mredParams)
		z[4] = MRedLazy(x[4]+t-a, y[4], t, mredParams)
		z[5] = MRedLazy(x[5]+t-a, y[5], t, mredParams)
		z[6] = MRedLazy(x[6]+t-a, y[6], t, mredParams)
		z[7] = MRedLazy(x[7]+t-a, y[7], t, mredParams)
	}
}
