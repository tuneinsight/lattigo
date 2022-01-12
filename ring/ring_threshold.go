package ring

import (
	"errors"
)

type Scalar []uint64

func (r *Ring) NewScalar() Scalar {
	return make(Scalar, len(r.Modulus))
}

func (r *Ring) NewScalarFromUInt64(v uint64) Scalar {
	s := make(Scalar, len(r.Modulus))
	for i, qi := range r.Modulus {
		s[i] = v % qi
	}
	return s
}

// Unit sets all coeff of the target polynomial to 1.
func (pol *Poly) Unit() {
	for i := range pol.Coeffs {
		ptmp := pol.Coeffs[i]
		for j := range ptmp {
			ptmp[j] = 1
		}
	}
}

// IsInvertible returns true iff poly p is invertible
//in the ring r. Expects p in NTT domain.
func (r *Ring) IsInvertible(p *Poly) (invertible bool) {
	invertible = true
	for m := range r.Modulus {
		for _, coeff := range p.Coeffs[m] {
			invertible = invertible && (coeff != 0)
		}
	}
	return
}

func (r *Ring) SubScalarCRT(s1, s2, sout Scalar) {
	for i, qi := range r.Modulus {
		if s1[i] > s2[i] {
			sout[i] = s2[i] + qi - s1[i]
		} else {
			sout[i] = s2[i] - s1[i]
		}
	}
}

func (r *Ring) ScalarMulCRT(s1, s2, sout Scalar) {
	for i, qi := range r.Modulus {
		sout[i] = MRedConstant(s1[i], s2[i], qi, r.MredParams[i])
	}
}

// InvMultPolyMontgomeryNTT computes the multiplicative inverse
// of a polynomial p1 (in NTT and Montgomery domain) and stores it in p2
// Returns a non-nil error iff p1 is not invertible.
func (r *Ring) InvMultPolyMontgomeryNTT(p1 *Poly, p2 *Poly) error {

	tmp := p1.CopyNew()
	p2.Unit()
	r.MForm(p2, p2)

	if !(r.IsInvertible(p1)) {
		return errors.New("invalid polynomial inversion (element has no multiplicative inverse)")
	}

	//EEA
	for i, qi := range r.Modulus {

		p2tmp, ptmp := p2.Coeffs[i], tmp.Coeffs[i]
		mredParam := r.MredParams[i]

		for k := qi - 2; k > 0; k >>= 1 {

			if k&1 == 1 {
				for j := 0; j < r.N; j++ {
					p2tmp[j] = MRedConstant(p2tmp[j], ptmp[j], qi, mredParam)
				}
			}

			for j := 0; j < r.N; j++ {
				ptmp[j] = MRedConstant(ptmp[j], ptmp[j], qi, mredParam)
			}
		}
	}
	return nil
}

// InverseCRT computes the modular inverse of a scalar a expressed in a CRT decomposition.
// The inversion is done in-place and assumes that a is in Montgomery form.
func (r *Ring) InverseCRT(a []uint64) {
	for i, qi := range r.Modulus {
		a[i] = ModexpMontgomery(a[i], int(qi-2), qi, r.MredParams[i], r.BredParams[i])
	}
}

// EvalPolMontgomeryNTT evaluate the polynomial pol at pk and writes the result in p3
func (r *Ring) EvalPolMontgomeryNTT(pol []*Poly, pk *Poly, p3 *Poly) {
	p3.Copy(pol[len(pol)-1])
	for i := len(pol) - 1; i > 0; i-- {
		r.MulCoeffsMontgomeryConstant(p3, pk, p3)
		r.AddNoMod(p3, pol[i-1], p3)
	}
	r.Reduce(p3, p3)
}

// EvalPolMontgomeryScalarNTT evaluate the polynomial pol at pk and writes the result in p3
func (r *Ring) EvalPolMontgomeryScalarNTT(pol []*Poly, pk uint64, p3 *Poly) {
	p3.Copy(pol[len(pol)-1])
	for i := len(pol) - 1; i > 0; i-- {
		//r.MulCoeffsMontgomeryConstant(p3, pk, p3)
		r.MulScalar(p3, pk, p3)
		r.AddNoMod(p3, pol[i-1], p3)
	}
	r.Reduce(p3, p3)
}
