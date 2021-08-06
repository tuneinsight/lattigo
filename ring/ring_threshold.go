package ring

import (
	"errors"
)

//Sets all coeff of the target polynomial to 1.
func (pol *Poly) Unit() {
	for i := range pol.Coeffs {
		ptmp := pol.Coeffs[i]
		for j := range ptmp {
			ptmp[j] = 1
		}
	}
}

//isInvertible returns true iff poly p is invertible
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

//InvMultPolyMontgomeryNTT computes the multiplicative inverse
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

// EvalPol evaluate the polynomial pol at pk and writes the result in p3
func (r *Ring) EvalPolMontgomeryNTT(pol []*Poly, pk *Poly, p3 *Poly) {
	p3.Copy(pol[len(pol)-1])
	for i := len(pol) - 1; i > 0; i-- {
		r.MulCoeffsMontgomeryConstant(p3, pk, p3)
		r.AddNoMod(p3, pol[i-1], p3)
	}
	r.Reduce(p3, p3)
}
