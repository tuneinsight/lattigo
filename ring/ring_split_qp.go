package ring

// PolyQP represents a polynomial in the ring of polynomial modulo Q*P.
// This type is simply the union type between two ring.Poly, each one
// containing the modulus Q and P coefficients of that polynomial.
// The modulus Q represent the ciphertext modulus and the modulus P
// the special primes for the RNS decomposition during homomorphic
// operations involving keys.
type PolyQP struct {
	Q, P *Poly
}

// Equals returns true if the receiver PolyQP is equal to the provided other PolyQP.
// This method checks for equality of its two sub-polynomials.
func (p *PolyQP) Equals(other *PolyQP) bool {
	return p == other || (p.P.Equals(other.P) && p.Q.Equals(other.Q))
}

// CopyValues copies the coefficients of p1 on the target polynomial.
// This method simply calls the CopyValues method for each of its sub-polynomials.
func (p *PolyQP) CopyValues(other *PolyQP) {
	p.Q.CopyValues(other.Q)
	p.P.CopyValues(other.P)
}

// CopyNew creates an exact copy of the target polynomial.
func (p *PolyQP) CopyNew() *PolyQP {
	if p == nil {
		return nil
	}
	return &PolyQP{Q: p.Q.CopyNew(), P: p.P.CopyNew()}
}

// SplitRingQP is a structure that implements the operation in the ring R_QP.
// This type is simply a union type between the two Ring types representing
// R_Q and R_P.
type SplitRingQP struct {
	RingQ, RingP *Ring
}

// NewPoly creates a new polynomial with all coefficients set to 0.
func (r *SplitRingQP) NewPoly() *PolyQP {
	return &PolyQP{r.RingQ.NewPoly(), r.RingP.NewPoly()}
}

// NewPolyLvl creates a new polynomial with all coefficients set to 0.
func (r *SplitRingQP) NewPolyLvl(levelQ, levelP int) *PolyQP {
	return &PolyQP{r.RingQ.NewPolyLvl(levelQ), r.RingP.NewPolyLvl(levelP)}
}

// Add adds p1 to p2 coefficient-wise and writes the result on p3.
func (r *SplitRingQP) Add(p1, p2, pOut *PolyQP) {
	r.RingQ.Add(p1.Q, p2.Q, pOut.Q)
	r.RingP.Add(p1.P, p2.P, pOut.P)
}

// Sub subtracts p2 to p1 coefficient-wise and writes the result on p3.
func (r *SplitRingQP) Sub(p1, p2, pOut *PolyQP) {
	r.RingQ.Sub(p1.Q, p2.Q, pOut.Q)
	r.RingP.Sub(p1.P, p2.P, pOut.P)
}

// NTT computes the NTT of p1 and returns the result on p2.
func (r *SplitRingQP) NTT(p, pOut *PolyQP) {
	r.RingQ.NTT(p.Q, pOut.Q)
	r.RingP.NTT(p.P, pOut.P)
}

// NTTLvl computes the NTT of p1 and returns the result on p2.
// The operation is performed at levelQ for in the ciphertext ring
// and at levelP for the extended ring.
func (r *SplitRingQP) NTTLvl(levelQ, levelP int, p, pOut *PolyQP) {
	r.RingQ.NTTLvl(levelQ, p.Q, pOut.Q)
	r.RingP.NTTLvl(levelP, p.P, pOut.P)
}

// InvNTT computes the inverse-NTT of p1 and returns the result on p2.
func (r *SplitRingQP) InvNTT(p, pOut *PolyQP) {
	r.RingQ.InvNTT(p.Q, pOut.Q)
	r.RingP.InvNTT(p.P, pOut.P)
}

// InvNTTLvl computes the inverse-NTT of p1 and returns the result on p2.
// The operation is performed at levelQ for in the ciphertext ring
// and at levelP for the extended ring.
func (r *SplitRingQP) InvNTTLvl(levelQ, levelP int, p, pOut *PolyQP) {
	r.RingQ.InvNTTLvl(levelQ, p.Q, pOut.Q)
	r.RingP.InvNTTLvl(levelP, p.P, pOut.P)
}

// NTTLazy computes the NTT of p1 and returns the result on p2.
// Output values are in the range [0, 2q-1]
func (r *SplitRingQP) NTTLazy(p, pOut *PolyQP) {
	r.RingQ.NTTLazy(p.Q, pOut.Q)
	r.RingP.NTTLazy(p.P, pOut.P)
}

// MForm switches p1 to the Montgomery domain and writes the result on p2.
func (r *SplitRingQP) MForm(p, pOut *PolyQP) {
	r.RingQ.MForm(p.Q, pOut.Q)
	r.RingP.MForm(p.P, pOut.P)
}

// MFormLvl switches p1 to the Montgomery domain and writes the result on p2.
func (r *SplitRingQP) MFormLvl(levelQ, levelP int, p, pOut *PolyQP) {
	r.RingQ.MFormLvl(levelQ, p.Q, pOut.Q)
	r.RingP.MFormLvl(levelP, p.P, pOut.P)
}

// InvMForm switches back p1 from the Montgomery domain to the conventional domain and writes the result on p2.
func (r *SplitRingQP) InvMForm(p, pOut *PolyQP) {
	r.RingQ.InvMForm(p.Q, pOut.Q)
	r.RingP.InvMForm(p.P, pOut.P)
}

// MulCoeffsMontgomery multiplies p1 by p2 coefficient-wise with a
// Montgomery modular reduction and returns the result on p3.
func (r *SplitRingQP) MulCoeffsMontgomery(p1, p2, p3 *PolyQP) {
	r.RingQ.MulCoeffsMontgomery(p1.Q, p2.Q, p3.Q)
	r.RingP.MulCoeffsMontgomery(p1.P, p2.P, p3.P)
}

// MulCoeffsMontgomeryLvl multiplies p1 by p2 coefficient-wise with a Montgomery
// modular reduction for the moduli from q_0 up to q_level and returns the result on p3.
func (r *SplitRingQP) MulCoeffsMontgomeryLvl(levelQ, levelP int, p1, p2, p3 *PolyQP) {
	r.RingQ.MulCoeffsMontgomeryLvl(levelQ, p1.Q, p2.Q, p3.Q)
	r.RingP.MulCoeffsMontgomeryLvl(levelP, p1.P, p2.P, p3.P)
}

// MulCoeffsMontgomeryConstant multiplies p1 by p2 coefficient-wise with a
// constant-time Montgomery modular reduction and writes the result on p3.
func (r *SplitRingQP) MulCoeffsMontgomeryConstant(p1, p2, p3 *PolyQP) {
	r.RingQ.MulCoeffsMontgomeryConstant(p1.Q, p2.Q, p3.Q)
	r.RingP.MulCoeffsMontgomeryConstant(p1.P, p2.P, p3.P)
}

// MulCoeffsMontgomeryAndSub multiplies p1 by p2 coefficient-wise with
// a Montgomery modular reduction and subtracts the result from p3.
func (r *SplitRingQP) MulCoeffsMontgomeryAndSub(p1, p2, p3 *PolyQP) {
	r.RingQ.MulCoeffsMontgomeryAndSub(p1.Q, p2.Q, p3.Q)
	r.RingP.MulCoeffsMontgomeryAndSub(p1.P, p2.P, p3.P)
}

// MulCoeffsMontgomeryAndAdd multiplies p1 by p2 coefficient-wise with a
// Montgomery modular reduction and adds the result to p3.
func (r *SplitRingQP) MulCoeffsMontgomeryAndAdd(p1, p2, p3 *PolyQP) {
	r.RingQ.MulCoeffsMontgomeryAndAdd(p1.Q, p2.Q, p3.Q)
	r.RingP.MulCoeffsMontgomeryAndAdd(p1.P, p2.P, p3.P)
}

// ExtendBasisSmallNormAndCenter extends a small-norm polynomial polQ in R_Q to a polynomial
// polQP in R_QP.
func (r *SplitRingQP) ExtendBasisSmallNormAndCenter(polQ *Poly, levelP int, polQP *PolyQP) {
	var coeff, Q, QHalf, sign uint64
	Q = r.RingQ.Modulus[0]
	QHalf = Q >> 1

	if polQ != polQP.Q {
		polQP.Q.Copy(polQ)
	}

	for j := 0; j < r.RingQ.N; j++ {

		coeff = polQ.Coeffs[0][j]

		sign = 1
		if coeff > QHalf {
			coeff = Q - coeff
			sign = 0
		}

		for i, pi := range r.RingP.Modulus[:levelP+1] {
			polQP.P.Coeffs[i][j] = (coeff * sign) | (pi-coeff)*(sign^1)
		}
	}
}

// Copy copies the input polyQP on the target polyQP.
func (p *PolyQP) Copy(polFrom *PolyQP) {
	p.Q.Copy(polFrom.Q)
	p.P.Copy(polFrom.P)
}

// GetDataLen returns the length in byte of the target PolyQP
func (p PolyQP) GetDataLen(WithMetadata bool) (dataLen int) {
	return p.Q.GetDataLen(WithMetadata) + p.P.GetDataLen(WithMetadata)
}

// WriteTo writes a polyQP on the inpute data.
func (p *PolyQP) WriteTo(data []byte) (pt int, err error) {
	var inc int
	if inc, err = p.Q.WriteTo(data[pt:]); err != nil {
		return
	}
	pt += inc

	if inc, err = p.P.WriteTo(data[pt:]); err != nil {
		return
	}
	pt += inc

	return
}

// DecodePolyNew decodes the input bytes on the target polyQP.
func (p *PolyQP) DecodePolyNew(data []byte) (pt int, err error) {
	p.Q = new(Poly)
	var inc int
	if inc, err = p.Q.DecodePolyNew(data[pt:]); err != nil {
		return
	}
	pt += inc

	p.P = new(Poly)
	if inc, err = p.P.DecodePolyNew(data[pt:]); err != nil {
		return
	}
	pt += inc

	return
}
