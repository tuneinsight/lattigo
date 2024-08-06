package ringqp

import (
	"github.com/tuneinsight/lattigo/v6/ring"
)

// Add adds p1 to p2 coefficient-wise and writes the result on p3.
func (r Ring) Add(p1, p2, p3 Poly) {
	if r.RingQ != nil {
		r.RingQ.Add(p1.Q, p2.Q, p3.Q)
	}
	if r.RingP != nil {
		r.RingP.Add(p1.P, p2.P, p3.P)
	}
}

// AddLazy adds p1 to p2 coefficient-wise and writes the result on p3 without modular reduction.
func (r Ring) AddLazy(p1, p2, p3 Poly) {
	if r.RingQ != nil {
		r.RingQ.AddLazy(p1.Q, p2.Q, p3.Q)
	}
	if r.RingP != nil {
		r.RingP.AddLazy(p1.P, p2.P, p3.P)
	}
}

// Sub subtracts p2 to p1 coefficient-wise and writes the result on p3.
func (r Ring) Sub(p1, p2, p3 Poly) {
	if r.RingQ != nil {
		r.RingQ.Sub(p1.Q, p2.Q, p3.Q)
	}
	if r.RingP != nil {
		r.RingP.Sub(p1.P, p2.P, p3.P)
	}
}

// Neg negates p1 coefficient-wise and writes the result on p2.
func (r Ring) Neg(p1, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.Neg(p1.Q, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.Neg(p1.P, p2.P)
	}
}

// NewRNSScalar creates a new Scalar value (i.e., a degree-0 polynomial) in the RingQP.
func (r Ring) NewRNSScalar() ring.RNSScalar {
	modlen := r.RingQ.ModuliChainLength()
	if r.RingP != nil {
		modlen += r.RingP.ModuliChainLength()
	}
	return make(ring.RNSScalar, modlen)
}

// NewRNSScalarFromUInt64 creates a new Scalar in the RingQP initialized with value v.
func (r Ring) NewRNSScalarFromUInt64(v uint64) ring.RNSScalar {
	var scalarQ, scalarP []uint64
	if r.RingQ != nil {
		scalarQ = r.RingQ.NewRNSScalarFromUInt64(v)
	}
	if r.RingP != nil {
		scalarP = r.RingP.NewRNSScalarFromUInt64(v)
	}
	return append(scalarQ, scalarP...)
}

// SubRNSScalar subtracts s2 to s1 and stores the result in sout.
func (r Ring) SubRNSScalar(s1, s2, sout ring.RNSScalar) {
	qlen := r.RingQ.ModuliChainLength()
	if r.RingQ != nil {
		r.RingQ.SubRNSScalar(s1[:qlen], s2[:qlen], sout[:qlen])
	}
	if r.RingP != nil {
		r.RingP.SubRNSScalar(s1[qlen:], s2[qlen:], sout[qlen:])

	}
}

// MulRNSScalar multiplies s1 and s2 and stores the result in sout.
func (r Ring) MulRNSScalar(s1, s2, sout ring.RNSScalar) {
	qlen := r.RingQ.ModuliChainLength()
	if r.RingQ != nil {
		r.RingQ.MulRNSScalar(s1[:qlen], s2[:qlen], sout[:qlen])
	}
	if r.RingP != nil {
		r.RingP.MulRNSScalar(s1[qlen:], s2[qlen:], sout[qlen:])
	}
}

// EvalPolyScalar evaluate the polynomial pol at pt and writes the result in p3
func (r Ring) EvalPolyScalar(pol []Poly, pt uint64, p3 Poly) {
	polQ, polP := make([]ring.Poly, len(pol)), make([]ring.Poly, len(pol))
	for i, coeff := range pol {
		polQ[i] = coeff.Q
		polP[i] = coeff.P
	}
	r.RingQ.EvalPolyScalar(polQ, pt, p3.Q)
	if r.RingP != nil {
		r.RingP.EvalPolyScalar(polP, pt, p3.P)
	}
}

// MulScalar multiplies p1 by scalar and returns the result in p2.
func (r Ring) MulScalar(p1 Poly, scalar uint64, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.MulScalar(p1.Q, scalar, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.MulScalar(p1.P, scalar, p2.P)
	}
}

// NTT computes the NTT of p1 and returns the result on p2.
func (r Ring) NTT(p1, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.NTT(p1.Q, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.NTT(p1.P, p2.P)
	}
}

// INTT computes the inverse-NTT of p1 and returns the result on p2.
func (r Ring) INTT(p1, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.INTT(p1.Q, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.INTT(p1.P, p2.P)
	}
}

// NTTLazy computes the NTT of p1 and returns the result on p2.
// Output values are in the range [0, 2q-1].
func (r Ring) NTTLazy(p1, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.NTTLazy(p1.Q, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.NTTLazy(p1.P, p2.P)
	}
}

// INTTLazy computes the inverse-NTT of p1 and returns the result on p2.
// Output values are in the range [0, 2q-1].
func (r Ring) INTTLazy(p1, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.INTTLazy(p1.Q, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.INTTLazy(p1.P, p2.P)
	}
}

// MForm switches p1 to the Montgomery domain and writes the result on p2.
func (r Ring) MForm(p1, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.MForm(p1.Q, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.MForm(p1.P, p2.P)
	}
}

// IMForm switches back p1 from the Montgomery domain to the conventional domain and writes the result on p2.
func (r Ring) IMForm(p1, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.IMForm(p1.Q, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.IMForm(p1.P, p2.P)
	}
}

// MulCoeffsMontgomery multiplies p1 by p2 coefficient-wise with a Montgomery modular reduction.
func (r Ring) MulCoeffsMontgomery(p1, p2, p3 Poly) {
	if r.RingQ != nil {
		r.RingQ.MulCoeffsMontgomery(p1.Q, p2.Q, p3.Q)
	}
	if r.RingP != nil {
		r.RingP.MulCoeffsMontgomery(p1.P, p2.P, p3.P)
	}
}

// MulCoeffsMontgomeryLazy multiplies p1 by p2 coefficient-wise with a constant-time Montgomery modular reduction.
// Result is within [0, 2q-1].
func (r Ring) MulCoeffsMontgomeryLazy(p1, p2, p3 Poly) {
	if r.RingQ != nil {
		r.RingQ.MulCoeffsMontgomeryLazy(p1.Q, p2.Q, p3.Q)
	}
	if r.RingP != nil {
		r.RingP.MulCoeffsMontgomeryLazy(p1.P, p2.P, p3.P)
	}
}

// MulCoeffsMontgomeryLazyThenAddLazy multiplies p1 by p2 coefficient-wise with a
// constant-time Montgomery modular reduction and adds the result on p3.
// Result is within [0, 2q-1]
func (r Ring) MulCoeffsMontgomeryLazyThenAddLazy(p1, p2, p3 Poly) {
	if r.RingQ != nil {
		r.RingQ.MulCoeffsMontgomeryLazyThenAddLazy(p1.Q, p2.Q, p3.Q)
	}
	if r.RingP != nil {
		r.RingP.MulCoeffsMontgomeryLazyThenAddLazy(p1.P, p2.P, p3.P)
	}
}

// MulCoeffsMontgomeryThenSub multiplies p1 by p2 coefficient-wise with
// a Montgomery modular reduction and subtracts the result from p3.
func (r Ring) MulCoeffsMontgomeryThenSub(p1, p2, p3 Poly) {
	if r.RingQ != nil {
		r.RingQ.MulCoeffsMontgomeryThenSub(p1.Q, p2.Q, p3.Q)
	}
	if r.RingP != nil {
		r.RingP.MulCoeffsMontgomeryThenSub(p1.P, p2.P, p3.P)
	}
}

// MulCoeffsMontgomeryLazyThenSubLazy multiplies p1 by p2 coefficient-wise with
// a Montgomery modular reduction and subtracts the result from p3.
func (r Ring) MulCoeffsMontgomeryLazyThenSubLazy(p1, p2, p3 Poly) {
	if r.RingQ != nil {
		r.RingQ.MulCoeffsMontgomeryLazyThenSubLazy(p1.Q, p2.Q, p3.Q)
	}
	if r.RingP != nil {
		r.RingP.MulCoeffsMontgomeryLazyThenSubLazy(p1.P, p2.P, p3.P)
	}
}

// MulCoeffsMontgomeryThenAdd multiplies p1 by p2 coefficient-wise with a
// Montgomery modular reduction and adds the result to p3.
func (r Ring) MulCoeffsMontgomeryThenAdd(p1, p2, p3 Poly) {
	if r.RingQ != nil {
		r.RingQ.MulCoeffsMontgomeryThenAdd(p1.Q, p2.Q, p3.Q)
	}
	if r.RingP != nil {
		r.RingP.MulCoeffsMontgomeryThenAdd(p1.P, p2.P, p3.P)
	}
}

// MulRNSScalarMontgomery multiplies p with a scalar value expressed in the CRT decomposition.
// It assumes the scalar decomposition to be in Montgomery form.
func (r Ring) MulRNSScalarMontgomery(p Poly, scalar []uint64, pOut Poly) {
	scalarQ, scalarP := scalar[:r.RingQ.ModuliChainLength()], scalar[r.RingQ.ModuliChainLength():]
	if r.RingQ != nil {
		r.RingQ.MulRNSScalarMontgomery(p.Q, scalarQ, pOut.Q)
	}
	if r.RingP != nil {
		r.RingP.MulRNSScalarMontgomery(p.P, scalarP, pOut.P)
	}
}

// Inverse computes the modular inverse of a scalar a expressed in a CRT decomposition.
// The inversion is done in-place and assumes that a is in Montgomery form.
func (r Ring) Inverse(scalar ring.RNSScalar) {
	scalarQ, scalarP := scalar[:r.RingQ.ModuliChainLength()], scalar[r.RingQ.ModuliChainLength():]
	if r.RingQ != nil {
		r.RingQ.Inverse(scalarQ)
	}
	if r.RingP != nil {
		r.RingP.Inverse(scalarP)
	}
}

// Reduce applies the modular reduction on the coefficients of p1 and returns the result on p2.
func (r Ring) Reduce(p1, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.Reduce(p1.Q, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.Reduce(p1.P, p2.P)
	}
}

// Automorphism applies the automorphism X^{i} -> X^{i*gen} on p1 and writes the result on p2.
// Method is not in place.
func (r Ring) Automorphism(p1 Poly, galEl uint64, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.Automorphism(p1.Q, galEl, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.Automorphism(p1.P, galEl, p2.P)
	}
}

// AutomorphismNTT applies the automorphism X^{i} -> X^{i*gen} on p1 and writes the result on p2.
// Method is not in place.
// Inputs are assumed to be in the NTT domain.
func (r Ring) AutomorphismNTT(p1 Poly, galEl uint64, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.AutomorphismNTT(p1.Q, galEl, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.AutomorphismNTT(p1.P, galEl, p2.P)
	}
}

// AutomorphismNTTWithIndex applies the automorphism X^{i} -> X^{i*gen} on p1 and writes the result on p2.
// Index of automorphism must be provided.
// Method is not in place.
func (r Ring) AutomorphismNTTWithIndex(p1 Poly, index []uint64, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.AutomorphismNTTWithIndex(p1.Q, index, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.AutomorphismNTTWithIndex(p1.P, index, p2.P)
	}
}

// AutomorphismNTTWithIndexThenAddLazy applies the automorphism X^{i} -> X^{i*gen} on p1 and adds the result on p2.
// Index of automorphism must be provided.
// Method is not in place.
func (r Ring) AutomorphismNTTWithIndexThenAddLazy(p1 Poly, index []uint64, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.AutomorphismNTTWithIndexThenAddLazy(p1.Q, index, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.AutomorphismNTTWithIndexThenAddLazy(p1.P, index, p2.P)
	}
}

// ExtendBasisSmallNormAndCenter extends a small-norm polynomial polQ in R_Q to a polynomial
// polQP in R_QP.
func (r Ring) ExtendBasisSmallNormAndCenter(polyInQ ring.Poly, levelP int, polyOutQ, polyOutP ring.Poly) {
	var coeff, Q, QHalf, sign uint64
	Q = r.RingQ.SubRings[0].Modulus
	QHalf = Q >> 1

	if !polyInQ.Equal(&polyOutQ) {
		polyOutQ.Copy(polyInQ)
	}

	P := r.RingP.ModuliChain()
	N := r.RingQ.N()

	for j := 0; j < N; j++ {

		coeff = polyInQ.Coeffs[0][j]

		sign = 1
		if coeff > QHalf {
			coeff = Q - coeff
			sign = 0
		}

		for i, pi := range P[:levelP+1] {
			polyOutP.Coeffs[i][j] = (coeff * sign) | (pi-coeff)*(sign^1)
		}
	}
}
