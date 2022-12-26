// Package ringqp is implements a wrapper for both the ringQ and ringP.
package ringqp

import (
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// Poly represents a polynomial in the ring of polynomial modulo Q*P.
// This type is simply the union type between two ring.Poly, each one
// containing the modulus Q and P coefficients of that polynomial.
// The modulus Q represent the ciphertext modulus and the modulus P
// the special primes for the RNS decomposition during homomorphic
// operations involving keys.
type Poly struct {
	Q, P *ring.Poly
}

// LevelQ returns the level of the polynomial modulo Q.
// Returns -1 if the modulus Q is absent.
func (p *Poly) LevelQ() int {
	if p.Q != nil {
		return p.Q.Level()
	}
	return -1
}

// LevelP returns the level of the polynomial modulo P.
// Returns -1 if the modulus P is absent.
func (p *Poly) LevelP() int {
	if p.P != nil {
		return p.P.Level()
	}
	return -1
}

// Equals returns true if the receiver Poly is equal to the provided other Poly.
// This method checks for equality of its two sub-polynomials.
func (p *Poly) Equals(other Poly) (v bool) {

	if p == &other {
		return true
	}

	v = true
	if p.Q != nil {
		v = p.Q.Equals(other.Q)
	}
	if p.P != nil {
		v = v && p.P.Equals(other.P)
	}
	return v
}

// Copy copies the coefficients of other on the target polynomial.
// This method simply calls the Copy method for each of its sub-polynomials.
func (p *Poly) Copy(other Poly) {
	if p.Q != nil {
		copy(p.Q.Buff, other.Q.Buff)
	}

	if p.P != nil {
		copy(p.P.Buff, other.P.Buff)
	}
}

// CopyLvl copies the values of p1 on p2.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func CopyLvl(levelQ, levelP int, p1, p2 Poly) {

	if p1.Q != nil && p2.Q != nil {
		ring.CopyLvl(levelQ, p1.Q, p2.Q)
	}

	if p1.P != nil && p2.P != nil {
		ring.CopyLvl(levelP, p1.P, p2.P)
	}
}

// CopyNew creates an exact copy of the target polynomial.
func (p *Poly) CopyNew() Poly {
	if p == nil {
		return Poly{}
	}

	var Q, P *ring.Poly
	if p.Q != nil {
		Q = p.Q.CopyNew()
	}

	if p.P != nil {
		P = p.P.CopyNew()
	}

	return Poly{Q, P}
}

// Ring is a structure that implements the operation in the ring R_QP.
// This type is simply a union type between the two Ring types representing
// R_Q and R_P.
type Ring struct {
	RingQ, RingP *ring.Ring
}

func (r *Ring) AtLevel(levelQ, levelP int) *Ring {

	var ringQ, ringP *ring.Ring

	if levelQ > -1 && r.RingQ != nil {
		ringQ = r.RingQ.AtLevel(levelQ)
	}

	if levelP > -1 && r.RingP != nil {
		ringP = r.RingP.AtLevel(levelP)
	}

	return &Ring{
		RingQ: ringQ,
		RingP: ringP,
	}
}

// NewPoly creates a new polynomial with all coefficients set to 0.
func (r *Ring) NewPoly() Poly {
	var Q, P *ring.Poly
	if r.RingQ != nil {
		Q = r.RingQ.NewPoly()
	}

	if r.RingP != nil {
		P = r.RingP.NewPoly()
	}
	return Poly{Q, P}
}

// AddLvl adds p1 to p2 coefficient-wise and writes the result on p3.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *Ring) Add(p1, p2, p3 Poly) {
	if r.RingQ != nil {
		r.RingQ.Add(p1.Q, p2.Q, p3.Q)
	}
	if r.RingP != nil {
		r.RingP.Add(p1.P, p2.P, p3.P)
	}
}

// AddNoModLvl adds p1 to p2 coefficient-wise and writes the result on p3 without modular reduction.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *Ring) AddNoMod(p1, p2, p3 Poly) {
	if r.RingQ != nil {
		r.RingQ.AddNoMod(p1.Q, p2.Q, p3.Q)
	}
	if r.RingP != nil {
		r.RingP.AddNoMod(p1.P, p2.P, p3.P)
	}
}

// SubLvl subtracts p2 to p1 coefficient-wise and writes the result on p3.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *Ring) Sub(p1, p2, p3 Poly) {
	if r.RingQ != nil {
		r.RingQ.Sub(p1.Q, p2.Q, p3.Q)
	}
	if r.RingP != nil {
		r.RingP.Sub(p1.P, p2.P, p3.P)
	}
}

// NegLvl negates p1 coefficient-wise and writes the result on p2.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *Ring) Neg(p1, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.Neg(p1.Q, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.Neg(p1.P, p2.P)
	}
}

// NewRNSScalar creates a new Scalar value (i.e., a degree-0 polynomial) in the RingQP.
func (r *Ring) NewRNSScalar() ring.RNSScalar {
	modlen := r.RingQ.NbModuli()
	if r.RingP != nil {
		modlen += r.RingP.NbModuli()
	}
	return make(ring.RNSScalar, modlen)
}

// NewRNSScalarFromUInt64 creates a new Scalar in the RingQP initialized with value v.
func (r *Ring) NewRNSScalarFromUInt64(v uint64) ring.RNSScalar {
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
func (r *Ring) SubRNSScalar(s1, s2, sout ring.RNSScalar) {
	qlen := r.RingQ.NbModuli()
	if r.RingQ != nil {
		r.RingQ.SubRNSScalar(s1[:qlen], s2[:qlen], sout[:qlen])
	}
	if r.RingP != nil {
		r.RingP.SubRNSScalar(s1[qlen:], s2[qlen:], sout[qlen:])

	}
}

// MulRNSScalar multiplies s1 and s2 and stores the result in sout.
func (r *Ring) MulRNSScalar(s1, s2, sout ring.RNSScalar) {
	qlen := r.RingQ.NbModuli()
	if r.RingQ != nil {
		r.RingQ.MulRNSScalar(s1[:qlen], s2[:qlen], sout[:qlen])
	}
	if r.RingP != nil {
		r.RingP.MulRNSScalar(s1[qlen:], s2[qlen:], sout[qlen:])
	}
}

// EvalPolyScalar evaluate the polynomial pol at pt and writes the result in p3
func (r *Ring) EvalPolyScalar(pol []Poly, pt uint64, p3 Poly) {
	polQ, polP := make([]*ring.Poly, len(pol)), make([]*ring.Poly, len(pol))
	for i, coeff := range pol {
		polQ[i] = coeff.Q
		polP[i] = coeff.P
	}
	r.RingQ.EvalPolyScalar(polQ, pt, p3.Q)
	if r.RingP != nil {
		r.RingP.EvalPolyScalar(polP, pt, p3.P)
	}
}

// MulScalarLvl multiplies p1 by scalar and returns the result in p2.
func (r *Ring) MulScalar(p1 Poly, scalar uint64, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.MulScalar(p1.Q, scalar, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.MulScalar(p1.P, scalar, p2.P)
	}
}

// NTTLvl computes the NTT of p1 and returns the result on p2.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *Ring) NTT(p1, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.NTT(p1.Q, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.NTT(p1.P, p2.P)
	}
}

// InvNTTLvl computes the inverse-NTT of p1 and returns the result on p2.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *Ring) InvNTT(p1, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.InvNTT(p1.Q, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.InvNTT(p1.P, p2.P)
	}
}

// NTTLazyLvl computes the NTT of p1 and returns the result on p2.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
// Output values are in the range [0, 2q-1].
func (r *Ring) NTTLazy(p1, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.NTTLazy(p1.Q, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.NTTLazy(p1.P, p2.P)
	}
}

// MFormLvl switches p1 to the Montgomery domain and writes the result on p2.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *Ring) MForm(p1, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.MForm(p1.Q, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.MForm(p1.P, p2.P)
	}
}

// InvMFormLvl switches back p1 from the Montgomery domain to the conventional domain and writes the result on p2.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *Ring) InvMForm(p1, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.InvMForm(p1.Q, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.InvMForm(p1.P, p2.P)
	}
}

// MulCoeffsMontgomeryLvl multiplies p1 by p2 coefficient-wise with a Montgomery modular reduction.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *Ring) MulCoeffsMontgomery(p1, p2, p3 Poly) {
	if r.RingQ != nil {
		r.RingQ.MulCoeffsMontgomery(p1.Q, p2.Q, p3.Q)
	}
	if r.RingP != nil {
		r.RingP.MulCoeffsMontgomery(p1.P, p2.P, p3.P)
	}
}

// MulCoeffsMontgomeryConstantLvl multiplies p1 by p2 coefficient-wise with a constant-time Montgomery modular reduction.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
// Result is within [0, 2q-1].
func (r *Ring) MulCoeffsMontgomeryConstant(p1, p2, p3 Poly) {
	if r.RingQ != nil {
		r.RingQ.MulCoeffsMontgomeryConstant(p1.Q, p2.Q, p3.Q)
	}
	if r.RingP != nil {
		r.RingP.MulCoeffsMontgomeryConstant(p1.P, p2.P, p3.P)
	}
}

// MulCoeffsMontgomeryConstantAndAddNoModLvl multiplies p1 by p2 coefficient-wise with a
// constant-time Montgomery modular reduction and adds the result on p3.
// Result is within [0, 2q-1]
func (r *Ring) MulCoeffsMontgomeryConstantAndAddNoMod(p1, p2, p3 Poly) {
	if r.RingQ != nil {
		r.RingQ.MulCoeffsMontgomeryConstantAndAddNoMod(p1.Q, p2.Q, p3.Q)
	}
	if r.RingP != nil {
		r.RingP.MulCoeffsMontgomeryConstantAndAddNoMod(p1.P, p2.P, p3.P)
	}
}

// MulCoeffsMontgomeryAndSubLvl multiplies p1 by p2 coefficient-wise with
// a Montgomery modular reduction and subtracts the result from p3.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *Ring) MulCoeffsMontgomeryAndSub(p1, p2, p3 Poly) {
	if r.RingQ != nil {
		r.RingQ.MulCoeffsMontgomeryAndSub(p1.Q, p2.Q, p3.Q)
	}
	if r.RingP != nil {
		r.RingP.MulCoeffsMontgomeryAndSub(p1.P, p2.P, p3.P)
	}
}

// MulCoeffsMontgomeryConstantAndSubNoModLvl multiplies p1 by p2 coefficient-wise with
// a Montgomery modular reduction and subtracts the result from p3.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *Ring) MulCoeffsMontgomeryConstantAndSubNoMod(p1, p2, p3 Poly) {
	if r.RingQ != nil {
		r.RingQ.MulCoeffsMontgomeryConstantAndSubNoMod(p1.Q, p2.Q, p3.Q)
	}
	if r.RingP != nil {
		r.RingP.MulCoeffsMontgomeryConstantAndSubNoMod(p1.P, p2.P, p3.P)
	}
}

// MulCoeffsMontgomeryAndAddLvl multiplies p1 by p2 coefficient-wise with a
// Montgomery modular reduction and adds the result to p3.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *Ring) MulCoeffsMontgomeryAndAdd(p1, p2, p3 Poly) {
	if r.RingQ != nil {
		r.RingQ.MulCoeffsMontgomeryAndAdd(p1.Q, p2.Q, p3.Q)
	}
	if r.RingP != nil {
		r.RingP.MulCoeffsMontgomeryAndAdd(p1.P, p2.P, p3.P)
	}
}

// MulRNSScalarMontgomery multiplies p with a scalar value expressed in the CRT decomposition.
// It assumes the scalar decomposition to be in Montgomery form.
func (r *Ring) MulRNSScalarMontgomery(p Poly, scalar []uint64, pOut Poly) {
	scalarQ, scalarP := scalar[:r.RingQ.NbModuli()], scalar[r.RingQ.NbModuli():]
	if r.RingQ != nil {
		r.RingQ.MulRNSScalarMontgomery(p.Q, scalarQ, pOut.Q)
	}
	if r.RingP != nil {
		r.RingP.MulRNSScalarMontgomery(p.P, scalarP, pOut.P)
	}
}

// Inverse computes the modular inverse of a scalar a expressed in a CRT decomposition.
// The inversion is done in-place and assumes that a is in Montgomery form.
func (r *Ring) Inverse(scalar ring.RNSScalar) {
	scalarQ, scalarP := scalar[:r.RingQ.NbModuli()], scalar[r.RingQ.NbModuli():]
	if r.RingQ != nil {
		r.RingQ.Inverse(scalarQ)
	}
	if r.RingP != nil {
		r.RingP.Inverse(scalarP)
	}
}

// ReduceLvl applies the modular reduction on the coefficients of p1 and returns the result on p2.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *Ring) Reduce(p1, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.Reduce(p1.Q, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.Reduce(p1.P, p2.P)
	}
}

// PermuteNTTWithIndexLvl applies the automorphism X^{5^j} on p1 and writes the result on p2.
// Index of automorphism must be provided.
// Method is not in place.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *Ring) PermuteNTTWithIndex(p1 Poly, index []uint64, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.PermuteNTTWithIndex(p1.Q, index, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.PermuteNTTWithIndex(p1.P, index, p2.P)
	}
}

// PermuteNTTWithIndexAndAddNoModLvl applies the automorphism X^{5^j} on p1 and adds the result on p2.
// Index of automorphism must be provided.
// Method is not in place.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *Ring) PermuteNTTWithIndexAndAddNoMod(p1 Poly, index []uint64, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.PermuteNTTWithIndexAndAddNoMod(p1.Q, index, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.PermuteNTTWithIndexAndAddNoMod(p1.P, index, p2.P)
	}
}

// ExtendBasisSmallNormAndCenter extends a small-norm polynomial polQ in R_Q to a polynomial
// polQP in R_QP.
func (r *Ring) ExtendBasisSmallNormAndCenter(polyInQ *ring.Poly, levelP int, polyOutQ, polyOutP *ring.Poly) {
	var coeff, Q, QHalf, sign uint64
	Q = r.RingQ.Tables[0].Modulus
	QHalf = Q >> 1

	if polyInQ != polyOutQ && polyOutQ != nil {
		polyOutQ.Copy(polyInQ)
	}

	P := r.RingP.Moduli()

	for j := 0; j < r.RingQ.N(); j++ {

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

// MarshalBinarySize64 returns the length in byte of the target Poly.
// Assumes that each coefficient uses 8 bytes.
func (p *Poly) MarshalBinarySize64() (dataLen int) {

	dataLen = 2

	if p.Q != nil {
		dataLen += p.Q.MarshalBinarySize64()
	}
	if p.P != nil {
		dataLen += p.P.MarshalBinarySize64()
	}

	return
}

// Encode64 writes a Poly on the input data.
// Encodes each coefficient on 8 bytes.
func (p *Poly) Encode64(data []byte) (pt int, err error) {
	var inc int

	if p.Q != nil {
		data[0] = 1
	}

	if p.P != nil {
		data[1] = 1
	}

	pt = 2

	if data[0] == 1 {
		if inc, err = p.Q.Encode64(data[pt:]); err != nil {
			return
		}
		pt += inc
	}

	if data[1] == 1 {
		if inc, err = p.P.Encode64(data[pt:]); err != nil {
			return
		}
		pt += inc
	}

	return
}

// Decode64 decodes the input bytes on the target Poly.
// Writes on pre-allocated coefficients.
// Assumes that each coefficient is encoded on 8 bytes.
func (p *Poly) Decode64(data []byte) (pt int, err error) {

	var inc int
	pt = 2

	if data[0] == 1 {

		if p.Q == nil {
			p.Q = new(ring.Poly)
		}

		if inc, err = p.Q.Decode64(data[pt:]); err != nil {
			return
		}
		pt += inc
	}

	if data[1] == 1 {

		if p.P == nil {
			p.P = new(ring.Poly)
		}

		if inc, err = p.P.Decode64(data[pt:]); err != nil {
			return
		}
		pt += inc
	}

	return
}

func (p *Poly) MarshalBinary() ([]byte, error) {
	b := make([]byte, p.MarshalBinarySize64())
	_, err := p.Encode64(b)
	return b, err
}

func (p *Poly) UnmarshalBinary(b []byte) error {
	_, err := p.Decode64(b)
	return err
}

// UniformSampler is a type for sampling polynomials in Ring.
type UniformSampler struct {
	samplerQ, samplerP *ring.UniformSampler
}

// NewUniformSampler instantiates a new UniformSampler from a given PRNG.
func NewUniformSampler(prng utils.PRNG, r Ring) (s UniformSampler) {
	if r.RingQ != nil {
		s.samplerQ = ring.NewUniformSampler(prng, r.RingQ)
	}

	if r.RingP != nil {
		s.samplerP = ring.NewUniformSampler(prng, r.RingP)
	}

	return s
}

func (s UniformSampler) AtLevel(levelQ, levelP int) UniformSampler {

	var samplerQ, samplerP *ring.UniformSampler

	if levelQ > -1 {
		samplerQ = s.samplerQ.AtLevel(levelQ)
	}

	if levelP > -1 {
		samplerP = s.samplerP.AtLevel(levelP)
	}

	return UniformSampler{
		samplerQ: samplerQ,
		samplerP: samplerP,
	}
}

// Read samples a new polynomial in Ring and stores it into p.
func (s UniformSampler) Read(p Poly) {
	if p.Q != nil && s.samplerQ != nil {
		s.samplerQ.Read(p.Q)
	}

	if p.P != nil && s.samplerP != nil {
		s.samplerP.Read(p.P)
	}
}

func (s UniformSampler) WithPRNG(prng utils.PRNG) UniformSampler {
	sp := UniformSampler{samplerQ: s.samplerQ.WithPRNG(prng)}
	if s.samplerP != nil {
		sp.samplerP = s.samplerP.WithPRNG(prng)
	}
	return sp
}
