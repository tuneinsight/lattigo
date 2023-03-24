// Package ringqp is implements a wrapper for both the ringQ and ringP.
package ringqp

import (
	"bufio"
	"io"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils/buffer"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
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

// AtLevel returns a shallow copy of the target ring configured to
// carry on operations at the specified levels.
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

// LevelQ returns the level at which the target
// ring operates for the modulus Q.
func (r *Ring) LevelQ() int {
	if r.RingQ != nil {
		return r.RingQ.Level()
	}

	return -1
}

// LevelP returns the level at which the target
// ring operates for the modulus P.
func (r *Ring) LevelP() int {
	if r.RingP != nil {
		return r.RingP.Level()
	}

	return -1
}

func (r *Ring) Equal(p1, p2 Poly) (v bool) {
	v = true
	if r.RingQ != nil {
		v = v && r.RingQ.Equal(p1.Q, p2.Q)
	}

	if r.RingP != nil {
		v = v && r.RingP.Equal(p1.P, p2.P)
	}

	return
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

// Add adds p1 to p2 coefficient-wise and writes the result on p3.
func (r *Ring) Add(p1, p2, p3 Poly) {
	if r.RingQ != nil {
		r.RingQ.Add(p1.Q, p2.Q, p3.Q)
	}
	if r.RingP != nil {
		r.RingP.Add(p1.P, p2.P, p3.P)
	}
}

// AddLazy adds p1 to p2 coefficient-wise and writes the result on p3 without modular reduction.
func (r *Ring) AddLazy(p1, p2, p3 Poly) {
	if r.RingQ != nil {
		r.RingQ.AddLazy(p1.Q, p2.Q, p3.Q)
	}
	if r.RingP != nil {
		r.RingP.AddLazy(p1.P, p2.P, p3.P)
	}
}

// Sub subtracts p2 to p1 coefficient-wise and writes the result on p3.
func (r *Ring) Sub(p1, p2, p3 Poly) {
	if r.RingQ != nil {
		r.RingQ.Sub(p1.Q, p2.Q, p3.Q)
	}
	if r.RingP != nil {
		r.RingP.Sub(p1.P, p2.P, p3.P)
	}
}

// Neg negates p1 coefficient-wise and writes the result on p2.
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
	modlen := r.RingQ.ModuliChainLength()
	if r.RingP != nil {
		modlen += r.RingP.ModuliChainLength()
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
	qlen := r.RingQ.ModuliChainLength()
	if r.RingQ != nil {
		r.RingQ.SubRNSScalar(s1[:qlen], s2[:qlen], sout[:qlen])
	}
	if r.RingP != nil {
		r.RingP.SubRNSScalar(s1[qlen:], s2[qlen:], sout[qlen:])

	}
}

// MulRNSScalar multiplies s1 and s2 and stores the result in sout.
func (r *Ring) MulRNSScalar(s1, s2, sout ring.RNSScalar) {
	qlen := r.RingQ.ModuliChainLength()
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

// MulScalar multiplies p1 by scalar and returns the result in p2.
func (r *Ring) MulScalar(p1 Poly, scalar uint64, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.MulScalar(p1.Q, scalar, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.MulScalar(p1.P, scalar, p2.P)
	}
}

// NTT computes the NTT of p1 and returns the result on p2.
func (r *Ring) NTT(p1, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.NTT(p1.Q, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.NTT(p1.P, p2.P)
	}
}

// INTT computes the inverse-NTT of p1 and returns the result on p2.
func (r *Ring) INTT(p1, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.INTT(p1.Q, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.INTT(p1.P, p2.P)
	}
}

// NTTLazy computes the NTT of p1 and returns the result on p2.
// Output values are in the range [0, 2q-1].
func (r *Ring) NTTLazy(p1, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.NTTLazy(p1.Q, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.NTTLazy(p1.P, p2.P)
	}
}

// INTTLazy computes the inverse-NTT of p1 and returns the result on p2.
// Output values are in the range [0, 2q-1].
func (r *Ring) INTTLazy(p1, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.INTTLazy(p1.Q, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.INTTLazy(p1.P, p2.P)
	}
}

// MForm switches p1 to the Montgomery domain and writes the result on p2.
func (r *Ring) MForm(p1, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.MForm(p1.Q, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.MForm(p1.P, p2.P)
	}
}

// IMForm switches back p1 from the Montgomery domain to the conventional domain and writes the result on p2.
func (r *Ring) IMForm(p1, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.IMForm(p1.Q, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.IMForm(p1.P, p2.P)
	}
}

// MulCoeffsMontgomery multiplies p1 by p2 coefficient-wise with a Montgomery modular reduction.
func (r *Ring) MulCoeffsMontgomery(p1, p2, p3 Poly) {
	if r.RingQ != nil {
		r.RingQ.MulCoeffsMontgomery(p1.Q, p2.Q, p3.Q)
	}
	if r.RingP != nil {
		r.RingP.MulCoeffsMontgomery(p1.P, p2.P, p3.P)
	}
}

// MulCoeffsMontgomeryLazy multiplies p1 by p2 coefficient-wise with a constant-time Montgomery modular reduction.
// Result is within [0, 2q-1].
func (r *Ring) MulCoeffsMontgomeryLazy(p1, p2, p3 Poly) {
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
func (r *Ring) MulCoeffsMontgomeryLazyThenAddLazy(p1, p2, p3 Poly) {
	if r.RingQ != nil {
		r.RingQ.MulCoeffsMontgomeryLazyThenAddLazy(p1.Q, p2.Q, p3.Q)
	}
	if r.RingP != nil {
		r.RingP.MulCoeffsMontgomeryLazyThenAddLazy(p1.P, p2.P, p3.P)
	}
}

// MulCoeffsMontgomeryThenSub multiplies p1 by p2 coefficient-wise with
// a Montgomery modular reduction and subtracts the result from p3.
func (r *Ring) MulCoeffsMontgomeryThenSub(p1, p2, p3 Poly) {
	if r.RingQ != nil {
		r.RingQ.MulCoeffsMontgomeryThenSub(p1.Q, p2.Q, p3.Q)
	}
	if r.RingP != nil {
		r.RingP.MulCoeffsMontgomeryThenSub(p1.P, p2.P, p3.P)
	}
}

// MulCoeffsMontgomeryLazyThenSubLazy multiplies p1 by p2 coefficient-wise with
// a Montgomery modular reduction and subtracts the result from p3.
func (r *Ring) MulCoeffsMontgomeryLazyThenSubLazy(p1, p2, p3 Poly) {
	if r.RingQ != nil {
		r.RingQ.MulCoeffsMontgomeryLazyThenSubLazy(p1.Q, p2.Q, p3.Q)
	}
	if r.RingP != nil {
		r.RingP.MulCoeffsMontgomeryLazyThenSubLazy(p1.P, p2.P, p3.P)
	}
}

// MulCoeffsMontgomeryThenAdd multiplies p1 by p2 coefficient-wise with a
// Montgomery modular reduction and adds the result to p3.
func (r *Ring) MulCoeffsMontgomeryThenAdd(p1, p2, p3 Poly) {
	if r.RingQ != nil {
		r.RingQ.MulCoeffsMontgomeryThenAdd(p1.Q, p2.Q, p3.Q)
	}
	if r.RingP != nil {
		r.RingP.MulCoeffsMontgomeryThenAdd(p1.P, p2.P, p3.P)
	}
}

// MulRNSScalarMontgomery multiplies p with a scalar value expressed in the CRT decomposition.
// It assumes the scalar decomposition to be in Montgomery form.
func (r *Ring) MulRNSScalarMontgomery(p Poly, scalar []uint64, pOut Poly) {
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
func (r *Ring) Inverse(scalar ring.RNSScalar) {
	scalarQ, scalarP := scalar[:r.RingQ.ModuliChainLength()], scalar[r.RingQ.ModuliChainLength():]
	if r.RingQ != nil {
		r.RingQ.Inverse(scalarQ)
	}
	if r.RingP != nil {
		r.RingP.Inverse(scalarP)
	}
}

// Reduce applies the modular reduction on the coefficients of p1 and returns the result on p2.
func (r *Ring) Reduce(p1, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.Reduce(p1.Q, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.Reduce(p1.P, p2.P)
	}
}

// Automorphism applies the automorphism X^{i} -> X^{i*gen} on p1 and writes the result on p2.
// Method is not in place.
func (r *Ring) Automorphism(p1 Poly, galEl uint64, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.Automorphism(p1.Q, galEl, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.Automorphism(p1.P, galEl, p2.P)
	}
}

// AutomorphismNTTWithIndex applies the automorphism X^{i} -> X^{i*gen} on p1 and writes the result on p2.
// Index of automorphism must be provided.
// Method is not in place.
func (r *Ring) AutomorphismNTTWithIndex(p1 Poly, index []uint64, p2 Poly) {
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
func (r *Ring) AutomorphismNTTWithIndexThenAddLazy(p1 Poly, index []uint64, p2 Poly) {
	if r.RingQ != nil {
		r.RingQ.AutomorphismNTTWithIndexThenAddLazy(p1.Q, index, p2.Q)
	}
	if r.RingP != nil {
		r.RingP.AutomorphismNTTWithIndexThenAddLazy(p1.P, index, p2.P)
	}
}

// ExtendBasisSmallNormAndCenter extends a small-norm polynomial polQ in R_Q to a polynomial
// polQP in R_QP.
func (r *Ring) ExtendBasisSmallNormAndCenter(polyInQ *ring.Poly, levelP int, polyOutQ, polyOutP *ring.Poly) {
	var coeff, Q, QHalf, sign uint64
	Q = r.RingQ.SubRings[0].Modulus
	QHalf = Q >> 1

	if polyInQ != polyOutQ && polyOutQ != nil {
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

// BinarySize returns the size in bytes that the object once marshalled into a binary form.
// Assumes that each coefficient takes 8 bytes.
func (p *Poly) BinarySize() (dataLen int) {

	dataLen = 2

	if p.Q != nil {
		dataLen += p.Q.BinarySize()
	}
	if p.P != nil {
		dataLen += p.P.BinarySize()
	}

	return
}

// WriteTo writes the object on an io.Writer.
// To ensure optimal efficiency and minimal allocations, the user is encouraged
// to provide a struct implementing the interface buffer.Writer, which defines
// a subset of the method of the bufio.Writer.
// If w is not compliant to the buffer.Writer interface, it will be wrapped in
// a new bufio.Writer.
// For additional information, see lattigo/utils/buffer/writer.go.
func (p *Poly) WriteTo(w io.Writer) (n int64, err error) {

	switch w := w.(type) {
	case buffer.Writer:

		if p.Q != nil {

			var inc int
			if inc, err = buffer.WriteUint8(w, 1); err != nil {
				return int64(n), err
			}

			n += int64(inc)

		} else {
			var inc int
			if inc, err = buffer.WriteUint8(w, 0); err != nil {
				return int64(n), err
			}

			n += int64(inc)
		}

		if p.P != nil {
			var inc int
			if inc, err = buffer.WriteUint8(w, 1); err != nil {
				return int64(n), err
			}

			n += int64(inc)
		} else {
			var inc int
			if inc, err = buffer.WriteUint8(w, 0); err != nil {
				return int64(n), err
			}

			n += int64(inc)
		}

		if p.Q != nil {
			var inc int64
			if inc, err = p.Q.WriteTo(w); err != nil {
				return n + inc, err
			}

			n += inc
		}

		if p.P != nil {
			var inc int64
			if inc, err = p.P.WriteTo(w); err != nil {
				return n + inc, err
			}

			n += inc
		}

		return n, w.Flush()

	default:
		return p.WriteTo(bufio.NewWriter(w))
	}
}

// ReadFrom reads on the object from an io.Writer.
// To ensure optimal efficiency and minimal allocations, the user is encouraged
// to provide a struct implementing the interface buffer.Reader, which defines
// a subset of the method of the bufio.Reader.
// If r is not compliant to the buffer.Reader interface, it will be wrapped in
// a new bufio.Reader.
// For additional information, see lattigo/utils/buffer/reader.go.
func (p *Poly) ReadFrom(r io.Reader) (n int64, err error) {
	switch r := r.(type) {
	case buffer.Reader:

		var hasQ, hasP uint8

		var inc int
		if inc, err = buffer.ReadUint8(r, &hasQ); err != nil {
			return n + int64(inc), err
		}

		n += int64(inc)

		if inc, err = buffer.ReadUint8(r, &hasP); err != nil {
			return n + int64(inc), err
		}

		n += int64(inc)

		if hasQ == 1 {

			if p.Q == nil {
				p.Q = new(ring.Poly)
			}

			var inc int64
			if inc, err = p.Q.ReadFrom(r); err != nil {
				return n + inc, err
			}

			n += inc
		}

		if hasP == 1 {

			if p.P == nil {
				p.P = new(ring.Poly)
			}

			var inc int64
			if inc, err = p.P.ReadFrom(r); err != nil {
				return n + inc, err
			}

			n += inc
		}

		return

	default:
		return p.ReadFrom(bufio.NewReader(r))
	}
}

// Read encodes the object into a binary form on a preallocated slice of bytes
// and returns the number of bytes written.
func (p *Poly) Read(data []byte) (n int, err error) {
	var inc int

	if p.Q != nil {
		data[0] = 1
	}

	if p.P != nil {
		data[1] = 1
	}

	n = 2

	if data[0] == 1 {
		if inc, err = p.Q.Read(data[n:]); err != nil {
			return
		}
		n += inc
	}

	if data[1] == 1 {
		if inc, err = p.P.Read(data[n:]); err != nil {
			return
		}
		n += inc
	}

	return
}

// Write decodes a slice of bytes generated by MarshalBinary or
// Read on the object and returns the number of bytes read.
func (p *Poly) Write(data []byte) (n int, err error) {

	var inc int
	n = 2

	if data[0] == 1 {

		if p.Q == nil {
			p.Q = new(ring.Poly)
		}

		if inc, err = p.Q.Write(data[n:]); err != nil {
			return
		}
		n += inc
	}

	if data[1] == 1 {

		if p.P == nil {
			p.P = new(ring.Poly)
		}

		if inc, err = p.P.Write(data[n:]); err != nil {
			return
		}
		n += inc
	}

	return
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (p *Poly) MarshalBinary() (data []byte, err error) {
	data = make([]byte, p.BinarySize())
	_, err = p.Read(data)
	return
}

// UnmarshalBinary decodes a slice of bytes generated by MarshalBinary
// or Read on the object.
func (p *Poly) UnmarshalBinary(data []byte) (err error) {
	_, err = p.Write(data)
	return err
}

// UniformSampler is a type for sampling polynomials in Ring.
type UniformSampler struct {
	samplerQ, samplerP *ring.UniformSampler
}

// NewUniformSampler instantiates a new UniformSampler from a given PRNG.
func NewUniformSampler(prng sampling.PRNG, r Ring) (s UniformSampler) {
	if r.RingQ != nil {
		s.samplerQ = ring.NewUniformSampler(prng, r.RingQ)
	}

	if r.RingP != nil {
		s.samplerP = ring.NewUniformSampler(prng, r.RingP)
	}

	return s
}

// AtLevel returns a shallow copy of the target sampler that operates at the specified levels.
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

func (s UniformSampler) WithPRNG(prng sampling.PRNG) UniformSampler {
	sp := UniformSampler{samplerQ: s.samplerQ.WithPRNG(prng)}
	if s.samplerP != nil {
		sp.samplerP = s.samplerP.WithPRNG(prng)
	}
	return sp
}
