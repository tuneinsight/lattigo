package ring

import (
	"math/big"
)

// RNSScalar represents a scalar value in the Ring (i.e., a degree-0 polynomial) in RNS form.
type RNSScalar []uint64

// NewRNSScalar creates a new Scalar value.
func (r *Ring) NewRNSScalar() RNSScalar {
	return make(RNSScalar, r.level+1)
}

// NewRNSScalarFromUInt64 creates a new Scalar initialized with value v.
func (r *Ring) NewRNSScalarFromUInt64(v uint64) (rns RNSScalar) {
	rns = make(RNSScalar, r.level+1)
	for i, s := range r.SubRings[:r.level+1] {
		rns[i] = v % s.Modulus
	}
	return rns
}

// NewRNSScalarFromBigint creates a new Scalar initialized with value v.
func (r *Ring) NewRNSScalarFromBigint(v *big.Int) (rns RNSScalar) {
	rns = make(RNSScalar, r.level+1)
	tmp0 := new(big.Int)
	tmp1 := new(big.Int)
	for i, s := range r.SubRings[:r.level+1] {
		rns[i] = tmp0.Mod(v, tmp1.SetUint64(s.Modulus)).Uint64()
	}
	return rns
}

// MFormRNSScalar switches an RNS scalar to the Montgomery domain.
// s2 = s1<<64 mod Q
func (r *Ring) MFormRNSScalar(s1, s2 RNSScalar) {
	for i, s := range r.SubRings[:r.level+1] {
		s2[i] = MForm(s1[i], s.Modulus, s.BRedConstant)
	}
}

// NegRNSScalar evaluates s2 = -s1.
func (r *Ring) NegRNSScalar(s1, s2 RNSScalar) {
	for i, s := range r.SubRings[:r.level+1] {
		s2[i] = s.Modulus - s1[i]
	}
}

// SubRNSScalar subtracts s2 to s1 and stores the result in sout.
func (r *Ring) SubRNSScalar(s1, s2, sout RNSScalar) {
	for i, s := range r.SubRings[:r.level+1] {
		if s2[i] > s1[i] {
			sout[i] = s1[i] + s.Modulus - s2[i]
		} else {
			sout[i] = s1[i] - s2[i]
		}
	}
}

// MulRNSScalar multiplies s1 and s2 and stores the result in sout.
// Multiplication is operated with Montgomery.
func (r *Ring) MulRNSScalar(s1, s2, sout RNSScalar) {
	for i, s := range r.SubRings[:r.level+1] {
		sout[i] = MRedLazy(s1[i], s2[i], s.Modulus, s.MRedConstant)
	}
}

// Inverse computes the modular inverse of a scalar a expressed in a CRT decomposition.
// The inversion is done in-place and assumes that a is in Montgomery form.
func (r *Ring) Inverse(a RNSScalar) {
	for i, s := range r.SubRings[:r.level+1] {
		a[i] = ModexpMontgomery(a[i], int(s.Modulus-2), s.Modulus, s.MRedConstant, s.BRedConstant)
	}
}
