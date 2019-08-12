package ring

import (
	"math/bits"
)

//============================
//=== MONTGOMERY REDUCTION ===
//============================

// MForm returns a*2^64 mod q.
func MForm(a, q uint64, u []uint64) (r uint64) {
	mhi, _ := bits.Mul64(a, u[1])
	r = -(a*u[0] + mhi) * q
	if r >= q {
		r -= q
	}
	return
}

// MFormConstant is identical to MForm, except that it runs in constant time
// and returns a value in [0, 2q-1]
func MFormConstant(a, q uint64, u []uint64) (r uint64) {
	mhi, _ := bits.Mul64(a, u[1])
	r = -(a*u[0] + mhi) * q
	return
}

// InvMForm returns a*(1/2^64) mod q.
func InvMForm(a, q, qInv uint64) (r uint64) {
	r, _ = bits.Mul64(a*qInv, q)
	r = q - r
	if r >= q {
		r -= q
	}
	return
}

// InvMFormConstan is indentical to InvMForm, except that it runs in constant time
// and returns a value in [0, 2q-1]
func InvMFormConstant(a, q, qInv uint64) (r uint64) {
	r, _ = bits.Mul64(a*qInv, q)
	r = q - r
	return
}

// MRedParams computes the parameter qInv = (q^-1) mod 2^64,
// required for MontgomeryReduce.
func MRedParams(q uint64) (qInv uint64) {
	var x uint64
	qInv = 1
	x = q
	for i := 0; i < 63; i++ {
		qInv *= x
		qInv &= 0xFFFFFFFFFFFFFFFF
		x *= x
		x &= 0xFFFFFFFFFFFFFFFF
	}
	return
}

// MRed operates a 64x64 bit multiplication with
// a montgomery reduction over a radix of 2^64.
func MRed(x, y, q, qInv uint64) (r uint64) {
	ahi, alo := bits.Mul64(x, y)
	R := alo * qInv
	H, _ := bits.Mul64(R, q)
	r = ahi - H + q
	if r >= q {
		r -= q
	}
	return
}

// MRedConstant is identical to MRed except it runs in
// constant time and returns a value in [0, 2q-1].
func MRedConstant(x, y, q, qInv uint64) (r uint64) {
	ahi, alo := bits.Mul64(x, y)
	R := alo * qInv
	H, _ := bits.Mul64(R, q)
	r = ahi - H + q
	return
}

//==========================
//=== BARRETT REDUCTION  ===
//==========================

// BRedParams computes the parameters required for the Barret reduction with
// a radix of 2^128.
func BRedParams(q uint64) (params []uint64) {
	bigR := new(Int).Lsh(NewUint(1), 128)
	bigR.Div(bigR, NewUint(q))

	// 2^radix // q
	mhi := new(Int).Rsh(bigR, 64).Uint64()
	mlo := bigR.Uint64()

	return []uint64{mhi, mlo}
}

// BRedAdd reduces a 64 bit integer by q.
// Assumes that x <= 64bits.
func BRedAdd(x, q uint64, u []uint64) (r uint64) {
	s0, _ := bits.Mul64(x, u[0])
	r = x - s0*q
	if r >= q {
		r -= q
	}
	return
}

// BRedAddConstant is indentical to BReAdd, except it runs
// in constant time and returns a value in [0, 2q-1]
func BRedAddConstant(x, q uint64, u []uint64) uint64 {
	s0, _ := bits.Mul64(x, u[0])
	return x - s0*q
}

// BRed operates a 64x64 bit multiplication with
// a barrett reduction.
func BRed(x, y, q uint64, u []uint64) (r uint64) {

	var lhi, mhi, mlo, s0, s1, carry uint64

	ahi, alo := bits.Mul64(x, y)

	// (alo*ulo)>>64

	lhi, _ = bits.Mul64(alo, u[1])

	// ((ahi*ulo + alo*uhi) + (alo*ulo))>>64

	mhi, mlo = bits.Mul64(alo, u[0])

	s0, carry = bits.Add64(mlo, lhi, 0)

	s1 = mhi + carry

	mhi, mlo = bits.Mul64(ahi, u[1])

	_, carry = bits.Add64(mlo, s0, 0)

	lhi = mhi + carry

	// (ahi*uhi) + (((ahi*ulo + alo*uhi) + (alo*ulo))>>64)

	s0 = ahi*u[0] + s1 + lhi

	r = alo - s0*q

	if r >= q {
		r -= q
	}

	return
}

// BRedConstant is indentical to BRed, except it runs
// in constant time and returns a value in [0, 2q-1]
func BRedConstant(x, y, q uint64, u []uint64) (r uint64) {

	var lhi, mhi, mlo, s0, s1, carry uint64

	ahi, alo := bits.Mul64(x, y)

	// alo*ulo

	lhi, _ = bits.Mul64(alo, u[1])

	// ahi*ulo + alo*uhi

	mhi, mlo = bits.Mul64(alo, u[0])

	s0, carry = bits.Add64(mlo, lhi, 0)

	s1 = mhi + carry

	mhi, mlo = bits.Mul64(ahi, u[1])

	_, carry = bits.Add64(mlo, s0, 0)

	lhi = mhi + carry

	// ahi*uhi

	s0 = ahi*u[0] + s1 + lhi

	r = alo - s0*q

	return
}

//===============================
//==== CONDITIONAL REDUCTION ====
//===============================

// CRed reduce returns a mod q, where,
// a is required to be in the range [0, 2q-1].
func CRed(a, q uint64) uint64 {
	if a >= q {
		return a - q
	}
	return a
}
