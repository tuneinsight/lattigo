package ring

import (
	"math/big"
	"math/bits"
)

// MForm switches a to the Montgomery domain by computing
// a*2^64 mod q.
func MForm(a, q uint64, u []uint64) (r uint64) {
	mhi, _ := bits.Mul64(a, u[1])
	r = -(a*u[0] + mhi) * q
	if r >= q {
		r -= q
	}
	return
}

// MForm switches a to the Montgomery domain by computing
// a*2^64 mod q in constant time.
// The result is between 0 and 2*q-1.
func MFormConstant(a, q uint64, u []uint64) (r uint64) {
	mhi, _ := bits.Mul64(a, u[1])
	r = -(a*u[0] + mhi) * q
	return
}

// InvMForm switches a from the Montgomery domain back to the
// standard domain by computing a*(1/2^64) mod q.
func InvMForm(a, q, qInv uint64) (r uint64) {
	r, _ = bits.Mul64(a*qInv, q)
	r = q - r
	if r >= q {
		r -= q
	}
	return
}

// InvMForm switches a from the Montgomery domain back to the
// standard domain by computing a*(1/2^64) mod q in constant time.
// The result is between 0 and 2*q-1.
func InvMFormConstant(a, q, qInv uint64) (r uint64) {
	r, _ = bits.Mul64(a*qInv, q)
	r = q - r
	return
}

// MRedParams computes the parameter qInv = (q^-1) mod 2^64,
// required for MRed.
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

// MRed computes x * y * (1/2^64) mod q.
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

// MRed computes x * y * (1/2^64) mod q in constant time.
// The result is between 0 and 2*q-1.
func MRedConstant(x, y, q, qInv uint64) (r uint64) {
	ahi, alo := bits.Mul64(x, y)
	R := alo * qInv
	H, _ := bits.Mul64(R, q)
	r = ahi - H + q
	return
}

// BRedParams computes the parameters for the BRed algorithm.
// Returns ((2^128)/q)/(2^64) and (2^128)/q mod 2^64.
func BRedParams(q uint64) (params []uint64) {
	bigR := new(big.Int).Lsh(NewUint(1), 128)
	bigR.Quo(bigR, NewUint(q))

	// 2^radix // q
	mhi := new(big.Int).Rsh(bigR, 64).Uint64()
	mlo := bigR.Uint64()

	return []uint64{mhi, mlo}
}

// BRedAdd computes a mod q.
func BRedAdd(x, q uint64, u []uint64) (r uint64) {
	s0, _ := bits.Mul64(x, u[0])
	r = x - s0*q
	if r >= q {
		r -= q
	}
	return
}

// BRedAdd computes a mod q in constant time.
// The result is between 0 and 2*q-1.
func BRedAddConstant(x, q uint64, u []uint64) uint64 {
	s0, _ := bits.Mul64(x, u[0])
	return x - s0*q
}

// BRed compute x*y mod q.
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

// BRed compute x*y mod q in constant time.
// The result is between 0 and 2*q-1.
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

// CRed reduce returns a mod q where a is between 0 and 2*q-1.
func CRed(a, q uint64) uint64 {
	if a >= q {
		return a - q
	}
	return a
}
