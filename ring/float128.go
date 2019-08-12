// ==== Initial package ====
// This package implements 128-bit ("double double") floating point using
// a pair of 64-bit hardware floating point values and standard hardware
// floating point operations. It is based directly on libqd by Yozo Hida,
// Xiaoye S. Li, David H. Bailey, Yves Renard and E. Jason Riedy. Source:
// http://crd.lbl.gov/~dhbailey/mpdist/qd-2.3.13.tar.gz
// ==== Current package ====
// It has been modified and adapted its required use in the BFV. Unused
// functions have been removed, and some functions simplified to save
// computation (full precision is not needed).
package ring

// TODO : implement e^n from https://golang.org/src/math/exp.go?s=368:395#L4 and https://golang.org/src/math/ldexp.go

import (
	"math"
)

// A Float128 represents a double-double floating point number with
// 106 bits of mantissa, or about 32 decimal digits. The zero value
// for a Float128 represents the value 0.
type Float128 [2]float64 // float128 represented by two float64s

//
// SET/GET
//

// Set 128-bit floating point object to 0.0
func Float128SetZero() (result Float128) {
	result[0] = 0.0
	result[1] = 0.0
	return
}

// Set 128-bit floating point object from an uint64.
// To be used only if the integer is 53 or less bits
// and if the computation will stay under 54 bits.
func Float128SetUint53(i uint64) (result Float128) {
	result[0] = float64(i)
	result[1] = 0.0
	return
}

// Set 128-bit floating point object from uint64
// Allows to import integers bigger than 53 bits and
// do computation with result of up to 64 bits.
func Float128SetUint64(i uint64) (result Float128) {
	result[0] = float64(i >> 12)
	result[1] = float64(i&0xfff) / float64(4096)
	return
}

// Set 128-bit floating point object from uint64
// Allows to import integers bigger than 53 bits and
// do computation with result of up to 64 bits.
func Float128SetInt64(i int64) (result Float128) {
	sign := i < 0
	if sign {
		i *= -1
	}
	result[0] = float64(uint64(i) >> 12)
	result[1] = float64(uint64(i)&0xfff) / float64(4096)

	if sign {
		result[0] *= -1
		result[1] *= -1
	}
	return
}

// Get the uint64 value from 128-bit floating point object
// rounded down, to be used if the result of the computation
// is under 53 bits.
func Float128ToUint53(f Float128) uint64 {
	return uint64(f[0])
}

// Reconstruct an uint64 integer from a Float128 that was imported with the Float128SetUint64() function.
// Isolates the integer part from the floatting point part, then adds the rounded floating point
// part to the integer part (this prevents occasional rounding errors when the second floating element is negative)
func Float128ToUint64(f Float128) uint64 {
	return uint64(f[0]*4096) + uint64(math.Round((f[0]*4096)-float64(uint64(f[0]*4096))+f[1]*4096))
}

//
// ADDITION
//

// Compute fl(a+b) and err(a+b).
func twoSum(a, b float64) (s, err float64) {
	s = a + b
	bb := s - a
	err = (a - (s - bb)) + (b - bb)
	return
}

// Compute fl(a+b) and err(a+b).  Assumes |a| >= |b|.
func quickTwoSum(a, b float64) (s, err float64) {
	s = a + b
	err = b - (s - a)
	return
}

// Compute D = D + D
func Float128Add(a, b Float128) (f Float128) {
	s1, s2 := twoSum(a[0], b[0])
	t1, t2 := twoSum(a[1], b[1])
	s2 += t1
	s1, s2 = quickTwoSum(s1, s2)
	s2 += t2
	f[0], f[1] = quickTwoSum(s1, s2)
	return
}

//
// SUBTRACTION
//

// Compute fl(a-b) and err(a-b).
func twoDiff(a, b float64) (s, err float64) {
	s = a - b
	bb := s - a
	err = (a - (s - bb)) - (b + bb)
	return
}

// Compute D = D - D
func Float128Sub(a, b Float128) (f Float128) {
	s1, s2 := twoDiff(a[0], b[0])
	t1, t2 := twoDiff(a[1], b[1])
	s2 += t1
	s1, s2 = quickTwoSum(s1, s2)
	s2 += t2
	f[0], f[1] = quickTwoSum(s1, s2)

	return
}

//
// MULTIPLICATION
//

// Compute high and lo words of a float64 value
func split(a float64) (hi, lo float64) {
	temp := 134217729.0 * a
	hi = temp - (temp - a)
	lo = a - hi
	return
}

// Compute fl(a*b) and err(a*b).
func twoProd(a, b float64) (p, err float64) {
	p = a * b
	aHi, aLo := split(a)
	bHi, bLo := split(b)
	err = ((aHi*bHi - p) + aHi*bLo + aLo*bHi) + aLo*bLo
	return
}

// Compute D = D * D
func Float128Mul(a, b Float128) (f Float128) {
	p1, p2 := twoProd(a[0], b[0])
	p2 += a[0]*b[1] + a[1]*b[0]
	f[0], f[1] = quickTwoSum(p1, p2)
	return
}

//
// DIVISION
//

// Compute D = D / D
func Float128Div(a, b Float128) Float128 {

	var q1, p1, p2, p3, p4, v1, v2, r, t0, t1 float64
	var f Float128

	//var q2, q3, p3, p4, v3, v4, r1, r2 float64
	//var f2, f3, t2 Float128

	q1 = a[0] / b[0] //0.3ns

	// MUL ~10ns
	p1, p2 = twoProd(q1, b[0])
	p2 += q1 * b[1]
	t0 = p1 + p2
	t1 = p2 - (t0 - p1)
	//

	// SUB ~10ns
	p3, p4 = twoDiff(a[0], t0)
	v1, v2 = twoDiff(a[1], t1)
	p4 += v1
	p3, p4 = quickTwoSum(p3, p4)
	p4 += v2

	r = (p3 + p4) / b[0] //r1

	// ======= Additional Precision ========

	// Not required for its current usage, so it is removed to improve the speed

	//r1[0], r1[1] = quickTwoSum(p3, p4)

	// MUL - SUB
	//q2 = r1[0] / b[0] //0.3ns
	//f2[0], f2[1] = q2, 0
	//p1, p2 = twoProd(f2[0], b[0])
	//p2 += f2[0]*b[1] + f2[1] * b[0]
	//t2[0], t2[1] = quickTwoSum(p1, p2)

	// SUB
	//r = Float128Sub(r, t2)
	//p3, p4 = twoDiff(r1[0], t2[0])
	//v3, v4 = twoDiff(r1[1], t2[1])
	//p4 += v3
	//p3, p4 = quickTwoSum(p3, p4)
	//p2 += v4
	//r2, _ = quickTwoSum(p3, p4)
	//

	//q3 = r2/b[0]
	//f3[0], f3[1] = q3, 0

	// ======= ================== ========

	// quickTwoSum(q1, r)
	f[0] = q1 + r
	f[1] = r - (f[0] - q1)

	//f = Float128Add(f, f3)

	return f
}
