package ring

import (
	"unsafe"
)

// Add evaluates p3 = p1 + p2 (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) Add(p1, p2, p3 []uint64) {

	N := len(p1)
	modulus := s.Modulus

	for j := 0; j < N; j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = CRed(x[0]+y[0], modulus)
		z[1] = CRed(x[1]+y[1], modulus)
		z[2] = CRed(x[2]+y[2], modulus)
		z[3] = CRed(x[3]+y[3], modulus)
		z[4] = CRed(x[4]+y[4], modulus)
		z[5] = CRed(x[5]+y[5], modulus)
		z[6] = CRed(x[6]+y[6], modulus)
		z[7] = CRed(x[7]+y[7], modulus)
	}
}

// AddLazy evaluates p3 = p1 + p2.
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) AddLazy(p1, p2, p3 []uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = x[0] + y[0]
		z[1] = x[1] + y[1]
		z[2] = x[2] + y[2]
		z[3] = x[3] + y[3]
		z[4] = x[4] + y[4]
		z[5] = x[5] + y[5]
		z[6] = x[6] + y[6]
		z[7] = x[7] + y[7]
	}
}

// Sub evaluates p3 = p1 - p2 (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) Sub(p1, p2, p3 []uint64) {

	modulus := s.Modulus
	N := len(p1)

	for j := 0; j < N; j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = CRed((x[0]+modulus)-y[0], modulus)
		z[1] = CRed((x[1]+modulus)-y[1], modulus)
		z[2] = CRed((x[2]+modulus)-y[2], modulus)
		z[3] = CRed((x[3]+modulus)-y[3], modulus)
		z[4] = CRed((x[4]+modulus)-y[4], modulus)
		z[5] = CRed((x[5]+modulus)-y[5], modulus)
		z[6] = CRed((x[6]+modulus)-y[6], modulus)
		z[7] = CRed((x[7]+modulus)-y[7], modulus)
	}
}

// SubLazy evaluates p3 = p1 - p2.
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) SubLazy(p1, p2, p3 []uint64) {

	N := len(p1)
	modulus := s.Modulus

	for j := 0; j < N; j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = x[0] + modulus - y[0]
		z[1] = x[1] + modulus - y[1]
		z[2] = x[2] + modulus - y[2]
		z[3] = x[3] + modulus - y[3]
		z[4] = x[4] + modulus - y[4]
		z[5] = x[5] + modulus - y[5]
		z[6] = x[6] + modulus - y[6]
		z[7] = x[7] + modulus - y[7]
	}
}

// Neg evaluates p2 = -p1 (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) Neg(p1, p2 []uint64) {

	N := len(p1)
	modulus := s.Modulus

	for j := 0; j < N; j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = modulus - x[0]
		z[1] = modulus - x[1]
		z[2] = modulus - x[2]
		z[3] = modulus - x[3]
		z[4] = modulus - x[4]
		z[5] = modulus - x[5]
		z[6] = modulus - x[6]
		z[7] = modulus - x[7]
	}
}

// Reduce evaluates p2 = p1 (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) Reduce(p1, p2 []uint64) {

	N := len(p1)
	modulus := s.Modulus
	brc := s.BRedConstant

	for j := 0; j < N; j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = BRedAdd(x[0], modulus, brc)
		z[1] = BRedAdd(x[1], modulus, brc)
		z[2] = BRedAdd(x[2], modulus, brc)
		z[3] = BRedAdd(x[3], modulus, brc)
		z[4] = BRedAdd(x[4], modulus, brc)
		z[5] = BRedAdd(x[5], modulus, brc)
		z[6] = BRedAdd(x[6], modulus, brc)
		z[7] = BRedAdd(x[7], modulus, brc)
	}
}

// ReduceLazy evaluates p2 = p1 (mod modulus) with p2 in range [0, 2*modulus-1].
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) ReduceLazy(p1, p2 []uint64) {

	N := len(p1)
	modulus := s.Modulus
	brc := s.BRedConstant

	for j := 0; j < N; j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = BRedAddLazy(x[0], modulus, brc)
		z[1] = BRedAddLazy(x[1], modulus, brc)
		z[2] = BRedAddLazy(x[2], modulus, brc)
		z[3] = BRedAddLazy(x[3], modulus, brc)
		z[4] = BRedAddLazy(x[4], modulus, brc)
		z[5] = BRedAddLazy(x[5], modulus, brc)
		z[6] = BRedAddLazy(x[6], modulus, brc)
		z[7] = BRedAddLazy(x[7], modulus, brc)
	}
}

// MulCoeffsLazy evaluates p3 = p1*p2.
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsLazy(p1, p2, p3 []uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = x[0] * y[0]
		z[1] = x[1] * y[1]
		z[2] = x[2] * y[2]
		z[3] = x[3] * y[3]
		z[4] = x[4] * y[4]
		z[5] = x[5] * y[5]
		z[6] = x[6] * y[6]
		z[7] = x[7] * y[7]
	}
}

// MulCoeffsLazyThenAddLazy evaluates p3 += p1*p2.
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsLazyThenAddLazy(p1, p2, p3 []uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] += x[0] * y[0]
		z[1] += x[1] * y[1]
		z[2] += x[2] * y[2]
		z[3] += x[3] * y[3]
		z[4] += x[4] * y[4]
		z[5] += x[5] * y[5]
		z[6] += x[6] * y[6]
		z[7] += x[7] * y[7]
	}
}

// MulCoeffs evaluates p3 = p1*p2 (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffs(p1, p2, p3 []uint64) {

	N := len(p1)
	modulus := s.Modulus
	brc := s.BRedConstant

	for j := 0; j < N; j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = BRed(x[0], y[0], modulus, brc)
		z[1] = BRed(x[1], y[1], modulus, brc)
		z[2] = BRed(x[2], y[2], modulus, brc)
		z[3] = BRed(x[3], y[3], modulus, brc)
		z[4] = BRed(x[4], y[4], modulus, brc)
		z[5] = BRed(x[5], y[5], modulus, brc)
		z[6] = BRed(x[6], y[6], modulus, brc)
		z[7] = BRed(x[7], y[7], modulus, brc)
	}
}

// MulCoeffsThenAdd evaluates p3 = p3 + (p1*p2) (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsThenAdd(p1, p2, p3 []uint64) {

	N := len(p1)
	modulus := s.Modulus
	brc := s.BRedConstant

	for j := 0; j < N; j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = CRed(z[0]+BRed(x[0], y[0], modulus, brc), modulus)
		z[1] = CRed(z[1]+BRed(x[1], y[1], modulus, brc), modulus)
		z[2] = CRed(z[2]+BRed(x[2], y[2], modulus, brc), modulus)
		z[3] = CRed(z[3]+BRed(x[3], y[3], modulus, brc), modulus)
		z[4] = CRed(z[4]+BRed(x[4], y[4], modulus, brc), modulus)
		z[5] = CRed(z[5]+BRed(x[5], y[5], modulus, brc), modulus)
		z[6] = CRed(z[6]+BRed(x[6], y[6], modulus, brc), modulus)
		z[7] = CRed(z[7]+BRed(x[7], y[7], modulus, brc), modulus)
	}
}

// MulCoeffsThenAddLazy evaluates p3 = p3 + p1*p2 (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsThenAddLazy(p1, p2, p3 []uint64) {

	N := len(p1)
	modulus := s.Modulus
	brc := s.BRedConstant

	for j := 0; j < N; j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] += BRed(x[0], y[0], modulus, brc)
		z[1] += BRed(x[1], y[1], modulus, brc)
		z[2] += BRed(x[2], y[2], modulus, brc)
		z[3] += BRed(x[3], y[3], modulus, brc)
		z[4] += BRed(x[4], y[4], modulus, brc)
		z[5] += BRed(x[5], y[5], modulus, brc)
		z[6] += BRed(x[6], y[6], modulus, brc)
		z[7] += BRed(x[7], y[7], modulus, brc)
	}
}

// MulCoeffsMontgomery evaluates p3 = p1*p2 (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsMontgomery(p1, p2, p3 []uint64) {

	N := len(p1)
	modulus := s.Modulus
	mrc := s.MRedConstant

	for j := 0; j < N; j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = MRed(x[0], y[0], modulus, mrc)
		z[1] = MRed(x[1], y[1], modulus, mrc)
		z[2] = MRed(x[2], y[2], modulus, mrc)
		z[3] = MRed(x[3], y[3], modulus, mrc)
		z[4] = MRed(x[4], y[4], modulus, mrc)
		z[5] = MRed(x[5], y[5], modulus, mrc)
		z[6] = MRed(x[6], y[6], modulus, mrc)
		z[7] = MRed(x[7], y[7], modulus, mrc)
	}
}

// MulCoeffsMontgomeryLazy evaluates p3 = p1*p2 (mod modulus) with p3 in range [0, 2*modulus-1].
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsMontgomeryLazy(p1, p2, p3 []uint64) {

	N := len(p1)
	modulus := s.Modulus
	mrc := s.MRedConstant

	for j := 0; j < N; j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = MRedLazy(x[0], y[0], modulus, mrc)
		z[1] = MRedLazy(x[1], y[1], modulus, mrc)
		z[2] = MRedLazy(x[2], y[2], modulus, mrc)
		z[3] = MRedLazy(x[3], y[3], modulus, mrc)
		z[4] = MRedLazy(x[4], y[4], modulus, mrc)
		z[5] = MRedLazy(x[5], y[5], modulus, mrc)
		z[6] = MRedLazy(x[6], y[6], modulus, mrc)
		z[7] = MRedLazy(x[7], y[7], modulus, mrc)
	}
}

// MulCoeffsMontgomeryThenAdd evaluates p3 = p3 + (p1*p2) (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsMontgomeryThenAdd(p1, p2, p3 []uint64) {

	N := len(p1)
	modulus := s.Modulus
	mrc := s.MRedConstant

	for j := 0; j < N; j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = CRed(z[0]+MRed(x[0], y[0], modulus, mrc), modulus)
		z[1] = CRed(z[1]+MRed(x[1], y[1], modulus, mrc), modulus)
		z[2] = CRed(z[2]+MRed(x[2], y[2], modulus, mrc), modulus)
		z[3] = CRed(z[3]+MRed(x[3], y[3], modulus, mrc), modulus)
		z[4] = CRed(z[4]+MRed(x[4], y[4], modulus, mrc), modulus)
		z[5] = CRed(z[5]+MRed(x[5], y[5], modulus, mrc), modulus)
		z[6] = CRed(z[6]+MRed(x[6], y[6], modulus, mrc), modulus)
		z[7] = CRed(z[7]+MRed(x[7], y[7], modulus, mrc), modulus)
	}
}

// MulCoeffsMontgomeryThenAddLazy evaluates p3 = p3 + (p1*p2 (mod modulus)).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsMontgomeryThenAddLazy(p1, p2, p3 []uint64) {

	N := len(p1)
	modulus := s.Modulus
	mrc := s.MRedConstant

	for j := 0; j < N; j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] += MRed(x[0], y[0], modulus, mrc)
		z[1] += MRed(x[1], y[1], modulus, mrc)
		z[2] += MRed(x[2], y[2], modulus, mrc)
		z[3] += MRed(x[3], y[3], modulus, mrc)
		z[4] += MRed(x[4], y[4], modulus, mrc)
		z[5] += MRed(x[5], y[5], modulus, mrc)
		z[6] += MRed(x[6], y[6], modulus, mrc)
		z[7] += MRed(x[7], y[7], modulus, mrc)
	}
}

// MulCoeffsMontgomeryLazyThenAddLazy evaluates p3 = p3 + p1*p2 (mod modulus) with p3 in range [0, 3modulus-2].
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsMontgomeryLazyThenAddLazy(p1, p2, p3 []uint64) {

	N := len(p1)
	modulus := s.Modulus
	mrc := s.MRedConstant

	for j := 0; j < N; j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] += MRedLazy(x[0], y[0], modulus, mrc)
		z[1] += MRedLazy(x[1], y[1], modulus, mrc)
		z[2] += MRedLazy(x[2], y[2], modulus, mrc)
		z[3] += MRedLazy(x[3], y[3], modulus, mrc)
		z[4] += MRedLazy(x[4], y[4], modulus, mrc)
		z[5] += MRedLazy(x[5], y[5], modulus, mrc)
		z[6] += MRedLazy(x[6], y[6], modulus, mrc)
		z[7] += MRedLazy(x[7], y[7], modulus, mrc)
	}
}

// MulCoeffsMontgomeryThenSub evaluates p3 = p3 - p1*p2 (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsMontgomeryThenSub(p1, p2, p3 []uint64) {

	N := len(p1)
	modulus := s.Modulus
	mrc := s.MRedConstant

	for j := 0; j < N; j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = CRed(z[0]+(modulus-MRed(x[0], y[0], modulus, mrc)), modulus)
		z[1] = CRed(z[1]+(modulus-MRed(x[1], y[1], modulus, mrc)), modulus)
		z[2] = CRed(z[2]+(modulus-MRed(x[2], y[2], modulus, mrc)), modulus)
		z[3] = CRed(z[3]+(modulus-MRed(x[3], y[3], modulus, mrc)), modulus)
		z[4] = CRed(z[4]+(modulus-MRed(x[4], y[4], modulus, mrc)), modulus)
		z[5] = CRed(z[5]+(modulus-MRed(x[5], y[5], modulus, mrc)), modulus)
		z[6] = CRed(z[6]+(modulus-MRed(x[6], y[6], modulus, mrc)), modulus)
		z[7] = CRed(z[7]+(modulus-MRed(x[7], y[7], modulus, mrc)), modulus)
	}
}

// MulCoeffsMontgomeryThenSubLazy evaluates p3 = p3 - p1*p2 (mod modulus) with p3 in range [0, 2*modulus-2].
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsMontgomeryThenSubLazy(p1, p2, p3 []uint64) {

	N := len(p1)
	modulus := s.Modulus
	mrc := s.MRedConstant

	for j := 0; j < N; j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] += (modulus - MRed(x[0], y[0], modulus, mrc))
		z[1] += (modulus - MRed(x[1], y[1], modulus, mrc))
		z[2] += (modulus - MRed(x[2], y[2], modulus, mrc))
		z[3] += (modulus - MRed(x[3], y[3], modulus, mrc))
		z[4] += (modulus - MRed(x[4], y[4], modulus, mrc))
		z[5] += (modulus - MRed(x[5], y[5], modulus, mrc))
		z[6] += (modulus - MRed(x[6], y[6], modulus, mrc))
		z[7] += (modulus - MRed(x[7], y[7], modulus, mrc))
	}
}

// MulCoeffsMontgomeryLazyThenSubLazy evaluates p3 = p3 - p1*p2 (mod modulus) with p3 in range [0, 3*modulus-2].
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsMontgomeryLazyThenSubLazy(p1, p2, p3 []uint64) {

	N := len(p1)
	modulus := s.Modulus
	mrc := s.MRedConstant
	twomodulus := modulus << 1

	for j := 0; j < N; j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] += twomodulus - MRedLazy(x[0], y[0], modulus, mrc)
		z[1] += twomodulus - MRedLazy(x[1], y[1], modulus, mrc)
		z[2] += twomodulus - MRedLazy(x[2], y[2], modulus, mrc)
		z[3] += twomodulus - MRedLazy(x[3], y[3], modulus, mrc)
		z[4] += twomodulus - MRedLazy(x[4], y[4], modulus, mrc)
		z[5] += twomodulus - MRedLazy(x[5], y[5], modulus, mrc)
		z[6] += twomodulus - MRedLazy(x[6], y[6], modulus, mrc)
		z[7] += twomodulus - MRedLazy(x[7], y[7], modulus, mrc)
	}
}

// MulCoeffsMontgomeryLazyThenNeg evaluates p3 = - p1*p2 (mod modulus) with p3 in range [0, 2*modulus-2].
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsMontgomeryLazyThenNeg(p1, p2, p3 []uint64) {

	N := len(p1)
	modulus := s.Modulus
	mrc := s.MRedConstant
	twomodulus := modulus << 1

	for j := 0; j < N; j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = twomodulus - MRedLazy(x[0], y[0], modulus, mrc)
		z[1] = twomodulus - MRedLazy(x[1], y[1], modulus, mrc)
		z[2] = twomodulus - MRedLazy(x[2], y[2], modulus, mrc)
		z[3] = twomodulus - MRedLazy(x[3], y[3], modulus, mrc)
		z[4] = twomodulus - MRedLazy(x[4], y[4], modulus, mrc)
		z[5] = twomodulus - MRedLazy(x[5], y[5], modulus, mrc)
		z[6] = twomodulus - MRedLazy(x[6], y[6], modulus, mrc)
		z[7] = twomodulus - MRedLazy(x[7], y[7], modulus, mrc)
	}
}

// MulCoeffsBarrettLazy evaluates p3 = p1*p2 (mod modulus) with p3 in range [0, 2*modulus-1].
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulCoeffsBarrettLazy(p1, p2, p3 []uint64) {

	N := len(p1)
	modulus := s.Modulus
	brc := s.BRedConstant

	for j := 0; j < N; j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = BRedLazy(x[0], y[0], modulus, brc)
		z[1] = BRedLazy(x[1], y[1], modulus, brc)
		z[2] = BRedLazy(x[2], y[2], modulus, brc)
		z[3] = BRedLazy(x[3], y[3], modulus, brc)
		z[4] = BRedLazy(x[4], y[4], modulus, brc)
		z[5] = BRedLazy(x[5], y[5], modulus, brc)
		z[6] = BRedLazy(x[6], y[6], modulus, brc)
		z[7] = BRedLazy(x[7], y[7], modulus, brc)
	}
}

// AddLazyThenMulScalarMontgomery evaluates p3 = (p1+p2)*scalarMont (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) AddLazyThenMulScalarMontgomery(p1, p2 []uint64, scalarMont uint64, p3 []uint64) {

	N := len(p1)
	modulus := s.Modulus
	mrc := s.MRedConstant

	for j := 0; j < N; j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = MRed(x[0]+y[0], scalarMont, modulus, mrc)
		z[1] = MRed(x[1]+y[1], scalarMont, modulus, mrc)
		z[2] = MRed(x[2]+y[2], scalarMont, modulus, mrc)
		z[3] = MRed(x[3]+y[3], scalarMont, modulus, mrc)
		z[4] = MRed(x[4]+y[4], scalarMont, modulus, mrc)
		z[5] = MRed(x[5]+y[5], scalarMont, modulus, mrc)
		z[6] = MRed(x[6]+y[6], scalarMont, modulus, mrc)
		z[7] = MRed(x[7]+y[7], scalarMont, modulus, mrc)
	}
}

// AddScalarLazyThenMulScalarMontgomery evaluates p3 = (scalarMont0+p2)*scalarMont1 (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) AddScalarLazyThenMulScalarMontgomery(p1 []uint64, scalar0, scalarMont1 uint64, p2 []uint64) {

	N := len(p1)
	modulus := s.Modulus
	mrc := s.MRedConstant

	for j := 0; j < N; j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = MRed(x[0]+scalar0, scalarMont1, modulus, mrc)
		z[1] = MRed(x[1]+scalar0, scalarMont1, modulus, mrc)
		z[2] = MRed(x[2]+scalar0, scalarMont1, modulus, mrc)
		z[3] = MRed(x[3]+scalar0, scalarMont1, modulus, mrc)
		z[4] = MRed(x[4]+scalar0, scalarMont1, modulus, mrc)
		z[5] = MRed(x[5]+scalar0, scalarMont1, modulus, mrc)
		z[6] = MRed(x[6]+scalar0, scalarMont1, modulus, mrc)
		z[7] = MRed(x[7]+scalar0, scalarMont1, modulus, mrc)
	}
}

// AddScalar evaluates p2 = p1 + scalar (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) AddScalar(p1 []uint64, scalar uint64, p2 []uint64) {

	N := len(p1)
	modulus := s.Modulus

	for j := 0; j < N; j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = CRed(x[0]+scalar, modulus)
		z[1] = CRed(x[1]+scalar, modulus)
		z[2] = CRed(x[2]+scalar, modulus)
		z[3] = CRed(x[3]+scalar, modulus)
		z[4] = CRed(x[4]+scalar, modulus)
		z[5] = CRed(x[5]+scalar, modulus)
		z[6] = CRed(x[6]+scalar, modulus)
		z[7] = CRed(x[7]+scalar, modulus)
	}
}

// AddScalarLazy evaluates p2 = p1 + scalar.
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) AddScalarLazy(p1 []uint64, scalar uint64, p2 []uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = x[0] + scalar
		z[1] = x[1] + scalar
		z[2] = x[2] + scalar
		z[3] = x[3] + scalar
		z[4] = x[4] + scalar
		z[5] = x[5] + scalar
		z[6] = x[6] + scalar
		z[7] = x[7] + scalar
	}
}

// AddScalarLazyThenNegTwoModulusLazy evaluates p2 = 2*modulus - p1 + scalar.
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) AddScalarLazyThenNegTwoModulusLazy(p1 []uint64, scalar uint64, p2 []uint64) {

	N := len(p1)
	modulus := s.Modulus
	twomodulus := modulus << 1

	for j := 0; j < N; j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = scalar + twomodulus - x[0]
		z[1] = scalar + twomodulus - x[1]
		z[2] = scalar + twomodulus - x[2]
		z[3] = scalar + twomodulus - x[3]
		z[4] = scalar + twomodulus - x[4]
		z[5] = scalar + twomodulus - x[5]
		z[6] = scalar + twomodulus - x[6]
		z[7] = scalar + twomodulus - x[7]
	}
}

// SubScalar evaluates p2 = p1 - scalar (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) SubScalar(p1 []uint64, scalar uint64, p2 []uint64) {

	N := len(p1)
	modulus := s.Modulus

	for j := 0; j < N; j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = CRed(x[0]+modulus-scalar, modulus)
		z[1] = CRed(x[1]+modulus-scalar, modulus)
		z[2] = CRed(x[2]+modulus-scalar, modulus)
		z[3] = CRed(x[3]+modulus-scalar, modulus)
		z[4] = CRed(x[4]+modulus-scalar, modulus)
		z[5] = CRed(x[5]+modulus-scalar, modulus)
		z[6] = CRed(x[6]+modulus-scalar, modulus)
		z[7] = CRed(x[7]+modulus-scalar, modulus)
	}
}

// MulScalarMontgomery evaluates p2 = p1*scalarMont (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulScalarMontgomery(p1 []uint64, scalarMont uint64, p2 []uint64) {

	N := len(p1)
	modulus := s.Modulus
	mrc := s.MRedConstant

	for j := 0; j < N; j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = MRed(x[0], scalarMont, modulus, mrc)
		z[1] = MRed(x[1], scalarMont, modulus, mrc)
		z[2] = MRed(x[2], scalarMont, modulus, mrc)
		z[3] = MRed(x[3], scalarMont, modulus, mrc)
		z[4] = MRed(x[4], scalarMont, modulus, mrc)
		z[5] = MRed(x[5], scalarMont, modulus, mrc)
		z[6] = MRed(x[6], scalarMont, modulus, mrc)
		z[7] = MRed(x[7], scalarMont, modulus, mrc)
	}
}

// MulScalarMontgomeryLazy evaluates p2 = p1*scalarMont (mod modulus) with p2 in range [0, 2*modulus-1].
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulScalarMontgomeryLazy(p1 []uint64, scalarMont uint64, p2 []uint64) {

	N := len(p1)
	modulus := s.Modulus
	mrc := s.MRedConstant

	for j := 0; j < N; j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = MRedLazy(x[0], scalarMont, modulus, mrc)
		z[1] = MRedLazy(x[1], scalarMont, modulus, mrc)
		z[2] = MRedLazy(x[2], scalarMont, modulus, mrc)
		z[3] = MRedLazy(x[3], scalarMont, modulus, mrc)
		z[4] = MRedLazy(x[4], scalarMont, modulus, mrc)
		z[5] = MRedLazy(x[5], scalarMont, modulus, mrc)
		z[6] = MRedLazy(x[6], scalarMont, modulus, mrc)
		z[7] = MRedLazy(x[7], scalarMont, modulus, mrc)
	}
}

// MulScalarMontgomeryThenAdd evaluates p2 = p2 + p1*scalarMont (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulScalarMontgomeryThenAdd(p1 []uint64, scalarMont uint64, p2 []uint64) {

	N := len(p1)
	modulus := s.Modulus
	mrc := s.MRedConstant

	for j := 0; j < N; j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = CRed(z[0]+MRed(x[0], scalarMont, modulus, mrc), modulus)
		z[1] = CRed(z[1]+MRed(x[1], scalarMont, modulus, mrc), modulus)
		z[2] = CRed(z[2]+MRed(x[2], scalarMont, modulus, mrc), modulus)
		z[3] = CRed(z[3]+MRed(x[3], scalarMont, modulus, mrc), modulus)
		z[4] = CRed(z[4]+MRed(x[4], scalarMont, modulus, mrc), modulus)
		z[5] = CRed(z[5]+MRed(x[5], scalarMont, modulus, mrc), modulus)
		z[6] = CRed(z[6]+MRed(x[6], scalarMont, modulus, mrc), modulus)
		z[7] = CRed(z[7]+MRed(x[7], scalarMont, modulus, mrc), modulus)
	}
}

// MulScalarMontgomeryThenAddScalar evaluates p2 = scalar + p1*scalarMont (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MulScalarMontgomeryThenAddScalar(p1 []uint64, scalar0, scalarMont1 uint64, p2 []uint64) {

	N := len(p1)
	modulus := s.Modulus
	mrc := s.MRedConstant

	for j := 0; j < N; j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = CRed(MRed(x[0], scalarMont1, modulus, mrc)+scalar0, modulus)
		z[1] = CRed(MRed(x[1], scalarMont1, modulus, mrc)+scalar0, modulus)
		z[2] = CRed(MRed(x[2], scalarMont1, modulus, mrc)+scalar0, modulus)
		z[3] = CRed(MRed(x[3], scalarMont1, modulus, mrc)+scalar0, modulus)
		z[4] = CRed(MRed(x[4], scalarMont1, modulus, mrc)+scalar0, modulus)
		z[5] = CRed(MRed(x[5], scalarMont1, modulus, mrc)+scalar0, modulus)
		z[6] = CRed(MRed(x[6], scalarMont1, modulus, mrc)+scalar0, modulus)
		z[7] = CRed(MRed(x[7], scalarMont1, modulus, mrc)+scalar0, modulus)
	}
}

// SubThenMulScalarMontgomeryTwoModulus evaluates p3 = (p1 + twomodulus - p2) * scalarMont (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) SubThenMulScalarMontgomeryTwoModulus(p1, p2 []uint64, scalarMont uint64, p3 []uint64) {

	N := len(p1)
	modulus := s.Modulus
	mrc := s.MRedConstant
	twomodulus := modulus << 1

	for j := 0; j < N; j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = MRed(twomodulus-y[0]+x[0], scalarMont, modulus, mrc)
		z[1] = MRed(twomodulus-y[1]+x[1], scalarMont, modulus, mrc)
		z[2] = MRed(twomodulus-y[2]+x[2], scalarMont, modulus, mrc)
		z[3] = MRed(twomodulus-y[3]+x[3], scalarMont, modulus, mrc)
		z[4] = MRed(twomodulus-y[4]+x[4], scalarMont, modulus, mrc)
		z[5] = MRed(twomodulus-y[5]+x[5], scalarMont, modulus, mrc)
		z[6] = MRed(twomodulus-y[6]+x[6], scalarMont, modulus, mrc)
		z[7] = MRed(twomodulus-y[7]+x[7], scalarMont, modulus, mrc)

	}
}

// MForm evaluates p2 = p1 * 2^64 (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MForm(p1, p2 []uint64) {

	N := len(p1)
	modulus := s.Modulus
	brc := s.BRedConstant

	for j := 0; j < N; j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = MForm(x[0], modulus, brc)
		z[1] = MForm(x[1], modulus, brc)
		z[2] = MForm(x[2], modulus, brc)
		z[3] = MForm(x[3], modulus, brc)
		z[4] = MForm(x[4], modulus, brc)
		z[5] = MForm(x[5], modulus, brc)
		z[6] = MForm(x[6], modulus, brc)
		z[7] = MForm(x[7], modulus, brc)
	}
}

// MFormLazy evaluates p2 = p1 * 2^64 (mod modulus) with p2 in the range [0, 2*modulus-1].
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) MFormLazy(p1, p2 []uint64) {

	N := len(p1)
	modulus := s.Modulus
	brc := s.BRedConstant

	for j := 0; j < N; j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = MFormLazy(x[0], modulus, brc)
		z[1] = MFormLazy(x[1], modulus, brc)
		z[2] = MFormLazy(x[2], modulus, brc)
		z[3] = MFormLazy(x[3], modulus, brc)
		z[4] = MFormLazy(x[4], modulus, brc)
		z[5] = MFormLazy(x[5], modulus, brc)
		z[6] = MFormLazy(x[6], modulus, brc)
		z[7] = MFormLazy(x[7], modulus, brc)
	}
}

// IMForm evaluates p2 = p1 * (2^64)^-1 (mod modulus).
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func (s *SubRing) IMForm(p1, p2 []uint64) {

	N := len(p1)
	modulus := s.Modulus
	mrc := s.MRedConstant

	for j := 0; j < N; j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = IMForm(x[0], modulus, mrc)
		z[1] = IMForm(x[1], modulus, mrc)
		z[2] = IMForm(x[2], modulus, mrc)
		z[3] = IMForm(x[3], modulus, mrc)
		z[4] = IMForm(x[4], modulus, mrc)
		z[5] = IMForm(x[5], modulus, mrc)
		z[6] = IMForm(x[6], modulus, mrc)
		z[7] = IMForm(x[7], modulus, mrc)
	}
}

// Zero sets all values of p1 to zero.
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func Zero(p1 []uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		z := (*[8]uint64)(unsafe.Pointer(&p1[j]))

		z[0] = 0
		z[1] = 0
		z[2] = 0
		z[3] = 0
		z[4] = 0
		z[5] = 0
		z[6] = 0
		z[7] = 0
	}
}

// Mask evaluates p2 = (p1>>w) & mask
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func Mask(p1 []uint64, w int, mask uint64, p2 []uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = (x[0] >> w) & mask
		z[1] = (x[1] >> w) & mask
		z[2] = (x[2] >> w) & mask
		z[3] = (x[3] >> w) & mask
		z[4] = (x[4] >> w) & mask
		z[5] = (x[5] >> w) & mask
		z[6] = (x[6] >> w) & mask
		z[7] = (x[7] >> w) & mask
	}
}
