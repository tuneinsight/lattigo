package ring

import (
	"math/big"
	"math/bits"
	"unsafe"
)

// AddVec returns p3 = p1 + p2 mod qi.
func AddVec(p1, p2, p3 []uint64, qi uint64) {
	for j := 0; j < len(p1); j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = CRed(x[0]+y[0], qi)
		z[1] = CRed(x[1]+y[1], qi)
		z[2] = CRed(x[2]+y[2], qi)
		z[3] = CRed(x[3]+y[3], qi)
		z[4] = CRed(x[4]+y[4], qi)
		z[5] = CRed(x[5]+y[5], qi)
		z[6] = CRed(x[6]+y[6], qi)
		z[7] = CRed(x[7]+y[7], qi)
	}
}

// AddVecNoMod returns p3 = p1 + p2.
func AddVecNoMod(p1, p2, p3 []uint64) {
	for j := 0; j < len(p1); j = j + 8 {
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

// SubVec returns p3 = p1 - p2 mod qi.
func SubVec(p1, p2, p3 []uint64, qi uint64) {
	for j := 0; j < len(p1); j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = CRed((x[0]+qi)-y[0], qi)
		z[1] = CRed((x[1]+qi)-y[1], qi)
		z[2] = CRed((x[2]+qi)-y[2], qi)
		z[3] = CRed((x[3]+qi)-y[3], qi)
		z[4] = CRed((x[4]+qi)-y[4], qi)
		z[5] = CRed((x[5]+qi)-y[5], qi)
		z[6] = CRed((x[6]+qi)-y[6], qi)
		z[7] = CRed((x[7]+qi)-y[7], qi)
	}
}

// SubVecNomod returns p3 = p1 + qi - p2.
func SubVecNomod(p1, p2, p3 []uint64, qi uint64) {
	for j := 0; j < len(p1); j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = x[0] + qi - y[0]
		z[1] = x[1] + qi - y[1]
		z[2] = x[2] + qi - y[2]
		z[3] = x[3] + qi - y[3]
		z[4] = x[4] + qi - y[4]
		z[5] = x[5] + qi - y[5]
		z[6] = x[6] + qi - y[6]
		z[7] = x[7] + qi - y[7]
	}
}

// NegVec returns p2 = -p1 mod qi.
func NegVec(p1, p2 []uint64, qi uint64) {
	for j := 0; j < len(p1); j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = qi - x[0]
		z[1] = qi - x[1]
		z[2] = qi - x[2]
		z[3] = qi - x[3]
		z[4] = qi - x[4]
		z[5] = qi - x[5]
		z[6] = qi - x[6]
		z[7] = qi - x[7]
	}
}

// ReduceVec returns p2 = p1 mod qi.
func ReduceVec(p1, p2 []uint64, qi uint64, bredParams []uint64) {
	for j := 0; j < len(p1); j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = BRedAdd(x[0], qi, bredParams)
		z[1] = BRedAdd(x[1], qi, bredParams)
		z[2] = BRedAdd(x[2], qi, bredParams)
		z[3] = BRedAdd(x[3], qi, bredParams)
		z[4] = BRedAdd(x[4], qi, bredParams)
		z[5] = BRedAdd(x[5], qi, bredParams)
		z[6] = BRedAdd(x[6], qi, bredParams)
		z[7] = BRedAdd(x[7], qi, bredParams)
	}
}

// ReduceConstantVec returns p2 = p1 mod qi with output coefficients range [0, 2qi-1].
func ReduceConstantVec(p1, p2 []uint64, qi uint64, bredParams []uint64) {
	for j := 0; j < len(p1); j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = BRedAddConstant(x[0], qi, bredParams)
		z[1] = BRedAddConstant(x[1], qi, bredParams)
		z[2] = BRedAddConstant(x[2], qi, bredParams)
		z[3] = BRedAddConstant(x[3], qi, bredParams)
		z[4] = BRedAddConstant(x[4], qi, bredParams)
		z[5] = BRedAddConstant(x[5], qi, bredParams)
		z[6] = BRedAddConstant(x[6], qi, bredParams)
		z[7] = BRedAddConstant(x[7], qi, bredParams)
	}
}

// ModVec returns p2 = p1 mod m.
func ModVec(p1, p2 []uint64, m uint64, bredParams []uint64) {
	for j := 0; j < len(p1); j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = BRedAdd(x[0], m, bredParams)
		z[1] = BRedAdd(x[1], m, bredParams)
		z[2] = BRedAdd(x[2], m, bredParams)
		z[3] = BRedAdd(x[3], m, bredParams)
		z[4] = BRedAdd(x[4], m, bredParams)
		z[5] = BRedAdd(x[5], m, bredParams)
		z[6] = BRedAdd(x[6], m, bredParams)
		z[7] = BRedAdd(x[7], m, bredParams)
	}
}

// MulCoeffsVec returns p3 = p1*p2 mod qi.
func MulCoeffsVec(p1, p2, p3 []uint64, qi uint64, bredParams []uint64) {
	for j := 0; j < len(p1); j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = BRed(x[0], y[0], qi, bredParams)
		z[1] = BRed(x[1], y[1], qi, bredParams)
		z[2] = BRed(x[2], y[2], qi, bredParams)
		z[3] = BRed(x[3], y[3], qi, bredParams)
		z[4] = BRed(x[4], y[4], qi, bredParams)
		z[5] = BRed(x[5], y[5], qi, bredParams)
		z[6] = BRed(x[6], y[6], qi, bredParams)
		z[7] = BRed(x[7], y[7], qi, bredParams)
	}
}

// MulCoeffsAndAddVec returns p3 = p3 + (p1*p2) mod qi.
func MulCoeffsAndAddVec(p1, p2, p3 []uint64, qi uint64, bredParams []uint64) {
	for j := 0; j < len(p1); j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = CRed(z[0]+BRed(x[0], y[0], qi, bredParams), qi)
		z[1] = CRed(z[1]+BRed(x[1], y[1], qi, bredParams), qi)
		z[2] = CRed(z[2]+BRed(x[2], y[2], qi, bredParams), qi)
		z[3] = CRed(z[3]+BRed(x[3], y[3], qi, bredParams), qi)
		z[4] = CRed(z[4]+BRed(x[4], y[4], qi, bredParams), qi)
		z[5] = CRed(z[5]+BRed(x[5], y[5], qi, bredParams), qi)
		z[6] = CRed(z[6]+BRed(x[6], y[6], qi, bredParams), qi)
		z[7] = CRed(z[7]+BRed(x[7], y[7], qi, bredParams), qi)
	}
}

// MulCoeffsAndAddNoModVec returns p3 = p3 + (p1*p2 mod qi).
func MulCoeffsAndAddNoModVec(p1, p2, p3 []uint64, qi uint64, bredParams []uint64) {
	for j := 0; j < len(p1); j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] += BRed(x[0], y[0], qi, bredParams)
		z[1] += BRed(x[1], y[1], qi, bredParams)
		z[2] += BRed(x[2], y[2], qi, bredParams)
		z[3] += BRed(x[3], y[3], qi, bredParams)
		z[4] += BRed(x[4], y[4], qi, bredParams)
		z[5] += BRed(x[5], y[5], qi, bredParams)
		z[6] += BRed(x[6], y[6], qi, bredParams)
		z[7] += BRed(x[7], y[7], qi, bredParams)
	}
}

// MulCoeffsMontgomeryVec returns p3 = p1*p2 mod qi.
func MulCoeffsMontgomeryVec(p1, p2, p3 []uint64, qi, mredParams uint64) {
	for j := 0; j < len(p1); j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = MRed(x[0], y[0], qi, mredParams)
		z[1] = MRed(x[1], y[1], qi, mredParams)
		z[2] = MRed(x[2], y[2], qi, mredParams)
		z[3] = MRed(x[3], y[3], qi, mredParams)
		z[4] = MRed(x[4], y[4], qi, mredParams)
		z[5] = MRed(x[5], y[5], qi, mredParams)
		z[6] = MRed(x[6], y[6], qi, mredParams)
		z[7] = MRed(x[7], y[7], qi, mredParams)
	}
}

// MulCoeffsMontgomeryConstantVec returns p3 = p1*p2 mod qi with output coefficients in range [0, 2qi-1].
func MulCoeffsMontgomeryConstantVec(p1, p2, p3 []uint64, qi, mredParams uint64) {
	for j := 0; j < len(p1); j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = MRedConstant(x[0], y[0], qi, mredParams)
		z[1] = MRedConstant(x[1], y[1], qi, mredParams)
		z[2] = MRedConstant(x[2], y[2], qi, mredParams)
		z[3] = MRedConstant(x[3], y[3], qi, mredParams)
		z[4] = MRedConstant(x[4], y[4], qi, mredParams)
		z[5] = MRedConstant(x[5], y[5], qi, mredParams)
		z[6] = MRedConstant(x[6], y[6], qi, mredParams)
		z[7] = MRedConstant(x[7], y[7], qi, mredParams)
	}
}

// MulCoeffsMontgomeryAndAddVec returns p3 = p3 + (p1*p2) mod qi.
func MulCoeffsMontgomeryAndAddVec(p1, p2, p3 []uint64, qi, mredParams uint64) {
	for j := 0; j < len(p1); j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = CRed(z[0]+MRed(x[0], y[0], qi, mredParams), qi)
		z[1] = CRed(z[1]+MRed(x[1], y[1], qi, mredParams), qi)
		z[2] = CRed(z[2]+MRed(x[2], y[2], qi, mredParams), qi)
		z[3] = CRed(z[3]+MRed(x[3], y[3], qi, mredParams), qi)
		z[4] = CRed(z[4]+MRed(x[4], y[4], qi, mredParams), qi)
		z[5] = CRed(z[5]+MRed(x[5], y[5], qi, mredParams), qi)
		z[6] = CRed(z[6]+MRed(x[6], y[6], qi, mredParams), qi)
		z[7] = CRed(z[7]+MRed(x[7], y[7], qi, mredParams), qi)
	}
}

// MulCoeffsMontgomeryAndAddNoModVec returns p3 = p3 + (p1*p2 mod qi).
func MulCoeffsMontgomeryAndAddNoModVec(p1, p2, p3 []uint64, qi, mredParams uint64) {
	for j := 0; j < len(p1); j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] += MRed(x[0], y[0], qi, mredParams)
		z[1] += MRed(x[1], y[1], qi, mredParams)
		z[2] += MRed(x[2], y[2], qi, mredParams)
		z[3] += MRed(x[3], y[3], qi, mredParams)
		z[4] += MRed(x[4], y[4], qi, mredParams)
		z[5] += MRed(x[5], y[5], qi, mredParams)
		z[6] += MRed(x[6], y[6], qi, mredParams)
		z[7] += MRed(x[7], y[7], qi, mredParams)
	}
}

// MulCoeffsMontgomeryConstantAndAddNoModVec returns p3 = p3 + p1*p2 mod qi with output coefficients in range [0, 3qi-2].
func MulCoeffsMontgomeryConstantAndAddNoModVec(p1, p2, p3 []uint64, qi, mredParams uint64) {
	for j := 0; j < len(p1); j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] += MRedConstant(x[0], y[0], qi, mredParams)
		z[1] += MRedConstant(x[1], y[1], qi, mredParams)
		z[2] += MRedConstant(x[2], y[2], qi, mredParams)
		z[3] += MRedConstant(x[3], y[3], qi, mredParams)
		z[4] += MRedConstant(x[4], y[4], qi, mredParams)
		z[5] += MRedConstant(x[5], y[5], qi, mredParams)
		z[6] += MRedConstant(x[6], y[6], qi, mredParams)
		z[7] += MRedConstant(x[7], y[7], qi, mredParams)
	}
}

// MulCoeffsMontgomeryAndSubVec returns p3 = p3 - p1*p2 mod qi.
func MulCoeffsMontgomeryAndSubVec(p1, p2, p3 []uint64, qi, mredParams uint64) {
	for j := 0; j < len(p1); j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = CRed(z[0]+(qi-MRed(x[0], y[0], qi, mredParams)), qi)
		z[1] = CRed(z[1]+(qi-MRed(x[1], y[1], qi, mredParams)), qi)
		z[2] = CRed(z[2]+(qi-MRed(x[2], y[2], qi, mredParams)), qi)
		z[3] = CRed(z[3]+(qi-MRed(x[3], y[3], qi, mredParams)), qi)
		z[4] = CRed(z[4]+(qi-MRed(x[4], y[4], qi, mredParams)), qi)
		z[5] = CRed(z[5]+(qi-MRed(x[5], y[5], qi, mredParams)), qi)
		z[6] = CRed(z[6]+(qi-MRed(x[6], y[6], qi, mredParams)), qi)
		z[7] = CRed(z[7]+(qi-MRed(x[7], y[7], qi, mredParams)), qi)
	}
}

// MulCoeffsMontgomeryAndSubNoMod returns p3 = p3 - p1*p2 mod qi with output coefficients in range [0, 2qi-2].
func MulCoeffsMontgomeryAndSubNoMod(p1, p2, p3 []uint64, qi, mredParams uint64) {
	for j := 0; j < len(p1); j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] += (qi - MRed(x[0], y[0], qi, mredParams))
		z[1] += (qi - MRed(x[1], y[1], qi, mredParams))
		z[2] += (qi - MRed(x[2], y[2], qi, mredParams))
		z[3] += (qi - MRed(x[3], y[3], qi, mredParams))
		z[4] += (qi - MRed(x[4], y[4], qi, mredParams))
		z[5] += (qi - MRed(x[5], y[5], qi, mredParams))
		z[6] += (qi - MRed(x[6], y[6], qi, mredParams))
		z[7] += (qi - MRed(x[7], y[7], qi, mredParams))
	}
}

// MulCoeffsMontgomeryConstantAndSubNoMod returns p3 = p3 - p1*p2 mod qi with output coefficients in range [0, 3qi-2].
func MulCoeffsMontgomeryConstantAndSubNoMod(p1, p2, p3 []uint64, qi, mredParams uint64) {
	twoQi := qi << 1
	for j := 0; j < len(p1); j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] += twoQi - MRedConstant(x[0], y[0], qi, mredParams)
		z[1] += twoQi - MRedConstant(x[1], y[1], qi, mredParams)
		z[2] += twoQi - MRedConstant(x[2], y[2], qi, mredParams)
		z[3] += twoQi - MRedConstant(x[3], y[3], qi, mredParams)
		z[4] += twoQi - MRedConstant(x[4], y[4], qi, mredParams)
		z[5] += twoQi - MRedConstant(x[5], y[5], qi, mredParams)
		z[6] += twoQi - MRedConstant(x[6], y[6], qi, mredParams)
		z[7] += twoQi - MRedConstant(x[7], y[7], qi, mredParams)
	}
}

// MulCoeffsMontgomeryConstantAndNeg returns p3 = - p1*p2 mod qi with output coefficients in range [0, 2qi-2].
func MulCoeffsMontgomeryConstantAndNeg(p1, p2, p3 []uint64, qi, mredParams uint64) {
	twoqi := qi << 1
	for j := 0; j < len(p1); j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = twoqi - MRedConstant(x[0], y[0], qi, mredParams)
		z[1] = twoqi - MRedConstant(x[1], y[1], qi, mredParams)
		z[2] = twoqi - MRedConstant(x[2], y[2], qi, mredParams)
		z[3] = twoqi - MRedConstant(x[3], y[3], qi, mredParams)
		z[4] = twoqi - MRedConstant(x[4], y[4], qi, mredParams)
		z[5] = twoqi - MRedConstant(x[5], y[5], qi, mredParams)
		z[6] = twoqi - MRedConstant(x[6], y[6], qi, mredParams)
		z[7] = twoqi - MRedConstant(x[7], y[7], qi, mredParams)
	}
}

// MulCoeffsConstantVec returns p3 = p1*p2 mod qi with output coefficients in range [0, 2qi-1].
func MulCoeffsConstantVec(p1, p2, p3 []uint64, qi uint64, bredParams []uint64) {
	for j := 0; j < len(p1); j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = BRedConstant(x[0], y[0], qi, bredParams)
		z[1] = BRedConstant(x[1], y[1], qi, bredParams)
		z[2] = BRedConstant(x[2], y[2], qi, bredParams)
		z[3] = BRedConstant(x[3], y[3], qi, bredParams)
		z[4] = BRedConstant(x[4], y[4], qi, bredParams)
		z[5] = BRedConstant(x[5], y[5], qi, bredParams)
		z[6] = BRedConstant(x[6], y[6], qi, bredParams)
		z[7] = BRedConstant(x[7], y[7], qi, bredParams)
	}
}

// AddVecNoModAndMulScalarMontgomeryVec returns p3 = (p1+p2)*scalarMont mod qi.
func AddVecNoModAndMulScalarMontgomeryVec(p1, p2, p3 []uint64, scalarMont, qi, mredParams uint64) {
	for j := 0; j < len(p1); j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = MRed(x[0]+y[0], scalarMont, qi, mredParams)
		z[1] = MRed(x[1]+y[1], scalarMont, qi, mredParams)
		z[2] = MRed(x[2]+y[2], scalarMont, qi, mredParams)
		z[3] = MRed(x[3]+y[3], scalarMont, qi, mredParams)
		z[4] = MRed(x[4]+y[4], scalarMont, qi, mredParams)
		z[5] = MRed(x[5]+y[5], scalarMont, qi, mredParams)
		z[6] = MRed(x[6]+y[6], scalarMont, qi, mredParams)
		z[7] = MRed(x[7]+y[7], scalarMont, qi, mredParams)
	}
}

// AddScalarVec returns p2 = p1 + scalar mod qi.
func AddScalarVec(p1, p2 []uint64, scalar, qi uint64) {
	for j := 0; j < len(p1); j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = CRed(x[0]+scalar, qi)
		z[1] = CRed(x[1]+scalar, qi)
		z[2] = CRed(x[2]+scalar, qi)
		z[3] = CRed(x[3]+scalar, qi)
		z[4] = CRed(x[4]+scalar, qi)
		z[5] = CRed(x[5]+scalar, qi)
		z[6] = CRed(x[6]+scalar, qi)
		z[7] = CRed(x[7]+scalar, qi)
	}
}

// AddScalarNoModVec returns p2 = p1 + scalar.
func AddScalarNoModVec(p1, p2 []uint64, scalar uint64) {
	for j := 0; j < len(p1); j = j + 8 {
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

// AddScalarNoModAndNegTwoQiNoModVec returns p2 = 2*qi - p1 + scalar.
func AddScalarNoModAndNegTwoQiNoModVec(p1, p2 []uint64, scalar, qi uint64) {
	twoqi := qi << 1
	for j := 0; j < len(p1); j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = scalar + twoqi - x[0]
		z[1] = scalar + twoqi - x[1]
		z[2] = scalar + twoqi - x[2]
		z[3] = scalar + twoqi - x[3]
		z[4] = scalar + twoqi - x[4]
		z[5] = scalar + twoqi - x[5]
		z[6] = scalar + twoqi - x[6]
		z[7] = scalar + twoqi - x[7]
	}
}

// SubScalarVec returns p2 = p1 - scalar mod qi.
func SubScalarVec(p1, p2 []uint64, scalar, qi uint64) {
	for j := 0; j < len(p1); j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = CRed(x[0]+qi-scalar, qi)
		z[1] = CRed(x[1]+qi-scalar, qi)
		z[2] = CRed(x[2]+qi-scalar, qi)
		z[3] = CRed(x[3]+qi-scalar, qi)
		z[4] = CRed(x[4]+qi-scalar, qi)
		z[5] = CRed(x[5]+qi-scalar, qi)
		z[6] = CRed(x[6]+qi-scalar, qi)
		z[7] = CRed(x[7]+qi-scalar, qi)
	}
}

// MulScalarMontgomeryVec returns p2 = p1*scalarMont mod qi.
func MulScalarMontgomeryVec(p1, p2 []uint64, scalarMont, qi, mredParams uint64) {
	for j := 0; j < len(p1); j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = MRed(x[0], scalarMont, qi, mredParams)
		z[1] = MRed(x[1], scalarMont, qi, mredParams)
		z[2] = MRed(x[2], scalarMont, qi, mredParams)
		z[3] = MRed(x[3], scalarMont, qi, mredParams)
		z[4] = MRed(x[4], scalarMont, qi, mredParams)
		z[5] = MRed(x[5], scalarMont, qi, mredParams)
		z[6] = MRed(x[6], scalarMont, qi, mredParams)
		z[7] = MRed(x[7], scalarMont, qi, mredParams)
	}
}

// MulScalarMontgomeryConstantVec returns p2 = p1*scalarMont mod qi with output coefficients in range [0, 2qi-1].
func MulScalarMontgomeryConstantVec(p1, p2 []uint64, scalarMont, qi, mredParams uint64) {
	for j := 0; j < len(p1); j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = MRedConstant(x[0], scalarMont, qi, mredParams)
		z[1] = MRedConstant(x[1], scalarMont, qi, mredParams)
		z[2] = MRedConstant(x[2], scalarMont, qi, mredParams)
		z[3] = MRedConstant(x[3], scalarMont, qi, mredParams)
		z[4] = MRedConstant(x[4], scalarMont, qi, mredParams)
		z[5] = MRedConstant(x[5], scalarMont, qi, mredParams)
		z[6] = MRedConstant(x[6], scalarMont, qi, mredParams)
		z[7] = MRedConstant(x[7], scalarMont, qi, mredParams)
	}
}

// MulScalarMontgomeryAndAddVec returns p2 = p2 + p1*scalarMont mod qi.
func MulScalarMontgomeryAndAddVec(p1, p2 []uint64, scalarMont, qi, mredParams uint64) {
	for j := 0; j < len(p1); j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = CRed(z[0]+MRed(x[0], scalarMont, qi, mredParams), qi)
		z[1] = CRed(z[1]+MRed(x[1], scalarMont, qi, mredParams), qi)
		z[2] = CRed(z[2]+MRed(x[2], scalarMont, qi, mredParams), qi)
		z[3] = CRed(z[3]+MRed(x[3], scalarMont, qi, mredParams), qi)
		z[4] = CRed(z[4]+MRed(x[4], scalarMont, qi, mredParams), qi)
		z[5] = CRed(z[5]+MRed(x[5], scalarMont, qi, mredParams), qi)
		z[6] = CRed(z[6]+MRed(x[6], scalarMont, qi, mredParams), qi)
		z[7] = CRed(z[7]+MRed(x[7], scalarMont, qi, mredParams), qi)
	}
}

// SubVecAndMulScalarMontgomeryTwoQiVec returns p3 = (p1 + twoqi - p2) * scalarMont mod qi.
func SubVecAndMulScalarMontgomeryTwoQiVec(p1, p2, p3 []uint64, scalarMont, qi, mredParams uint64) {
	twoqi := qi << 1
	for j := 0; j < len(p1); j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = MRed(twoqi-y[0]+x[0], scalarMont, qi, mredParams)
		z[1] = MRed(twoqi-y[1]+x[1], scalarMont, qi, mredParams)
		z[2] = MRed(twoqi-y[2]+x[2], scalarMont, qi, mredParams)
		z[3] = MRed(twoqi-y[3]+x[3], scalarMont, qi, mredParams)
		z[4] = MRed(twoqi-y[4]+x[4], scalarMont, qi, mredParams)
		z[5] = MRed(twoqi-y[5]+x[5], scalarMont, qi, mredParams)
		z[6] = MRed(twoqi-y[6]+x[6], scalarMont, qi, mredParams)
		z[7] = MRed(twoqi-y[7]+x[7], scalarMont, qi, mredParams)

	}
}

// MFormVec returns p2 = p1 * 2^64 mod qi.
func MFormVec(p1, p2 []uint64, qi uint64, bredParams []uint64) {
	for j := 0; j < len(p1); j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = MForm(x[0], qi, bredParams)
		z[1] = MForm(x[1], qi, bredParams)
		z[2] = MForm(x[2], qi, bredParams)
		z[3] = MForm(x[3], qi, bredParams)
		z[4] = MForm(x[4], qi, bredParams)
		z[5] = MForm(x[5], qi, bredParams)
		z[6] = MForm(x[6], qi, bredParams)
		z[7] = MForm(x[7], qi, bredParams)
	}
}

// MFormConstantVec returns p2 = p1 * 2^64 mod qi with result in the range [0, 2q-1]
func MFormConstantVec(p1, p2 []uint64, qi uint64, bredParams []uint64) {
	for j := 0; j < len(p1); j = j + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = MFormConstant(x[0], qi, bredParams)
		z[1] = MFormConstant(x[1], qi, bredParams)
		z[2] = MFormConstant(x[2], qi, bredParams)
		z[3] = MFormConstant(x[3], qi, bredParams)
		z[4] = MFormConstant(x[4], qi, bredParams)
		z[5] = MFormConstant(x[5], qi, bredParams)
		z[6] = MFormConstant(x[6], qi, bredParams)
		z[7] = MFormConstant(x[7], qi, bredParams)
	}
}

// InvMFormVec returns p2 = p1 * (2^64)^-1 mod qi.
func InvMFormVec(p1, p2 []uint64, qi, mredParams uint64) {
	for j := 0; j < len(p1); j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = InvMForm(x[0], qi, mredParams)
		z[1] = InvMForm(x[1], qi, mredParams)
		z[2] = InvMForm(x[2], qi, mredParams)
		z[3] = InvMForm(x[3], qi, mredParams)
		z[4] = InvMForm(x[4], qi, mredParams)
		z[5] = InvMForm(x[5], qi, mredParams)
		z[6] = InvMForm(x[6], qi, mredParams)
		z[7] = InvMForm(x[7], qi, mredParams)
	}
}

// MulByPow2Vec returns p2 = p1 * 2^pow2 mod qi.
func MulByPow2Vec(p1, p2 []uint64, pow2 int, qi, mredParams uint64) {
	for j := 0; j < len(p1); j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = PowerOf2(x[0], pow2, qi, mredParams)
		z[1] = PowerOf2(x[1], pow2, qi, mredParams)
		z[2] = PowerOf2(x[2], pow2, qi, mredParams)
		z[3] = PowerOf2(x[3], pow2, qi, mredParams)
		z[4] = PowerOf2(x[4], pow2, qi, mredParams)
		z[5] = PowerOf2(x[5], pow2, qi, mredParams)
		z[6] = PowerOf2(x[6], pow2, qi, mredParams)
		z[7] = PowerOf2(x[7], pow2, qi, mredParams)
	}
}

// ScaleUpVec takes a Poly pIn in ringT, scales its coefficients up by (Q/T) mod Q, and writes the result in a
// Poly pOut in ringQ.
func ScaleUpVec(ringQ, ringT *Ring, rescaleParams, tmp []uint64, pIn, pOut *Poly) {

	qModTmontgomery := MForm(new(big.Int).Mod(ringQ.ModulusBigint, ringT.ModulusBigint).Uint64(), ringT.Modulus[0], ringT.BredParams[0])

	t := ringT.Modulus[0]
	tHalf := t >> 1
	tInv := ringT.MredParams[0]

	// (x * Q + T/2) mod T
	for i := 0; i < ringQ.N; i = i + 8 {
		x := (*[8]uint64)(unsafe.Pointer(&pIn.Coeffs[0][i]))
		z := (*[8]uint64)(unsafe.Pointer(&tmp[i]))

		z[0] = CRed(MRed(x[0], qModTmontgomery, t, tInv)+tHalf, t)
		z[1] = CRed(MRed(x[1], qModTmontgomery, t, tInv)+tHalf, t)
		z[2] = CRed(MRed(x[2], qModTmontgomery, t, tInv)+tHalf, t)
		z[3] = CRed(MRed(x[3], qModTmontgomery, t, tInv)+tHalf, t)
		z[4] = CRed(MRed(x[4], qModTmontgomery, t, tInv)+tHalf, t)
		z[5] = CRed(MRed(x[5], qModTmontgomery, t, tInv)+tHalf, t)
		z[6] = CRed(MRed(x[6], qModTmontgomery, t, tInv)+tHalf, t)
		z[7] = CRed(MRed(x[7], qModTmontgomery, t, tInv)+tHalf, t)
	}

	// (x * T^-1 - T/2) mod Qi
	for i := 0; i < len(pOut.Coeffs); i++ {
		p0tmp := tmp
		p1tmp := pOut.Coeffs[i]
		qi := ringQ.Modulus[i]
		bredParams := ringQ.BredParams[i]
		mredParams := ringQ.MredParams[i]
		rescaleParams := qi - rescaleParams[i]

		tHalfNegQi := qi - BRedAdd(tHalf, qi, bredParams)

		for j := 0; j < ringQ.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p0tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))

			z[0] = MRed(x[0]+tHalfNegQi, rescaleParams, qi, mredParams)
			z[1] = MRed(x[1]+tHalfNegQi, rescaleParams, qi, mredParams)
			z[2] = MRed(x[2]+tHalfNegQi, rescaleParams, qi, mredParams)
			z[3] = MRed(x[3]+tHalfNegQi, rescaleParams, qi, mredParams)
			z[4] = MRed(x[4]+tHalfNegQi, rescaleParams, qi, mredParams)
			z[5] = MRed(x[5]+tHalfNegQi, rescaleParams, qi, mredParams)
			z[6] = MRed(x[6]+tHalfNegQi, rescaleParams, qi, mredParams)
			z[7] = MRed(x[7]+tHalfNegQi, rescaleParams, qi, mredParams)
		}
	}
}

// SpecialFFTVec performs the CKKS special FFT transform in place.
func SpecialFFTVec(values []complex128, N, M int, rotGroup []int, roots []complex128) {

	SliceBitReverseInPlaceComplex128(values, N)

	logN := int(bits.Len64(uint64(N))) - 1
	logM := int(bits.Len64(uint64(M))) - 1

	for loglen := 1; loglen <= logN; loglen++ {
		len := 1 << loglen
		lenh := len >> 1
		lenq := len << 2
		logGap := logM - 2 - loglen
		mask := lenq - 1

		if lenh > 8 {
			for i := 0; i < N; i += len {

				for j, k := 0, i; j < lenh; j, k = j+8, k+8 {

					u := (*[8]complex128)(unsafe.Pointer(&values[k]))
					v := (*[8]complex128)(unsafe.Pointer(&values[k+lenh]))
					w := (*[8]int)(unsafe.Pointer(&rotGroup[j]))

					v[0] *= roots[(w[0]&mask)<<logGap]
					v[1] *= roots[(w[1]&mask)<<logGap]
					v[2] *= roots[(w[2]&mask)<<logGap]
					v[3] *= roots[(w[3]&mask)<<logGap]
					v[4] *= roots[(w[4]&mask)<<logGap]
					v[5] *= roots[(w[5]&mask)<<logGap]
					v[6] *= roots[(w[6]&mask)<<logGap]
					v[7] *= roots[(w[7]&mask)<<logGap]

					u[0], v[0] = u[0]+v[0], u[0]-v[0]
					u[1], v[1] = u[1]+v[1], u[1]-v[1]
					u[2], v[2] = u[2]+v[2], u[2]-v[2]
					u[3], v[3] = u[3]+v[3], u[3]-v[3]
					u[4], v[4] = u[4]+v[4], u[4]-v[4]
					u[5], v[5] = u[5]+v[5], u[5]-v[5]
					u[6], v[6] = u[6]+v[6], u[6]-v[6]
					u[7], v[7] = u[7]+v[7], u[7]-v[7]
				}
			}
		} else if lenh == 8 {

			psi0 := roots[(rotGroup[0]&mask)<<logGap]
			psi1 := roots[(rotGroup[1]&mask)<<logGap]
			psi2 := roots[(rotGroup[2]&mask)<<logGap]
			psi3 := roots[(rotGroup[3]&mask)<<logGap]
			psi4 := roots[(rotGroup[4]&mask)<<logGap]
			psi5 := roots[(rotGroup[5]&mask)<<logGap]
			psi6 := roots[(rotGroup[6]&mask)<<logGap]
			psi7 := roots[(rotGroup[7]&mask)<<logGap]

			for i := 0; i < N; i += 16 {

				u := (*[16]complex128)(unsafe.Pointer(&values[i]))

				u[8] *= psi0
				u[9] *= psi1
				u[10] *= psi2
				u[11] *= psi3
				u[12] *= psi4
				u[13] *= psi5
				u[14] *= psi6
				u[15] *= psi7

				u[0], u[8] = u[0]+u[8], u[0]-u[8]
				u[1], u[9] = u[1]+u[9], u[1]-u[9]
				u[2], u[10] = u[2]+u[10], u[2]-u[10]
				u[3], u[11] = u[3]+u[11], u[3]-u[11]
				u[4], u[12] = u[4]+u[12], u[4]-u[12]
				u[5], u[13] = u[5]+u[13], u[5]-u[13]
				u[6], u[14] = u[6]+u[14], u[6]-u[14]
				u[7], u[15] = u[7]+u[15], u[7]-u[15]

			}
		} else if lenh == 4 {

			psi0 := roots[(rotGroup[0]&mask)<<logGap]
			psi1 := roots[(rotGroup[1]&mask)<<logGap]
			psi2 := roots[(rotGroup[2]&mask)<<logGap]
			psi3 := roots[(rotGroup[3]&mask)<<logGap]

			for i := 0; i < N; i += 16 {

				u := (*[16]complex128)(unsafe.Pointer(&values[i]))

				u[4] *= psi0
				u[5] *= psi1
				u[6] *= psi2
				u[7] *= psi3
				u[12] *= psi0
				u[13] *= psi1
				u[14] *= psi2
				u[15] *= psi3

				u[0], u[4] = u[0]+u[4], u[0]-u[4]
				u[1], u[5] = u[1]+u[5], u[1]-u[5]
				u[2], u[6] = u[2]+u[6], u[2]-u[6]
				u[3], u[7] = u[3]+u[7], u[3]-u[7]
				u[8], u[12] = u[8]+u[12], u[8]-u[12]
				u[9], u[13] = u[9]+u[13], u[9]-u[13]
				u[10], u[14] = u[10]+u[14], u[10]-u[14]
				u[11], u[15] = u[11]+u[15], u[11]-u[15]
			}
		} else if lenh == 2 {

			psi0 := roots[(rotGroup[0]&mask)<<logGap]
			psi1 := roots[(rotGroup[1]&mask)<<logGap]

			for i := 0; i < N; i += 16 {

				u := (*[16]complex128)(unsafe.Pointer(&values[i]))

				u[2] *= psi0
				u[3] *= psi1
				u[6] *= psi0
				u[7] *= psi1
				u[10] *= psi0
				u[11] *= psi1
				u[14] *= psi0
				u[15] *= psi1

				u[0], u[2] = u[0]+u[2], u[0]-u[2]
				u[1], u[3] = u[1]+u[3], u[1]-u[3]
				u[4], u[6] = u[4]+u[6], u[4]-u[6]
				u[5], u[7] = u[5]+u[7], u[5]-u[7]
				u[8], u[10] = u[8]+u[10], u[8]-u[10]
				u[9], u[11] = u[9]+u[11], u[9]-u[11]
				u[12], u[14] = u[12]+u[14], u[12]-u[14]
				u[13], u[15] = u[13]+u[15], u[13]-u[15]
			}
		} else if lenh == 1 {

			psi0 := roots[(rotGroup[0]&mask)<<logGap]

			for i := 0; i < N; i += 16 {

				u := (*[16]complex128)(unsafe.Pointer(&values[i]))

				u[1] *= psi0
				u[3] *= psi0
				u[5] *= psi0
				u[7] *= psi0
				u[9] *= psi0
				u[11] *= psi0
				u[13] *= psi0
				u[15] *= psi0

				u[0], u[1] = u[0]+u[1], u[0]-u[1]
				u[2], u[3] = u[2]+u[3], u[2]-u[3]
				u[4], u[5] = u[4]+u[5], u[4]-u[5]
				u[6], u[7] = u[6]+u[7], u[6]-u[7]
				u[8], u[9] = u[8]+u[9], u[8]-u[9]
				u[10], u[11] = u[10]+u[11], u[10]-u[11]
				u[12], u[13] = u[12]+u[13], u[12]-u[13]
				u[14], u[15] = u[14]+u[15], u[14]-u[15]
			}
		}
	}
}

// SpecialInvFFTVec performs the CKKS special inverse FFT transform in place.
func SpecialInvFFTVec(values []complex128, N, M int, rotGroup []int, roots []complex128) {

	var lenh, lenq, mask, logGap int

	logN := int(bits.Len64(uint64(N))) - 1
	logM := int(bits.Len64(uint64(M))) - 1

	for loglen := logN; loglen > 0; loglen-- {
		len := 1 << loglen
		lenh = len >> 1
		lenq = len << 2
		logGap = logM - 2 - loglen
		mask = lenq - 1

		if lenh > 8 {
			for i := 0; i < N; i += len {
				for j, k := 0, i; j < lenh; j, k = j+8, k+8 {

					u := (*[8]complex128)(unsafe.Pointer(&values[k]))
					v := (*[8]complex128)(unsafe.Pointer(&values[k+lenh]))
					w := (*[8]int)(unsafe.Pointer(&rotGroup[j]))

					u[0], v[0] = u[0]+v[0], (u[0]-v[0])*roots[(lenq-(w[0]&mask))<<logGap]
					u[1], v[1] = u[1]+v[1], (u[1]-v[1])*roots[(lenq-(w[1]&mask))<<logGap]
					u[2], v[2] = u[2]+v[2], (u[2]-v[2])*roots[(lenq-(w[2]&mask))<<logGap]
					u[3], v[3] = u[3]+v[3], (u[3]-v[3])*roots[(lenq-(w[3]&mask))<<logGap]
					u[4], v[4] = u[4]+v[4], (u[4]-v[4])*roots[(lenq-(w[4]&mask))<<logGap]
					u[5], v[5] = u[5]+v[5], (u[5]-v[5])*roots[(lenq-(w[5]&mask))<<logGap]
					u[6], v[6] = u[6]+v[6], (u[6]-v[6])*roots[(lenq-(w[6]&mask))<<logGap]
					u[7], v[7] = u[7]+v[7], (u[7]-v[7])*roots[(lenq-(w[7]&mask))<<logGap]
				}
			}
		} else if lenh == 8 {

			psi0 := roots[(lenq-(rotGroup[0]&mask))<<logGap]
			psi1 := roots[(lenq-(rotGroup[1]&mask))<<logGap]
			psi2 := roots[(lenq-(rotGroup[2]&mask))<<logGap]
			psi3 := roots[(lenq-(rotGroup[3]&mask))<<logGap]
			psi4 := roots[(lenq-(rotGroup[4]&mask))<<logGap]
			psi5 := roots[(lenq-(rotGroup[5]&mask))<<logGap]
			psi6 := roots[(lenq-(rotGroup[6]&mask))<<logGap]
			psi7 := roots[(lenq-(rotGroup[7]&mask))<<logGap]

			for i := 0; i < N; i += 16 {

				u := (*[16]complex128)(unsafe.Pointer(&values[i]))

				u[0], u[8] = u[0]+u[8], (u[0]-u[8])*psi0
				u[1], u[9] = u[1]+u[9], (u[1]-u[9])*psi1
				u[2], u[10] = u[2]+u[10], (u[2]-u[10])*psi2
				u[3], u[11] = u[3]+u[11], (u[3]-u[11])*psi3
				u[4], u[12] = u[4]+u[12], (u[4]-u[12])*psi4
				u[5], u[13] = u[5]+u[13], (u[5]-u[13])*psi5
				u[6], u[14] = u[6]+u[14], (u[6]-u[14])*psi6
				u[7], u[15] = u[7]+u[15], (u[7]-u[15])*psi7
			}

		} else if lenh == 4 {

			psi0 := roots[(lenq-(rotGroup[0]&mask))<<logGap]
			psi1 := roots[(lenq-(rotGroup[1]&mask))<<logGap]
			psi2 := roots[(lenq-(rotGroup[2]&mask))<<logGap]
			psi3 := roots[(lenq-(rotGroup[3]&mask))<<logGap]

			for i := 0; i < N; i += 16 {

				u := (*[16]complex128)(unsafe.Pointer(&values[i]))

				u[0], u[4] = u[0]+u[4], (u[0]-u[4])*psi0
				u[1], u[5] = u[1]+u[5], (u[1]-u[5])*psi1
				u[2], u[6] = u[2]+u[6], (u[2]-u[6])*psi2
				u[3], u[7] = u[3]+u[7], (u[3]-u[7])*psi3
				u[8], u[12] = u[8]+u[12], (u[8]-u[12])*psi0
				u[9], u[13] = u[9]+u[13], (u[9]-u[13])*psi1
				u[10], u[14] = u[10]+u[14], (u[10]-u[14])*psi2
				u[11], u[15] = u[11]+u[15], (u[11]-u[15])*psi3
			}
		} else if lenh == 2 {

			psi0 := roots[(lenq-(rotGroup[0]&mask))<<logGap]
			psi1 := roots[(lenq-(rotGroup[1]&mask))<<logGap]

			for i := 0; i < N; i += 16 {

				u := (*[16]complex128)(unsafe.Pointer(&values[i]))

				u[0], u[2] = u[0]+u[2], (u[0]-u[2])*psi0
				u[1], u[3] = u[1]+u[3], (u[1]-u[3])*psi1
				u[4], u[6] = u[4]+u[6], (u[4]-u[6])*psi0
				u[5], u[7] = u[5]+u[7], (u[5]-u[7])*psi1
				u[8], u[10] = u[8]+u[10], (u[8]-u[10])*psi0
				u[9], u[11] = u[9]+u[11], (u[9]-u[11])*psi1
				u[12], u[14] = u[12]+u[14], (u[12]-u[14])*psi0
				u[13], u[15] = u[13]+u[15], (u[13]-u[15])*psi1
			}
		} else if lenh == 1 {

			psi0 := roots[(lenq-(rotGroup[0]&mask))<<logGap]

			for i := 0; i < N; i += 16 {

				u := (*[16]complex128)(unsafe.Pointer(&values[i]))

				u[0], u[1] = u[0]+u[1], (u[0]-u[1])*psi0
				u[2], u[3] = u[2]+u[3], (u[2]-u[3])*psi0
				u[4], u[5] = u[4]+u[5], (u[4]-u[5])*psi0
				u[6], u[7] = u[6]+u[7], (u[6]-u[7])*psi0
				u[8], u[9] = u[8]+u[9], (u[8]-u[9])*psi0
				u[10], u[11] = u[10]+u[11], (u[10]-u[11])*psi0
				u[12], u[13] = u[12]+u[13], (u[12]-u[13])*psi0
				u[14], u[15] = u[14]+u[15], (u[14]-u[15])*psi0
			}
		}
	}

	DivideComplex128SliceVec(values, complex(float64(N), 0))

	SliceBitReverseInPlaceComplex128(values, N)
}

// DivideComplex128SliceVec divides the entries in values by scaleVal in place.
func DivideComplex128SliceVec(values []complex128, scaleVal complex128) {
	lenValues := len(values)
	for i := 0; i < lenValues; i = i + 8 {

		v := (*[8]complex128)(unsafe.Pointer(&values[i]))

		v[0] /= scaleVal
		v[1] /= scaleVal
		v[2] /= scaleVal
		v[3] /= scaleVal
		v[4] /= scaleVal
		v[5] /= scaleVal
		v[6] /= scaleVal
		v[7] /= scaleVal
	}
}
