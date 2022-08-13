package ring

import (
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

// AddScalarNoModAndMulScalarMontgomeryVec returns p3 = (scalarMont0+p2)*scalarMont1 mod qi.
func AddScalarNoModAndMulScalarMontgomeryVec(p1, p2 []uint64, scalar0, scalarMont1, qi, mredParams uint64) {
	for j := 0; j < len(p1); j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = MRed(x[0]+scalar0, scalarMont1, qi, mredParams)
		z[1] = MRed(x[1]+scalar0, scalarMont1, qi, mredParams)
		z[2] = MRed(x[2]+scalar0, scalarMont1, qi, mredParams)
		z[3] = MRed(x[3]+scalar0, scalarMont1, qi, mredParams)
		z[4] = MRed(x[4]+scalar0, scalarMont1, qi, mredParams)
		z[5] = MRed(x[5]+scalar0, scalarMont1, qi, mredParams)
		z[6] = MRed(x[6]+scalar0, scalarMont1, qi, mredParams)
		z[7] = MRed(x[7]+scalar0, scalarMont1, qi, mredParams)
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

// MulScalarMontgomeryAndAddScalarVec returns p2 = scalar + p1*scalarMont mod qi.
func MulScalarMontgomeryAndAddScalarVec(p1, p2 []uint64, scalar0, scalarMont1, qi, mredParams uint64) {
	for j := 0; j < len(p1); j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = CRed(MRed(x[0], scalarMont1, qi, mredParams)+scalar0, qi)
		z[1] = CRed(MRed(x[1], scalarMont1, qi, mredParams)+scalar0, qi)
		z[2] = CRed(MRed(x[2], scalarMont1, qi, mredParams)+scalar0, qi)
		z[3] = CRed(MRed(x[3], scalarMont1, qi, mredParams)+scalar0, qi)
		z[4] = CRed(MRed(x[4], scalarMont1, qi, mredParams)+scalar0, qi)
		z[5] = CRed(MRed(x[5], scalarMont1, qi, mredParams)+scalar0, qi)
		z[6] = CRed(MRed(x[6], scalarMont1, qi, mredParams)+scalar0, qi)
		z[7] = CRed(MRed(x[7], scalarMont1, qi, mredParams)+scalar0, qi)
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

// ZeroVec sets all values of p1 to zero.
func ZeroVec(p1 []uint64) {
	for j := 0; j < len(p1); j = j + 8 {

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
