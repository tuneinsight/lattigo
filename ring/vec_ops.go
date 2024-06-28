package ring

import (
	"unsafe"
)

func addvec(p1, p2, p3 []uint64, modulus uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p3)%8 */
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

func addlazyvec(p1, p2, p3 []uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p3)%8 */
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

func subvec(p1, p2, p3 []uint64, modulus uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p3)%8 */
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

func sublazyvec(p1, p2, p3 []uint64, modulus uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p3)%8 */
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

func negvec(p1, p2 []uint64, modulus uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
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

func reducevec(p1, p2 []uint64, modulus uint64, bredconstant [2]uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = BRedAdd(x[0], modulus, bredconstant)
		z[1] = BRedAdd(x[1], modulus, bredconstant)
		z[2] = BRedAdd(x[2], modulus, bredconstant)
		z[3] = BRedAdd(x[3], modulus, bredconstant)
		z[4] = BRedAdd(x[4], modulus, bredconstant)
		z[5] = BRedAdd(x[5], modulus, bredconstant)
		z[6] = BRedAdd(x[6], modulus, bredconstant)
		z[7] = BRedAdd(x[7], modulus, bredconstant)
	}
}

func reducelazyvec(p1, p2 []uint64, modulus uint64, bredconstant [2]uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = BRedAddLazy(x[0], modulus, bredconstant)
		z[1] = BRedAddLazy(x[1], modulus, bredconstant)
		z[2] = BRedAddLazy(x[2], modulus, bredconstant)
		z[3] = BRedAddLazy(x[3], modulus, bredconstant)
		z[4] = BRedAddLazy(x[4], modulus, bredconstant)
		z[5] = BRedAddLazy(x[5], modulus, bredconstant)
		z[6] = BRedAddLazy(x[6], modulus, bredconstant)
		z[7] = BRedAddLazy(x[7], modulus, bredconstant)
	}
}

func mulcoeffslazyvec(p1, p2, p3 []uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p3)%8 */
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

func mulcoeffslazythenaddlazyvec(p1, p2, p3 []uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p3)%8 */
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

func mulcoeffsbarrettvec(p1, p2, p3 []uint64, modulus uint64, bredconstant [2]uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p3)%8 */
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = BRed(x[0], y[0], modulus, bredconstant)
		z[1] = BRed(x[1], y[1], modulus, bredconstant)
		z[2] = BRed(x[2], y[2], modulus, bredconstant)
		z[3] = BRed(x[3], y[3], modulus, bredconstant)
		z[4] = BRed(x[4], y[4], modulus, bredconstant)
		z[5] = BRed(x[5], y[5], modulus, bredconstant)
		z[6] = BRed(x[6], y[6], modulus, bredconstant)
		z[7] = BRed(x[7], y[7], modulus, bredconstant)
	}
}

func mulcoeffsbarrettlazyvec(p1, p2, p3 []uint64, modulus uint64, bredconstant [2]uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p3)%8 */
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = BRedLazy(x[0], y[0], modulus, bredconstant)
		z[1] = BRedLazy(x[1], y[1], modulus, bredconstant)
		z[2] = BRedLazy(x[2], y[2], modulus, bredconstant)
		z[3] = BRedLazy(x[3], y[3], modulus, bredconstant)
		z[4] = BRedLazy(x[4], y[4], modulus, bredconstant)
		z[5] = BRedLazy(x[5], y[5], modulus, bredconstant)
		z[6] = BRedLazy(x[6], y[6], modulus, bredconstant)
		z[7] = BRedLazy(x[7], y[7], modulus, bredconstant)
	}
}

func mulcoeffsthenaddvec(p1, p2, p3 []uint64, modulus uint64, bredconstant [2]uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p3)%8 */
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = CRed(z[0]+BRed(x[0], y[0], modulus, bredconstant), modulus)
		z[1] = CRed(z[1]+BRed(x[1], y[1], modulus, bredconstant), modulus)
		z[2] = CRed(z[2]+BRed(x[2], y[2], modulus, bredconstant), modulus)
		z[3] = CRed(z[3]+BRed(x[3], y[3], modulus, bredconstant), modulus)
		z[4] = CRed(z[4]+BRed(x[4], y[4], modulus, bredconstant), modulus)
		z[5] = CRed(z[5]+BRed(x[5], y[5], modulus, bredconstant), modulus)
		z[6] = CRed(z[6]+BRed(x[6], y[6], modulus, bredconstant), modulus)
		z[7] = CRed(z[7]+BRed(x[7], y[7], modulus, bredconstant), modulus)
	}
}

func mulcoeffsbarrettthenaddlazyvec(p1, p2, p3 []uint64, modulus uint64, bredconstant [2]uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p3)%8 */
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] += BRed(x[0], y[0], modulus, bredconstant)
		z[1] += BRed(x[1], y[1], modulus, bredconstant)
		z[2] += BRed(x[2], y[2], modulus, bredconstant)
		z[3] += BRed(x[3], y[3], modulus, bredconstant)
		z[4] += BRed(x[4], y[4], modulus, bredconstant)
		z[5] += BRed(x[5], y[5], modulus, bredconstant)
		z[6] += BRed(x[6], y[6], modulus, bredconstant)
		z[7] += BRed(x[7], y[7], modulus, bredconstant)
	}
}

func mulcoeffsmontgomeryvec(p1, p2, p3 []uint64, modulus, mredconstant uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p3)%8 */
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = MRed(x[0], y[0], modulus, mredconstant)
		z[1] = MRed(x[1], y[1], modulus, mredconstant)
		z[2] = MRed(x[2], y[2], modulus, mredconstant)
		z[3] = MRed(x[3], y[3], modulus, mredconstant)
		z[4] = MRed(x[4], y[4], modulus, mredconstant)
		z[5] = MRed(x[5], y[5], modulus, mredconstant)
		z[6] = MRed(x[6], y[6], modulus, mredconstant)
		z[7] = MRed(x[7], y[7], modulus, mredconstant)
	}
}

func mulcoeffsmontgomerylazyvec(p1, p2, p3 []uint64, modulus, mredconstant uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p3)%8 */
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = MRedLazy(x[0], y[0], modulus, mredconstant)
		z[1] = MRedLazy(x[1], y[1], modulus, mredconstant)
		z[2] = MRedLazy(x[2], y[2], modulus, mredconstant)
		z[3] = MRedLazy(x[3], y[3], modulus, mredconstant)
		z[4] = MRedLazy(x[4], y[4], modulus, mredconstant)
		z[5] = MRedLazy(x[5], y[5], modulus, mredconstant)
		z[6] = MRedLazy(x[6], y[6], modulus, mredconstant)
		z[7] = MRedLazy(x[7], y[7], modulus, mredconstant)
	}
}

func mulcoeffsmontgomerythenaddvec(p1, p2, p3 []uint64, modulus, mredconstant uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p3)%8 */
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = CRed(z[0]+MRed(x[0], y[0], modulus, mredconstant), modulus)
		z[1] = CRed(z[1]+MRed(x[1], y[1], modulus, mredconstant), modulus)
		z[2] = CRed(z[2]+MRed(x[2], y[2], modulus, mredconstant), modulus)
		z[3] = CRed(z[3]+MRed(x[3], y[3], modulus, mredconstant), modulus)
		z[4] = CRed(z[4]+MRed(x[4], y[4], modulus, mredconstant), modulus)
		z[5] = CRed(z[5]+MRed(x[5], y[5], modulus, mredconstant), modulus)
		z[6] = CRed(z[6]+MRed(x[6], y[6], modulus, mredconstant), modulus)
		z[7] = CRed(z[7]+MRed(x[7], y[7], modulus, mredconstant), modulus)
	}
}

func mulcoeffsmontgomerythenaddlazyvec(p1, p2, p3 []uint64, modulus, mredconstant uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p3)%8 */
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] += MRed(x[0], y[0], modulus, mredconstant)
		z[1] += MRed(x[1], y[1], modulus, mredconstant)
		z[2] += MRed(x[2], y[2], modulus, mredconstant)
		z[3] += MRed(x[3], y[3], modulus, mredconstant)
		z[4] += MRed(x[4], y[4], modulus, mredconstant)
		z[5] += MRed(x[5], y[5], modulus, mredconstant)
		z[6] += MRed(x[6], y[6], modulus, mredconstant)
		z[7] += MRed(x[7], y[7], modulus, mredconstant)
	}
}

func mulcoeffsmontgomerylazythenaddlazyvec(p1, p2, p3 []uint64, modulus, mredconstant uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p3)%8 */
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] += MRedLazy(x[0], y[0], modulus, mredconstant)
		z[1] += MRedLazy(x[1], y[1], modulus, mredconstant)
		z[2] += MRedLazy(x[2], y[2], modulus, mredconstant)
		z[3] += MRedLazy(x[3], y[3], modulus, mredconstant)
		z[4] += MRedLazy(x[4], y[4], modulus, mredconstant)
		z[5] += MRedLazy(x[5], y[5], modulus, mredconstant)
		z[6] += MRedLazy(x[6], y[6], modulus, mredconstant)
		z[7] += MRedLazy(x[7], y[7], modulus, mredconstant)
	}
}

func mulcoeffsmontgomerythensubvec(p1, p2, p3 []uint64, modulus, mredconstant uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p3)%8 */
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = CRed(z[0]+(modulus-MRed(x[0], y[0], modulus, mredconstant)), modulus)
		z[1] = CRed(z[1]+(modulus-MRed(x[1], y[1], modulus, mredconstant)), modulus)
		z[2] = CRed(z[2]+(modulus-MRed(x[2], y[2], modulus, mredconstant)), modulus)
		z[3] = CRed(z[3]+(modulus-MRed(x[3], y[3], modulus, mredconstant)), modulus)
		z[4] = CRed(z[4]+(modulus-MRed(x[4], y[4], modulus, mredconstant)), modulus)
		z[5] = CRed(z[5]+(modulus-MRed(x[5], y[5], modulus, mredconstant)), modulus)
		z[6] = CRed(z[6]+(modulus-MRed(x[6], y[6], modulus, mredconstant)), modulus)
		z[7] = CRed(z[7]+(modulus-MRed(x[7], y[7], modulus, mredconstant)), modulus)
	}
}

func mulcoeffsmontgomerythensublazyvec(p1, p2, p3 []uint64, modulus, mredconstant uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p3)%8 */
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] += (modulus - MRed(x[0], y[0], modulus, mredconstant))
		z[1] += (modulus - MRed(x[1], y[1], modulus, mredconstant))
		z[2] += (modulus - MRed(x[2], y[2], modulus, mredconstant))
		z[3] += (modulus - MRed(x[3], y[3], modulus, mredconstant))
		z[4] += (modulus - MRed(x[4], y[4], modulus, mredconstant))
		z[5] += (modulus - MRed(x[5], y[5], modulus, mredconstant))
		z[6] += (modulus - MRed(x[6], y[6], modulus, mredconstant))
		z[7] += (modulus - MRed(x[7], y[7], modulus, mredconstant))
	}
}

func mulcoeffsmontgomerylazythensublazyvec(p1, p2, p3 []uint64, modulus, mredconstant uint64) {

	N := len(p1)
	twomodulus := modulus << 1

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p3)%8 */
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] += twomodulus - MRedLazy(x[0], y[0], modulus, mredconstant)
		z[1] += twomodulus - MRedLazy(x[1], y[1], modulus, mredconstant)
		z[2] += twomodulus - MRedLazy(x[2], y[2], modulus, mredconstant)
		z[3] += twomodulus - MRedLazy(x[3], y[3], modulus, mredconstant)
		z[4] += twomodulus - MRedLazy(x[4], y[4], modulus, mredconstant)
		z[5] += twomodulus - MRedLazy(x[5], y[5], modulus, mredconstant)
		z[6] += twomodulus - MRedLazy(x[6], y[6], modulus, mredconstant)
		z[7] += twomodulus - MRedLazy(x[7], y[7], modulus, mredconstant)
	}
}

func mulcoeffsmontgomerylazythenNegvec(p1, p2, p3 []uint64, modulus, mredconstant uint64) {

	N := len(p1)
	twomodulus := modulus << 1

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p3)%8 */
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = twomodulus - MRedLazy(x[0], y[0], modulus, mredconstant)
		z[1] = twomodulus - MRedLazy(x[1], y[1], modulus, mredconstant)
		z[2] = twomodulus - MRedLazy(x[2], y[2], modulus, mredconstant)
		z[3] = twomodulus - MRedLazy(x[3], y[3], modulus, mredconstant)
		z[4] = twomodulus - MRedLazy(x[4], y[4], modulus, mredconstant)
		z[5] = twomodulus - MRedLazy(x[5], y[5], modulus, mredconstant)
		z[6] = twomodulus - MRedLazy(x[6], y[6], modulus, mredconstant)
		z[7] = twomodulus - MRedLazy(x[7], y[7], modulus, mredconstant)
	}
}

func addlazythenmulscalarmontgomeryvec(p1, p2 []uint64, scalarMont uint64, p3 []uint64, modulus, mredconstant uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p3)%8 */
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = MRed(x[0]+y[0], scalarMont, modulus, mredconstant)
		z[1] = MRed(x[1]+y[1], scalarMont, modulus, mredconstant)
		z[2] = MRed(x[2]+y[2], scalarMont, modulus, mredconstant)
		z[3] = MRed(x[3]+y[3], scalarMont, modulus, mredconstant)
		z[4] = MRed(x[4]+y[4], scalarMont, modulus, mredconstant)
		z[5] = MRed(x[5]+y[5], scalarMont, modulus, mredconstant)
		z[6] = MRed(x[6]+y[6], scalarMont, modulus, mredconstant)
		z[7] = MRed(x[7]+y[7], scalarMont, modulus, mredconstant)
	}
}

func addscalarlazythenmulscalarmontgomeryvec(p1 []uint64, scalar0, scalarMont1 uint64, p2 []uint64, modulus, mredconstant uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = MRed(x[0]+scalar0, scalarMont1, modulus, mredconstant)
		z[1] = MRed(x[1]+scalar0, scalarMont1, modulus, mredconstant)
		z[2] = MRed(x[2]+scalar0, scalarMont1, modulus, mredconstant)
		z[3] = MRed(x[3]+scalar0, scalarMont1, modulus, mredconstant)
		z[4] = MRed(x[4]+scalar0, scalarMont1, modulus, mredconstant)
		z[5] = MRed(x[5]+scalar0, scalarMont1, modulus, mredconstant)
		z[6] = MRed(x[6]+scalar0, scalarMont1, modulus, mredconstant)
		z[7] = MRed(x[7]+scalar0, scalarMont1, modulus, mredconstant)
	}
}

func addscalarvec(p1 []uint64, scalar uint64, p2 []uint64, modulus uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
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

func addscalarlazyvec(p1 []uint64, scalar uint64, p2 []uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
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

func addscalarlazythenNegTwoModuluslazyvec(p1 []uint64, scalar uint64, p2 []uint64, modulus uint64) {

	N := len(p1)
	twomodulus := modulus << 1

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
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

func subscalarvec(p1 []uint64, scalar uint64, p2 []uint64, modulus uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
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

func mulscalarmontgomeryvec(p1 []uint64, scalarMont uint64, p2 []uint64, modulus, mredconstant uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = MRed(x[0], scalarMont, modulus, mredconstant)
		z[1] = MRed(x[1], scalarMont, modulus, mredconstant)
		z[2] = MRed(x[2], scalarMont, modulus, mredconstant)
		z[3] = MRed(x[3], scalarMont, modulus, mredconstant)
		z[4] = MRed(x[4], scalarMont, modulus, mredconstant)
		z[5] = MRed(x[5], scalarMont, modulus, mredconstant)
		z[6] = MRed(x[6], scalarMont, modulus, mredconstant)
		z[7] = MRed(x[7], scalarMont, modulus, mredconstant)
	}
}

func mulscalarmontgomerylazyvec(p1 []uint64, scalarMont uint64, p2 []uint64, modulus, mredconstant uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = MRedLazy(x[0], scalarMont, modulus, mredconstant)
		z[1] = MRedLazy(x[1], scalarMont, modulus, mredconstant)
		z[2] = MRedLazy(x[2], scalarMont, modulus, mredconstant)
		z[3] = MRedLazy(x[3], scalarMont, modulus, mredconstant)
		z[4] = MRedLazy(x[4], scalarMont, modulus, mredconstant)
		z[5] = MRedLazy(x[5], scalarMont, modulus, mredconstant)
		z[6] = MRedLazy(x[6], scalarMont, modulus, mredconstant)
		z[7] = MRedLazy(x[7], scalarMont, modulus, mredconstant)
	}
}

func mulscalarmontgomerythenaddvec(p1 []uint64, scalarMont uint64, p2 []uint64, modulus, mredconstant uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = CRed(z[0]+MRed(x[0], scalarMont, modulus, mredconstant), modulus)
		z[1] = CRed(z[1]+MRed(x[1], scalarMont, modulus, mredconstant), modulus)
		z[2] = CRed(z[2]+MRed(x[2], scalarMont, modulus, mredconstant), modulus)
		z[3] = CRed(z[3]+MRed(x[3], scalarMont, modulus, mredconstant), modulus)
		z[4] = CRed(z[4]+MRed(x[4], scalarMont, modulus, mredconstant), modulus)
		z[5] = CRed(z[5]+MRed(x[5], scalarMont, modulus, mredconstant), modulus)
		z[6] = CRed(z[6]+MRed(x[6], scalarMont, modulus, mredconstant), modulus)
		z[7] = CRed(z[7]+MRed(x[7], scalarMont, modulus, mredconstant), modulus)
	}
}

func mulscalarmontgomerythenaddscalarvec(p1 []uint64, scalar0, scalarMont1 uint64, p2 []uint64, modulus, mredconstant uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = CRed(MRed(x[0], scalarMont1, modulus, mredconstant)+scalar0, modulus)
		z[1] = CRed(MRed(x[1], scalarMont1, modulus, mredconstant)+scalar0, modulus)
		z[2] = CRed(MRed(x[2], scalarMont1, modulus, mredconstant)+scalar0, modulus)
		z[3] = CRed(MRed(x[3], scalarMont1, modulus, mredconstant)+scalar0, modulus)
		z[4] = CRed(MRed(x[4], scalarMont1, modulus, mredconstant)+scalar0, modulus)
		z[5] = CRed(MRed(x[5], scalarMont1, modulus, mredconstant)+scalar0, modulus)
		z[6] = CRed(MRed(x[6], scalarMont1, modulus, mredconstant)+scalar0, modulus)
		z[7] = CRed(MRed(x[7], scalarMont1, modulus, mredconstant)+scalar0, modulus)
	}
}

func subthenmulscalarmontgomeryTwoModulusvec(p1, p2 []uint64, scalarMont uint64, p3 []uint64, modulus, mredconstant uint64) {

	N := len(p1)
	twomodulus := modulus << 1

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		y := (*[8]uint64)(unsafe.Pointer(&p2[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p3)%8 */
		z := (*[8]uint64)(unsafe.Pointer(&p3[j]))

		z[0] = MRed(twomodulus-y[0]+x[0], scalarMont, modulus, mredconstant)
		z[1] = MRed(twomodulus-y[1]+x[1], scalarMont, modulus, mredconstant)
		z[2] = MRed(twomodulus-y[2]+x[2], scalarMont, modulus, mredconstant)
		z[3] = MRed(twomodulus-y[3]+x[3], scalarMont, modulus, mredconstant)
		z[4] = MRed(twomodulus-y[4]+x[4], scalarMont, modulus, mredconstant)
		z[5] = MRed(twomodulus-y[5]+x[5], scalarMont, modulus, mredconstant)
		z[6] = MRed(twomodulus-y[6]+x[6], scalarMont, modulus, mredconstant)
		z[7] = MRed(twomodulus-y[7]+x[7], scalarMont, modulus, mredconstant)

	}
}

func mformvec(p1, p2 []uint64, modulus uint64, bredconstant [2]uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = MForm(x[0], modulus, bredconstant)
		z[1] = MForm(x[1], modulus, bredconstant)
		z[2] = MForm(x[2], modulus, bredconstant)
		z[3] = MForm(x[3], modulus, bredconstant)
		z[4] = MForm(x[4], modulus, bredconstant)
		z[5] = MForm(x[5], modulus, bredconstant)
		z[6] = MForm(x[6], modulus, bredconstant)
		z[7] = MForm(x[7], modulus, bredconstant)
	}
}

func mformlazyvec(p1, p2 []uint64, modulus uint64, bredconstant [2]uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = MFormLazy(x[0], modulus, bredconstant)
		z[1] = MFormLazy(x[1], modulus, bredconstant)
		z[2] = MFormLazy(x[2], modulus, bredconstant)
		z[3] = MFormLazy(x[3], modulus, bredconstant)
		z[4] = MFormLazy(x[4], modulus, bredconstant)
		z[5] = MFormLazy(x[5], modulus, bredconstant)
		z[6] = MFormLazy(x[6], modulus, bredconstant)
		z[7] = MFormLazy(x[7], modulus, bredconstant)
	}
}

func imformvec(p1, p2 []uint64, modulus, mredconstant uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
		z := (*[8]uint64)(unsafe.Pointer(&p2[j]))

		z[0] = IMForm(x[0], modulus, mredconstant)
		z[1] = IMForm(x[1], modulus, mredconstant)
		z[2] = IMForm(x[2], modulus, mredconstant)
		z[3] = IMForm(x[3], modulus, mredconstant)
		z[4] = IMForm(x[4], modulus, mredconstant)
		z[5] = IMForm(x[5], modulus, mredconstant)
		z[6] = IMForm(x[6], modulus, mredconstant)
		z[7] = IMForm(x[7], modulus, mredconstant)
	}
}

// ZeroVec sets all values of p1 to zero.
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func ZeroVec(p1 []uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
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

// MaskVec evaluates p2 = vec(p1>>w) & mask
// Iteration is done with respect to len(p1).
// All input must have a size which is a multiple of 8.
func MaskVec(p1 []uint64, w int, mask uint64, p2 []uint64) {

	N := len(p1)

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p1)%8 */
		x := (*[8]uint64)(unsafe.Pointer(&p1[j]))
		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(p2)%8 */
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
