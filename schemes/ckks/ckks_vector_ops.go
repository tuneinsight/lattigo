package ckks

import (
	"fmt"
	"math/big"
	"math/bits"
	"unsafe"

	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

const (
	minVecLenForLoopUnrolling = 16
)

// SpecialIFFTDouble performs the CKKS special inverse FFT transform in place.
func SpecialIFFTDouble(values []complex128, N, M int, rotGroup []int, roots []complex128) {

	// Sanity check
	if len(values) < N || len(rotGroup) < N || len(roots) < M+1 {
		panic(fmt.Sprintf("invalid call of SpecialIFFTDouble: len(values)=%d or len(rotGroup)=%d < N=%d or len(roots)=%d < M+1=%d", len(values), len(rotGroup), N, len(roots), M))
	}

	logN := int(bits.Len64(uint64(N))) - 1
	logM := int(bits.Len64(uint64(M))) - 1
	for loglen := logN; loglen > 0; loglen-- {
		len := 1 << loglen
		lenh := len >> 1
		lenq := len << 2
		logGap := logM - 2 - loglen
		mask := lenq - 1
		for i := 0; i < N; i += len {
			for j, k := 0, i; j < lenh; j, k = j+1, k+1 {
				values[k], values[k+lenh] = values[k]+values[k+lenh], (values[k]-values[k+lenh])*roots[(lenq-(rotGroup[j]&mask))<<logGap]
			}
		}
	}

	for i := 0; i < N; i++ {
		values[i] /= complex(float64(N), 0)
	}

	utils.BitReverseInPlaceSlice(values, N)
}

// SpecialFFTDouble performs the CKKS special FFT transform in place.
func SpecialFFTDouble(values []complex128, N, M int, rotGroup []int, roots []complex128) {

	// Sanity check
	if len(values) < N || len(rotGroup) < N || len(roots) < M+1 {
		panic(fmt.Sprintf("invalid call of SpecialFFTDouble: len(values)=%d or len(rotGroup)=%d < N=%d or len(roots)=%d < M+1=%d", len(values), len(rotGroup), N, len(roots), M))
	}

	utils.BitReverseInPlaceSlice(values, N)
	logN := int(bits.Len64(uint64(N))) - 1
	logM := int(bits.Len64(uint64(M))) - 1
	for loglen := 1; loglen <= logN; loglen++ {
		len := 1 << loglen
		lenh := len >> 1
		lenq := len << 2
		logGap := logM - 2 - loglen
		mask := lenq - 1
		for i := 0; i < N; i += len {
			for j, k := 0, i; j < lenh; j, k = j+1, k+1 {
				values[k+lenh] *= roots[(rotGroup[j]&mask)<<logGap]
				values[k], values[k+lenh] = values[k]+values[k+lenh], values[k]-values[k+lenh]
			}
		}
	}
}

// SpecialFFTArbitrary evaluates the decoding matrix on a slice of ring.Complex values.
func SpecialFFTArbitrary(values []*bignum.Complex, N, M int, rotGroup []int, roots []*bignum.Complex) {

	// Sanity check
	if len(values) < N || len(rotGroup) < N || len(roots) < M+1 {
		panic(fmt.Sprintf("invalid call of SpecialFFTArbitrary: len(values)=%d or len(rotGroup)=%d < N=%d or len(roots)=%d < M+1=%d", len(values), len(rotGroup), N, len(roots), M))
	}

	u := &bignum.Complex{new(big.Float), new(big.Float)}
	v := &bignum.Complex{new(big.Float), new(big.Float)}

	utils.BitReverseInPlaceSlice(values, N)

	cMul := bignum.NewComplexMultiplier()

	logN := int(bits.Len64(uint64(N))) - 1
	logM := int(bits.Len64(uint64(M))) - 1
	for loglen := 1; loglen <= logN; loglen++ {
		len := 1 << loglen
		lenh := len >> 1
		lenq := len << 2
		logGap := logM - 2 - loglen
		mask := lenq - 1
		for i := 0; i < N; i += len {
			for j, k := 0, i; j < lenh; j, k = j+1, k+1 {
				u.Set(values[i+j])
				v.Set(values[i+j+lenh])
				cMul.Mul(v, roots[(rotGroup[j]&mask)<<logGap], v)
				values[i+j].Add(u, v)
				values[i+j+lenh].Sub(u, v)
			}
		}
	}
}

// SpecialIFFTArbitrary evaluates the encoding matrix on a slice of ring.Complex values.
func SpecialIFFTArbitrary(values []*bignum.Complex, N, M int, rotGroup []int, roots []*bignum.Complex) {

	// Sanity check
	if len(values) < N || len(rotGroup) < N || len(roots) < M+1 {
		panic(fmt.Sprintf("invalid call of SpecialIFFTArbitrary: len(values)=%d or len(rotGroup)=%d < N=%d or len(roots)=%d < M+1=%d", len(values), len(rotGroup), N, len(roots), M))
	}

	u := &bignum.Complex{new(big.Float), new(big.Float)}
	v := &bignum.Complex{new(big.Float), new(big.Float)}

	cMul := bignum.NewComplexMultiplier()

	logN := int(bits.Len64(uint64(N))) - 1
	logM := int(bits.Len64(uint64(M))) - 1
	for loglen := logN; loglen > 0; loglen-- {
		len := 1 << loglen
		lenh := len >> 1
		lenq := len << 2
		logGap := logM - 2 - loglen
		mask := lenq - 1
		for i := 0; i < N; i += len {
			for j, k := 0, i; j < lenh; j, k = j+1, k+1 {
				u.Add(values[i+j], values[i+j+lenh])
				v.Sub(values[i+j], values[i+j+lenh])
				cMul.Mul(v, roots[(lenq-(rotGroup[j]&mask))<<logGap], v)
				values[i+j].Set(u)
				values[i+j+lenh].Set(v)

			}
		}
	}

	NBig := new(big.Float).SetInt64(int64(N))
	for i := range values {
		values[i][0].Quo(values[i][0], NBig)
		values[i][1].Quo(values[i][1], NBig)
	}

	utils.BitReverseInPlaceSlice(values, N)
}

// SpecialFFTDoubleUL8 performs the CKKS special FFT transform in place with unrolled loops of size 8.
func SpecialFFTDoubleUL8(values []complex128, N, M int, rotGroup []int, roots []complex128) {

	// Sanity check
	if len(values) < minVecLenForLoopUnrolling {
		panic(fmt.Sprintf("unsafe call of SpecialFFTDoubleUL8: len(values)=%d < %d", len(values), minVecLenForLoopUnrolling))
	}

	// Sanity check
	if len(values) < N || len(rotGroup) < N || len(roots) < M+1 {
		panic(fmt.Sprintf("invalid call of SpecialFFTDoubleUL8: len(values)=%d or len(rotGroup)=%d < N=%d or len(roots)=%d < M+1=%d", len(values), len(rotGroup), N, len(roots), M))
	}

	utils.BitReverseInPlaceSlice(values, N)

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

					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(values)%8 != 0 */
					u := (*[8]complex128)(unsafe.Pointer(&values[k]))
					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(values)%8 != 0 */
					v := (*[8]complex128)(unsafe.Pointer(&values[k+lenh]))
					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(values)%8 != 0 */
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

				/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(values)%8 != 0 */
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

				/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(values)%8 != 0 */
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

				/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(values)%8 != 0 */
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

				/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(values)%8 != 0 */
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

// SpecialiFFTDoubleUnrolled8 performs the CKKS special inverse FFT transform in place with unrolled loops of size 8.
func SpecialiFFTDoubleUnrolled8(values []complex128, N, M int, rotGroup []int, roots []complex128) {

	// Sanity check
	if len(values) < minVecLenForLoopUnrolling {
		panic(fmt.Sprintf("unsafe call of SpecialiFFTDoubleUnrolled8: len(values)=%d < %d", len(values), minVecLenForLoopUnrolling))
	}

	// Sanity check
	if len(values) < N || len(rotGroup) < N || len(roots) < M+1 {
		panic(fmt.Sprintf("invalid call of SpecialiFFTDoubleUnrolled8: len(values)=%d or len(rotGroup)=%d < N=%d or len(roots)=%d < M+1=%d", len(values), len(rotGroup), N, len(roots), M))
	}

	logN := int(bits.Len64(uint64(N))) - 1
	logM := int(bits.Len64(uint64(M))) - 1

	for loglen := logN; loglen > 0; loglen-- {

		len := 1 << loglen
		lenh := len >> 1
		lenq := len << 2
		logGap := logM - 2 - loglen
		mask := lenq - 1
		if lenh > 8 {
			for i := 0; i < N; i += len {
				for j, k := 0, i; j < lenh; j, k = j+8, k+8 {

					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(values)%8 != 0 */
					u := (*[8]complex128)(unsafe.Pointer(&values[k]))
					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(values)%8 != 0 */
					v := (*[8]complex128)(unsafe.Pointer(&values[k+lenh]))
					/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(rotGroup)%8 != 0 */
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

				/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(rotGroup)%8 != 0  */
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

				/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(rotGroup)%8 != 0  */
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

				/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(rotGroup)%8 != 0  */
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

				/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(rotGroup)%8 != 0  */
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

	divideComplex128SliceUnrolled8(values, complex(float64(N), 0))

	utils.BitReverseInPlaceSlice(values, N)
}

// divideComplex128SliceUnrolled8 divides the entries in values by scaleVal in place.
func divideComplex128SliceUnrolled8(values []complex128, scaleVal complex128) {
	lenValues := len(values)
	for i := 0; i < lenValues; i = i + 8 {

		/* #nosec G103 -- behavior and consequences well understood */
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
