package buffer

import (
	"encoding/binary"
	"fmt"
	"unsafe"
)

// ReadAsUint64 reads an uint64 from r and stores the result into c with pointer type casting into type T.
func ReadAsUint64[T any](r Reader, c *T) (n int64, err error) {
	/* #nosec G103 -- behavior and consequences well understood, pointer type cast */
	return ReadUint64(r, (*uint64)(unsafe.Pointer(c)))
}

// ReadAsUint32 reads an uint32 from r and stores the result into c with pointer type casting into type T.
func ReadAsUint32[T any](r Reader, c *T) (n int64, err error) {
	/* #nosec G103 -- behavior and consequences well understood, pointer type cast */
	return ReadUint32(r, (*uint32)(unsafe.Pointer(c)))
}

// ReadAsUint16 reads an uint16 from r and stores the result into c with pointer type casting into type T.
func ReadAsUint16[T any](r Reader, c *T) (n int64, err error) {
	/* #nosec G103 -- behavior and consequences well understood, pointer type cast */
	return ReadUint16(r, (*uint16)(unsafe.Pointer(c)))
}

// ReadAsUint8 reads an uint8 from r and stores the result into c with pointer type casting into type T.
func ReadAsUint8[T any](r Reader, c *T) (n int64, err error) {
	/* #nosec G103 -- behavior and consequences well understood, pointer type cast */
	return ReadUint8(r, (*uint8)(unsafe.Pointer(c)))
}

// ReadAsUint64Slice reads a slice of uint64 from r and stores the result into c with pointer type casting into type T.
func ReadAsUint64Slice[T any](r Reader, c []T) (n int64, err error) {
	/* #nosec G103 -- behavior and consequences well understood, pointer type cast */
	return ReadUint64Slice(r, *(*[]uint64)(unsafe.Pointer(&c)))
}

// ReadAsUint32Slice reads a slice of uint32 from r and stores the result into c with pointer type casting into type T.
func ReadAsUint32Slice[T any](r Reader, c []T) (n int64, err error) {
	/* #nosec G103 -- behavior and consequences well understood, pointer type cast */
	return ReadUint32Slice(r, *(*[]uint32)(unsafe.Pointer(&c)))
}

// ReadAsUint16Slice reads a slice of uint16 from r and stores the result into c with pointer type casting into type T.
func ReadAsUint16Slice[T any](r Reader, c []T) (n int64, err error) {
	/* #nosec G103 -- behavior and consequences well understood, pointer type cast */
	return ReadUint16Slice(r, *(*[]uint16)(unsafe.Pointer(&c)))
}

// ReadAsUint8Slice reads a slice of uint8 from r and stores the result into c with pointer type casting into type T.
func ReadAsUint8Slice[T any](r Reader, c []T) (n int64, err error) {
	/* #nosec G103 -- behavior and consequences well understood, pointer type cast */
	return ReadUint8Slice(r, *(*[]uint8)(unsafe.Pointer(&c)))
}

// Read reads a slice of bytes from r and copies it on c.
func Read(r Reader, c []byte) (n int64, err error) {
	slice, err := r.Peek(len(c))
	if err != nil {
		return int64(len(slice)), err
	}
	copy(c, slice)
	nint, err := r.Discard(len(c))
	return int64(nint), err
}

// ReadUint8 reads a byte from r and stores the result into *c.
func ReadUint8(r Reader, c *uint8) (n int64, err error) {

	if c == nil {
		return 0, fmt.Errorf("cannot ReadUint8: c is nil")
	}

	slice, err := r.Peek(1)
	if err != nil {
		return int64(len(slice)), err
	}

	// Reads one byte
	*c = uint8(slice[0])

	nint, err := r.Discard(1)

	return int64(nint), err
}

// ReadUint8Slice reads a slice of byte from r and stores the result into c.
func ReadUint8Slice(r Reader, c []uint8) (n int64, err error) {
	nint, err := r.Read(c)
	return int64(nint), err
}

// ReadUint16 reads a uint16 from r and stores the result into *c.
func ReadUint16(r Reader, c *uint16) (n int64, err error) {

	if c == nil {
		return 0, fmt.Errorf("cannot ReadUint16: c is nil")
	}

	slice, err := r.Peek(2)
	if err != nil {
		return int64(len(slice)), err
	}

	// Reads one byte
	*c = binary.LittleEndian.Uint16(slice)

	nint, err := r.Discard(2)

	return int64(nint), err
}

// ReadUint16Slice reads a slice of uint16 from r and stores the result into c.
func ReadUint16Slice(r Reader, c []uint16) (n int64, err error) {

	// c is empty, return
	if len(c) == 0 {
		return
	}

	var slice []byte

	size := r.Size()
	if len(c)<<1 < size {
		size = len(c) << 1
	}

	// Then returns the written bytes
	if slice, err = r.Peek(size); err != nil {
		return int64(len(slice)), err
	}

	buffered := len(slice) >> 1

	// If the slice to write on is equal or smaller than the amount peaked
	if N := len(c); N <= buffered {

		for i, j := 0, 0; i < N; i, j = i+1, j+2 {
			c[i] = binary.LittleEndian.Uint16(slice[j:])
		}

		nint, err := r.Discard(N << 1) // Discards what was read

		return int64(nint), err
	}

	// Decodes the maximum
	for i, j := 0, 0; i < buffered; i, j = i+1, j+2 {
		c[i] = binary.LittleEndian.Uint16(slice[j:])
	}

	// Discard what was peeked
	var inc int
	if inc, err = r.Discard(len(slice)); err != nil {
		return n + int64(inc), err
	}

	n += int64(inc)

	// Recurses on the remaining slice to fill
	var inc64 int64
	inc64, err = ReadUint16Slice(r, c[buffered:])

	return n + inc64, err
}

// ReadUint32 reads a uint32 from r and stores the result into *c.
func ReadUint32(r Reader, c *uint32) (n int64, err error) {

	if c == nil {
		return 0, fmt.Errorf("cannot ReadUint32: c is nil")
	}

	slice, err := r.Peek(4)
	if err != nil {
		return int64(len(slice)), err
	}

	// Reads one byte
	*c = binary.LittleEndian.Uint32(slice)

	nint, err := r.Discard(4)

	return int64(nint), err
}

// ReadUint32Slice reads a slice of uint32 from r and stores the result into c.
func ReadUint32Slice(r Reader, c []uint32) (n int64, err error) {

	// c is empty, return
	if len(c) == 0 {
		return
	}

	var slice []byte

	// Avoid EOF
	size := r.Size()
	if len(c)<<2 < size {
		size = len(c) << 2
	}

	// Then returns the written bytes
	if slice, err = r.Peek(size); err != nil {
		return int64(len(slice)), err
	}

	buffered := len(slice) >> 2

	// If the slice to write on is equal or smaller than the amount peaked
	if N := len(c); N <= buffered {

		for i, j := 0, 0; i < N; i, j = i+1, j+4 {
			c[i] = binary.LittleEndian.Uint32(slice[j:])
		}

		nint, err := r.Discard(N << 2) // Discards what was read

		return int64(nint), err
	}

	// Decodes the maximum
	for i, j := 0, 0; i < buffered; i, j = i+1, j+4 {
		c[i] = binary.LittleEndian.Uint32(slice[j:])
	}

	// Discard what was peeked
	var inc int
	if inc, err = r.Discard(len(slice)); err != nil {
		return n + int64(inc), err
	}

	n += int64(inc)

	// Recurses on the remaining slice to fill
	var inc64 int64
	inc64, err = ReadUint32Slice(r, c[buffered:])

	return n + inc64, err
}

// ReadUint64 reads a uint64 from r and stores the result into c.
func ReadUint64(r Reader, c *uint64) (n int64, err error) {

	if c == nil {
		return 0, fmt.Errorf("cannot ReadUint64: c is nil")
	}

	bytes, err := r.Peek(8)
	if err != nil {
		return int64(len(bytes)), err
	}

	// Reads one byte
	*c = binary.LittleEndian.Uint64(bytes)

	nint, err := r.Discard(8)

	return int64(nint), err
}

// ReadUint64Slice reads a slice of uint64 from r and stores the result into c.
func ReadUint64Slice(r Reader, c []uint64) (n int64, err error) {

	// c is empty, return
	if len(c) == 0 {
		return
	}

	var slice []byte

	// Avoid EOF
	size := r.Size()
	if len(c)<<3 < size {
		size = len(c) << 3
	}

	// Then returns the written bytes
	if slice, err = r.Peek(size); err != nil {
		return int64(len(slice)), err
	}

	buffered := len(slice) >> 3

	// If the slice to write on is equal or smaller than the amount peaked
	if N := len(c); N <= buffered {

		for i, j := 0, 0; i < N; i, j = i+1, j+8 {
			c[i] = binary.LittleEndian.Uint64(slice[j:])
		}

		nint, err := r.Discard(N << 3) // Discards what was read

		return int64(nint), err
	}

	// Decodes the maximum
	for i, j := 0, 0; i < buffered; i, j = i+1, j+8 {
		c[i] = binary.LittleEndian.Uint64(slice[j:])
	}

	// Discard what was peeked
	var inc int
	if inc, err = r.Discard(len(slice)); err != nil {
		return n + int64(inc), err
	}

	n += int64(inc)

	// Recurses on the remaining slice to fill
	var inc64 int64
	inc64, err = ReadUint64Slice(r, c[buffered:])

	return n + inc64, err
}
