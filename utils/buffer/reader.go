package buffer

import (
	"encoding/binary"
	"fmt"

	"github.com/tuneinsight/lattigo/v4/utils"
)

// ReadInt reads an int values from r and stores the result into *c.
func ReadInt(r Reader, c *int) (n int, err error) {

	if c == nil {
		return 0, fmt.Errorf("cannot ReadInt: c is nil")
	}

	return ReadUint64(r, utils.PointyIntToPointUint64(c))
}

// ReadUint8 reads a byte from r and stores the result into *c.
func ReadUint8(r Reader, c *uint8) (n int, err error) {

	if c == nil {
		return 0, fmt.Errorf("cannot ReadUint8: c is nil")
	}

	slice, err := r.Peek(1)
	if err != nil {
		return len(slice), err
	}

	// Reads one byte
	*c = uint8(slice[0])

	return r.Discard(1)
}

// ReadUint8Slice reads a slice of byte from r and stores the result into c.
func ReadUint8Slice(r Reader, c []uint8) (n int, err error) {
	return r.Read(c)
}

// ReadUint16 reads a uint16 from r and stores the result into *c.
func ReadUint16(r Reader, c *uint16) (n int, err error) {

	if c == nil {
		return 0, fmt.Errorf("cannot ReadUint16: c is nil")
	}

	slice, err := r.Peek(2)
	if err != nil {
		return len(slice), err
	}

	// Reads one byte
	*c = binary.LittleEndian.Uint16(slice)

	return r.Discard(2)
}

// ReadUint16Slice reads a slice of uint16 from r and stores the result into c.
func ReadUint16Slice(r Reader, c []uint16) (n int, err error) {

	// c is empty, return
	if len(c) == 0 {
		return
	}

	var slice []byte

	size := r.Size()
	if len(c)<<1 < size {
		size = len(c) << 1
	}

	// Then returns the writen bytes
	if slice, err = r.Peek(size); err != nil {
		return len(slice), err
	}

	buffered := len(slice) >> 1

	// If the slice to write on is equal or smaller than the amount peaked
	if N := len(c); N <= buffered {

		for i, j := 0, 0; i < N; i, j = i+1, j+2 {
			c[i] = binary.LittleEndian.Uint16(slice[j:])
		}

		return r.Discard(N << 1) // Discards what was read
	}

	// Decodes the maximum
	for i, j := 0, 0; i < buffered; i, j = i+1, j+2 {
		c[i] = binary.LittleEndian.Uint16(slice[j:])
	}

	// Discard what was peeked
	var inc int
	if inc, err = r.Discard(len(slice)); err != nil {
		return n + inc, err
	}

	n += inc

	// Recurses on the remaining slice to fill
	if inc, err = ReadUint16Slice(r, c[buffered:]); err != nil {
		return n + inc, err
	}

	return n + inc, nil
}

// ReadUint32 reads a uint32 from r and stores the result into *c.
func ReadUint32(r Reader, c *uint32) (n int, err error) {

	if c == nil {
		return 0, fmt.Errorf("cannot ReadUint32: c is nil")
	}

	slice, err := r.Peek(4)
	if err != nil {
		return len(slice), err
	}

	// Reads one byte
	*c = binary.LittleEndian.Uint32(slice)

	return r.Discard(4)
}

// ReadUint32Slice reads a slice of uint32 from r and stores the result into c.
func ReadUint32Slice(r Reader, c []uint32) (n int, err error) {

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

	// Then returns the writen bytes
	if slice, err = r.Peek(size); err != nil {
		return len(slice), err
	}

	buffered := len(slice) >> 2

	// If the slice to write on is equal or smaller than the amount peaked
	if N := len(c); N <= buffered {

		for i, j := 0, 0; i < N; i, j = i+1, j+4 {
			c[i] = binary.LittleEndian.Uint32(slice[j:])
		}

		return r.Discard(N << 2) // Discards what was read
	}

	// Decodes the maximum
	for i, j := 0, 0; i < buffered; i, j = i+1, j+4 {
		c[i] = binary.LittleEndian.Uint32(slice[j:])
	}

	// Discard what was peeked
	var inc int
	if inc, err = r.Discard(len(slice)); err != nil {
		return n + inc, err
	}

	n += inc

	// Recurses on the remaining slice to fill
	if inc, err = ReadUint32Slice(r, c[buffered:]); err != nil {
		return n + inc, err
	}

	return n + inc, nil
}

// ReadUint64 reads a uint64 from r and stores the result into c.
func ReadUint64(r Reader, c *uint64) (n int, err error) {

	if c == nil {
		return 0, fmt.Errorf("cannot ReadUint64: c is nil")
	}

	bytes, err := r.Peek(8)
	if err != nil {
		return len(bytes), err
	}

	// Reads one byte
	*c = binary.LittleEndian.Uint64(bytes)

	return r.Discard(8)
}

// ReadUint64Slice reads a slice of uint64 from r and stores the result into c.
func ReadUint64Slice(r Reader, c []uint64) (n int, err error) {

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

	// Then returns the writen bytes
	if slice, err = r.Peek(size); err != nil {
		return
	}

	buffered := len(slice) >> 3

	// If the slice to write on is equal or smaller than the amount peaked
	if N := len(c); N <= buffered {

		for i, j := 0, 0; i < N; i, j = i+1, j+8 {
			c[i] = binary.LittleEndian.Uint64(slice[j:])
		}

		return r.Discard(N << 3) // Discards what was read
	}

	// Decodes the maximum
	for i, j := 0, 0; i < buffered; i, j = i+1, j+8 {
		c[i] = binary.LittleEndian.Uint64(slice[j:])
	}

	// Discard what was peeked
	var inc int
	if inc, err = r.Discard(len(slice)); err != nil {
		return n + inc, err
	}

	n += inc

	// Recurses on the remaining slice to fill
	if inc, err = ReadUint64Slice(r, c[buffered:]); err != nil {
		return n + inc, err
	}

	return n + inc, nil
}
