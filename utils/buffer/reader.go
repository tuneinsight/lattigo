package buffer

import (
	"encoding/binary"
	"fmt"
	"io"

	"github.com/tuneinsight/lattigo/v4/utils"
)

// Reader defines a interface comprising of the minimum subset
// of methods defined by the type bufio.Reader necessary to run
// the functions defined in this file.
// See the documentation of bufio.Reader: https://pkg.go.dev/bufio.
type Reader interface {
	io.Reader
	Size() int
	Peek(n int) ([]byte, error)
	Discard(n int) (discarded int, err error)
}

func ReadInt(r Reader, c *int) (n int, err error) {

	if c == nil {
		return 0, fmt.Errorf("cannot ReadInt: c is nil")
	}

	return ReadUint64(r, utils.PointyIntToPointUint64(c))
}

func ReadUint8(r Reader, c *uint8) (n int, err error) {

	if c == nil {
		return 0, fmt.Errorf("cannot ReadUint8: c is nil")
	}

	var bb = [1]byte{}

	if n, err = r.Read(bb[:]); err != nil {
		return
	}

	// Reads one byte
	*c = uint8(bb[0])

	return n, nil
}

func ReadUint8Slice(r Reader, c []uint8) (n int, err error) {
	return r.Read(c)
}

func ReadUint16(r Reader, c *uint16) (n int, err error) {

	if c == nil {
		return 0, fmt.Errorf("cannot ReadUint16: c is nil")
	}

	var bb = [2]byte{}

	if n, err = r.Read(bb[:]); err != nil {
		return
	}

	// Reads one byte
	*c = binary.LittleEndian.Uint16(bb[:])

	return n, nil
}

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
		fmt.Println(err)
		return
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

func ReadUint32(r Reader, c *uint32) (n int, err error) {

	if c == nil {
		return 0, fmt.Errorf("cannot ReadUint32: c is nil")
	}

	var bb = [4]byte{}

	if n, err = r.Read(bb[:]); err != nil {
		return
	}

	// Reads one byte
	*c = binary.LittleEndian.Uint32(bb[:])

	return n, nil
}

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
		fmt.Println(err)
		return
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

func ReadUint64(r Reader, c *uint64) (n int, err error) {

	if c == nil {
		return 0, fmt.Errorf("cannot ReadUint64: c is nil")
	}

	var bb = [8]byte{}

	if n, err = r.Read(bb[:]); err != nil {
		return
	}

	// Reads one byte
	*c = binary.LittleEndian.Uint64(bb[:])

	return n, nil
}

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
