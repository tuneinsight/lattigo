package utils

import (
	"encoding/binary"
	"io"
)

type Reader struct {
	io.Reader
	buff []byte
}

func NewReader(r io.Reader) *Reader {
	return &Reader{
		Reader: r,
		buff:   make([]byte, 1<<12),
	}
}

func (r *Reader) Read(p []byte) (n int, err error) {
	return r.Reader.Read(p)
}

func (r *Reader) ReadUint8(c *uint8) (n int, err error) {

	if n, err = r.Reader.Read(r.buff[:1]); err != nil {
		return
	}

	// Reads one byte
	*c = uint8(r.buff[0])

	return n, nil
}

func (r *Reader) ReadUint8Slice(c []uint8) (n int, err error) {

	buff := r.buff

	if n, err = r.Reader.Read(r.buff); err != nil {
		return
	}

	available := len(buff)

	// If the slice to write on is smaller than the available buffer
	if len(c) < available {

		// Copy the maximum on c
		copy(c, buff)

		return len(c), nil
	}

	// Copy the maximum on c
	copy(c, buff)

	// Updates the number of bytes read
	n += available

	// Recurses on the remaining slice to fill
	var inc int
	if inc, err = r.ReadUint8Slice(c[available:]); err != nil {
		return n + inc, err
	}

	return n + inc, nil
}

func (r *Reader) ReadUint16(c *uint16) (n int, err error) {

	if n, err = r.Reader.Read(r.buff[:2]); err != nil {
		return
	}

	// Reads one byte
	*c = binary.LittleEndian.Uint16(r.buff[:2])

	return n, nil
}

func (r *Reader) ReadUint16Slice(c []uint16) (n int, err error) {

	if len(c) == 0 {
		return
	}

	buff := r.buff

	if n, err = r.Reader.Read(r.buff); err != nil {
		return
	}

	available := len(buff) >> 1

	// If the slice to write on is smaller than the available buffer
	if N := len(c); N < available {

		for i, j := 0, 0; i < N; i, j = i+1, j+2 {
			c[i] = binary.LittleEndian.Uint16(buff[j:])
		}

		return n, nil
	}

	// Writes the maximum on c
	for i, j := 0, 0; i < available; i, j = i+1, j+2 {
		c[i] = binary.LittleEndian.Uint16(buff[j:])
	}

	// Recurses on the remaining slice to fill
	var inc int
	if inc, err = r.ReadUint16Slice(c[available:]); err != nil {
		return n + inc, err
	}

	return n + inc, nil
}

func (r *Reader) ReadUint32(c *uint32) (n int, err error) {

	if n, err = r.Reader.Read(r.buff[:4]); err != nil {
		return
	}

	*c = binary.LittleEndian.Uint32(r.buff[:4])

	return n, nil
}

func (r *Reader) ReadUint32Slice(c []uint32) (n int, err error) {

	if len(c) == 0 {
		return
	}

	buff := r.buff

	if n, err = r.Reader.Read(r.buff); err != nil {
		return
	}

	available := len(buff) >> 2

	// If the slice to write on is smaller than the available buffer
	if N := len(c); N < available {

		for i, j := 0, 0; i < N; i, j = i+1, j+4 {
			c[i] = binary.LittleEndian.Uint32(buff[j:])
		}

		return n, nil
	}

	// Writes the maximum on c
	for i, j := 0, 0; i < available; i, j = i+1, j+4 {
		c[i] = binary.LittleEndian.Uint32(buff[j:])
	}

	// Recurses on the remaining slice to fill
	var inc int
	if inc, err = r.ReadUint32Slice(c[available:]); err != nil {
		return n + inc, err
	}

	return n + inc, nil
}

func (r *Reader) ReadUint64(c *uint64) (n int, err error) {

	if n, err = r.Reader.Read(r.buff[:8]); err != nil {
		return
	}

	// Reads one byte
	*c = binary.LittleEndian.Uint64(r.buff[:8])

	return n, nil
}

func (r *Reader) ReadUint64Slice(c []uint64) (n int, err error) {

	if len(c) == 0 {
		return
	}

	buff := r.buff

	if n, err = r.Reader.Read(r.buff); err != nil {
		return
	}

	available := len(buff) >> 3

	// If the slice to write on is smaller than the available buffer
	if N := len(c); N < available {

		for i, j := 0, 0; i < N; i, j = i+1, j+8 {
			c[i] = binary.LittleEndian.Uint64(buff[j:])
		}

		return n, nil
	}

	// Writes the maximum on c
	for i, j := 0, 0; i < available; i, j = i+1, j+8 {
		c[i] = binary.LittleEndian.Uint64(buff[j:])
	}

	// Recurses on the remaining slice to fill
	var inc int
	if inc, err = r.ReadUint64Slice(c[available:]); err != nil {
		return n + inc, err
	}

	return n + inc, nil
}
