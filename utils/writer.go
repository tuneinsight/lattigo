package utils

import (
	"encoding/binary"
	"io"
)

type Writer struct {
	io.Writer
	buff []byte
	n    int
}

func NewWriter(w io.Writer) *Writer {
	return &Writer{
		Writer: w,
		buff:   make([]byte, 1<<10), //1KB of buffer
		n:      0,
	}
}

func (w *Writer) Flush() (err error) {
	if _, err = w.Writer.Write(w.buff[:w.n]); err != nil {
		return
	}

	w.n = 0

	return
}

func (w *Writer) Write(p []byte) (n int, err error) {
	return w.Writer.Write(p)
}

func (w *Writer) WriteUint8(c uint8) (n int, err error) {

	if len(w.buff[w.n:]) < 1 {
		if err = w.Flush(); err != nil {
			return
		}
	}

	w.buff[w.n] = c

	w.n++

	return 1, nil
}

func (w *Writer) WriteUint8Slice(c []uint8) (n int, err error) {
	return w.Write(c)
}

func (w *Writer) WriteUint16(c uint16) (n int, err error) {

	if len(w.buff[w.n:]) < 2 {
		if err = w.Flush(); err != nil {
			return
		}
	}

	binary.LittleEndian.PutUint16(w.buff[w.n:], c)

	w.n += 2

	return 2, nil
}

func (w *Writer) WriteUint16Slice(c []uint16) (n int, err error) {

	buff := w.buff[w.n:]

	// Remaining available space in the internal buffer
	available := len(buff) >> 1

	if len(c) < available { // If there is enough space in the available buffer

		N := len(c)

		for i, j := 0, 0; i < N; i, j = i+1, j+2 {
			binary.LittleEndian.PutUint16(buff[j:], c[i])
		}

		w.n += N << 1

		return N << 1, nil
	}

	// First fills the space
	for i, j := 0, 0; i < available; i, j = i+1, j+2 {
		binary.LittleEndian.PutUint16(buff[j:], c[i])
	}

	w.n += available << 1 // Updates pointer

	n += available << 1 // Updates number of bytes written

	// Flushes
	if err = w.Flush(); err != nil {
		return n, err
	}

	// Then recurses on itself with the remaining slice
	var inc int
	if inc, err = w.WriteUint16Slice(c[available:]); err != nil {
		return n + inc, err
	}

	return n + inc, nil
}

func (w *Writer) WriteUint32(c uint32) (n int, err error) {

	if len(w.buff[w.n:]) < 4 {
		if err = w.Flush(); err != nil {
			return
		}
	}

	binary.LittleEndian.PutUint32(w.buff[w.n:], c)

	w.n += 4

	return 4, nil
}

func (w *Writer) WriteUint32Slice(c []uint32) (n int, err error) {

	buff := w.buff[w.n:]

	// Remaining available space in the internal buffer
	available := len(buff) >> 2

	if len(c) < available { // If there is enough space in the available buffer

		N := len(c)

		for i, j := 0, 0; i < N; i, j = i+1, j+4 {
			binary.LittleEndian.PutUint32(buff[j:], c[i])
		}

		w.n += N << 2

		return N << 2, nil
	}

	// First fills the space
	for i, j := 0, 0; i < available; i, j = i+1, j+4 {
		binary.LittleEndian.PutUint32(buff[j:], c[i])
	}

	w.n += available << 2 // Updates pointer

	n += available << 2 // Updates number of bytes written

	// Flushes
	if err = w.Flush(); err != nil {
		return n, err
	}

	// Then recurses on itself with the remaining slice
	var inc int
	if inc, err = w.WriteUint32Slice(c[available:]); err != nil {
		return n + inc, err
	}

	return n + inc, nil
}

func (w *Writer) WriteUint64(c uint64) (n int, err error) {

	if len(w.buff[w.n:]) < 8 {
		if err = w.Flush(); err != nil {
			return
		}
	}

	binary.LittleEndian.PutUint64(w.buff[w.n:], c)

	w.n += 8

	return 8, nil
}

func (w *Writer) WriteUint64Slice(c []uint64) (n int, err error) {

	buff := w.buff[w.n:]

	// Remaining available space in the internal buffer
	available := len(buff) >> 3

	if len(c) < available { // If there is enough space in the available buffer

		N := len(c)

		for i, j := 0, 0; i < N; i, j = i+1, j+8 {
			binary.LittleEndian.PutUint64(buff[j:], c[i])
		}

		w.n += N << 3

		return N << 3, nil
	}

	// First fills the space
	for i, j := 0, 0; i < available; i, j = i+1, j+8 {
		binary.LittleEndian.PutUint64(buff[j:], c[i])
	}

	w.n += available << 3 // Updates pointer

	n += available << 3 // Updates number of bytes written

	// Flushes
	if err = w.Flush(); err != nil {
		return n, err
	}

	// Then recurses on itself with the remaining slice
	var inc int
	if inc, err = w.WriteUint64Slice(c[available:]); err != nil {
		return n + inc, err
	}

	return n + inc, nil
}
