package buffer

import (
	"encoding/binary"
	"io"
)

// Writer defines a interface comprising of the minimum subset
// of methods defined by the type bufio.Writer necessary to run
// the functions defined in this file.
// See the documentation of bufio.Writer: https://pkg.go.dev/bufio.
type Writer interface {
	io.Writer
	Flush() (err error)
	AvailableBuffer() []byte
	Available() int
}

func WriteInt(w Writer, c int) (n int, err error) {
	return WriteUint64(w, uint64(c))
}

func WriteUint8(w Writer, c uint8) (n int, err error) {
	return w.Write([]byte{c})
}

func WriteUint8Slice(w Writer, c []uint8) (n int, err error) {
	return w.Write(c)
}

func WriteUint16(w Writer, c uint16) (n int, err error) {

	buf := w.AvailableBuffer()

	if w.Available()>>1 == 0 {
		if err = w.Flush(); err != nil {
			return
		}
	}

	var bb = [2]byte{}
	binary.LittleEndian.PutUint16(bb[:], c)
	buf = append(buf, bb[:]...)

	return w.Write(buf)
}

func WriteUint16Slice(w Writer, c []uint16) (n int, err error) {

	if len(c) == 0 {
		return
	}

	buf := w.AvailableBuffer()

	// Remaining available space in the internal buffer
	available := w.Available() >> 1

	if available == 0 {
		if err = w.Flush(); err != nil {
			return
		}

		available = w.Available() >> 1
	}

	var bb = [2]byte{}

	if N := len(c); N <= available { // If there is enough space in the available buffer

		for i := 0; i < N; i++ {
			binary.LittleEndian.PutUint16(bb[:], c[i])
			buf = append(buf, bb[:]...)
		}

		return w.Write(buf)
	}

	// First fills the space
	for i := 0; i < available; i++ {
		binary.LittleEndian.PutUint16(bb[:], c[i])
		buf = append(buf, bb[:]...)
	}

	var inc int
	if inc, err = w.Write(buf); err != nil {
		return n + inc, err
	}

	n += inc

	// Flushes
	if err = w.Flush(); err != nil {
		return n, err
	}

	// Then recurses on itself with the remaining slice
	if inc, err = WriteUint16Slice(w, c[available:]); err != nil {
		return n + inc, err
	}

	return n + inc, nil
}

func WriteUint32(w Writer, c uint32) (n int, err error) {

	buf := w.AvailableBuffer()

	if w.Available()>>2 == 0 {
		if err = w.Flush(); err != nil {
			return
		}
	}

	var bb = [4]byte{}
	binary.LittleEndian.PutUint32(bb[:], c)
	buf = append(buf, bb[:]...)

	return w.Write(buf)
}

func WriteUint32Slice(w Writer, c []uint32) (n int, err error) {

	if len(c) == 0 {
		return
	}

	buf := w.AvailableBuffer()

	// Remaining available space in the internal buffer
	available := w.Available() >> 2

	if available == 0 {
		if err = w.Flush(); err != nil {
			return
		}

		available = w.Available() >> 2
	}

	var bb = [4]byte{}

	if N := len(c); N <= available { // If there is enough space in the available buffer

		for i := 0; i < N; i++ {
			binary.LittleEndian.PutUint32(bb[:], c[i])
			buf = append(buf, bb[:]...)
		}

		return w.Write(buf)
	}

	// First fills the space
	for i := 0; i < available; i++ {
		binary.LittleEndian.PutUint32(bb[:], c[i])
		buf = append(buf, bb[:]...)
	}

	var inc int
	if inc, err = w.Write(buf); err != nil {
		return n + inc, err
	}

	n += inc

	// Flushes
	if err = w.Flush(); err != nil {
		return n, err
	}

	// Then recurses on itself with the remaining slice
	if inc, err = WriteUint32Slice(w, c[available:]); err != nil {
		return n + inc, err
	}

	return n + inc, nil
}

func WriteUint64(w Writer, c uint64) (n int, err error) {

	buf := w.AvailableBuffer()

	if w.Available()>>3 == 0 {
		if err = w.Flush(); err != nil {
			return
		}
	}

	var bb = [8]byte{}
	binary.LittleEndian.PutUint64(bb[:], c)
	buf = append(buf, bb[:]...)

	return w.Write(buf)
}

func WriteUint64Slice(w Writer, c []uint64) (n int, err error) {

	if len(c) == 0 {
		return
	}

	buf := w.AvailableBuffer()

	// Remaining available space in the internal buffer
	available := w.Available() >> 3

	if available == 0 {
		if err = w.Flush(); err != nil {
			return
		}

		available = w.Available() >> 3
	}

	var bb = [8]byte{}

	if N := len(c); N <= available { // If there is enough space in the available buffer

		for i := 0; i < N; i++ {
			binary.LittleEndian.PutUint64(bb[:], c[i])
			buf = append(buf, bb[:]...)
		}

		return w.Write(buf)
	}

	// First fills the space
	for i := 0; i < available; i++ {
		binary.LittleEndian.PutUint64(bb[:], c[i])
		buf = append(buf, bb[:]...)
	}

	var inc int
	if inc, err = w.Write(buf); err != nil {
		return n + inc, err
	}

	n += inc

	// Flushes
	if err = w.Flush(); err != nil {
		return n, err
	}

	// Then recurses on itself with the remaining slice
	if inc, err = WriteUint64Slice(w, c[available:]); err != nil {
		return n + inc, err
	}

	return n + inc, nil
}
