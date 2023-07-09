package buffer

import (
	"encoding/binary"
	"fmt"
)

// WriteInt writes an int c to w.
func WriteInt(w Writer, c int) (n int, err error) {
	return WriteUint64(w, uint64(c))
}

// WriteUint8 writes a byte c to w.
func WriteUint8(w Writer, c uint8) (n int, err error) {
	return w.Write([]byte{c})
}

// WriteUint8Slice writes a slice of bytes c to w.
func WriteUint8Slice(w Writer, c []uint8) (n int, err error) {
	return w.Write(c)
}

// WriteUint16 writes a uint16 c to w.
func WriteUint16(w Writer, c uint16) (n int, err error) {

	buf := w.AvailableBuffer()

	if w.Available()>>1 == 0 {
		if err = w.Flush(); err != nil {
			return
		}

		if w.Available()>>1 == 0 {
			return 0, fmt.Errorf("cannot WriteUint16: available buffer/2 is zero even after flush")
		}
	}

	binary.LittleEndian.PutUint16(buf[:2], c)

	return w.Write(buf[:2])
}

// WriteUint16Slice writes a slice of uint16 c to w.
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

		if available == 0 {
			return 0, fmt.Errorf("cannot WriteUint16Slice: available buffer/2 is zero even after flush")
		}
	}

	if N := len(c); N <= available { // If there is enough space in the available buffer
		buf = buf[:N<<1]
		for i := 0; i < N; i++ {
			binary.LittleEndian.PutUint16(buf[i<<2:(i<<2)+2], c[i])
		}

		return w.Write(buf)
	}

	// First fills the space
	for i := 0; i < available; i++ {
		buf = buf[:available<<1]
		binary.LittleEndian.PutUint16(buf[i<<1:(i<<1)+2], c[i])
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

// WriteUint32 writes a uint32 c into w.
func WriteUint32(w Writer, c uint32) (n int, err error) {

	buf := w.AvailableBuffer()

	if w.Available()>>2 == 0 {
		if err = w.Flush(); err != nil {
			return
		}

		if w.Available()>>2 == 0 {
			return 0, fmt.Errorf("cannot WriteUint32: available buffer/4 is zero even after flush")
		}
	}

	buf = buf[:4]
	binary.LittleEndian.PutUint32(buf, c)
	return w.Write(buf)
}

// WriteUint32Slice writes a slice of uint32 c into w.
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

		if available == 0 {
			return 0, fmt.Errorf("cannot WriteUint32Slice: available buffer/4 is zero even after flush")
		}
	}

	if N := len(c); N <= available { // If there is enough space in the available buffer
		buf = buf[:N<<2]
		for i := 0; i < N; i++ {
			binary.LittleEndian.PutUint32(buf[i<<2:(i<<2)+4], c[i])
		}
		return w.Write(buf)
	}

	// First fills the space
	buf = buf[:available<<2]
	for i := 0; i < available; i++ {
		binary.LittleEndian.PutUint32(buf[i<<2:(i<<2)+4], c[i])
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

// WriteUint64 writes a uint64 c into w.
func WriteUint64(w Writer, c uint64) (n int, err error) {

	buf := w.AvailableBuffer()

	if w.Available()>>3 == 0 {
		if err = w.Flush(); err != nil {
			return
		}

		if w.Available()>>3 == 0 {
			return 0, fmt.Errorf("cannot WriteUint64: available buffer/8 is zero even after flush")
		}
	}

	binary.LittleEndian.PutUint64(buf[:8], c)

	return w.Write(buf[:8])
}

// WriteUint64Slice writes a slice of uint64 into w.
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

		if available == 0 {
			return 0, fmt.Errorf("cannot WriteUint64Slice: available buffer/8 is zero even after flush")
		}
	}

	if N := len(c); N <= available { // If there is enough space in the available buffer
		buf = buf[:N<<3]
		for i := 0; i < N; i++ {
			binary.LittleEndian.PutUint64(buf[i<<3:(i<<3)+8], c[i])
		}
		return w.Write(buf)
	}

	// First fills the space
	buf = buf[:available<<3]
	for i := 0; i < available; i++ {
		binary.LittleEndian.PutUint64(buf[i<<3:(i<<3)+8], c[i])
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
