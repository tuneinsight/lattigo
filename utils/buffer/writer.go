package buffer

import (
	"encoding/binary"
	"fmt"
	"io"
)

const (
	DefaultWriterBufferSize = 1024
)

// Writer implements buffering for an io.Writer object.
// If an error occurs writing to a Writer, no more data will be accepted and all subsequent writes, and Flush, will return the error.
// After all data has been written, the client should call the Flush method to guarantee all data has been forwarded to the underlying io.Writer.
type Writer struct {
	io.Writer
	buff []byte
	n    int
	err  error
}

// NewWriter returns a new Writer whose buffer has the default size DefaultWriterBufferSize.
// If the argument io.Writer is already a Writer with large enough buffer size, it returns the underlying Writer.
func NewWriter(w io.Writer) *Writer {
	return NewWriterSize(w, DefaultWriterBufferSize)
}

// NewWriterSize returns a new Writer whose buffer has the specified size.
func NewWriterSize(w io.Writer, size int) *Writer {

	switch w := w.(type) {
	case *Writer:
		if w.Size() >= size {
			return w
		}
	}

	return &Writer{
		Writer: w,
		buff:   make([]byte, size),
		n:      0,
	}
}

// Available returns how many bytes are unused in the buffer.
func (w *Writer) Available() int {
	return len(w.buff[w.n:])
}

// AvailableBuffer returns an empty buffer with b.Available() capacity.
// This buffer is intended to be appended to and passed to an immediately succeeding Write call.
// The buffer is only valid until the next write operation on b.
func (w *Writer) AvailableBuffer() []byte {
	return make([]byte, w.Available())
}

// Size returns the size of the underlying buffer in bytes.
func (w *Writer) Size() int {
	return len(w.buff)
}

// Buffered returns the number of bytes that have been written into the current buffer.
func (w *Writer) Buffered() int {
	return w.n
}

// Flush writes any buffered data to the underlying io.Writer.
func (w *Writer) Flush() (err error) {

	if w.err != nil {
		return fmt.Errorf("cannot flush: previous error: %w", w.err)
	}

	if _, err = w.Writer.Write(w.buff[:w.n]); err != nil {
		w.err = err
		return fmt.Errorf("cannot flush: %w", err)
	}

	w.n = 0

	return
}

// Reset discards any unflushed buffered data, clears any error, and resets b to write its output to w.
// Calling Reset on the zero value of Writer initializes the internal buffer to the default size.
func (w *Writer) Reset() {
	w.err = nil
	buff := w.buff
	for i := range buff {
		buff[i] = 0
	}
	w.n = 0
}

// Write flushes the internal buffer on the io.Writer and writes p directly on the underlying io.Writer.
// It returns the number of bytes written.
func (w *Writer) Write(p []byte) (n int, err error) {

	if w.err != nil {
		return n, fmt.Errorf("cannot Write: previous error: %w", w.err)
	}

	// First we flush because we bypass the internal buffer
	if err = w.Flush(); err != nil {
		w.err = err
		return
	}

	return w.Writer.Write(p)
}

// WriteUint8 writes a single uint8.
func (w *Writer) WriteUint8(c uint8) (n int, err error) {

	if w.err != nil {
		return n, fmt.Errorf("cannot WriteUint8: previous error: %w", w.err)
	}

	if len(w.buff[w.n:]) < 1 {
		if err = w.Flush(); err != nil {
			w.err = err
			return n, fmt.Errorf("cannot WriteUint8: %w", err)
		}
	}

	w.buff[w.n] = c

	w.n++

	return 1, nil
}

// WriteUint8Slice writes a slice of uint8.
func (w *Writer) WriteUint8Slice(c []uint8) (n int, err error) {

	if w.err != nil {
		return n, fmt.Errorf("cannot WriteUint8Slice: previous error: %w", w.err)
	}

	return w.Write(c)
}

// WriteUint16 writes a single uint16.
func (w *Writer) WriteUint16(c uint16) (n int, err error) {

	if w.err != nil {
		return n, fmt.Errorf("cannot WriteUint16: previous error: %w", w.err)
	}

	if len(w.buff[w.n:]) < 2 {
		if err = w.Flush(); err != nil {
			w.err = err
			return n, fmt.Errorf("cannot WriteUint16: %w", err)
		}
	}

	binary.LittleEndian.PutUint16(w.buff[w.n:], c)

	w.n += 2

	return 2, nil
}

// WriteUint16Slice writes a slice of uint16.
func (w *Writer) WriteUint16Slice(c []uint16) (n int, err error) {

	if w.err != nil {
		return n, fmt.Errorf("cannot WriteUint16Slice: previous error: %w", w.err)
	}

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
		w.err = err
		return n + inc, err
	}

	return n + inc, nil
}

// WriteUint32 writes a single uint32.
func (w *Writer) WriteUint32(c uint32) (n int, err error) {

	if w.err != nil {
		return n, fmt.Errorf("cannot WriteUint32: previous error: %w", w.err)
	}

	if len(w.buff[w.n:]) < 4 {
		if err = w.Flush(); err != nil {
			w.err = err
			return n, fmt.Errorf("cannot WriteUint32: %w", err)
		}
	}

	binary.LittleEndian.PutUint32(w.buff[w.n:], c)

	w.n += 4

	return 4, nil
}

// WriteUint32Slice writes a slice of uint32.
func (w *Writer) WriteUint32Slice(c []uint32) (n int, err error) {

	if w.err != nil {
		return n, fmt.Errorf("cannot WriteUint32Slice: previous error: %w", w.err)
	}

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
		w.err = err
		return n + inc, err
	}

	return n + inc, nil
}

// WriteUint64 writes a single uint64.
func (w *Writer) WriteUint64(c uint64) (n int, err error) {

	if w.err != nil {
		return n, fmt.Errorf("cannot WriteUint64: previous error: %w", w.err)
	}

	if len(w.buff[w.n:]) < 8 {
		if err = w.Flush(); err != nil {
			w.err = err
			return n, fmt.Errorf("cannot WriteUint64: %w", err)
		}
	}

	binary.LittleEndian.PutUint64(w.buff[w.n:], c)

	w.n += 8

	return 8, nil
}

// WriteUint64Slice writes a slice of uint64.
func (w *Writer) WriteUint64Slice(c []uint64) (n int, err error) {

	if w.err != nil {
		return n, fmt.Errorf("cannot WriteUint64Slice: previous error: %w", w.err)
	}

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
		w.err = err
		return n + inc, err
	}

	return n + inc, nil
}
