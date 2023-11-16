// Package buffer implement methods for efficiently writing and reading values
// to and from io.Writer and io.Reader that also expose their internal buffers.
package buffer

import (
	"fmt"
	"io"
)

// Writer is an interface for writers that expose their internal
// buffers.
// This interface is notably implemented by the bufio.Writer type
// (see https://pkg.go.dev/bufio#Writer) and by the Buffer type.
type Writer interface {
	io.Writer
	Flush() (err error)
	AvailableBuffer() []byte
	Available() int
}

// Reader is an interface for readers that expose their internal
// buffers.
// This interface is notably implemented by the bufio.Reader type
// (see https://pkg.go.dev/bufio#Reader) and by the Buffer type.
type Reader interface {
	io.Reader
	Size() int
	Peek(n int) ([]byte, error)
	Discard(n int) (discarded int, err error)
}

// Buffer is a simple []byte-based buffer that complies to the
// Writer and Reader interfaces. This type assumes that its
// backing slice has a fixed size and won't attempt to extend
// it. Instead, writes beyond capacity will result in an error.
type Buffer struct {
	buf []byte
	n   int
	off int
}

// NewBuffer creates a new Buffer struct with buff as a backing
// []byte. The read and write offset are initialized at buff[0].
// Hence, writing new data will overwrite the content of buff.
func NewBuffer(buff []byte) *Buffer {
	b := new(Buffer)
	b.buf = buff
	return b
}

// NewBufferSize creates a new Buffer with size capacity.
func NewBufferSize(size int) *Buffer {
	b := new(Buffer)
	b.buf = make([]byte, size)
	return b
}

// Write writes p into b. It returns the number of bytes written
// and an error if attempting to write passed the initial capacity
// of the buffer. Note that the case where p shares the same backing
// memory as b is optimized.
func (b *Buffer) Write(p []byte) (n int, err error) {
	if len(p)+b.n > cap(b.buf) {
		return 0, fmt.Errorf("buffer too small")
	}
	inc := copy(b.buf[b.n:], p) // This is optimized if &b.buf[b.n:][0] == &p[0]
	b.n += inc
	return inc, nil
}

// Flush doesn't do anything on this slice-based buffer.
func (b *Buffer) Flush() (err error) {
	return nil
}

// AvailableBuffer returns an empty buffer with b.Available() capacity, to be
// directly appended to and passed to a Write call. The buffer is only valid
// until the next write operation on b.
func (b *Buffer) AvailableBuffer() []byte {
	return b.buf[b.n:][:0]
}

// Available returns the number of bytes available for writes on the buffer.
func (b *Buffer) Available() int {
	return len(b.buf) - b.n
}

// Bytes returns the backing slice.
func (b *Buffer) Bytes() []byte {
	return b.buf
}

// Reset re-initializes the read and write offsets of b.
func (b *Buffer) Reset() {
	b.n = 0
	b.off = 0
}

// Read reads len(p) bytes from the read offset of b into p. It returns the
// number n of bytes read and an error if n < len(p).
func (b *Buffer) Read(p []byte) (n int, err error) {
	n = copy(p, b.buf[b.off:])
	b.off += n
	if n < len(p) {
		return n, io.EOF
	}
	return n, nil
}

// Size returns the size of the buffer available for read.
func (b *Buffer) Size() int {
	return len(b.buf) - b.off
}

// Peek returns the next n bytes without advancing the read offset, directly
// as a reslice of the internal buffer. It returns an error if the number of
// returned bytes is smaller than n.
func (b *Buffer) Peek(n int) ([]byte, error) {
	if b.off+n > len(b.buf) {
		return b.buf[b.off:], io.EOF
	}
	return b.buf[b.off : b.off+n], nil
}

// Discard skips the next n bytes, returning the number of bytes discarded. If
// Discard skips fewer than n bytes, it also returns an error.
func (b *Buffer) Discard(n int) (discarded int, err error) {
	remain := len(b.buf) - b.off
	if n > remain {
		b.off = len(b.buf)
		return remain, io.EOF
	}
	b.off += n
	return n, nil
}
