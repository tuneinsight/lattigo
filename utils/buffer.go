// Package containing helper structures and function
package utils

// Buffer is a simple wrapper around a []byte to facilitate efficient marshaling of lattigo's objects
type Buffer struct {
	buf []byte
}

// NewBuffer creates a new buffer from the provided backing slice
func NewBuffer(s []byte) *Buffer {
	return &Buffer{s}
}

func (b *Buffer) WriteUint8(c byte) {
	b.buf = append(b.buf, c)
}

func (b *Buffer) WriteUint64(v uint64) {
	b.buf = append(b.buf, byte(v>>56),
		byte(v>>48),
		byte(v>>40),
		byte(v>>32),
		byte(v>>24),
		byte(v>>16),
		byte(v>>8),
		byte(v))
}

func (b *Buffer) WriteUint64Slice(s []uint64) {
	for _, v := range s {
		b.WriteUint64(v)
	}
}

func (b *Buffer) WriteUint8Slice(s []uint8) {
	for _, v := range s {
		b.WriteUint8(v)
	}
}

func (b *Buffer) ReadUint8() byte {
	v := b.buf[0]
	b.buf = b.buf[1:]
	return v
}

func (b *Buffer) ReadUint64() uint64 {
	v := b.buf[:8]
	b.buf = b.buf[8:]
	return uint64(v[7]) | uint64(v[6])<<8 | uint64(v[5])<<16 | uint64(v[4])<<24 |
		uint64(v[3])<<32 | uint64(v[2])<<40 | uint64(v[1])<<48 | uint64(v[0])<<56
}

func (b *Buffer) ReadUint64Slice(rec []uint64) {
	for i := range rec {
		rec[i] = b.ReadUint64()
	}
}

func (b *Buffer) ReadUint8Slice(rec []uint8) {
	for i := range rec {
		rec[i] = b.ReadUint8()
	}
}

func (b *Buffer) Bytes() []byte {
	return b.buf
}
