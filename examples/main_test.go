package main_test

import (
	"fmt"
	"testing"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

func Benchmark(b *testing.B) {

	LogN := 15
	Qi := []uint64{0x1fffffffffe00001, 0x1fffffffffc80001, 0x1fffffffffb40001, 0x1fffffffff500001,
		0x1fffffffff380001, 0x1fffffffff000001, 0x1ffffffffef00001, 0x1ffffffffee80001,
		0x1ffffffffeb40001, 0x1ffffffffe780001, 0x1ffffffffe600001, 0x1ffffffffe4c0001}

	r, err := ring.NewRing(1<<LogN, Qi)

	if err != nil {
		b.Fatal(err)
	}

	prng, err := utils.NewPRNG()
	if err != nil {
		b.Fatal(err)
	}

	sampler := ring.NewUniformSampler(prng, r)

	pol := sampler.ReadNew()

	b.Run("Read([]byte)", func(b *testing.B) {
		data := make([]byte, pol.MarshalBinarySize())
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if _, err = pol.Read(data); err != nil {
				b.Fatal(err)
			}
		}
	})

	b.Run("WriteTo(utils.Writer)", func(b *testing.B) {
		writer := NewWriter(pol.MarshalBinarySize())
		w := utils.NewWriter(writer)
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if _, err = pol.WriteTo(w); err != nil {
				b.Fatal(err)
			}

			if err = w.Flush(); err != nil {
				b.Fatal(err)
			}

			writer.n = 0
		}
	})

	b.Run("Write([]byte)", func(b *testing.B) {

		data := make([]byte, pol.MarshalBinarySize())

		if _, err = pol.Read(data); err != nil {
			b.Fatal(err)
		}

		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if _, err = pol.Write(data); err != nil {
				b.Fatal(err)
			}
		}
	})

	b.Run("ReadFrom(utils.Reader)", func(b *testing.B) {

		writer := NewWriter(pol.MarshalBinarySize())

		w := utils.NewWriter(writer)

		if _, err = pol.WriteTo(w); err != nil {
			b.Fatal(err)
		}

		if err = w.Flush(); err != nil {
			b.Fatal(err)
		}

		reader := NewReader(writer.buff)

		r := utils.NewReader(reader)

		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if _, err = pol.ReadFrom(r); err != nil {
				b.Fatal(err)
			}

			reader.n = 0
		}
	})

}

type Reader struct {
	buff []byte
	n    int
}

func NewReader(b []byte) *Reader {
	return &Reader{
		buff: b,
		n:    0,
	}
}

func (r *Reader) Read(b []byte) (n int, err error) {
	if len(b) > len(r.buff[r.n:]) {
		return 0, fmt.Errorf("cannot read: len(b)=%d > %d", len(b), len(r.buff[r.n:]))
	}

	copy(b, r.buff[r.n:])

	r.n += len(b)

	return len(b), nil
}

type Writer struct {
	buff []byte
	n    int
}

func NewWriter(size int) *Writer {
	return &Writer{
		buff: make([]byte, size),
		n:    0,
	}
}

func (w *Writer) Write(b []byte) (n int, err error) {

	if len(b) > len(w.buff[w.n:]) {
		return 0, fmt.Errorf("cannot write len(b)=%d > %d", len(b), len(w.buff[w.n:]))
	}

	copy(w.buff[w.n:], b)

	w.n += len(b)

	return len(b), nil
}
