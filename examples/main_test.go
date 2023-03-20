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

	data := make([]byte, pol.MarshalBinarySize64())

	b.Run("Encode/Native", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			if _, err = pol.Encode64(data); err != nil {
				b.Fatal(err)
			}

		}
	})

	fmt.Println(data[:8])

	writer := NewWriter(len(data))

	w := utils.NewWriter(writer)

	b.Run("Encode/utils.Writer", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			if _, err = pol.Write(w); err != nil {
				b.Fatal(err)
			}
			writer.n = 0
		}
	})

	fmt.Println(data[:8])
	fmt.Println(writer.buff[:8])
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
