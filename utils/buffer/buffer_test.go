package buffer

import (
	"encoding/binary"
	"errors"
	"io"
	"testing"
)

// TestReadUintSliceTruncated checks that reading a slice of 16-, 32- or
// 64-bit words from a reader that holds fewer bytes than requested returns
// an ordinary error instead of recursing indefinitely.
func TestReadUintSliceTruncated(t *testing.T) {

	const n = 1024

	// read dispatches on the word size and returns the number of bytes
	// consumed, as reported by the ReadUint*Slice functions.
	read := func(wordSize int, r Reader) (int64, error) {
		switch wordSize {
		case 2:
			return ReadUint16Slice(r, make([]uint16, n))
		case 4:
			return ReadUint32Slice(r, make([]uint32, n))
		default:
			return ReadUint64Slice(r, make([]uint64, n))
		}
	}

	testCases := []struct {
		name     string
		wordSize int
		dataLen  int
	}{
		{"Uint16/EmptyInput", 2, 0},
		{"Uint16/PartialWord", 2, 1},
		{"Uint16/OneByteShort", 2, 2*n - 1},
		{"Uint32/EmptyInput", 4, 0},
		{"Uint32/PartialWord", 4, 3},
		{"Uint32/OneByteShort", 4, 4*n - 1},
		{"Uint64/EmptyInput", 8, 0},
		{"Uint64/PartialWord", 8, 7},
		{"Uint64/OneByteShort", 8, 8*n - 1},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			b := NewBuffer(make([]byte, tc.dataLen))

			nRead, err := read(tc.wordSize, b)
			if !errors.Is(err, io.ErrUnexpectedEOF) {
				t.Fatalf("expected io.ErrUnexpectedEOF, got %v", err)
			}
			if nRead != int64(tc.dataLen) {
				t.Fatalf("expected %d bytes consumed, got %d", tc.dataLen, nRead)
			}
		})
	}
}

// TestReadUintSliceComplete checks that complete and empty inputs are
// decoded as before.
func TestReadUintSliceComplete(t *testing.T) {

	const n = 1024

	t.Run("Uint64/Complete", func(t *testing.T) {
		data := make([]byte, 8*n)
		for i := range data {
			data[i] = byte(i)
		}

		c := make([]uint64, n)
		nRead, err := ReadUint64Slice(NewBuffer(data), c)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if nRead != int64(len(data)) {
			t.Fatalf("expected %d bytes read, got %d", len(data), nRead)
		}
		for i := range c {
			if expected := binary.LittleEndian.Uint64(data[8*i:]); c[i] != expected {
				t.Fatalf("element %d: expected %d, got %d", i, expected, c[i])
			}
		}
	})

	t.Run("EmptySlice", func(t *testing.T) {
		b := NewBuffer(make([]byte, 8))
		if nRead, err := ReadUint64Slice(b, nil); err != nil || nRead != 0 {
			t.Fatalf("expected (0, nil), got (%d, %v)", nRead, err)
		}
		if b.Size() != 8 {
			t.Fatalf("expected reader to be untouched, %d bytes left", b.Size())
		}
	})
}
