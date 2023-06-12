package rlwe

import (
	"bufio"
	"fmt"
	"io"

	"github.com/tuneinsight/lattigo/v4/utils/buffer"
	"github.com/tuneinsight/lattigo/v4/utils/structs"
)

// EvaluationKeySet is an interface implementing methods
// to load the RelinearizationKey and GaloisKeys in the Evaluator.
// Implementations of this interface must be safe for concurrent use.
type EvaluationKeySet interface {

	// GetGaloisKey retrieves the Galois key for the automorphism X^{i} -> X^{i*galEl}.
	GetGaloisKey(galEl uint64) (evk *GaloisKey, err error)

	// GetGaloisKeysList returns the list of all the Galois elements
	// for which a Galois key exists in the object.
	GetGaloisKeysList() (galEls []uint64)

	// GetRelinearizationKey retrieves the RelinearizationKey.
	GetRelinearizationKey() (evk *RelinearizationKey, err error)
}

// MemEvaluationKeySet is a basic in-memory implementation of the EvaluationKeySet interface.
type MemEvaluationKeySet struct {
	Rlk *RelinearizationKey
	Gks structs.Map[uint64, GaloisKey]
}

// NewMemEvaluationKeySet returns a new EvaluationKeySet with the provided RelinearizationKey and GaloisKeys.
func NewMemEvaluationKeySet(relinKey *RelinearizationKey, galoisKeys ...*GaloisKey) (eks *MemEvaluationKeySet) {
	eks = &MemEvaluationKeySet{Gks: map[uint64]*GaloisKey{}}
	eks.Rlk = relinKey
	for _, k := range galoisKeys {
		eks.Gks[k.GaloisElement] = k
	}
	return eks
}

// GetGaloisKey retrieves the Galois key for the automorphism X^{i} -> X^{i*galEl}.
func (evk *MemEvaluationKeySet) GetGaloisKey(galEl uint64) (gk *GaloisKey, err error) {
	var ok bool
	if gk, ok = evk.Gks[galEl]; !ok {
		return nil, fmt.Errorf("GaloiKey[%d] is nil", galEl)
	}

	return
}

// GetGaloisKeysList returns the list of all the Galois elements
// for which a Galois key exists in the object.
func (evk *MemEvaluationKeySet) GetGaloisKeysList() (galEls []uint64) {

	if evk == nil || evk.Gks == nil {
		return []uint64{}
	}

	galEls = make([]uint64, len(evk.Gks))

	var i int
	for galEl := range evk.Gks {
		galEls[i] = galEl
		i++
	}

	return
}

// GetRelinearizationKey retrieves the RelinearizationKey.
func (evk *MemEvaluationKeySet) GetRelinearizationKey() (rk *RelinearizationKey, err error) {
	if evk.Rlk != nil {
		return evk.Rlk, nil
	}

	return nil, fmt.Errorf("RelinearizationKey is nil")
}

func (evk *MemEvaluationKeySet) BinarySize() (size int) {

	size++
	if evk.Rlk != nil {
		size += evk.Rlk.BinarySize()
	}

	size++
	if evk.Gks != nil {
		size += evk.Gks.BinarySize()
	}

	return
}

// WriteTo writes the object on an io.Writer. It implements the io.WriterTo
// interface, and will write exactly object.BinarySize() bytes on w.
//
// Unless w implements the buffer.Writer interface (see lattigo/utils/buffer/writer.go),
// it will be wrapped into a bufio.Writer. Since this requires allocations, it
// is preferable to pass a buffer.Writer directly:
//
//   - When writing multiple times to a io.Writer, it is preferable to first wrap the
//     io.Writer in a pre-allocated bufio.Writer.
//   - When writing to a pre-allocated var b []byte, it is preferable to pass
//     buffer.NewBuffer(b) as w (see lattigo/utils/buffer/buffer.go).
func (evk *MemEvaluationKeySet) WriteTo(w io.Writer) (int64, error) {
	switch w := w.(type) {
	case buffer.Writer:

		var inc int
		var n, inc64 int64
		var err error

		if evk.Rlk != nil {
			if inc, err = buffer.WriteUint8(w, 1); err != nil {
				return int64(inc), err
			}

			n += int64(inc)

			if inc64, err = evk.Rlk.WriteTo(w); err != nil {
				return n + inc64, err
			}

			n += inc64

		} else {
			if inc, err = buffer.WriteUint8(w, 0); err != nil {
				return int64(inc), err
			}
			n += int64(inc)
		}

		if evk.Gks != nil {
			if inc, err = buffer.WriteUint8(w, 1); err != nil {
				return int64(inc), err
			}

			n += int64(inc)

			if inc64, err = evk.Gks.WriteTo(w); err != nil {
				return n + inc64, err
			}

			n += inc64

		} else {
			if inc, err = buffer.WriteUint8(w, 0); err != nil {
				return int64(inc), err
			}
			n += int64(inc)
		}

		return n, w.Flush()

	default:
		return evk.WriteTo(bufio.NewWriter(w))
	}
}

// ReadFrom reads on the object from an io.Writer. It implements the
// io.ReaderFrom interface.
//
// Unless r implements the buffer.Reader interface (see see lattigo/utils/buffer/reader.go),
// it will be wrapped into a bufio.Reader. Since this requires allocation, it
// is preferable to pass a buffer.Reader directly:
//
//   - When reading multiple values from a io.Reader, it is preferable to first
//     first wrap io.Reader in a pre-allocated bufio.Reader.
//   - When reading from a var b []byte, it is preferable to pass a buffer.NewBuffer(b)
//     as w (see lattigo/utils/buffer/buffer.go).
func (evk *MemEvaluationKeySet) ReadFrom(r io.Reader) (n int64, err error) {
	switch r := r.(type) {
	case buffer.Reader:
		var inc int
		var n, inc64 int64
		var err error

		var hasKey uint8

		if inc, err = buffer.ReadUint8(r, &hasKey); err != nil {
			return int64(inc), err
		}

		n += int64(inc)

		if hasKey == 1 {

			if evk.Rlk == nil {
				evk.Rlk = new(RelinearizationKey)
			}

			if inc64, err = evk.Rlk.ReadFrom(r); err != nil {
				return n + inc64, err
			}

			n += inc64
		}

		if inc, err = buffer.ReadUint8(r, &hasKey); err != nil {
			return int64(inc), err
		}

		n += int64(inc)

		if hasKey == 1 {

			if evk.Gks == nil {
				evk.Gks = structs.Map[uint64, GaloisKey]{}
			}

			if inc64, err = evk.Gks.ReadFrom(r); err != nil {
				return n + inc64, err
			}

			n += inc64
		}

		return n, nil

	default:
		return evk.ReadFrom(bufio.NewReader(r))
	}
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (evk *MemEvaluationKeySet) MarshalBinary() (p []byte, err error) {
	buf := buffer.NewBufferSize(evk.BinarySize())
	_, err = evk.WriteTo(buf)
	return buf.Bytes(), err
}

// UnmarshalBinary decodes a slice of bytes generated by
// MarshalBinary or WriteTo on the object.
func (evk *MemEvaluationKeySet) UnmarshalBinary(p []byte) (err error) {
	_, err = evk.ReadFrom(buffer.NewBuffer(p))
	return
}
