package rlwe

import (
	"bufio"
	"fmt"
	"io"

	"github.com/tuneinsight/lattigo/v4/utils/buffer"
	"github.com/tuneinsight/lattigo/v4/utils/structs"
)

// EvaluationKeySetInterface is an interface implementing methods
// to load the RelinearizationKey and GaloisKeys in the Evaluator.
// This interface must support concurrent calls on the methods
// GetGaloisKey and GetRelinearizationKey.
type EvaluationKeySetInterface interface {

	// GetGaloisKey retrieves the Galois key for the automorphism X^{i} -> X^{i*galEl}.
	GetGaloisKey(galEl uint64) (evk *GaloisKey, err error)

	// GetGaloisKeysList returns the list of all the Galois elements
	// for which a Galois key exists in the object.
	GetGaloisKeysList() (galEls []uint64)

	// GetRelinearizationKey retrieves the RelinearizationKey.
	GetRelinearizationKey() (evk *RelinearizationKey, err error)
}

// EvaluationKeySet is a generic struct that complies to the EvaluationKeySetInterface interface.
// This interface can be re-implemented by users to suit application specific requirement.
type EvaluationKeySet struct {
	*RelinearizationKey
	GaloisKeys structs.Map[uint64, GaloisKey]
}

// NewEvaluationKeySet returns a new EvaluationKeySet with nil RelinearizationKey and empty GaloisKeys map.
func NewEvaluationKeySet() (evk *EvaluationKeySet) {
	return &EvaluationKeySet{
		RelinearizationKey: nil,
		GaloisKeys:         map[uint64]*GaloisKey{},
	}
}

// GetGaloisKey retrieves the Galois key for the automorphism X^{i} -> X^{i*galEl}.
func (evk *EvaluationKeySet) GetGaloisKey(galEl uint64) (gk *GaloisKey, err error) {
	var ok bool
	if gk, ok = evk.GaloisKeys[galEl]; !ok {
		return nil, fmt.Errorf("GaloiKey[%d] is nil", galEl)
	}

	return
}

// GetGaloisKeysList returns the list of all the Galois elements
// for which a Galois key exists in the object.
func (evk *EvaluationKeySet) GetGaloisKeysList() (galEls []uint64) {

	if evk == nil || evk.GaloisKeys == nil {
		return []uint64{}
	}

	galEls = make([]uint64, len(evk.GaloisKeys))

	var i int
	for galEl := range evk.GaloisKeys {
		galEls[i] = galEl
		i++
	}

	return
}

// GetRelinearizationKey retrieves the RelinearizationKey.
func (evk *EvaluationKeySet) GetRelinearizationKey() (rk *RelinearizationKey, err error) {
	if evk.RelinearizationKey != nil {
		return evk.RelinearizationKey, nil
	}

	return nil, fmt.Errorf("RelinearizationKey is nil")
}

func (evk *EvaluationKeySet) BinarySize() (size int) {

	size++
	if evk.RelinearizationKey != nil {
		size += evk.RelinearizationKey.BinarySize()
	}

	size++
	if evk.GaloisKeys != nil {
		size += evk.GaloisKeys.BinarySize()
	}

	return
}

func (evk *EvaluationKeySet) MarshalBinary() (p []byte, err error) {
	p = make([]byte, evk.BinarySize())
	_, err = evk.Read(p)
	return
}

func (evk *EvaluationKeySet) Read(p []byte) (n int, err error) {
	var inc int
	if evk.RelinearizationKey != nil {
		p[n] = 1
		n++

		if inc, err = evk.RelinearizationKey.Read(p[n:]); err != nil {
			return n + inc, err
		}

		n += inc

	} else {
		n++
	}

	if evk.GaloisKeys != nil {
		p[n] = 1
		n++

		if inc, err = evk.GaloisKeys.Read(p[n:]); err != nil {

			return n + inc, err
		}

		n += inc

	} else {
		n++
	}

	return
}

func (evk *EvaluationKeySet) WriteTo(w io.Writer) (int64, error) {
	switch w := w.(type) {
	case buffer.Writer:

		var inc int
		var n, inc64 int64
		var err error

		if evk.RelinearizationKey != nil {
			if inc, err = buffer.WriteUint8(w, 1); err != nil {
				return int64(inc), err
			}

			n += int64(inc)

			if inc64, err = evk.RelinearizationKey.WriteTo(w); err != nil {
				return n + inc64, err
			}

			n += inc64

		} else {
			if inc, err = buffer.WriteUint8(w, 0); err != nil {
				return int64(inc), err
			}
			n += int64(inc)
		}

		if evk.GaloisKeys != nil {
			if inc, err = buffer.WriteUint8(w, 1); err != nil {
				return int64(inc), err
			}

			n += int64(inc)

			if inc64, err = evk.GaloisKeys.WriteTo(w); err != nil {
				return n + inc64, err
			}

			n += inc64

		} else {
			if inc, err = buffer.WriteUint8(w, 0); err != nil {
				return int64(inc), err
			}
			n += int64(inc)
		}

		return n, nil

	default:
		return evk.WriteTo(bufio.NewWriter(w))
	}
}

func (evk *EvaluationKeySet) UnmarshalBinary(p []byte) (err error) {
	_, err = evk.Write(p)
	return
}

func (evk *EvaluationKeySet) Write(p []byte) (n int, err error) {
	var inc int
	if p[n] == 1 {
		n++

		if evk.RelinearizationKey == nil {
			evk.RelinearizationKey = new(RelinearizationKey)
		}

		if inc, err = evk.RelinearizationKey.Write(p[n:]); err != nil {
			return n + inc, err
		}

		n += inc

	} else {
		n++
	}

	if p[n] == 1 {
		n++

		if evk.GaloisKeys == nil {
			evk.GaloisKeys = structs.Map[uint64, GaloisKey]{}
		}

		if inc, err = evk.GaloisKeys.Write(p[n:]); err != nil {
			return n + inc, err
		}

		n += inc

	} else {
		n++
	}

	return
}

func (evk *EvaluationKeySet) ReadFrom(r io.Reader) (n int64, err error) {
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

			if evk.RelinearizationKey == nil {
				evk.RelinearizationKey = new(RelinearizationKey)
			}

			if inc64, err = evk.RelinearizationKey.ReadFrom(r); err != nil {
				return n + inc64, err
			}

			n += inc64
		}

		if inc, err = buffer.ReadUint8(r, &hasKey); err != nil {
			return int64(inc), err
		}

		n += int64(inc)

		if hasKey == 1 {

			if evk.GaloisKeys == nil {
				evk.GaloisKeys = structs.Map[uint64, GaloisKey]{}
			}

			if inc64, err = evk.GaloisKeys.ReadFrom(r); err != nil {
				return n + inc64, err
			}

			n += inc64
		}

		return n, nil

	default:
		return evk.ReadFrom(bufio.NewReader(r))
	}
}
