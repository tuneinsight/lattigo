package rlwe

import (
	"encoding/json"
	"fmt"
	"io"
	"math/big"

	"github.com/google/go-cmp/cmp"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// MetaData is a struct storing metadata.
type MetaData struct {
	PlaintextMetaData
	CiphertextMetaData
}

// CopyNew returns a copy of the target.
func (m MetaData) CopyNew() *MetaData {
	return &m
}

func (m *MetaData) Equal(other *MetaData) (res bool) {
	return m.PlaintextMetaData.Equal(&other.PlaintextMetaData) && m.CiphertextMetaData.Equal(&other.CiphertextMetaData)
}

// BinarySize returns the size in bytes that the object once marshalled into a binary form.
func (m MetaData) BinarySize() int {
	return 44 + m.PlaintextMetaData.BinarySize() + m.CiphertextMetaData.BinarySize()
}

// WriteTo writes the object on an [io.Writer]. It implements the [io.WriterTo]
// interface, and will write exactly object.BinarySize() bytes on w.
func (m MetaData) WriteTo(w io.Writer) (int64, error) {
	if p, err := m.MarshalBinary(); err != nil {
		return 0, err
	} else {
		if n, err := w.Write(p); err != nil {
			return int64(n), err
		} else {
			return int64(n), nil
		}
	}
}

// ReadFrom reads on the object from an [io.Writer]. It implements the
// [io.ReaderFrom] interface.
//
// Unless r implements the [buffer.Reader] interface (see see lattigo/utils/buffer/reader.go),
// it will be wrapped into a [bufio.Reader]. Since this requires allocation, it
// is preferable to pass a [buffer.Reader] directly:
//
//   - When reading multiple values from a [io.Reader], it is preferable to first
//     first wrap [io.Reader] in a pre-allocated [bufio.Reader].
//   - When reading from a var b []byte, it is preferable to pass a buffer.NewBuffer(b)
//     as w (see lattigo/utils/buffer/buffer.go).
func (m *MetaData) ReadFrom(r io.Reader) (int64, error) {
	p := make([]byte, m.BinarySize())
	if n, err := r.Read(p); err != nil {
		return int64(n), err
	} else {
		return int64(n), m.UnmarshalBinary(p)
	}
}

func (m MetaData) MarshalJSON() (p []byte, err error) {
	aux := &struct {
		PlaintextMetaData  PlaintextMetaData
		CiphertextMetaData CiphertextMetaData
	}{
		PlaintextMetaData:  m.PlaintextMetaData,
		CiphertextMetaData: m.CiphertextMetaData,
	}

	return json.Marshal(aux)
}

func (m MetaData) MarshalBinary() (p []byte, err error) {
	return m.MarshalJSON()
}

func (m *MetaData) UnmarshalJSON(p []byte) (err error) {
	aux := &struct {
		PlaintextMetaData  PlaintextMetaData
		CiphertextMetaData CiphertextMetaData
	}{
		PlaintextMetaData:  m.PlaintextMetaData,
		CiphertextMetaData: m.CiphertextMetaData,
	}

	if err = json.Unmarshal(p, aux); err != nil {
		return
	}

	m.PlaintextMetaData = aux.PlaintextMetaData
	m.CiphertextMetaData = aux.CiphertextMetaData
	return
}

func (m *MetaData) UnmarshalBinary(p []byte) (err error) {
	return m.UnmarshalJSON(p)
}

// PlaintextMetaData is a struct storing metadata related to the plaintext.
type PlaintextMetaData struct {
	// Scale is the scaling factor of the plaintext.
	Scale Scale

	// LogDimensions is the Log2 of the 2D plaintext matrix dimensions.
	LogDimensions ring.Dimensions

	// IsBatched is a flag indicating if the underlying plaintext is encoded
	// in such a way that product in R[X]/(X^N+1) acts as a point-wise multiplication
	// in the plaintext space.
	IsBatched bool

	// IsBitReversed is a flag indicating if the underlying plaintext is
	// bit-reversed. This can be true for both batch and non-batched plaintexts.
	IsBitReversed bool
}

// Slots returns the total number of slots that the plaintext holds.
func (m PlaintextMetaData) Slots() int {
	return 1 << m.LogSlots()
}

// LogSlots returns the log2 of the total number of slots that the plaintext holds.
func (m PlaintextMetaData) LogSlots() int {
	return m.LogDimensions.Cols + m.LogDimensions.Rows
}

// LogScale returns log2(scale).
func (m PlaintextMetaData) LogScale() float64 {
	ln := bignum.Log(&m.Scale.Value)
	ln.Quo(ln, bignum.Log2(ln.Prec()))
	log2, _ := ln.Float64()
	return log2
}

func (m *PlaintextMetaData) Equal(other *PlaintextMetaData) (res bool) {
	res = cmp.Equal(&m.Scale, &other.Scale)
	res = res && m.IsBatched == other.IsBatched
	res = res && m.IsBitReversed == other.IsBitReversed
	res = res && m.LogDimensions == other.LogDimensions
	return
}

// BinarySize returns the size in bytes that the object once marshalled into a binary form.
func (m PlaintextMetaData) BinarySize() int {
	return 84 + m.Scale.BinarySize()
}

// WriteTo writes the object on an [io.Writer]. It implements the [io.WriterTo]
// interface, and will write exactly object.BinarySize() bytes on w.
//
// Unless w implements the [buffer.Writer] interface (see lattigo/utils/buffer/writer.go),
// it will be wrapped into a [bufio.Writer]. Since this requires allocations, it
// is preferable to pass a [buffer.Writer] directly:
//
//   - When writing multiple times to a [io.Writer], it is preferable to first wrap the
//     io.Writer in a pre-allocated [bufio.Writer].
//   - When writing to a pre-allocated var b []byte, it is preferable to pass
//     buffer.NewBuffer(b) as w (see lattigo/utils/buffer/buffer.go).
func (m PlaintextMetaData) WriteTo(w io.Writer) (int64, error) {
	if p, err := m.MarshalBinary(); err != nil {
		return 0, err
	} else {
		if n, err := w.Write(p); err != nil {
			return int64(n), err
		} else {
			return int64(n), nil
		}
	}
}

// ReadFrom reads on the object from an [io.Writer]. It implements the
// [io.ReaderFrom] interface.
//
// Unless r implements the [buffer.Reader] interface (see see lattigo/utils/buffer/reader.go),
// it will be wrapped into a [bufio.Reader]. Since this requires allocation, it
// is preferable to pass a [buffer.Reader] directly:
//
//   - When reading multiple values from a [io.Reader], it is preferable to first
//     first wrap [io.Reader] in a pre-allocated [bufio.Reader].
//   - When reading from a var b []byte, it is preferable to pass a buffer.NewBuffer(b)
//     as w (see lattigo/utils/buffer/buffer.go).
func (m *PlaintextMetaData) ReadFrom(r io.Reader) (int64, error) {
	p := make([]byte, m.BinarySize())
	if n, err := r.Read(p); err != nil {
		return int64(n), err
	} else {
		return int64(n), m.UnmarshalBinary(p)
	}
}

func (m PlaintextMetaData) MarshalJSON() (p []byte, err error) {

	var IsBatched uint8
	if m.IsBatched {
		IsBatched = 1
	}

	var IsBitReversed uint8
	if m.IsBitReversed {
		IsBitReversed = 1
	}

	aux := &struct {
		Scale         Scale
		IsBatched     string
		IsBitReversed string
		LogDimensions [2]string
	}{
		Scale:         m.Scale,
		IsBatched:     fmt.Sprintf("0x%02x", IsBatched),
		IsBitReversed: fmt.Sprintf("0x%02x", IsBitReversed),
		LogDimensions: [2]string{fmt.Sprintf("0x%02x", uint8(m.LogDimensions.Rows)), fmt.Sprintf("0x%02x", uint8(m.LogDimensions.Cols))},
	}

	p, err = json.Marshal(aux)

	return
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (m PlaintextMetaData) MarshalBinary() (p []byte, err error) {
	return m.MarshalJSON()
}

func (m *PlaintextMetaData) UnmarshalJSON(p []byte) (err error) {
	aux := &struct {
		Scale         Scale
		IsBatched     string
		IsBitReversed string
		LogDimensions [2]string
	}{}

	if err = json.Unmarshal(p, aux); err != nil {
		return
	}

	m.Scale = aux.Scale

	if y, err := hexconv(aux.IsBatched); err != nil {
		return err
	} else if y == 1 {
		m.IsBatched = true
	} else {
		m.IsBatched = false
	}

	if y, err := hexconv(aux.IsBitReversed); err != nil {
		return err
	} else if y == 1 {
		m.IsBitReversed = true
	} else {
		m.IsBitReversed = false
	}

	logRows, err := hexconv(aux.LogDimensions[0])

	if err != nil {
		return err
	}

	logCols, err := hexconv(aux.LogDimensions[1])

	if err != nil {
		return err
	}

	m.LogDimensions = ring.Dimensions{Rows: int(int8(logRows)), Cols: int(int8(logCols))}

	return
}

// UnmarshalBinary decodes a slice of bytes generated by
// [PlaintextMetaData.MarshalBinary] or [PlaintextMetaData.WriteTo] on the object.
func (m *PlaintextMetaData) UnmarshalBinary(p []byte) (err error) {
	return m.UnmarshalJSON(p)
}

// CiphertextMetaData is a struct storing metadata related to the ciphertext.
type CiphertextMetaData struct {
	// IsNTT is a flag indicating if the ciphertext is in the NTT domain.
	IsNTT bool
	// IsMontgomery is a flag indicating if the ciphertext is in the Montgomery domain.
	IsMontgomery bool
}

// Equal returns true if two MetaData structs are identical.
func (m *CiphertextMetaData) Equal(other *CiphertextMetaData) (res bool) {
	res = m.IsNTT == other.IsNTT
	res = res && m.IsMontgomery == other.IsMontgomery
	return
}

// BinarySize returns the size in bytes that the object once marshalled into a binary form.
func (m *CiphertextMetaData) BinarySize() int {
	return 38
}

// WriteTo writes the object on an [io.Writer]. It implements the [io.WriterTo]
// interface, and will write exactly object.BinarySize() bytes on w.
//
// Unless w implements the [buffer.Writer] interface (see lattigo/utils/buffer/writer.go),
// it will be wrapped into a [bufio.Writer]. Since this requires allocations, it
// is preferable to pass a [buffer.Writer] directly:
//
//   - When writing multiple times to a [io.Writer], it is preferable to first wrap the
//     io.Writer in a pre-allocated [bufio.Writer].
//   - When writing to a pre-allocated var b []byte, it is preferable to pass
//     buffer.NewBuffer(b) as w (see lattigo/utils/buffer/buffer.go).
func (m *CiphertextMetaData) WriteTo(w io.Writer) (int64, error) {
	if p, err := m.MarshalBinary(); err != nil {
		return 0, err
	} else {
		if n, err := w.Write(p); err != nil {
			return int64(n), err
		} else {
			return int64(n), nil
		}
	}
}

// ReadFrom reads on the object from an [io.Writer]. It implements the
// [io.ReaderFrom] interface.
//
// Unless r implements the [buffer.Reader] interface (see see lattigo/utils/buffer/reader.go),
// it will be wrapped into a [bufio.Reader]. Since this requires allocation, it
// is preferable to pass a [buffer.Reader] directly:
//
//   - When reading multiple values from a [io.Reader], it is preferable to first
//     first wrap [io.Reader] in a pre-allocated [bufio.Reader].
//   - When reading from a var b []byte, it is preferable to pass a buffer.NewBuffer(b)
//     as w (see lattigo/utils/buffer/buffer.go).
func (m *CiphertextMetaData) ReadFrom(r io.Reader) (int64, error) {
	p := make([]byte, m.BinarySize())
	if n, err := r.Read(p); err != nil {
		return int64(n), err
	} else {
		return int64(n), m.UnmarshalBinary(p)
	}
}

func (m CiphertextMetaData) MarshalJSON() (p []byte, err error) {
	var IsNTT, IsMontgomery uint8

	if m.IsNTT {
		IsNTT = 1
	}

	if m.IsMontgomery {
		IsMontgomery = 1
	}

	aux := &struct {
		IsNTT        string
		IsMontgomery string
	}{
		IsNTT:        fmt.Sprintf("0x%02x", IsNTT),
		IsMontgomery: fmt.Sprintf("0x%02x", IsMontgomery),
	}

	return json.Marshal(aux)
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (m CiphertextMetaData) MarshalBinary() (p []byte, err error) {
	return m.MarshalJSON()
}

func (m *CiphertextMetaData) UnmarshalJSON(p []byte) (err error) {
	aux := &struct {
		IsNTT        string
		IsMontgomery string
	}{}

	if err = json.Unmarshal(p, aux); err != nil {
		return
	}

	if y, err := hexconv(aux.IsNTT); err != nil {
		return err
	} else if y == 1 {
		m.IsNTT = true
	}

	if y, err := hexconv(aux.IsMontgomery); err != nil {
		return err
	} else if y == 1 {
		m.IsMontgomery = true
	}

	return
}

// UnmarshalBinary decodes a slice of bytes generated by
// [CiphertextMetaData.MarshalBinary] or [CiphertextMetaData.WriteTo] on the object.
func (m *CiphertextMetaData) UnmarshalBinary(p []byte) (err error) {
	return m.UnmarshalJSON(p)
}

func hexconv(x string) (uint64, error) {
	yBig, err := new(big.Int).SetString(x, 0)
	if !err {
		return 0, fmt.Errorf("hexconv: unsuccessful SetString")
	}
	return yBig.Uint64(), nil
}
