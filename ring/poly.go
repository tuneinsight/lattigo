package ring

import (
	"bufio"
	"fmt"
	"io"

	"github.com/tuneinsight/lattigo/v4/utils/buffer"
)

// Poly is the structure that contains the coefficients of a polynomial.
type Poly struct {
	Coeffs [][]uint64 // Dimension-2 slice of coefficients (re-slice of Buff)
	Buff   []uint64   // Dimension-1 slice of coefficient
}

// NewPoly creates a new polynomial with N coefficients set to zero and Level+1 moduli.
func NewPoly(N, Level int) (pol *Poly) {
	pol = new(Poly)

	pol.Buff = make([]uint64, N*(Level+1))
	pol.Coeffs = make([][]uint64, Level+1)
	for i := 0; i < Level+1; i++ {
		pol.Coeffs[i] = pol.Buff[i*N : (i+1)*N]
	}

	return
}

// Resize resizes the level of the target polynomial to the provided level.
// If the provided level is larger than the current level, then allocates zero
// coefficients, otherwise dereferences the coefficients above the provided level.
func (pol *Poly) Resize(level int) {
	N := pol.N()
	if pol.Level() > level {
		pol.Buff = pol.Buff[:N*(level+1)]
		pol.Coeffs = pol.Coeffs[:level+1]
	} else if level > pol.Level() {
		pol.Buff = append(pol.Buff, make([]uint64, N*(level-pol.Level()))...)
		pol.Coeffs = append(pol.Coeffs, make([][]uint64, level-pol.Level())...)
		for i := 0; i < level+1; i++ {
			pol.Coeffs[i] = pol.Buff[i*N : (i+1)*N]
		}
	}
}

// N returns the number of coefficients of the polynomial, which equals the degree of the Ring cyclotomic polynomial.
func (pol *Poly) N() int {
	return len(pol.Coeffs[0])
}

// Level returns the current number of moduli minus 1.
func (pol *Poly) Level() int {
	return len(pol.Coeffs) - 1
}

// Zero sets all coefficients of the target polynomial to 0.
func (pol *Poly) Zero() {
	ZeroVec(pol.Buff)
}

// CopyNew creates an exact copy of the target polynomial.
func (pol *Poly) CopyNew() (p1 *Poly) {
	p1 = NewPoly(pol.N(), pol.Level())
	copy(p1.Buff, pol.Buff)
	return
}

// Copy copies the coefficients of p0 on p1 within the given Ring. It requires p1 to be at least as big p0.
// Expects the degree of both polynomials to be identical.
func Copy(p0, p1 *Poly) {
	copy(p1.Buff, p0.Buff)
}

// CopyLvl copies the coefficients of p0 on p1 within the given Ring.
// Copies for up to level+1 moduli.
// Expects the degree of both polynomials to be identical.
func CopyLvl(level int, p0, p1 *Poly) {
	copy(p1.Buff[:p1.N()*(level+1)], p0.Buff)
}

// CopyValues copies the coefficients of p1 on the target polynomial.
// Onyl copies minLevel(pol, p1) levels.
// Expects the degree of both polynomials to be identical.
func (pol *Poly) CopyValues(p1 *Poly) {
	if pol != p1 {
		copy(pol.Buff, p1.Buff)
	}
}

// Copy copies the coefficients of p1 on the target polynomial.
// Onyl copies minLevel(pol, p1) levels.
func (pol *Poly) Copy(p1 *Poly) {
	pol.CopyValues(p1)
}

// Equal returns true if the receiver Poly is equal to the provided other Poly.
// This function checks for strict equality between the polynomial coefficients
// (i.e., it does not consider congruence as equality within the ring like
// `Ring.Equal` does).
func (pol Poly) Equal(other *Poly) bool {

	if &pol == other {
		return true
	}

	if &pol != nil && other != nil && len(pol.Buff) == len(other.Buff) {
		for i := range pol.Buff {
			if other.Buff[i] != pol.Buff[i] {
				return false
			}
		}
		return true
	}
	return false
}

// polyBinarySize returns the size in bytes of the Poly object.
func polyBinarySize(N, Level int) (size int) {
	return 16 + N*(Level+1)<<3
}

// BinarySize returns the serialized size of the object in bytes.
func (pol *Poly) BinarySize() (size int) {
	return polyBinarySize(pol.N(), pol.Level())
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
func (pol *Poly) WriteTo(w io.Writer) (int64, error) {

	switch w := w.(type) {
	case buffer.Writer:

		var err error

		var n, inc int

		if n, err = buffer.WriteInt(w, pol.N()); err != nil {
			return int64(n), err
		}

		if inc, err = buffer.WriteInt(w, pol.Level()); err != nil {
			return int64(n + inc), err
		}

		n += inc

		if inc, err = buffer.WriteUint64Slice(w, pol.Buff); err != nil {
			return int64(n + inc), err
		}

		return int64(n + inc), w.Flush()

	default:
		return pol.WriteTo(bufio.NewWriter(w))
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
func (pol *Poly) ReadFrom(r io.Reader) (int64, error) {

	switch r := r.(type) {
	case buffer.Reader:
		var err error

		var n, inc int

		var N int
		if n, err = buffer.ReadInt(r, &N); err != nil {
			return int64(n), fmt.Errorf("cannot ReadFrom: N: %w", err)
		}

		n += inc

		if N <= 0 {
			return int64(n), fmt.Errorf("error ReadFrom: N cannot be 0 or negative")
		}

		var Level int
		if inc, err = buffer.ReadInt(r, &Level); err != nil {
			return int64(n + inc), fmt.Errorf("cannot ReadFrom: Level: %w", err)
		}

		n += inc

		if Level < 0 {
			return int64(n), fmt.Errorf("invalid encoding: Level cannot be negative")
		}

		if pol.Buff == nil || len(pol.Buff) != N*(Level+1) {
			pol.Buff = make([]uint64, N*int(Level+1))
		}

		if inc, err = buffer.ReadUint64Slice(r, pol.Buff); err != nil {
			return int64(n + inc), fmt.Errorf("cannot ReadFrom: pol.Buff: %w", err)
		}

		n += inc

		// Reslice
		if len(pol.Coeffs) != Level+1 {
			pol.Coeffs = make([][]uint64, Level+1)
		}

		for i := 0; i < Level+1; i++ {
			pol.Coeffs[i] = pol.Buff[i*N : (i+1)*N]
		}

		return int64(n), nil

	default:
		return pol.ReadFrom(bufio.NewReader(r))
	}
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (pol *Poly) MarshalBinary() (p []byte, err error) {
	buf := buffer.NewBufferSize(pol.BinarySize())
	_, err = pol.WriteTo(buf)
	return buf.Bytes(), err
}

// UnmarshalBinary decodes a slice of bytes generated by
// MarshalBinary or WriteTo on the object.
func (pol *Poly) UnmarshalBinary(p []byte) (err error) {
	if _, err = pol.ReadFrom(buffer.NewBuffer(p)); err != nil {
		return
	}
	return
}
