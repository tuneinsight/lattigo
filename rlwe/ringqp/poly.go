package ringqp

import (
	"bufio"
	"io"

	"github.com/google/go-cmp/cmp"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils/buffer"
)

// Poly represents a polynomial in the ring of polynomial modulo Q*P.
// This type is simply the union type between two ring.Poly, each one
// containing the modulus Q and P coefficients of that polynomial.
// The modulus Q represent the ciphertext modulus and the modulus P
// the special primes for the RNS decomposition during homomorphic
// operations involving keys.
type Poly struct {
	Q, P *ring.Poly
}

// NewPoly creates a new polynomial at the given levels.
// If levelQ or levelP are negative, the corresponding polynomial will be nil.
func NewPoly(N, levelQ, levelP int) *Poly {
	var Q, P *ring.Poly

	if levelQ >= 0 {
		Q = ring.NewPoly(N, levelQ)
	}

	if levelP >= 0 {
		P = ring.NewPoly(N, levelP)
	}

	return &Poly{Q, P}
}

// LevelQ returns the level of the polynomial modulo Q.
// Returns -1 if the modulus Q is absent.
func (p *Poly) LevelQ() int {
	if p.Q != nil {
		return p.Q.Level()
	}
	return -1
}

// LevelP returns the level of the polynomial modulo P.
// Returns -1 if the modulus P is absent.
func (p *Poly) LevelP() int {
	if p.P != nil {
		return p.P.Level()
	}
	return -1
}

// Equal returns true if the receiver Poly is equal to the provided other Poly.
func (p *Poly) Equal(other *Poly) (v bool) {
	return cmp.Equal(p.Q, other.Q) && cmp.Equal(p.P, other.P)
}

// Copy copies the coefficients of other on the target polynomial.
// This method simply calls the Copy method for each of its sub-polynomials.
func (p *Poly) Copy(other *Poly) {
	if p.Q != nil {
		copy(p.Q.Buff, other.Q.Buff)
	}

	if p.P != nil {
		copy(p.P.Buff, other.P.Buff)
	}
}

// CopyLvl copies the values of p1 on p2.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func CopyLvl(levelQ, levelP int, p1, p2 *Poly) {

	if p1.Q != nil && p2.Q != nil {
		ring.CopyLvl(levelQ, p1.Q, p2.Q)
	}

	if p1.P != nil && p2.P != nil {
		ring.CopyLvl(levelP, p1.P, p2.P)
	}
}

// CopyNew creates an exact copy of the target polynomial.
func (p *Poly) CopyNew() *Poly {
	if p == nil {
		return nil
	}

	var Q, P *ring.Poly
	if p.Q != nil {
		Q = p.Q.CopyNew()
	}

	if p.P != nil {
		P = p.P.CopyNew()
	}

	return &Poly{Q, P}
}

// Resize resizes the levels of the target polynomial to the provided levels.
// If the provided level is larger than the current level, then allocates zero
// coefficients, otherwise dereferences the coefficients above the provided level.
// Nil polynmials are unafected.
func (p *Poly) Resize(levelQ, levelP int) {
	if p.Q != nil {
		p.Q.Resize(levelQ)
	}

	if p.P != nil {
		p.P.Resize(levelP)
	}
}

// BinarySize returns the size in bytes of the object
// when encoded using Encode.
// Assumes that each coefficient takes 8 bytes.
func (p *Poly) BinarySize() (dataLen int) {

	dataLen = 1

	if p.Q != nil {
		dataLen += p.Q.BinarySize()
	}
	if p.P != nil {
		dataLen += p.P.BinarySize()
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
func (p *Poly) WriteTo(w io.Writer) (n int64, err error) {

	switch w := w.(type) {
	case buffer.Writer:

		var hasQP byte
		if p.Q != nil {
			hasQP = hasQP | 2
		}
		if p.P != nil {
			hasQP = hasQP | 1
		}

		var inc int
		if inc, err = buffer.WriteUint8(w, hasQP); err != nil {
			return int64(n), err
		}

		n += int64(inc)

		if p.Q != nil {
			var inc int64
			if inc, err = p.Q.WriteTo(w); err != nil {
				return n + inc, err
			}

			n += inc
		}

		if p.P != nil {
			var inc int64
			if inc, err = p.P.WriteTo(w); err != nil {
				return n + inc, err
			}

			n += inc
		}

		return n, w.Flush()

	default:
		return p.WriteTo(bufio.NewWriter(w))
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
func (p *Poly) ReadFrom(r io.Reader) (n int64, err error) {
	switch r := r.(type) {
	case buffer.Reader:

		var hasQP byte
		var inc int
		if inc, err = buffer.ReadUint8(r, &hasQP); err != nil {
			return n + int64(inc), err
		}

		n += int64(inc)

		if hasQP&2 == 2 {

			if p.Q == nil {
				p.Q = new(ring.Poly)
			}

			var inc64 int64
			if inc64, err = p.Q.ReadFrom(r); err != nil {
				return n + inc64, err
			}

			n += inc64
		}

		if hasQP&1 == 1 {

			if p.P == nil {
				p.P = new(ring.Poly)
			}

			var inc int64
			if inc, err = p.P.ReadFrom(r); err != nil {
				return n + inc, err
			}

			n += inc
		}

		return

	default:
		return p.ReadFrom(bufio.NewReader(r))
	}
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (p *Poly) MarshalBinary() (data []byte, err error) {
	buf := buffer.NewBufferSize(p.BinarySize())
	_, err = p.WriteTo(buf)
	return buf.Bytes(), err
}

// UnmarshalBinary decodes a slice of bytes generated by
// MarshalBinary or WriteTo on the object.
func (p *Poly) UnmarshalBinary(data []byte) (err error) {
	_, err = p.ReadFrom(buffer.NewBuffer(data))
	return err
}

// Encode encodes the object into a binary form on a preallocated slice of bytes
// and returns the number of bytes written.
func (p *Poly) Encode(data []byte) (n int, err error) {
	var inc int

	if p.Q != nil {
		data[0] = 1
	}

	if p.P != nil {
		data[1] = 1
	}

	n = 2

	if data[0] == 1 {
		if inc, err = p.Q.Encode(data[n:]); err != nil {
			return
		}
		n += inc
	}

	if data[1] == 1 {
		if inc, err = p.P.Encode(data[n:]); err != nil {
			return
		}
		n += inc
	}

	return
}

// Decode decodes a slice of bytes generated by Encode
// on the object and returns the number of bytes read.
func (p *Poly) Decode(data []byte) (n int, err error) {

	var inc int
	n = 2

	if data[0] == 1 {

		if p.Q == nil {
			p.Q = new(ring.Poly)
		}

		if inc, err = p.Q.Decode(data[n:]); err != nil {
			return
		}
		n += inc
	}

	if data[1] == 1 {

		if p.P == nil {
			p.P = new(ring.Poly)
		}

		if inc, err = p.P.Decode(data[n:]); err != nil {
			return
		}
		n += inc
	}

	return
}
