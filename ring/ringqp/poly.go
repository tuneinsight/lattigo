package ringqp

import (
	"bufio"
	"io"

	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils"
	"github.com/tuneinsight/lattigo/v5/utils/buffer"
)

// Poly represents a polynomial in the ring of polynomial modulo Q*P.
// This type is simply the union type between two ring.Poly, each one
// containing the modulus Q and P coefficients of that polynomial.
// The modulus Q represent the ciphertext modulus and the modulus P
// the special primes for the RNS decomposition during homomorphic
// operations involving keys.
type Poly struct {
	Q, P ring.Poly
}

// NewPoly creates a new polynomial at the given levels.
// If levelQ or levelP are negative, the corresponding polynomial will be nil.
func NewPoly(N, levelQ, levelP int) Poly {
	var Q, P ring.Poly

	if levelQ >= 0 {
		Q = ring.NewPoly(N, levelQ)
	}

	if levelP >= 0 {
		P = ring.NewPoly(N, levelP)
	}

	return Poly{Q, P}
}

// LevelQ returns the level of the polynomial modulo Q.
// Returns -1 if the modulus Q is absent.
func (p Poly) LevelQ() int {
	return p.Q.Level()
}

// LevelP returns the level of the polynomial modulo P.
// Returns -1 if the modulus P is absent.
func (p Poly) LevelP() int {
	return p.P.Level()
}

// Equal returns true if the receiver Poly is equal to the provided other Poly.
func (p Poly) Equal(other *Poly) (v bool) {
	return p.Q.Equal(&other.Q) && p.P.Equal(&other.P)
}

// Copy copies the coefficients of other on the target polynomial.
// This method simply calls the Copy method for each of its sub-polynomials.
func (p *Poly) Copy(other Poly) {
	p.CopyLvl(utils.Min(p.LevelQ(), other.LevelQ()), utils.Min(p.LevelP(), other.LevelP()), other)
}

// CopyLvl copies the values of other on the target polynomial.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (p *Poly) CopyLvl(levelQ, levelP int, other Poly) {

	if p.Q.Level() != -1 && other.Q.Level() != -1 {
		p.Q.CopyLvl(levelQ, other.Q)
	}

	if p.P.Level() != -1 && other.P.Level() != -1 {
		p.P.CopyLvl(levelP, other.P)
	}
}

// CopyNew creates an exact copy of the target polynomial.
func (p Poly) CopyNew() *Poly {
	return &Poly{*p.Q.CopyNew(), *p.P.CopyNew()}
}

// Resize resizes the levels of the target polynomial to the provided levels.
// If the provided level is larger than the current level, then allocates zero
// coefficients, otherwise dereferences the coefficients above the provided level.
// Nil polynomials are unaffected.
func (p *Poly) Resize(levelQ, levelP int) {
	p.Q.Resize(levelQ)
	p.P.Resize(levelP)
}

// BinarySize returns the serialized size of the object in bytes.
// It assumes that each coefficient takes 8 bytes.
func (p Poly) BinarySize() (dataLen int) {
	return p.Q.BinarySize() + p.P.BinarySize()
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
func (p Poly) WriteTo(w io.Writer) (n int64, err error) {

	switch w := w.(type) {
	case buffer.Writer:

		var inc int64

		if inc, err = p.Q.WriteTo(w); err != nil {
			return n + inc, err
		}

		n += inc

		if inc, err = p.P.WriteTo(w); err != nil {
			return n + inc, err
		}

		n += inc

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

		var inc int64

		if inc, err = p.Q.ReadFrom(r); err != nil {
			return n + inc, err
		}

		n += inc

		if inc, err = p.P.ReadFrom(r); err != nil {
			return n + inc, err
		}

		n += inc

		return

	default:
		return p.ReadFrom(bufio.NewReader(r))
	}
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (p Poly) MarshalBinary() (data []byte, err error) {
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
