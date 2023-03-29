package rlwe

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"io"

	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum/polynomial"
	"github.com/tuneinsight/lattigo/v4/utils/buffer"
)

// PowerBasis is a struct storing powers of a ciphertext.
type PowerBasis struct {
	polynomial.Basis
	Value map[int]*Ciphertext
}

// NewPowerBasis creates a new PowerBasis. It takes as input a ciphertext
// and a basistype. The struct treats the input ciphertext as a monomial X and
// can be used to generates power of this monomial X^{n} in the given BasisType.
func NewPowerBasis(ct *Ciphertext, basis polynomial.Basis) (p *PowerBasis) {
	p = new(PowerBasis)
	p.Value = make(map[int]*Ciphertext)
	p.Value[1] = ct.CopyNew()
	p.Basis = basis
	return
}

func (p *PowerBasis) BinarySize() (size int) {
	size = 5 // Type & #Ct
	for _, ct := range p.Value {
		size += 4 + ct.BinarySize()
	}

	return
}

// MarshalBinary encodes the target on a slice of bytes.
func (p *PowerBasis) MarshalBinary() (data []byte, err error) {
	data = make([]byte, p.BinarySize())
	_, err = p.Read(data)
	return
}

func (p *PowerBasis) WriteTo(w io.Writer) (n int64, err error) {

	switch w := w.(type) {
	case buffer.Writer:

		var inc1 int

		if inc1, err = buffer.WriteUint8(w, uint8(p.Basis)); err != nil {
			return n + int64(inc1), err
		}

		n += int64(inc1)

		if inc1, err = buffer.WriteUint32(w, uint32(len(p.Value))); err != nil {
			return n + int64(inc1), err
		}

		n += int64(inc1)

		for _, key := range utils.GetSortedKeys(p.Value) {

			ct := p.Value[key]

			if inc1, err = buffer.WriteUint32(w, uint32(key)); err != nil {
				return n + int64(inc1), err
			}

			n += int64(inc1)

			var inc2 int64
			if inc2, err = ct.WriteTo(w); err != nil {
				return n + inc2, err
			}

			n += inc2
		}

		return

	default:
		return p.WriteTo(bufio.NewWriter(w))
	}
}

func (p *PowerBasis) Read(data []byte) (n int, err error) {

	if len(data) < p.BinarySize() {
		return n, fmt.Errorf("cannot Read: len(data)=%d < %d", len(data), p.BinarySize())
	}

	data[n] = uint8(p.Basis)
	n++

	binary.LittleEndian.PutUint32(data[n:], uint32(len(p.Value)))
	n += 4

	for _, key := range utils.GetSortedKeys(p.Value) {

		ct := p.Value[key]

		binary.LittleEndian.PutUint32(data[n:], uint32(key))
		n += 4

		var inc int
		if inc, err = ct.Read(data[n:]); err != nil {
			return n + inc, err
		}

		n += inc
	}

	return
}

// UnmarshalBinary decodes a slice of bytes on the target.
func (p *PowerBasis) UnmarshalBinary(data []byte) (err error) {
	_, err = p.Write(data)
	return
}

func (p *PowerBasis) ReadFrom(r io.Reader) (n int64, err error) {
	switch r := r.(type) {
	case buffer.Reader:
		var inc1 int

		var Basis uint8

		if inc1, err = buffer.ReadUint8(r, &Basis); err != nil {
			return n + int64(inc1), err
		}

		n += int64(inc1)

		p.Basis = polynomial.Basis(Basis)

		var nbCts uint32
		if inc1, err = buffer.ReadUint32(r, &nbCts); err != nil {
			return n + int64(inc1), err
		}

		n += int64(inc1)

		p.Value = make(map[int]*Ciphertext)

		for i := 0; i < int(nbCts); i++ {

			var key uint32

			if inc1, err = buffer.ReadUint32(r, &key); err != nil {
				return n + int64(inc1), err
			}

			n += int64(inc1)

			if p.Value[int(key)] == nil {
				p.Value[int(key)] = new(Ciphertext)
			}

			var inc2 int64
			if inc2, err = p.Value[int(key)].ReadFrom(r); err != nil {
				return n + inc2, err
			}

			n += inc2
		}

		return

	default:
		return p.ReadFrom(bufio.NewReader(r))
	}
}

func (p *PowerBasis) Write(data []byte) (n int, err error) {

	p.Basis = polynomial.Basis(data[n])
	n++

	nbCts := int(binary.LittleEndian.Uint32(data[n:]))
	n += 4

	p.Value = make(map[int]*Ciphertext)

	for i := 0; i < nbCts; i++ {

		idx := int(binary.LittleEndian.Uint32(data[n:]))
		n += 4

		if p.Value[idx] == nil {
			p.Value[idx] = new(Ciphertext)
		}

		var inc int
		if inc, err = p.Value[idx].Write(data[n:]); err != nil {
			return n + inc, err
		}

		n += inc
	}

	return
}
