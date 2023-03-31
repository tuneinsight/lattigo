package ring

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"io"

	"github.com/tuneinsight/lattigo/v4/utils/buffer"
)

type PolyMatrix []*PolyVector

func NewPolyMatrix(N, Level, rows, cols int) *PolyMatrix {
	m := make([]*PolyVector, rows)

	for i := range m {
		m[i] = NewPolyVector(N, Level, cols)
	}

	pm := PolyMatrix(m)

	return &pm
}

func (pm *PolyMatrix) Set(polys [][]*Poly) {

	m := PolyMatrix(make([]*PolyVector, len(polys)))
	for i := range m {
		m[i] = new(PolyVector)
		m[i].Set(polys[i])
	}

	*pm = m
}

func (pm *PolyMatrix) Get() [][]*Poly {
	m := *pm
	polys := make([][]*Poly, len(m))
	for i := range polys {
		polys[i] = m[i].Get()
	}
	return polys
}

func (pm *PolyMatrix) BinarySize() (size int) {
	size += 8
	for _, m := range *pm {
		size += m.BinarySize()
	}
	return
}

func (pm *PolyMatrix) MarshalBinary() (p []byte, err error) {
	p = make([]byte, pm.BinarySize())
	_, err = pm.Read(p)
	return
}

func (pm *PolyMatrix) Read(b []byte) (n int, err error) {

	m := *pm

	binary.LittleEndian.PutUint64(b[n:], uint64(len(m)))
	n += 8

	var inc int
	for i := range m {
		if inc, err = m[i].Read(b[n:]); err != nil {
			return n + inc, err
		}

		n += inc
	}

	return
}

func (pm *PolyMatrix) WriteTo(w io.Writer) (int64, error) {
	switch w := w.(type) {
	case buffer.Writer:

		var err error
		var n int64

		m := *pm

		var inc int
		if inc, err = buffer.WriteInt(w, len(m)); err != nil {
			return int64(inc), err
		}

		n += int64(inc)

		for i := range m {
			var inc int64
			if inc, err = m[i].WriteTo(w); err != nil {
				return n + inc, err
			}

			n += inc
		}

		return n, nil

	default:
		return pm.WriteTo(bufio.NewWriter(w))
	}
}

func (pm *PolyMatrix) UnmarhsalBinary(p []byte) (err error) {
	_, err = pm.Write(p)
	return
}

func (pm *PolyMatrix) Write(p []byte) (n int, err error) {
	size := int(binary.LittleEndian.Uint64(p[n:]))
	n += 8

	if len(*pm) != size {
		*pm = make([]*PolyVector, size)
	}

	m := *pm

	var inc int
	for i := range m {
		if m[i] == nil {
			m[i] = new(PolyVector)
		}

		if inc, err = m[i].Write(p[n:]); err != nil {
			return n + inc, err
		}

		n += inc
	}

	return
}

func (pm *PolyMatrix) ReadFrom(r io.Reader) (int64, error) {
	switch r := r.(type) {
	case buffer.Reader:

		var err error
		var size, n int

		if n, err = buffer.ReadInt(r, &size); err != nil {
			return int64(n), fmt.Errorf("cannot ReadFrom: size: %w", err)
		}

		if len(*pm) != size {
			*pm = make([]*PolyVector, size)
		}

		m := *pm

		for i := range m {

			if m[i] == nil {
				m[i] = new(PolyVector)
			}

			var inc int64
			if inc, err = m[i].ReadFrom(r); err != nil {
				return int64(n) + inc, err
			}

			n += int(inc)
		}

		return int64(n), nil

	default:
		return pm.ReadFrom(bufio.NewReader(r))
	}
}
