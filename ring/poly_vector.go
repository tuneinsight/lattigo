package ring

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"io"

	"github.com/tuneinsight/lattigo/v4/utils/buffer"
)

type PolyVector []*Poly

func NewPolyVector(N, Level, size int) *PolyVector {
	v := make([]*Poly, size)

	for i := range v {
		v[i] = NewPoly(N, Level)
	}

	pv := PolyVector(v)

	return &pv
}

func (pv *PolyVector) Set(polys []*Poly) {
	*pv = PolyVector(polys)
}

func (pv *PolyVector) Get() []*Poly {
	return []*Poly(*pv)
}

func (pv *PolyVector) N() int {
	return (*pv)[0].N()
}

func (pv *PolyVector) Level() int {
	return (*pv)[0].Level()
}

func (pv *PolyVector) Resize(level, size int) {
	N := pv.N()

	v := *pv

	for i := range v {
		v[i].Resize(level)
	}

	if len(v) > level {
		v = v[:level+1]
	} else {
		for i := len(v); i < level+1; i++ {
			v = append(v, NewPoly(N, level))
		}
	}

	*pv = v
}

func (pv *PolyVector) BinarySize() (size int) {
	size += 8
	for _, v := range *pv {
		size += v.BinarySize()
	}
	return
}

func (pv *PolyVector) MarshalBinary() (p []byte, err error) {
	p = make([]byte, pv.BinarySize())
	_, err = pv.Read(p)
	return
}

func (pv *PolyVector) Read(b []byte) (n int, err error) {

	v := *pv

	binary.LittleEndian.PutUint64(b[n:], uint64(len(v)))
	n += 8

	var inc int
	for i := range v {
		if inc, err = v[i].Read(b[n:]); err != nil {
			return n + inc, err
		}

		n += inc
	}

	return
}

func (pv *PolyVector) WriteTo(w io.Writer) (int64, error) {
	switch w := w.(type) {
	case buffer.Writer:

		var err error
		var n int64

		v := *pv

		var inc int
		if inc, err = buffer.WriteInt(w, len(v)); err != nil {
			return int64(inc), err
		}

		n += int64(inc)

		for i := range v {
			var inc int64
			if inc, err = v[i].WriteTo(w); err != nil {
				return n + inc, err
			}

			n += inc
		}

		return n, nil

	default:
		return pv.WriteTo(bufio.NewWriter(w))
	}
}

func (pv *PolyVector) UnmarhsalBinary(p []byte) (err error) {
	_, err = pv.Write(p)
	return
}

func (pv *PolyVector) Write(p []byte) (n int, err error) {
	size := int(binary.LittleEndian.Uint64(p[n:]))
	n += 8

	if len(*pv) != size {
		*pv = make([]*Poly, size)
	}

	v := *pv

	var inc int
	for i := range v {
		if v[i] == nil {
			v[i] = new(Poly)
		}

		if inc, err = v[i].Write(p[n:]); err != nil {
			return n + inc, err
		}

		n += inc
	}

	return
}

func (pv *PolyVector) ReadFrom(r io.Reader) (int64, error) {
	switch r := r.(type) {
	case buffer.Reader:

		var err error
		var size, n int

		if n, err = buffer.ReadInt(r, &size); err != nil {
			return int64(n), fmt.Errorf("cannot ReadFrom: size: %w", err)
		}

		if len(*pv) != size {
			*pv = make([]*Poly, size)
		}

		v := *pv

		for i := range v {

			if v[i] == nil {
				v[i] = new(Poly)
			}

			var inc int64
			if inc, err = v[i].ReadFrom(r); err != nil {
				return int64(n) + inc, err
			}

			n += int(inc)
		}

		return int64(n), nil

	default:
		return pv.ReadFrom(bufio.NewReader(r))
	}
}
