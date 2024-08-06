package structs

import (
	"bufio"
	"fmt"
	"io"

	"github.com/tuneinsight/lattigo/v6/utils/buffer"
)

// Vector is a struct wrapping a slice of components of type T.
// T can be:
//   - uint, uint64, uint32, uint16, uint8/byte, int, int64, int32, int16, int8, float64, float32.
//   - Or any object that implements CopyNewer, CopyNewer, io.WriterTo or io.ReaderFrom depending on
//     the method called.
type Vector[T any] []T

// CopyNew returns a deep copy of the object.
// If T is a struct, this method requires that T implements CopyNewer.
func (v Vector[T]) CopyNew() (vcpy Vector[T]) {

	var t T
	switch any(t).(type) {
	case uint, uint64, uint32, uint16, uint8, int, int64, int32, int16, int8, float64, float32:
		vcpy = Vector[T](make([]T, len(v)))
		copy(vcpy, v)
	default:
		if _, isCopiable := any(t).(CopyNewer[T]); !isCopiable {
			panic(fmt.Errorf("vector component of type %T does not comply to %T", t, new(CopyNewer[T])))
		}

		vcpy = Vector[T](make([]T, len(v)))
		for i := range v {
			vcpy[i] = *any(&v[i]).(CopyNewer[T]).CopyNew()
		}
	}

	return
}

// BinarySize returns the serialized size of the object in bytes.
// If T is a struct, this method requires that T implements BinarySizer.
func (v Vector[T]) BinarySize() (size int) {

	var t T
	switch any(t).(type) {
	case uint, uint64, int, int64, float64:
		return 8 + len(v)*8
	case uint32, int32, float32:
		return 8 + len(v)*4
	case uint16, int16:
		return 8 + len(v)*2
	case uint8, int8:
		return 8 + len(v)*1
	default:
		if _, isSizable := any(t).(BinarySizer); !isSizable {
			panic(fmt.Errorf("vector component of type %T does not comply to %T", t, new(BinarySizer)))
		}

		size += 8
		for i := range v {
			size += any(&v[i]).(BinarySizer).BinarySize()
		}
	}

	return
}

// WriteTo writes the object on an io.Writer. It implements the io.WriterTo
// interface, and will write exactly object.BinarySize() bytes on w.
//
// If T is a struct, this method requires that T implements io.WriterTo.
//
// Unless w implements the buffer.Writer interface (see lattigo/utils/buffer/writer.go),
// it will be wrapped into a bufio.Writer. Since this requires allocations, it
// is preferable to pass a buffer.Writer directly:
//
//   - When writing multiple times to a io.Writer, it is preferable to first wrap the
//     io.Writer in a pre-allocated bufio.Writer.
//   - When writing to a pre-allocated var b []byte, it is preferable to pass
//     buffer.NewBuffer(b) as w (see lattigo/utils/buffer/buffer.go).
func (v Vector[T]) WriteTo(w io.Writer) (n int64, err error) {

	switch w := w.(type) {
	case buffer.Writer:

		var inc int64
		if inc, err = buffer.WriteAsUint64[int](w, len(v)); err != nil {
			return inc, fmt.Errorf("buffer.WriteAsUint64[int]: %w", err)
		}

		n += inc

		var t T
		switch t := any(t).(type) {
		case uint, uint64, int, int64, float64:

			if inc, err = buffer.WriteAsUint64Slice[T](w, v); err != nil {
				return n + inc, fmt.Errorf("buffer.WriteAsUint64Slice[%T]: %w", t, err)
			}

			n += inc

		case uint32, int32, float32:

			if inc, err = buffer.WriteAsUint32Slice[T](w, v); err != nil {
				return n + inc, fmt.Errorf("buffer.WriteAsUint32Slice[%T]: %w", t, err)
			}

			n += inc

		case uint16, int16:

			if inc, err = buffer.WriteAsUint16Slice[T](w, v); err != nil {
				return n + inc, fmt.Errorf("buffer.WriteAsUint16Slice[%T]: %w", t, err)
			}

			n += inc

		case uint8, int8:

			if inc, err = buffer.WriteAsUint8Slice[T](w, v); err != nil {
				return n + inc, fmt.Errorf("buffer.WriteAsUint8Slice[%T]: %w", t, err)
			}

			n += inc

		default:

			if _, isWritable := any(new(T)).(io.WriterTo); !isWritable {
				return 0, fmt.Errorf("vector component of type %T does not comply to %T", t, new(io.WriterTo))
			}

			for i := range v {
				if inc, err = any(&v[i]).(io.WriterTo).WriteTo(w); err != nil {
					return n + inc, fmt.Errorf("%T.WriteTo: %w", t, err)
				}
				n += inc
			}
		}

		return n, w.Flush()

	default:
		return v.WriteTo(bufio.NewWriter(w))
	}
}

// ReadFrom reads on the object from an io.Writer. It implements the
// io.ReaderFrom interface.
//
// If T is a struct, this method requires that T implements io.ReaderFrom.
//
// Unless r implements the buffer.Reader interface (see lattigo/utils/buffer/reader.go),
// it will be wrapped into a bufio.Reader. Since this requires allocation, it
// is preferable to pass a buffer.Reader directly:
//
//   - When reading multiple values from a io.Reader, it is preferable to first
//     first wrap io.Reader in a pre-allocated bufio.Reader.
//   - When reading from a var b []byte, it is preferable to pass a buffer.NewBuffer(b)
//     as w (see lattigo/utils/buffer/buffer.go).
func (v *Vector[T]) ReadFrom(r io.Reader) (n int64, err error) {

	switch r := r.(type) {
	case buffer.Reader:

		var inc int64

		var size int

		if inc, err = buffer.ReadAsUint64[int](r, &size); err != nil {
			return inc, fmt.Errorf("buffer.ReadAsUint64[int]: %w", err)
		}

		n += inc

		if cap(*v) < size {
			*v = make([]T, size)
		}

		*v = (*v)[:size]

		var t T
		switch any(t).(type) {
		case uint, uint64, int, int64, float64:

			if inc, err = buffer.ReadAsUint64Slice[T](r, *v); err != nil {
				return n + inc, fmt.Errorf("buffer.ReadAsUint64Slice[%T]: %w", t, err)
			}

			n += inc

		case uint32, int32, float32:

			if inc, err = buffer.ReadAsUint32Slice[T](r, *v); err != nil {
				return n + inc, fmt.Errorf("buffer.ReadAsUint32Slice[%T]: %w", t, err)
			}

			n += inc

		case uint16, int16:

			if inc, err = buffer.ReadAsUint16Slice[T](r, *v); err != nil {
				return n + inc, fmt.Errorf("buffer.ReadAsUint16Slice[%T]: %w", t, err)
			}

			n += inc

		case uint8, int8:

			if inc, err = buffer.ReadAsUint8Slice[T](r, *v); err != nil {
				return n + inc, fmt.Errorf("buffer.ReadAsUint8Slice[%T]: %w", t, err)
			}

			n += inc
		default:

			if _, isReadable := any(new(T)).(io.ReaderFrom); !isReadable {
				return 0, fmt.Errorf("vector component of type %T does not comply to %T", t, new(io.ReaderFrom))
			}

			for i := range *v {
				if inc, err = any(&(*v)[i]).(io.ReaderFrom).ReadFrom(r); err != nil {
					var t T
					return n + inc, fmt.Errorf("%T.ReadFrom: %w", t, err)
				}
				n += inc
			}
		}

		return n, nil

	default:
		return v.ReadFrom(bufio.NewReader(r))
	}
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
// If T is a struct, this method requires that T implements io.WriterTo.
func (v Vector[T]) MarshalBinary() (p []byte, err error) {
	buf := buffer.NewBufferSize(v.BinarySize())
	_, err = v.WriteTo(buf)
	return buf.Bytes(), err
}

// UnmarshalBinary decodes a slice of bytes generated by
// MarshalBinary or WriteTo on the object.
// If T is a struct, this method requires that T implements io.ReaderFrom.
func (v *Vector[T]) UnmarshalBinary(p []byte) (err error) {
	_, err = v.ReadFrom(buffer.NewBuffer(p))
	return
}

// Equal performs a deep equal.
// If T is a struct, this method requires that T implements Equatable.
func (v Vector[T]) Equal(other Vector[T]) (isEqual bool) {

	var t T
	switch any(t).(type) {
	case uint, uint64, int, int64, float64:
		return buffer.EqualAsUint64Slice([]T(v), []T(other))
	case uint32, int32, float32:
		return buffer.EqualAsUint32Slice([]T(v), []T(other))
	case uint16, int16:
		return buffer.EqualAsUint16Slice([]T(v), []T(other))
	case uint8, int8:
		return buffer.EqualAsUint8Slice([]T(v), []T(other))
	default:

		if _, isEquatable := any(t).(Equatable[T]); !isEquatable {
			panic(fmt.Errorf("vector component of type %T does not comply to %T", t, new(Equatable[T])))
		}

		for i, v := range v {
			/* #nosec G601 -- Implicit memory aliasing in for loop acknowledged */
			if !any(&v).(Equatable[T]).Equal(&other[i]) {
				return false
			}
		}
		return true
	}
}
