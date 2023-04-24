package structs

import (
	"encoding"
	"fmt"
	"io"
)

type Codec[V any] struct{}

type CopyNewer[V any] interface {
	CopyNew() *V
}

func (c *Codec[V]) CopynewWrapper(T interface{}) (*V, error) {

	copyer, ok := T.(CopyNewer[V])

	if !ok {
		return nil, fmt.Errorf("cannot CopyNew: type T=%T does not implement CopyNew", T)
	}

	return copyer.CopyNew(), nil
}

type BinarySizer interface {
	BinarySize() int
}

func (c *Codec[V]) BinarySizeWrapper(T interface{}) (size int, err error) {
	binarysizer, ok := T.(BinarySizer)

	if !ok {
		return 0, fmt.Errorf("cannot MarshalBinary: type T=%T does not implement BinarySizer", T)
	}

	return binarysizer.BinarySize(), nil
}

func (c *Codec[V]) MarshalBinaryWrapper(T interface{}) (p []byte, err error) {
	binarymarshaler, ok := T.(encoding.BinaryMarshaler)

	if !ok {
		return nil, fmt.Errorf("cannot MarshalBinary: type T=%T does not implement encoding.BinaryMarshaler", T)
	}

	return binarymarshaler.MarshalBinary()
}

func (c *Codec[V]) UnmarshalBinaryWrapper(p []byte, T interface{}) (err error) {
	binaryunmarshaler, ok := T.(encoding.BinaryUnmarshaler)

	if !ok {
		return fmt.Errorf("cannot UnmarshalBinary: type T=%T does not implement encoding.UnmarshalBinary", T)
	}

	return binaryunmarshaler.UnmarshalBinary(p)
}

type Encoder interface {
	Encode(p []byte) (n int, err error)
}

func (c *Codec[V]) EncodeWrapper(p []byte, T interface{}) (n int, err error) {
	encoder, ok := T.(Encoder)

	if !ok {
		return 0, fmt.Errorf("cannot Encode: type T=%T does not implement Encoder", T)
	}

	return encoder.Encode(p)
}

type Decoder interface {
	Decode(p []byte) (n int, err error)
}

func (c *Codec[V]) DecodeWrapper(p []byte, T interface{}) (n int, err error) {
	decoder, ok := T.(Decoder)

	if !ok {
		return 0, fmt.Errorf("cannot Decode: type T=%T does not implement Decoder", T)
	}

	return decoder.Decode(p)
}

func (c *Codec[V]) WriteToWrapper(w io.Writer, T interface{}) (n int64, err error) {
	writerto, ok := T.(io.WriterTo)

	if !ok {
		return 0, fmt.Errorf("cannot WriteTo: type T=%T does not implement io.WriterTo", T)
	}

	return writerto.WriteTo(w)
}

func (c *Codec[V]) ReadFromWrapper(r io.Reader, T interface{}) (n int64, err error) {
	readerfrom, ok := T.(io.ReaderFrom)

	if !ok {
		return 0, fmt.Errorf("cannot ReadFrom: type T=%T does not implement io.ReaderFrom", T)
	}

	return readerfrom.ReadFrom(r)
}
