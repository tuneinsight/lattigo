package structs

import (
	"encoding"
	"fmt"
	"io"
)

type BinarySizer interface {
	BinarySize() int
}

type Codec struct{}

var codec = Codec{}

func (c *Codec) BinarySizeWrapper(T interface{}) (size int, err error) {
	binarysizer, ok := T.(BinarySizer)

	if !ok {
		return 0, fmt.Errorf("cannot MarshalBinary: type T=%T does not implement BinarySizer", T)
	}

	return binarysizer.BinarySize(), nil
}

func (c *Codec) MarshalBinaryWrapper(T interface{}) (p []byte, err error) {
	binarymarshaler, ok := T.(encoding.BinaryMarshaler)

	if !ok {
		return nil, fmt.Errorf("cannot MarshalBinary: type T=%T does not implement encoding.BinaryMarshaler", T)
	}

	return binarymarshaler.MarshalBinary()
}

func (c *Codec) UnmarshalBinaryWrapper(p []byte, T interface{}) (err error) {
	binaryunmarshaler, ok := T.(encoding.BinaryUnmarshaler)

	if !ok {
		return fmt.Errorf("cannot UnmarshalBinary: type T=%T does not implement encoding.UnmarshalBinary", T)
	}

	return binaryunmarshaler.UnmarshalBinary(p)
}

func (c *Codec) ReadWrapper(p []byte, T interface{}) (n int, err error) {
	reader, ok := T.(io.Reader)

	if !ok {
		return 0, fmt.Errorf("cannot Read: type T=%T does not implement io.Reader", T)
	}

	return reader.Read(p)
}

func (c *Codec) WriteWrapper(p []byte, T interface{}) (n int, err error) {
	writer, ok := T.(io.Writer)

	if !ok {
		return 0, fmt.Errorf("cannot Read: type T=%T does not implement io.Writer", T)
	}

	return writer.Write(p)
}

func (c *Codec) WriteToWrapper(w io.Writer, T interface{}) (n int64, err error) {
	writerto, ok := T.(io.WriterTo)

	if !ok {
		return 0, fmt.Errorf("cannot Read: type T=%T does not implement io.WriterTo", T)
	}

	return writerto.WriteTo(w)
}

func (c *Codec) ReadFromWrapper(r io.Reader, T interface{}) (n int64, err error) {
	readerfrom, ok := T.(io.ReaderFrom)

	if !ok {
		return 0, fmt.Errorf("cannot Read: type T=%T does not implement io.ReaderFrom", T)
	}

	return readerfrom.ReadFrom(r)
}
