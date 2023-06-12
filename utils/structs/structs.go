// Package structs implements helpers to generalize vectors and matrices of structs, as well as their serialization.
package structs

import (
	"encoding"
	"io"
)

type CopyNewer[V any] interface {
	CopyNew() *V
}

type BinarySizer interface {
	BinarySize() int
}

// BinarySerializer is a testing interface for byte encoding and decoding.
type BinarySerializer interface {
	io.WriterTo
	io.ReaderFrom
	encoding.BinaryMarshaler
	encoding.BinaryUnmarshaler
}
