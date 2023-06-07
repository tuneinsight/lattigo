// Package structs implements helpers to generalize vectors and matrices of structs, as well as their serialization.
package structs

type CopyNewer[V any] interface {
	CopyNew() *V
}

type BinarySizer interface {
	BinarySize() int
}

type Encoder interface {
	Encode(p []byte) (n int, err error)
}

type Decoder interface {
	Decode(p []byte) (n int, err error)
}
