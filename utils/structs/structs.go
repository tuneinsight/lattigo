// Package structs implements helpers to generalize vectors and matrices of structs, as well as their serialization.
package structs

type Equatable[T any] interface {
	Equal(*T) bool
}

type CopyNewer[V any] interface {
	CopyNew() *V
}

type BinarySizer interface {
	BinarySize() int
}
