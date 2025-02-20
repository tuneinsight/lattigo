package structs

import "sync"

// BufferPool is an interface for all pools of buffers.
type BufferPool[T any] interface {
	Get() T
	Put(T)
}

// SyncPool is a wrapper around [sync.Pool] (it avoids doing type conversion after Get()).
type SyncPool[T any] struct {
	pool *sync.Pool
}

// NewSyncPool creates a new SyncPool.
// The input function f is the function that is used to create new objects if none is available in the pool.
func NewSyncPool[T any](f func() T) *SyncPool[T] {
	pool := &sync.Pool{
		New: func() any {
			return f()
		},
	}
	return &SyncPool[T]{pool: pool}
}

// Get returns a new object of type T from the pool.
func (spool *SyncPool[T]) Get() T {
	return spool.pool.Get().(T)
}

// Put returns the buff to the pool.
func (spool *SyncPool[T]) Put(buff T) {
	spool.pool.Put(buff)
}

// BuffFromUintPool represents a pool of objects built on []uint64 backing arrays.
// It implements the [BufferPool] interface.
type BuffFromUintPool[T any] struct {
	createObject  func() T
	recycleObject func(T)
}

// NewBuffFromUintPool returns a new BuffFromUintPool structure.
// The create (resp. recycle) function are meant to use an underlying
// pool of []uint64 to build (resp. recycle) an object of type T.
func NewBuffFromUintPool[T any](create func() T, recycle func(T)) *BuffFromUintPool[T] {
	return &BuffFromUintPool[T]{
		createObject:  create,
		recycleObject: recycle,
	}
}

// Get returns a new object of type T built from a []uint64 backing array obtained from a pool.
func (bu *BuffFromUintPool[T]) Get() T {
	return bu.createObject()
}

// Put recycle an object of type T. I.e. it returns the []uint64 backing arrays of obj to their pool.
func (bu *BuffFromUintPool[T]) Put(obj T) {
	bu.recycleObject(obj)
}
