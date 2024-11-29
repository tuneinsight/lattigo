package structs

import "sync"

type BufferPool[T any] interface {
	Get() T
	Put(T)
}
type SyncPool[T any] struct {
	pool *sync.Pool
}

func NewSyncPool[T any](f func() T) *SyncPool[T] {
	pool := &sync.Pool{
		New: func() any {
			return f()
		},
	}
	return &SyncPool[T]{pool: pool}
}

func (spool *SyncPool[T]) Get() T {
	return spool.pool.Get().(T)
}

func (spool *SyncPool[T]) Put(buff T) {
	spool.pool.Put(buff)
}
