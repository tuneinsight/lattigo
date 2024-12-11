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

type BuffFromUintPool[T any] struct {
	uintPool      BufferPool[*[]uint64] // Pool that must contain *[]uint64 objects
	createObject  func(BufferPool[*[]uint64]) T
	recycleObject func(BufferPool[*[]uint64], T)
}

func NewBuffFromUintPool[T any](pool BufferPool[*[]uint64], create func(BufferPool[*[]uint64]) T, recycle func(BufferPool[*[]uint64], T)) *BuffFromUintPool[T] {
	return &BuffFromUintPool[T]{
		uintPool:      pool,
		createObject:  create,
		recycleObject: recycle,
	}
}

func (bu *BuffFromUintPool[T]) Get() T {
	return bu.createObject(bu.uintPool)
}

func (bu *BuffFromUintPool[T]) Put(obj T) {
	bu.recycleObject(bu.uintPool, obj)
}

type FreeList[T any] struct {
	pool      chan T
	newObject func() T
	capacity  int
}

func NewFreeList[T any](capacity int, f func() T) *FreeList[T] {
	pool := make(chan T, capacity)
	for i := 0; i < capacity; i++ {
		pool <- f()
	}
	return &FreeList[T]{
		pool:      pool,
		newObject: f,
		capacity:  capacity,
	}
}

func (fl *FreeList[T]) Get() T {
	var obj T

	select {
	case obj = <-fl.pool:
	default:
		obj = fl.newObject()
	}
	return obj
}

func (fl *FreeList[T]) Put(obj T) {
	select {
	case fl.pool <- obj:
	default:
	}
}
