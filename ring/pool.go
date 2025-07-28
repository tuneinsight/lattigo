package ring

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/utils/structs"
)

// BufferPool represents a pool of polys that can be used (concurrently) to instantiate temporary polynomials.
type BufferPool struct {
	level      int
	bufferPool structs.BufferPool[*[]uint64]
}

// NewPool returns a new pool given a ring, and optionally a pool to draw the backing arrays from.
func NewPool(ring *Ring, pools ...structs.BufferPool[*[]uint64]) *BufferPool {
	newPool := &BufferPool{level: ring.level}
	switch lenPool := len(pools); lenPool {
	case 0:
		newPool.bufferPool = structs.NewSyncPoolUint64(ring.N())
	case 1:
		newPool.bufferPool = pools[0]
	default:
		panic(fmt.Errorf("the method takes at most 1 argument but %d were given", lenPool))
	}

	return newPool
}

// GetLevel returns the level of the polynomials obtained from the pool.
func (p BufferPool) GetLevel() int {
	return p.level
}

// AtLevel returns a new pool from which polynomials at the given level can be drawn.
func (p BufferPool) AtLevel(level int) *BufferPool {
	return &BufferPool{level, p.bufferPool}
}

// GetBuffUintArray returns a new  []uint64 slice obtained from a pool.
// After use, the slice should be recycled using the [BufferPool.RecycleBuffUintArray] method.
func (p BufferPool) GetBuffUintArray() *[]uint64 {
	return p.bufferPool.Get()
}

// RecycleBuffUintArray takes a reference to a []uint64 slice and puts it back in the pool.
func (p BufferPool) RecycleBuffUintArray(arr *[]uint64) {
	p.bufferPool.Put(arr)
}

// GetBuffPoly returns a new [Poly], built from backing []uint64 arrays obtained from a pool.
// After use, the [Poly] should be recycled using the [BufferPool.RecycleBuffPoly] method.
func (p BufferPool) GetBuffPoly() *Poly {
	return NewPolyFromUintPool(p.bufferPool, p.level)
}

// RecycleBuffPoly takes a reference to a [Poly] and recycles its backing []uint64 arrays
// (i.e. they are returned to a pool). The input [Poly] must not be used after calling this method.
func (p BufferPool) RecycleBuffPoly(pol *Poly) {
	RecyclePolyInUintPool(p.bufferPool, pol)
}
