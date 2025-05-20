package ring

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/utils/structs"
)

type Pool struct {
	level      int
	bufferPool structs.BufferPool[*[]uint64]
}

func NewPool(ring *Ring, pools ...structs.BufferPool[*[]uint64]) (*Pool, error) {
	newPool := &Pool{level: ring.level}
	switch lenPool := len(pools); lenPool {
	case 0:
		newPool.bufferPool = ring.bufferPool
	case 1:
		newPool.bufferPool = pools[0]
	default:
		return nil, fmt.Errorf("the method takes at most 2 arguments but %d were given", lenPool+1)
	}
	return newPool, nil
}

func (p Pool) AtLevel(level int) *Pool {
	return &Pool{level, p.bufferPool}
}

// GetBuffPoly returns a new [Poly], built from backing []uint64 arrays obtained from a pool.
// After use, the [Poly] should be recycled using the [Ring.RecycleBuffPoly] method.
func (p Pool) GetBuffPoly() *Poly {
	return NewPolyFromUintPool(p.bufferPool, p.level)
}

// RecycleBuffPoly takes a reference to a [Poly] and recycles its backing []uint64 arrays
// (i.e. they are returned to a pool). The input [Poly] must not be used after calling this method.
func (p Pool) RecycleBuffPoly(pol *Poly) {
	RecyclePolyInUintPool(p.bufferPool, pol)
}
