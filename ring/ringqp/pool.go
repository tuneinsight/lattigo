package ringqp

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils/structs"
)

// BufferPool represents a pool of polys that can be used (concurrently) to instantiate temporary polynomials in RingQP.
type BufferPool struct {
	*ring.BufferPool
	PoolP *ring.BufferPool
}

// NewPool returns a new pool given a RingQP, and optionally a pool to draw the backing arrays from.
func NewPool(ringqp *Ring, pools ...structs.BufferPool[*[]uint64]) *BufferPool {
	// If no backing pool is given, we create one here.
	switch lenPool := len(pools); lenPool {
	case 0:
		pools = append(pools, structs.NewSyncPoolUint64(ringqp.N()))
	case 1:
	default:
		panic(fmt.Errorf("the method takes at most 1 argument but %d were given", lenPool))
	}

	var poolQ, poolP *ring.BufferPool

	if ringqp.RingQ != nil {
		poolQ = ring.NewPool(ringqp.RingQ, pools...)
	}
	if ringqp.RingP != nil {
		poolP = ring.NewPool(ringqp.RingP, pools...)
	}

	return &BufferPool{poolQ, poolP}
}

// AtLevel returns a new pool from which polynomials at the given levels can be drawn.
// The method accepts up to two arguments:
// Zero level: the polyomials are returned at level 0.
// One level: the polynomials in RingQ (resp. RingP) are returned at the given level (resp. level 0).
// Two levels: the polynomials in RingQ (resp. RingP) are returned at levels[0] (resp. levels[1]).
func (p BufferPool) AtLevel(levels ...int) *BufferPool {
	var levelQ, levelP int
	switch nbParams := len(levels); nbParams {
	case 0:
	case 1:
		levelQ = levels[0]
	case 2:
		levelQ = levels[0]
		levelP = levels[1]
	default:
		panic(fmt.Errorf("atlevel takes 2 parameters at most"))
	}

	var poolQ, poolP *ring.BufferPool
	if p.BufferPool != nil {
		poolQ = p.BufferPool.AtLevel(levelQ)
	}
	if p.PoolP != nil {
		poolP = p.PoolP.AtLevel(levelP)
	}
	return &BufferPool{poolQ, poolP}
}

// GetBuffPolyQP returns a new [Poly], built from backing []uint64 arrays obtained from a pool.
// After use, the [Poly] should be recycled using the [BufferPool.RecycleBuffPolyQP] method.
func (p BufferPool) GetBuffPolyQP() *Poly {
	var Q, P ring.Poly
	if p.BufferPool != nil {
		buffQ := p.GetBuffPoly()
		Q = *buffQ
	}
	if p.PoolP != nil {
		buffP := p.PoolP.GetBuffPoly()
		P = *buffP
	}
	return &Poly{Q, P}
}

// RecycleBuffPolyQP takes a reference to a [Poly] and recycles its backing []uint64 arrays
// (i.e. they are returned to a pool). The input [Poly] must not be used after calling this method.
func (p BufferPool) RecycleBuffPolyQP(poly *Poly) {
	if p.BufferPool != nil {
		p.RecycleBuffPoly(&poly.Q)
	}
	if p.PoolP != nil {
		p.PoolP.RecycleBuffPoly(&poly.P)
	}
}
