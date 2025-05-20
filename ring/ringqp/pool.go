package ringqp

import (
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils/structs"
)

type Pool struct {
	poolQ *ring.Pool
	poolP *ring.Pool
}

func NewPool(ringqp *Ring, pools ...structs.BufferPool[*[]uint64]) (*Pool, error) {
	var poolQ, poolP *ring.Pool
	var err error

	if ringqp.RingQ != nil {
		poolQ, err = ring.NewPool(ringqp.RingQ, pools...)
		if err != nil {
			return nil, err
		}
	}
	if ringqp.RingP != nil {
		poolP, err = ring.NewPool(ringqp.RingP, pools...)
		if err != nil {
			return nil, err
		}
	}

	return &Pool{poolQ, poolP}, nil
}

func (p Pool) AtLevel(levelQ, levelP int) *Pool {
	var poolQ, poolP *ring.Pool
	if p.poolQ != nil {
		poolQ = p.poolQ.AtLevel(levelQ)
	}
	if p.poolP != nil {
		poolP = p.poolP.AtLevel(levelQ)
	}
	return &Pool{poolQ, poolP}
}

// GetBuffPolyQP returns a new [Poly], built from backing []uint64 arrays obtained from a pool.
// After use, the [Poly] should be recycled using the [Ring.RecycleBuffPolyQP] method.
func (p Pool) GetBuffPolyQP() *Poly {
	var Q, P ring.Poly
	if p.poolQ != nil {
		buffQ := p.poolQ.GetBuffPoly()
		Q = *buffQ
	}
	if p.poolP != nil {
		buffP := p.poolP.GetBuffPoly()
		P = *buffP
	}
	return &Poly{Q, P}
}

// RecycleBuffPolyQP takes a reference to a [Poly] and recycles its backing []uint64 arrays
// (i.e. they are returned to a pool). The input [Poly] must not be used after calling this method.
func (p Pool) RecycleBuffPolyQP(poly *Poly) {
	if p.poolQ != nil {
		p.poolQ.RecycleBuffPoly(&poly.Q)
	}
	if p.poolP != nil {
		p.poolP.RecycleBuffPoly(&poly.P)
	}
}
