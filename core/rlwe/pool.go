package rlwe

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
	"github.com/tuneinsight/lattigo/v6/utils/structs"
)

// Pool represents a pool of different objects (plaintexts, ciphertexts, polys) that can be used to instantiate temporary buffers.
type Pool struct {
	*ringqp.Pool
}

// NewPool returns a new pool given a RingQP, and optionally a pool to draw the backing arrays from.
func NewPool(rqp *ringqp.Ring, pools ...structs.BufferPool[*[]uint64]) *Pool {
	// If no backing pool is given, we create one here.
	switch lenPool := len(pools); lenPool {
	case 0:
		pools = append(pools, structs.NewSyncPoolUint64(rqp.N()))
	case 1:
	default:
		panic(fmt.Errorf("the method takes at most 1 argument but %d were given", lenPool))
	}

	ringqpPool := ringqp.NewPool(rqp, pools...)

	return &Pool{ringqpPool}
}

// GetBuffCt returns a ciphertext that can be used as a buffer for intermediate computations.
// After use, the ciphertext should be recycled with [Pool.RecycleBuffCt].
// The optional dimensions specify the degree and level of the ciphertext (default to 2, eval.params.ringQ.Level()).
func (pool *Pool) GetBuffCt(dimensions ...int) *Ciphertext {
	degree := 2
	level := pool.GetLevel()
	switch nbParams := len(dimensions); nbParams {
	case 0:
	case 1:
		degree = dimensions[0]
	case 2:
		degree = dimensions[0]
		level = dimensions[1]
	default:
		panic(fmt.Errorf("getbuffct takes 2 parameters at most"))
	}

	poolQ := pool.AtLevel(level, 0)
	polys := make([]ring.Poly, degree+1)
	for i := range polys {
		polys[i] = *poolQ.GetBuffPoly()
	}

	ct, err := NewCiphertextAtLevelFromPoly(level, polys)

	// sanity check: should not happen
	if err != nil {
		panic(fmt.Errorf("cannot create new ciphertext: %w", err))
	}
	return ct
}

// RecycleBuffCt recycles a temporary ciphertext (i.e. returns its backing uint64 arrays to the pool).
// The input ciphertext must not be used after calling this method.
func (pool *Pool) RecycleBuffCt(ct *Ciphertext) {
	for i := range ct.Value {
		pool.RecycleBuffPoly(&ct.Value[i])
	}
	ct = nil
}

// GetBuffPt returns a plaintext that can be used as a buffer for intermediate computations.
// After use, the plaintext should be recycled with [Pool.RecycleBuffPt].
// The optional argument specifies the level of the returned plaintext (default to eval.params.ringQ.Level()).
func (pool *Pool) GetBuffPt(level ...int) *Plaintext {
	lvl := pool.GetLevel()
	switch nbParams := len(level); nbParams {
	case 0:
	case 1:
		lvl = level[0]
	default:
		panic(fmt.Errorf("getbuffpt takes 1 parameter at most but %d were given", nbParams))
	}

	poly := pool.AtLevel(lvl, 0).GetBuffPoly()

	pt, err := NewPlaintextAtLevelFromPoly(lvl, *poly)

	// sanity check: should not happen
	if err != nil {
		panic(fmt.Errorf("cannot create new plaintext: %w", err))
	}
	return pt
}

// RecycleBuffPt recycles a temporary plaintext (i.e. returns its backing uint64 arrays to the pool).
// The input plaintext must not be used after calling this method.
func (pool *Pool) RecycleBuffPt(pt *Plaintext) {
	pool.RecycleBuffPoly(&pt.Value)
}
