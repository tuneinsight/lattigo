// Package ringqp is implements a wrapper for both the ringQ and ringP.
package ringqp

import (
	"github.com/tuneinsight/lattigo/v4/ring"
)

// Ring is a structure that implements the operation in the ring R_QP.
// This type is simply a union type between the two Ring types representing
// R_Q and R_P.
type Ring struct {
	RingQ, RingP *ring.Ring
}

// AtLevel returns a shallow copy of the target ring configured to
// carry on operations at the specified levels.
func (r Ring) AtLevel(levelQ, levelP int) Ring {

	var ringQ, ringP *ring.Ring

	if levelQ > -1 && r.RingQ != nil {
		ringQ = r.RingQ.AtLevel(levelQ)
	}

	if levelP > -1 && r.RingP != nil {
		ringP = r.RingP.AtLevel(levelP)
	}

	return Ring{
		RingQ: ringQ,
		RingP: ringP,
	}
}

// LevelQ returns the level at which the target
// ring operates for the modulus Q.
func (r Ring) LevelQ() int {
	if r.RingQ != nil {
		return r.RingQ.Level()
	}

	return -1
}

// LevelP returns the level at which the target
// ring operates for the modulus P.
func (r Ring) LevelP() int {
	if r.RingP != nil {
		return r.RingP.Level()
	}

	return -1
}

func (r Ring) Equal(p1, p2 Poly) (v bool) {
	v = true
	if r.RingQ != nil {
		v = v && r.RingQ.Equal(p1.Q, p2.Q)
	}

	if r.RingP != nil {
		v = v && r.RingP.Equal(p1.P, p2.P)
	}

	return
}

// NewPoly creates a new polynomial with all coefficients set to 0.
func (r Ring) NewPoly() Poly {
	var Q, P ring.Poly
	if r.RingQ != nil {
		Q = r.RingQ.NewPoly()
	}

	if r.RingP != nil {
		P = r.RingP.NewPoly()
	}
	return Poly{Q, P}
}
