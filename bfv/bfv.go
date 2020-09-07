// Package bfv implements a RNS-accelerated Fan-Vercauteren version of Brakerski's scale invariant homomorphic encryption scheme. It provides modular arithmetic over the integers.
package bfv

import (
	"github.com/ldsec/lattigo/ring"
)

// GaloisGen is an integer of order N/2 modulo M and that spans Z_M with the integer -1. The j-th ring automorphism takes the root zeta to zeta^(5j).
const GaloisGen uint64 = 5

// bfvContext is a struct which contains all the elements required to instantiate the BFV Scheme. This includes the parameters (polynomial degree, plaintext modulus, ciphertext modulus,
// polynomial contexts and other parameters required for the homomorphic operations).
type bfvContext struct {
	params *Parameters

	// Polynomial degree
	n uint64

	// Polynomial contexts
	ringT    *ring.Ring // Plaintext modulus
	ringQ    *ring.Ring // Ciphertext modulus
	ringQMul *ring.Ring // Ciphertext extended modulus (for multiplication)
	ringP    *ring.Ring // Keys additional modulus
	ringQP   *ring.Ring // Keys modulus

	galElRotRow      uint64   // Rows rotation generator
	galElRotColLeft  []uint64 // Columns right rotations generators
	galElRotColRight []uint64 // Columsn left rotations generators
}

func newBFVContext(params *Parameters) (context *bfvContext) {

	context = new(bfvContext)
	var err error

	N := uint64(1 << params.logN)

	context.n = N

	if context.ringT, err = ring.NewRing(N, []uint64{params.t}); err != nil {
		panic(err)
	}

	if context.ringQ, err = ring.NewRing(N, params.qi); err != nil {
		panic(err)
	}

	if context.ringQMul, err = ring.NewRing(N, params.qiMul); err != nil {
		panic(err)
	}

	if len(params.pi) != 0 {
		if context.ringP, err = ring.NewRing(N, params.pi); err != nil {
			panic(err)
		}
	}

	if context.ringQP, err = ring.NewRing(N, append(params.qi, params.pi...)); err != nil {
		panic(err)
	}

	context.galElRotColLeft = ring.GenGaloisParams(N, GaloisGen)
	context.galElRotColRight = ring.GenGaloisParams(N, ring.ModExp(GaloisGen, 2*N-1, 2*N))
	context.galElRotRow = 2*N - 1
	return
}
