//Package ckks implements a RNS-accelerated version of the Homomorphic Encryption for Arithmetic for Approximate Numbers
//(HEAAN, a.k.a. CKKS) scheme. It provides approximate arithmetic over the complex numbers.
package ckks

import (
	"github.com/ldsec/lattigo/ring"
	"math/big"
)

// GaloisGen is an integer of order N/2 modulo M and that spans Z_M with the integer -1. The j-th ring automorphism takes the root zeta to zeta^(5j).
// Any other integer or order N/2 modulo M and congruent with 1 modulo 4 could be used instead.
const GaloisGen uint64 = 5

// Context is a struct that contains all the elements required to instantiate the CKKS Scheme. This includes the parameters (polynomial degree, ciphertext modulus,
// polynomial contexts and other parameters required for the homomorphic operations).
type Context struct {

	// Context parameters
	logN     uint64
	scale    float64
	n        uint64
	maxSlots uint64

	// Number of available levels
	levels uint64

	bigintChain []*big.Int

	// Contexts
	contextQ  *ring.Context
	contextP  *ring.Context
	contextQP *ring.Context

	// Rotation params
	galElConjugate   uint64
	galElRotColLeft  []uint64
	galElRotColRight []uint64
}

// NewContext creates a new Context with the given parameters. It returns an error if one of the parameters would not ensure the
// correctness of the scheme (but it does not check for security).
func newContext(params *Parameters) (ckkscontext *Context) {

	var err error

	ckkscontext = new(Context)

	ckkscontext.logN = params.logN
	ckkscontext.n = params.N()
	ckkscontext.maxSlots = params.MaxSlots()
	ckkscontext.scale = params.scale

	ckkscontext.levels = params.QiCount()

	n := ckkscontext.n

	ckkscontext.bigintChain = genBigIntChain(params.qi)

	if ckkscontext.contextQ, err = ring.NewContextWithParams(n, params.qi); err != nil {
		panic(err)
	}

	if params.PiCount() != 0 {
		if ckkscontext.contextP, err = ring.NewContextWithParams(n, params.pi); err != nil {
			panic(err)
		}
	}

	if ckkscontext.contextQP, err = ring.NewContextWithParams(n, append(params.qi, params.pi...)); err != nil {
		panic(err)
	}

	ckkscontext.galElRotColLeft = ring.GenGaloisParams(n, GaloisGen)
	ckkscontext.galElRotColRight = ring.GenGaloisParams(n, ring.ModExp(GaloisGen, 2*n-1, 2*n))
	ckkscontext.galElConjugate = 2*n - 1

	return ckkscontext

}
