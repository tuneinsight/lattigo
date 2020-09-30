package main

import (
	"fmt"
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
	"math"
	"math/big"
	"sort"
)

var h = uint64(192)
var logN = uint64(16)
var trials = int(1000)

func main() {

	var err error

	var prng utils.PRNG
	if prng, err = utils.NewPRNG(); err != nil {
		panic(err)
	}

	qiQ := []uint64{0xffa0001}
	qiQh := []uint64{0xffa0001, 0xfffffffff840001}

	var ringQ, ringQh *ring.Ring
	if ringQ, err = ring.NewRing(1<<logN, qiQ); err != nil {
		panic(err)
	}
	if ringQh, err = ring.NewRing(1<<logN, qiQh); err != nil {
		panic(err)
	}

	coeffsBigint := make([]*big.Int, ringQh.N)
	C := make(map[uint64]uint64)

	ternarySampler := ring.NewTernarySamplerSparse(prng, ringQh, h, true)
	uniformSampler := ring.NewUniformSampler(prng, ringQ)
	c0 := ringQh.NewPoly()
	c1 := ringQh.NewPoly()
	sk := ringQh.NewPoly()

	for j := 0; j < trials; j++ {

		sk.Zero()

		// Sample sk
		ternarySampler.Read(sk)
		ringQh.NTT(sk, sk)

		// Sample c1
		uniformSampler.Read(c1)
		modUp(c1, ringQ, ringQh)
		ringQh.NTT(c1, c1)

		// c0 = -sk * c1
		ringQ.MulCoeffsMontgomery(sk, c1, c0)
		ringQ.Neg(c0, c0)

		// Extends the basis from q0 to q0*q1
		ringQ.InvNTT(c0, c0)
		modUp(c0, ringQ, ringQh)
		ringQh.NTT(c0, c0)

		// c0 = [[-sk*c1]_{q0} + sk*c1]_{q0*q1}
		ringQh.MulCoeffsMontgomeryAndAdd(sk, c1, c0)
		ringQh.InvNTT(c0, c0)

		// Decodes
		ringQh.PolyToBigint(c0, coeffsBigint)

		Qh := ringQh.ModulusBigint
		qHalf := new(big.Int)
		qHalf.Set(Qh)
		qHalf.Rsh(qHalf, 1)
		var sign int
		var c uint64
		q := qiQ[0]

		// Centers around q0*q1, divides by q0, puts in the map result
		for i := range coeffsBigint {
			sign = coeffsBigint[i].Cmp(qHalf)
			if sign == 1 || sign == 0 {
				coeffsBigint[i].Sub(coeffsBigint[i], Qh)
			}
			c = coeffsBigint[i].Abs(coeffsBigint[i]).Uint64() / q
			C[c]++
		}
	}

	keys := make([]int, len(C))
	var i int
	for k := range C {
		keys[i] = int(k)
		i++
	}

	sort.Ints(keys)

	fmt.Printf("[\n")
	for _, k := range keys {
		fmt.Printf("(%d,%6f),\n", k, math.Log2(float64(C[uint64(k)])/float64(trials*(1<<logN))))
	}
	fmt.Printf("]")

}

func modUp(c *ring.Poly, ringQ, ringQh *ring.Ring) {

	//Centers the values around Q0 and extends the basis from Q0 to QL
	Q := ringQ.Modulus[0]
	bredparams := ringQh.GetBredParams()

	var coeff, qi uint64

	for j := uint64(0); j < ringQh.N; j++ {

		coeff = c.Coeffs[0][j]

		for i := 1; i < len(ringQh.Modulus); i++ {

			qi = ringQh.Modulus[i]

			if coeff > (Q >> 1) {
				c.Coeffs[i][j] = qi - ring.BRedAdd(Q-coeff, qi, bredparams[i])
			} else {
				c.Coeffs[i][j] = ring.BRedAdd(coeff, qi, bredparams[i])
			}
		}
	}
}
