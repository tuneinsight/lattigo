package main

import (
	"flag"
	"fmt"
	"math/big"
	"time"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// Vectorized oblivious evaluation is a two-party protocol for the function f(x) = ax + b where a sender
// inputs the field elements a, b, and a receiver inputs x and learns f(x).
//
// This example implements the secret-key and passively-secure vOLE protocol of Figure 5 of
// "Efficient Protocols for Oblivious Linear Function Evaluation from Ring-LWE" by Baum C. et al.
// https://eprint.iacr.org/2020/970
//
//
// We work in the ring R = Z_{Q}/(X^N + 1), for Q the product of k NTT-friendly primes and N a power of two.
// We define Q > P > M with P|Q and M|Q and M|P, i.e. P and M are themselves factors of Q.
// Thus, we also have R_{M} -> R_{P} -> R_{Q}. A polynomial mod M implies that this polynomial belongs to R_{M};
// the same follows for mod P and mod Q.
//
// Given two parties Alice and Bob, with their respective secret inputs v and u (mod M), the goal of this protocol is
// for Alice and Bob to construct the values alpha, beta (mod M) such that alpha + beta = u*v (mod M).
//
// 1. Setup phase
//		(a) Alice and Bob sample each a low-norm secret polynomial: skAlice and skBob (mod Q) respectively and set
//			sigmaAlice + sigmaBob = skAlice * skBob (mod Q)
//		(b) Alice and Bob sample two public random polynomials: a and a' (mod Q)
// 2. First Message. Given u (mod M), Bob's input
//		(a) Bob samples a noise polynomial eBob (mod Q) from a small-variance truncated discrete Gaussian distribution
// 			and sends c = (Q/P) * u + a*skBob + eBob to Alice (mod P)
//		(b) Alice computes rhoAlice = (skAlice * c - a * sigmaAlice) * (P/Q) (mod P)
//		(c) Bob computes rhoBob = -(a * sigmaBob) * (P/Q) (mod P)
//
// At this step, it should hold that u * skAlice = rhoAlice + rhoBob (mod P)
//
// 3. Second Message. Given v, Alice's input
//		(a) Alice samples a noise polynomial eAlice (mod P) from a small-variance truncated discrete Gaussian distribution
//			and sends d = (P/M) * v + (a'*skAlice + eAlice) (mod P) to Bob
//		(b) Bob outputs beta = (u * d - a' * rhoBob) * (M/P) (mod M)
//		(c) Alice outputs alpha = -(a' * rhoAlice) * (M/P) (mod M)
//
// It should hold that u * v = alpha + beta (mod M)

var flagShort = flag.Bool("short", false, "run the example on the 3 first parameters sets only.")

type parameters struct {
	logN   int    // Ring degree
	logQ   [2]int // logQ[0] = #Primes, logQ[1] = Primes bit-size
	n      int    // Number of vOLE to run
	plevel int    // Maximum level of the modulus P (level 0 is the lowest)
	mlevel int    // Maximum level of the modulus M (level 0 is the lowest)
}

var benchParameters = []parameters{
	{logN: 13, logQ: [2]int{4, 60}, n: 1, plevel: 2, mlevel: 0},
	{logN: 14, logQ: [2]int{6, 60}, n: 1, plevel: 3, mlevel: 0},
	{logN: 14, logQ: [2]int{8, 60}, n: 1, plevel: 5, mlevel: 1},
	{logN: 13, logQ: [2]int{4, 60}, n: 8, plevel: 2, mlevel: 0},
	{logN: 14, logQ: [2]int{6, 60}, n: 8, plevel: 3, mlevel: 0},
	{logN: 14, logQ: [2]int{8, 60}, n: 8, plevel: 5, mlevel: 1},
	{logN: 13, logQ: [2]int{4, 60}, n: 64, plevel: 2, mlevel: 0},
	{logN: 14, logQ: [2]int{6, 60}, n: 64, plevel: 3, mlevel: 0},
	{logN: 14, logQ: [2]int{8, 60}, n: 64, plevel: 5, mlevel: 1},
	{logN: 13, logQ: [2]int{4, 60}, n: 128, plevel: 2, mlevel: 0},
	{logN: 14, logQ: [2]int{6, 60}, n: 128, plevel: 3, mlevel: 0},
	{logN: 14, logQ: [2]int{8, 60}, n: 128, plevel: 5, mlevel: 1},
}

type vOLErings struct {
	ringQ *ring.Ring // Z_{Q}/(X^N+1)
	qDivP *big.Int   // Q/P
	pDivM *big.Int   // P/M
}

func newvOLErings(params parameters) *vOLErings {

	if params.mlevel >= params.plevel {
		panic("mlevel must be strictly smaller than plevel")
	}

	var err error

	N := 1 << params.logN

	rings := new(vOLErings)

	// Generate logQ[0] NTT-friendly primes each close to 2^logQ[1]
	primes := ring.GenerateNTTPrimes(params.logQ[1], 2*N, params.logQ[0])

	if rings.ringQ, err = ring.NewRing(N, primes); err != nil {
		panic(err)
	}

	rings.qDivP = ring.NewInt(1)
	for _, qi := range primes[params.plevel+1:] {
		rings.qDivP.Mul(rings.qDivP, ring.NewUint(qi))
	}

	rings.pDivM = ring.NewInt(1)
	for _, qi := range primes[params.mlevel+1 : params.plevel+1] {
		rings.pDivM.Mul(rings.pDivM, ring.NewUint(qi))
	}

	return rings
}

type lowNormSampler struct {
	baseRing *ring.Ring
	coeffs   []*big.Int
}

func newLowNormSampler(baseRing *ring.Ring) (lns *lowNormSampler) {
	lns = new(lowNormSampler)
	lns.baseRing = baseRing
	lns.coeffs = make([]*big.Int, baseRing.N)
	return
}

// Samples a uniform polynomial in Z_{norm}/(X^N + 1)
func (lns *lowNormSampler) newPolyLowNorm(norm []uint64) (pol *ring.Poly) {
	bound := ring.NewUint(0)
	for _, qi := range norm {
		bound.Add(bound, ring.NewUint(qi))
	}

	pol = lns.baseRing.NewPoly()

	for i := range lns.coeffs {
		lns.coeffs[i] = ring.RandInt(bound)
	}

	lns.baseRing.SetCoefficientsBigint(lns.coeffs, pol)

	return
}

func main() {

	flag.Parse()

	paramsSets := benchParameters
	if *flagShort {
		paramsSets = paramsSets[:6] // skips the benchmarks for n=64 and 128 in -short mode
	}

	for _, param := range paramsSets {

		// Crypto Context Initialization

		volerings := newvOLErings(param)

		ringQ := volerings.ringQ

		n := param.n
		qlevel := len(ringQ.Modulus) - 1
		plevel := param.plevel
		mlevel := param.mlevel

		fmt.Printf("Params : n=%d logN=%d qlevel=%d plevel=%d mlevel=%d\n", n, param.logN, qlevel, plevel, mlevel)

		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}

		ternarySamplerMontgomeryQ := ring.NewTernarySampler(prng, ringQ, 1.0/3.0, true)
		gaussianSamplerQ := ring.NewGaussianSampler(prng, ringQ, 3.2, 19)
		uniformSamplerQ := ring.NewUniformSampler(prng, ringQ)
		lowNormUniformQ := newLowNormSampler(ringQ)

		var elapsed, TotalTime, AliceTime, BobTime time.Duration
		start := time.Now()

		// ********* 1. SETUP *********

		// NTT(MForm(skBob))
		skBob := ternarySamplerMontgomeryQ.ReadNew()
		ringQ.NTT(skBob, skBob)

		// NTT(MForm(skAlice))
		skAlice := ternarySamplerMontgomeryQ.ReadNew()
		ringQ.NTT(skAlice, skAlice)

		// NTT(MForm(sigmaAlicelice))
		sigmaAlice := uniformSamplerQ.ReadNew()

		// NTT(MForm(sigmaBobob))
		sigmaBob := ringQ.NewPoly()

		// NTT(MForm(skBob * skAlice))
		ringQ.MulCoeffsMontgomery(skAlice, skBob, sigmaBob)

		// NTT(MForm(sigmaBob)) = NTT(MForm(ska_a * skAlice) - MForm(sigmaAlice))
		ringQ.Sub(sigmaBob, sigmaAlice, sigmaBob)

		a := make([]*ring.Poly, n)
		aprime := make([]*ring.Poly, n)

		// Sample common random poly vectors
		// NTT(a) in Z_Q
		// NTT(MForm(a')) in Z_P
		for i := range a {
			a[i] = uniformSamplerQ.ReadNew()
			aprime[i] = uniformSamplerQ.ReadLvlNew(plevel)
			ringQ.MFormLvl(plevel, aprime[i], aprime[i])
		}

		elapsed = time.Since(start)
		fmt.Printf("Setup : %s\n", elapsed)

		// Generate inputs and allocate memory
		start = time.Now()
		u := make([]*ring.Poly, n)
		v := make([]*ring.Poly, n)
		c := make([]*ring.Poly, n)
		d := make([]*ring.Poly, n)
		rhoAlice := make([]*ring.Poly, n)
		rhoBob := make([]*ring.Poly, n)
		alpha := make([]*ring.Poly, n)
		beta := make([]*ring.Poly, n)
		tmp := ringQ.NewPoly()

		for i := 0; i < n; i++ {

			// Generate uniform inputs u and v (in R_m) in the time domain
			u[i] = lowNormUniformQ.newPolyLowNorm(ringQ.Modulus[:mlevel+1])

			ringQ.NTT(u[i], u[i])

			v[i] = lowNormUniformQ.newPolyLowNorm(ringQ.Modulus[:mlevel+1])

			c[i] = ringQ.NewPoly()
			d[i] = ringQ.NewPolyLvl(plevel)
			rhoAlice[i] = ringQ.NewPoly()
			rhoBob[i] = ringQ.NewPoly()
			alpha[i] = ringQ.NewPolyLvl(plevel)
			beta[i] = ringQ.NewPolyLvl(plevel)
		}

		elapsed = time.Since(start)
		fmt.Printf("Gen Inputs : %s\n", elapsed)

		// ********* 2. FIRST MESSAGE *********

		// First Message starts (Bob)
		start = time.Now()
		for i := 0; i < n; i++ {

			// c = e
			gaussianSamplerQ.ReadAndAddLvl(qlevel, c[i])

			// c = NTT(e)
			ringQ.NTT(c[i], c[i])

			// tmp = NTT(u * Q/P)
			ringQ.MulScalarBigint(u[i], volerings.qDivP, tmp)

			// c = NTT(u * (Q/P) + e)
			ringQ.Add(c[i], tmp, c[i])

			// c = NTT(u * (Q/P) + e) + NTT(a) * MForm(NTT(skAlice))
			// c = NTT(u * (Q/P) + e + a*skAlice)
			ringQ.MulCoeffsMontgomeryAndAdd(a[i], skBob, c[i])
		}

		elapsed = time.Since(start)
		fmt.Printf("2.(a) : %s\n", elapsed)
		TotalTime = elapsed
		BobTime = elapsed

		// Alice
		start = time.Now()
		for i := 0; i < n; i++ {

			// rhoAlice = NTT(MForm(skBob)) * NTT(u * (Q/P) + e + a*skAlice)
			// rhoAlice = NTT(skBob * (u * (Q/P) + e + a*skAlice))
			ringQ.MulCoeffsMontgomery(skAlice, c[i], rhoAlice[i])

			// rhoAlice = NTT(skBob * (u * (Q/P) + e + a*skAlice)) - NTT(a) * NTT(MForm(sigmaAlice))
			// rhoAlice = NTT(skBob * (u * (Q/P) + e + a*skAlice) - a*sigmaAlice)
			ringQ.MulCoeffsMontgomeryAndSub(a[i], sigmaAlice, rhoAlice[i])

			// rhoAlice = NTT(skBob * u)  + NTT(a*skBob*skAlice - a*sigmaAlice) * (P/Q)
			ringQ.DivRoundByLastModulusManyNTT(rhoAlice[i], rhoAlice[i], qlevel-plevel)
		}

		elapsed = time.Since(start)
		fmt.Printf("2.(c) : %s\n", elapsed)
		TotalTime += elapsed
		AliceTime = elapsed

		// Bob
		start = time.Now()
		for i := 0; i < n; i++ {

			// rhoBob = NTT(a) * NTT(MForm(sigmaBob))
			// rhoBob = NTT(a * sigmaBob)
			ringQ.MulCoeffsMontgomery(a[i], sigmaBob, rhoBob[i])

			// rhoBob = NTT(a * sigmaBob * (P/Q))
			ringQ.DivRoundByLastModulusManyNTT(rhoBob[i], rhoBob[i], qlevel-plevel)

			// rhoBob = NTT(-a * sigmaBob) * (P/Q)
			ringQ.NegLvl(plevel, rhoBob[i], rhoBob[i])
		}

		elapsed = time.Since(start)
		fmt.Printf("2.(d) : %s\n", elapsed)
		TotalTime += elapsed
		BobTime += elapsed

		fmt.Printf("Total time First Message : %s\n", TotalTime)
		fmt.Printf("Alice time First Message : %s\n", AliceTime)
		fmt.Printf("Bob   time First Message : %s\n", BobTime)

		// ********* VERIFY CORRECTNESS *********

		checkMessage1a := ringQ.NewPolyLvl(plevel)
		checkMessage1b := ringQ.NewPolyLvl(plevel)

		nerrors := 0
		for i := 0; i < n; i++ {

			// checkMessage1a = NTT(u) * NTT(MForm(skBob))

			ringQ.MulCoeffsMontgomeryLvl(plevel, u[i], skAlice, checkMessage1a)

			// rhoAlice = NTT(skBob * u)  + NTT(a*skBob*skAlice - a*sigmaAlice) * P/Q
			// rhoBob = NTT(-a * sigmaBob) * P/Q

			// NTT(a*skBob*skAlice - a*sigmaAlice) * (P/Q) = NTT(a*skBob*skAlice) - NTT(a) * NTT(MForm(ska_a * skAlice - sigmaBob)
			// = NTT(a*skBob*skAlice - a * (skBob*skAlice - sigmaBob))
			// = NTT(a*skBob*skAlice - a*skBob*sk-b + a*sigmaBob)
			// = NTT(a*sigmaBob) * P/Q

			// -> rhoAlice + rhoBob = NTT(skBob * u) + NTT(-a * sigmaBob) * (P/Q) + NTT(a*sigmaBob) * P/Q
			//					    = NTT(skBob * u)

			ringQ.AddLvl(plevel, rhoAlice[i], rhoBob[i], checkMessage1b)
			if !ringQ.EqualLvl(plevel, checkMessage1a, checkMessage1b) {
				nerrors++
			}
		}

		fmt.Printf("Errors First Message : %d\n", nerrors)

		// ********* 3. SECOND MESSAGE *********

		// Second Message, Alice
		start = time.Now()
		ringQ.InvMFormLvl(plevel, skAlice, skAlice)
		for i := 0; i < n; i++ {
			// d = v * P/M
			ringQ.MulScalarBigintLvl(plevel, v[i], volerings.pDivM, d[i])

			// d = v * (P/M) + e
			gaussianSamplerQ.ReadAndAddLvl(plevel, d[i])

			// d = NTT(v * (P/M) + e)
			ringQ.NTTLvl(plevel, d[i], d[i])

			// d = NTT(v * (P/M) + e) + NTT(MForm(a'))*NTT(skAlice)
			//   = NTT(v * (P/M) + e + a'*skAlice)
			ringQ.MulCoeffsMontgomeryAndAddLvl(plevel, aprime[i], skAlice, d[i])
		}

		elapsed = time.Since(start)
		fmt.Printf("3.(a) : %s\n", elapsed)
		TotalTime = elapsed
		AliceTime = elapsed

		// Bob
		start = time.Now()
		for i := 0; i < n; i++ {

			// beta = NTT(v * (P/M) + e + a'*skAlice) * NTT(MForm(u))
			//      = NTT(u * (v * (P/M) + e + a'*skAlice))

			ringQ.MFormLvl(plevel, u[i], u[i])
			ringQ.MulCoeffsMontgomeryLvl(plevel, u[i], d[i], beta[i])

			// beta = NTT(u * (v * (P/M) + e + a'*skAlice)) - NTT(MForm(a')) * NTT(-a * sigmaBob) * (P/Q)
			// 		= NTT(u * (v * (P/M) + e + a'*skAlice) - a'*-a*sigmaBob * (P/Q))
			ringQ.MulCoeffsMontgomeryAndSubLvl(plevel, aprime[i], rhoBob[i], beta[i])

			// beta = u*(v * (P/M) + e + a'*skAlice) * u -  a'*-a*sigmaBob * (P/Q)
			ringQ.InvNTTLvl(plevel, beta[i], beta[i])

			// beta = (u*(v * (P/M) + e + a'*skAlice) * u -  a'*-a*sigmaBob * (P/Q)) * (M/P)
			// 		= (M/P) * u * (v * (P/M) + e + a'*skAlice) - a'*-a*sigmaBob * (M/Q)
			//		= u*v + a'*skAlice*u*(M/P) + a'*a*sigmaBob * (M/Q)
			ringQ.DivRoundByLastModulusMany(beta[i], beta[i], plevel-mlevel)
		}

		elapsed = time.Since(start)
		fmt.Printf("3.(c) : %s\n", elapsed)
		TotalTime += elapsed
		BobTime = elapsed

		// Alice
		start = time.Now()
		for i := 0; i < n; i++ {
			// alpha = NTT(MForm(a')) * (NTT(skAlice * u)  + NTT(a*skBob*skAlice - a*sigmaAlice) * (P/Q))
			// 		 = NTT(a'*skAlice*u + (a'*a*skBob*skAlice - a'*a*sigmaAlice) * (P/Q))

			ringQ.MulCoeffsMontgomeryLvl(plevel, aprime[i], rhoAlice[i], alpha[i])
			ringQ.InvNTTLvl(plevel, alpha[i], alpha[i])

			// alpha = (a'*skAlice*u + (a'*a*skBob*skAlice - a'*a*sigmaAlice) * (P/Q)) * (M/P)
			//		 = a'*skAlice*u*(M/P) + a'*a*skBob*skAlice*(M/Q) - a'*a*sigmaAlice*(M/Q)
			ringQ.DivRoundByLastModulusMany(alpha[i], alpha[i], plevel-mlevel)

			// alpha = - a'*skAlice*u*(M/P) - a'*a*skBob*skAlice*(M/Q) + a'*a*sigmaAlice*(M/Q)
			// 	 	 = - a'*skAlice*u*(M/P) - a'*a*skBob*skAlice*(M/Q) + a'*a*(skBob*skAlice - sigmaBob)*(M/Q)
			//		 = - a'*skAlice*u*(M/P) - a'*a*skBob*skAlice*(M/Q) + a'*a*skBob*skAlice*(M/Q) - a'*a*sigmaBob*(M/Q)
			// 		 = - a'*skAlice*u*(M/P) - a'*a*sigmaBob*(M/Q)
			ringQ.NegLvl(mlevel, alpha[i], alpha[i])
		}

		elapsed = time.Since(start)
		fmt.Printf("3.(d) : %s\n", elapsed)
		TotalTime += elapsed
		AliceTime += elapsed

		fmt.Printf("Total time Second Message : %s\n", TotalTime)
		fmt.Printf("Alice time Second Message : %s\n", AliceTime)
		fmt.Printf("Bob   time Second Message : %s\n", BobTime)

		// ********* VERIFY CORRECTNESS *********
		nerrors = 0
		for i := 0; i < n; i++ {

			// alpha =     - a'*skAlice*u*(M/P) - a'*a*sigmaBob*(M/Q)
			// beta  = u*v + a'*skAlice*u*(M/P) + a'*a*sigmaBob*(M/Q)

			// u * v = alpha + beta

			ringQ.NTTLvl(mlevel, v[i], checkMessage1b)
			ringQ.MulCoeffsMontgomeryLvl(mlevel, u[i], checkMessage1b, checkMessage1a)

			ringQ.InvNTTLvl(mlevel, checkMessage1a, checkMessage1a)

			ringQ.AddLvl(mlevel, alpha[i], beta[i], checkMessage1b)

			if !ringQ.EqualLvl(mlevel, checkMessage1a, checkMessage1b) {
				nerrors++
			}
		}

		fmt.Printf("Errors on Second Message : %d\n", nerrors)
		fmt.Printf("\n")
	}

}
