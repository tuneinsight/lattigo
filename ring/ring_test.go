package ring

import (
	"bytes"
	"fmt"
	"math"
	"math/big"
	"testing"

	"github.com/tuneinsight/lattigo/v6/utils/buffer"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
	"github.com/tuneinsight/lattigo/v6/utils/structs"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

var T = uint64(0x3ee0001)
var DefaultSigma = 3.2
var DefaultBound = 6.0 * DefaultSigma

func testString(opname string, ringQ *Ring) string {
	return fmt.Sprintf("%s/N=%d/limbs=%d", opname, ringQ.N(), ringQ.ModuliChainLength())
}

type testParams struct {
	ringQ           *Ring
	ringP           *Ring
	prng            sampling.PRNG
	uniformSamplerQ *UniformSampler
	uniformSamplerP *UniformSampler
}

func genTestParams(defaultParams Parameters) (tc *testParams, err error) {

	tc = new(testParams)

	if tc.ringQ, err = NewRing(1<<defaultParams.logN, defaultParams.qi); err != nil {
		return nil, err
	}
	if tc.ringP, err = NewRing(1<<defaultParams.logN, defaultParams.pi); err != nil {
		return nil, err
	}
	if tc.prng, err = sampling.NewPRNG(); err != nil {
		return nil, err
	}
	tc.uniformSamplerQ = NewUniformSampler(tc.prng, tc.ringQ)
	tc.uniformSamplerP = NewUniformSampler(tc.prng, tc.ringP)
	return
}

func TestRing(t *testing.T) {

	var err error

	var defaultParams = testParameters[:] // the default test

	testNewRing(t)
	testShift(t)

	for _, defaultParam := range defaultParams[:] {

		var tc *testParams
		if tc, err = genTestParams(defaultParam); err != nil {
			t.Fatal(err)
		}

		testNTTConjugateInvariant(tc, t)
		testPRNG(tc, t)
		testGenerateNTTPrimes(tc, t)
		testDivFloorByLastModulusMany(tc, t)
		testDivRoundByLastModulusMany(tc, t)
		testMarshalBinary(tc, t)
		testWriterAndReader(tc, t)
		testSampler(tc, t)
		testModularReduction(tc, t)
		testMForm(tc, t)
		testMulScalarBigint(tc, t)
		testExtendBasis(tc, t)
		testMultByMonomial(tc, t)

	}
}

func testNTTConjugateInvariant(tc *testParams, t *testing.T) {

	t.Run(testString("NTTConjugateInvariant", tc.ringQ), func(t *testing.T) {

		ringQ := tc.ringQ
		Q := ringQ.ModuliChain()
		N := ringQ.N()
		ringQ2N, _ := NewRing(N<<1, Q)
		ringQConjugateInvariant, _ := NewRingFromType(N, Q, ConjugateInvariant)

		sampler := NewUniformSampler(tc.prng, ringQ)
		p1 := sampler.ReadNew()
		p2 := ringQ2N.NewPoly()

		for i, qi := range Q {
			copy(p2.Coeffs[i], p1.Coeffs[i])
			for j := 1; j < N; j++ {
				p2.Coeffs[i][N*2-j] = qi - p2.Coeffs[i][j]
			}
		}

		ringQ2N.NTT(p2, p2)
		ringQ2N.MForm(p2, p2)
		ringQ2N.MulCoeffsMontgomery(p2, p2, p2)
		ringQ2N.IMForm(p2, p2)
		ringQ2N.INTT(p2, p2)

		p1tmp := ringQ2N.NewPoly()

		ringQConjugateInvariant.NTT(p1, p1tmp)
		ringQConjugateInvariant.MForm(p1tmp, p1tmp)
		ringQConjugateInvariant.MulCoeffsMontgomery(p1tmp, p1tmp, p1tmp)
		ringQConjugateInvariant.IMForm(p1tmp, p1tmp)
		ringQConjugateInvariant.INTT(p1tmp, p1)

		for j := range Q {
			for i := 0; i < N; i++ {
				require.Equal(t, p1.Coeffs[j][i], p2.Coeffs[j][i])
			}
		}
	})
}

func testNewRing(t *testing.T) {
	t.Run("NewRing", func(t *testing.T) {
		r, err := NewRing(0, nil)
		require.Nil(t, r)
		require.Error(t, err)

		r, err = NewRing(0, []uint64{})
		require.Nil(t, r)
		require.Error(t, err)

		r, err = NewRing(4, []uint64{})
		require.Nil(t, r)
		require.Error(t, err)

		r, err = NewRing(8, []uint64{})
		require.Nil(t, r)
		require.Error(t, err)

		r, err = NewRing(16, []uint64{7}) // Passing non NTT-enabling coeff modulus
		require.NotNil(t, r)              // Should still return a Ring instance
		require.Error(t, err)             // Should also return an error due to non NTT

		r, err = NewRing(16, []uint64{4}) // Passing non prime moduli
		require.NotNil(t, r)              // Should still return a Ring instance
		require.Error(t, err)             // Should also return an error due to non NTT

		r, err = NewRing(16, []uint64{97, 7}) // Passing a NTT-enabling and a non NTT-enabling coeff modulus
		require.NotNil(t, r)                  // Should still return a Ring instance
		require.Error(t, err)                 // Should also return an error due to non NTT

		r, err = NewRing(16, []uint64{97, 97}) // Passing non CRT-enabling coeff modulus
		require.Nil(t, r)                      // Should not return a Ring instance
		require.Error(t, err)

		r, err = NewRing(16, []uint64{97}) // Passing NTT-enabling coeff modulus
		require.NotNil(t, r)
		require.NoError(t, err)

	})
}

func testPRNG(tc *testParams, t *testing.T) {

	t.Run(testString("PRNG", tc.ringQ), func(t *testing.T) {

		var err error

		var prng1, prng2 sampling.PRNG

		if prng1, err = sampling.NewKeyedPRNG(nil); err != nil {
			t.Fatal(err)
		}

		if prng2, err = sampling.NewKeyedPRNG(nil); err != nil {
			t.Fatal(err)
		}

		crsGenerator1 := NewUniformSampler(prng1, tc.ringQ)
		crsGenerator2 := NewUniformSampler(prng2, tc.ringQ)

		p0 := crsGenerator1.ReadNew()
		p1 := crsGenerator2.ReadNew()

		require.True(t, tc.ringQ.Equal(p0, p1))
	})

}

func testGenerateNTTPrimes(tc *testParams, t *testing.T) {

	t.Run(testString("GenerateNTTPrimes", tc.ringQ), func(t *testing.T) {

		NthRoot := tc.ringQ.NthRoot()

		g := NewNTTFriendlyPrimesGenerator(55, NthRoot)

		primes, err := g.NextAlternatingPrimes(tc.ringQ.ModuliChainLength())

		require.NoError(t, err)

		// Checks that all returned are unique pair-wise
		// primes with an Nth-primitive root.
		list := map[uint64]bool{}
		for _, q := range primes {
			require.Equal(t, q&uint64(NthRoot-1), uint64(1))
			require.True(t, IsPrime(q), q)
			_, ok := list[q]
			require.False(t, ok)
			list[q] = true
		}

		upstreamPrimes, err := g.NextUpstreamPrimes(tc.ringQ.ModuliChainLength())
		require.NoError(t, err)
		for i := range upstreamPrimes {
			if i == 0 {
				require.True(t, IsPrime(upstreamPrimes[i]))
			} else {
				require.True(t, IsPrime(upstreamPrimes[i]) && upstreamPrimes[i] > upstreamPrimes[i-1])
			}

		}

		downstreamPrimes, err := g.NextDownstreamPrimes(tc.ringQ.ModuliChainLength())
		require.NoError(t, err)
		for i := range downstreamPrimes {
			if i == 0 {
				require.True(t, IsPrime(downstreamPrimes[i]))
			} else {
				require.True(t, IsPrime(downstreamPrimes[i]) && downstreamPrimes[i] < downstreamPrimes[i-1])
			}
		}
	})
}

func testDivFloorByLastModulusMany(tc *testParams, t *testing.T) {

	t.Run(testString("DivFloorByLastModulusMany", tc.ringQ), func(t *testing.T) {

		prng, _ := sampling.NewPRNG()

		N := tc.ringQ.N()

		level := tc.ringQ.Level()

		ringQ := tc.ringQ.AtLevel(level)

		coeffs := make([]*big.Int, N)
		for i := 0; i < N; i++ {
			coeffs[i] = bignum.RandInt(prng, tc.ringQ.ModulusAtLevel[level])
			coeffs[i].Quo(coeffs[i], bignum.NewInt(10))
		}

		nbRescales := level

		coeffsWant := make([]*big.Int, N)
		for i := range coeffs {
			coeffsWant[i] = new(big.Int).Set(coeffs[i])
			for j := 0; j < nbRescales; j++ {
				coeffsWant[i].Quo(coeffsWant[i], bignum.NewInt(tc.ringQ.SubRings[level-j].Modulus))
			}
		}

		polTest0 := tc.ringQ.NewPoly()
		polTest1 := tc.ringQ.NewPoly()
		polWant := tc.ringQ.NewPoly()
		buff := tc.ringQ.NewPoly()

		ringQ.SetCoefficientsBigint(coeffs, polTest0)
		ringQ.SetCoefficientsBigint(coeffsWant, polWant)
		ringQ.DivFloorByLastModulusMany(nbRescales, polTest0, buff, polTest1)

		for i := 0; i < N; i++ {
			for j := 0; j < polTest0.Level()-nbRescales+1; j++ {
				require.Equalf(t, polWant.Coeffs[j][i], polTest1.Coeffs[j][i], "coeff %v Qi%v = %s", i, j, coeffs[i].String())
			}
		}
	})
}

func testDivRoundByLastModulusMany(tc *testParams, t *testing.T) {

	t.Run(testString("bignum.DivRoundByLastModulusMany", tc.ringQ), func(t *testing.T) {

		prng, _ := sampling.NewPRNG()

		N := tc.ringQ.N()

		level := tc.ringQ.Level()

		ringQ := tc.ringQ.AtLevel(level)

		coeffs := make([]*big.Int, N)
		for i := 0; i < N; i++ {
			coeffs[i] = bignum.RandInt(prng, tc.ringQ.ModulusAtLevel[level])
			coeffs[i].Quo(coeffs[i], bignum.NewInt(10))
		}

		nbRescals := level

		coeffsWant := make([]*big.Int, N)
		for i := range coeffs {
			coeffsWant[i] = new(big.Int).Set(coeffs[i])
			for j := 0; j < nbRescals; j++ {
				bignum.DivRound(coeffsWant[i], bignum.NewInt(tc.ringQ.SubRings[level-j].Modulus), coeffsWant[i])
			}
		}

		polTest0 := tc.ringQ.NewPoly()
		polTest1 := tc.ringQ.NewPoly()
		polWant := tc.ringQ.NewPoly()
		buff := tc.ringQ.NewPoly()

		ringQ.SetCoefficientsBigint(coeffs, polTest0)
		ringQ.SetCoefficientsBigint(coeffsWant, polWant)

		ringQ.DivRoundByLastModulusMany(nbRescals, polTest0, buff, polTest1)

		for i := 0; i < N; i++ {
			for j := 0; j < polTest0.Level()-nbRescals+1; j++ {
				require.Equalf(t, polWant.Coeffs[j][i], polTest1.Coeffs[j][i], "coeff %v Qi%v = %s", i, j, coeffs[i].String())
			}
		}
	})
}

func testMarshalBinary(tc *testParams, t *testing.T) {

	t.Run(testString("MarshalBinary/Ring", tc.ringQ), func(t *testing.T) {

		var err error

		var data []byte
		if data, err = tc.ringQ.MarshalBinary(); err != nil {
			t.Fatal(err)
		}

		ringQTest := new(Ring)
		if err = ringQTest.UnmarshalBinary(data); err != nil {
			t.Fatal(err)
		}

		require.Equal(t, ringQTest, tc.ringQ)
	})

	t.Run(testString("MarshalBinary/Poly", tc.ringQ), func(t *testing.T) {
		poly := tc.uniformSamplerQ.ReadNew()
		buffer.RequireSerializerCorrect(t, &poly)
	})

	t.Run(testString("structs/PolyVector", tc.ringQ), func(t *testing.T) {

		polys := make([]Poly, 4)

		for i := range polys {
			polys[i] = tc.uniformSamplerQ.ReadNew()
		}

		v := structs.Vector[Poly](polys)

		buffer.RequireSerializerCorrect(t, &v)
	})

	t.Run(testString("structs/PolyMatrix", tc.ringQ), func(t *testing.T) {

		polys := make([][]Poly, 4)

		for i := range polys {
			polys[i] = make([]Poly, 4)

			for j := range polys {
				polys[i][j] = tc.uniformSamplerQ.ReadNew()
			}
		}

		m := structs.Matrix[Poly](polys)

		buffer.RequireSerializerCorrect(t, &m)
	})

	t.Run(testString("structs/PolyMap", tc.ringQ), func(t *testing.T) {

		m := make(structs.Map[int, Poly], 4)

		for i := 0; i < 4; i++ {
			p := tc.uniformSamplerQ.ReadNew()
			m[i] = &p
		}

		buffer.RequireSerializerCorrect(t, &m)
	})
}

func testWriterAndReader(tc *testParams, t *testing.T) {

	t.Run(testString("WriterAndReader/Poly", tc.ringQ), func(t *testing.T) {

		p := tc.uniformSamplerQ.ReadNew()

		data := make([]byte, 0, p.BinarySize())

		buf := bytes.NewBuffer(data) // Complient to io.Writer and io.Reader

		if n, err := p.WriteTo(buf); err != nil {
			t.Fatal(err)
		} else {
			if int(n) != p.BinarySize() {
				t.Fatal()
			}
		}

		if data2, err := p.MarshalBinary(); err != nil {
			t.Fatal(err)
		} else {
			if !bytes.Equal(buf.Bytes(), data2) {
				t.Fatal()
			}
		}

		pTest := new(Poly)
		if n, err := pTest.ReadFrom(buf); err != nil {
			t.Fatal(err)
		} else {
			if int(n) != p.BinarySize() {
				t.Fatal()
			}
		}

		for i := range tc.ringQ.SubRings {
			require.Equal(t, p.Coeffs[i][:tc.ringQ.N()], pTest.Coeffs[i][:tc.ringQ.N()])
		}
	})
}

func testSampler(tc *testParams, t *testing.T) {

	N := tc.ringQ.N()

	t.Run(testString("Sampler/Uniform", tc.ringQ), func(t *testing.T) {
		pol := tc.ringQ.NewPoly()
		tc.uniformSamplerQ.Read(pol)

		for i, qi := range tc.ringQ.ModuliChain() {
			coeffs := pol.Coeffs[i]
			for j := 0; j < N; j++ {
				require.False(t, coeffs[j] > qi)
			}
		}
	})

	t.Run(testString("Sampler/Gaussian/SmallSigma", tc.ringQ), func(t *testing.T) {

		dist := DiscreteGaussian{Sigma: DefaultSigma, Bound: DefaultBound}

		sampler, err := NewSampler(tc.prng, tc.ringQ, dist, false)
		require.NoError(t, err)

		noiseBound := uint64(dist.Bound)

		pol := sampler.ReadNew()

		for i := 0; i < N; i++ {
			for j, s := range tc.ringQ.SubRings {
				require.False(t, noiseBound < pol.Coeffs[j][i] && pol.Coeffs[j][i] < (s.Modulus-noiseBound))
			}
		}
	})

	t.Run(testString("Sampler/Gaussian/LargeSigma", tc.ringQ), func(t *testing.T) {

		dist := DiscreteGaussian{Sigma: 1e21, Bound: 1e25}

		sampler, err := NewSampler(tc.prng, tc.ringQ, dist, false)
		require.NoError(t, err)

		pol := sampler.ReadNew()

		require.InDelta(t, math.Log2(1e21), tc.ringQ.Log2OfStandardDeviation(pol), 1)
	})

	for _, p := range []float64{.5, 1. / 3., 128. / 65536.} {
		t.Run(testString(fmt.Sprintf("Sampler/Ternary/p=%1.2f", p), tc.ringQ), func(t *testing.T) {

			sampler, err := NewSampler(tc.prng, tc.ringQ, Ternary{P: p}, false)
			require.NoError(t, err)

			pol := sampler.ReadNew()

			for i, s := range tc.ringQ.SubRings {
				minOne := s.Modulus - 1
				for _, c := range pol.Coeffs[i] {
					require.True(t, c == 0 || c == minOne || c == 1)
				}
			}
		})
	}

	for _, h := range []int{64, 96, 128, 256} {
		t.Run(testString(fmt.Sprintf("Sampler/Ternary/hw=%d", h), tc.ringQ), func(t *testing.T) {

			sampler, err := NewSampler(tc.prng, tc.ringQ, Ternary{H: h}, false)
			require.NoError(t, err)

			checkPoly := func(pol Poly) {
				for i := range tc.ringQ.SubRings {
					hw := 0
					for _, c := range pol.Coeffs[i] {
						if c != 0 {
							hw++
						}
					}

					require.True(t, hw == h)
				}
			}

			pol := sampler.ReadNew()

			checkPoly(pol)

			sampler.Read(pol)

			checkPoly(pol)
		})
	}
}

func testModularReduction(tc *testParams, t *testing.T) {

	t.Run(testString("ModularReduction/BRed", tc.ringQ), func(t *testing.T) {

		var x, y uint64
		var bigQ, result *big.Int

		for j, q := range tc.ringQ.ModuliChain() {

			bigQ = bignum.NewInt(q)

			brc := tc.ringQ.SubRings[j].BRedConstant

			x = 1
			y = 1

			result = bignum.NewInt(x)
			result.Mul(result, bignum.NewInt(y))
			result.Mod(result, bigQ)

			require.Equalf(t, BRed(x, y, q, brc), result.Uint64(), "x = %v, y=%v", x, y)

			x = 1
			y = q - 1

			result = bignum.NewInt(x)
			result.Mul(result, bignum.NewInt(y))
			result.Mod(result, bigQ)

			require.Equalf(t, BRed(x, y, q, brc), result.Uint64(), "x = %v, y=%v", x, y)

			x = 1
			y = 0xFFFFFFFFFFFFFFFF

			result = bignum.NewInt(x)
			result.Mul(result, bignum.NewInt(y))
			result.Mod(result, bigQ)

			require.Equalf(t, BRed(x, y, q, brc), result.Uint64(), "x = %v, y=%v", x, y)

			x = q - 1
			y = q - 1

			result = bignum.NewInt(x)
			result.Mul(result, bignum.NewInt(y))
			result.Mod(result, bigQ)

			require.Equalf(t, BRed(x, y, q, brc), result.Uint64(), "x = %v, y=%v", x, y)

			x = q - 1
			y = 0xFFFFFFFFFFFFFFFF

			result = bignum.NewInt(x)
			result.Mul(result, bignum.NewInt(y))
			result.Mod(result, bigQ)

			require.Equalf(t, BRed(x, y, q, brc), result.Uint64(), "x = %v, y=%v", x, y)

			x = 0xFFFFFFFFFFFFFFFF
			y = 0xFFFFFFFFFFFFFFFF

			result = bignum.NewInt(x)
			result.Mul(result, bignum.NewInt(y))
			result.Mod(result, bigQ)

			require.Equalf(t, BRed(x, y, q, brc), result.Uint64(), "x = %v, y=%v", x, y)
		}
	})

	t.Run(testString("ModularReduction/MRed", tc.ringQ), func(t *testing.T) {

		var x, y uint64
		var bigQ, result *big.Int

		for j, q := range tc.ringQ.ModuliChain() {

			bigQ = bignum.NewInt(q)

			brc := tc.ringQ.SubRings[j].BRedConstant
			mrc := tc.ringQ.SubRings[j].MRedConstant

			x = 1
			y = 1

			result = bignum.NewInt(x)
			result.Mul(result, bignum.NewInt(y))
			result.Mod(result, bigQ)

			require.Equalf(t, MRed(x, MForm(y, q, brc), q, mrc), result.Uint64(), "x = %v, y=%v", x, y)

			x = 1
			y = q - 1

			result = bignum.NewInt(x)
			result.Mul(result, bignum.NewInt(y))
			result.Mod(result, bigQ)

			require.Equalf(t, MRed(x, MForm(y, q, brc), q, mrc), result.Uint64(), "x = %v, y=%v", x, y)

			x = 1
			y = 0xFFFFFFFFFFFFFFFF

			result = bignum.NewInt(x)
			result.Mul(result, bignum.NewInt(y))
			result.Mod(result, bigQ)

			require.Equalf(t, MRed(x, MForm(y, q, brc), q, mrc), result.Uint64(), "x = %v, y=%v", x, y)

			x = q - 1
			y = q - 1

			result = bignum.NewInt(x)
			result.Mul(result, bignum.NewInt(y))
			result.Mod(result, bigQ)

			require.Equalf(t, MRed(x, MForm(y, q, brc), q, mrc), result.Uint64(), "x = %v, y=%v", x, y)

			x = q - 1
			y = 0xFFFFFFFFFFFFFFFF

			result = bignum.NewInt(x)
			result.Mul(result, bignum.NewInt(y))
			result.Mod(result, bigQ)

			require.Equalf(t, MRed(x, MForm(y, q, brc), q, mrc), result.Uint64(), "x = %v, y=%v", x, y)

			x = 0xFFFFFFFFFFFFFFFF
			y = 0xFFFFFFFFFFFFFFFF

			result = bignum.NewInt(x)
			result.Mul(result, bignum.NewInt(y))
			result.Mod(result, bigQ)

			require.Equalf(t, MRed(x, MForm(y, q, brc), q, mrc), result.Uint64(), "x = %v, y=%v", x, y)
		}
	})
}

func testMForm(tc *testParams, t *testing.T) {

	t.Run(testString("MForm", tc.ringQ), func(t *testing.T) {

		polWant := tc.uniformSamplerQ.ReadNew()
		polTest := tc.ringQ.NewPoly()

		tc.ringQ.MForm(polWant, polTest)
		tc.ringQ.IMForm(polTest, polTest)

		require.True(t, tc.ringQ.Equal(polWant, polTest))
	})
}

func testMulScalarBigint(tc *testParams, t *testing.T) {

	t.Run(testString("MulScalarBigint", tc.ringQ), func(t *testing.T) {

		polWant := tc.uniformSamplerQ.ReadNew()
		polTest := *polWant.CopyNew()

		rand1 := RandUniform(tc.prng, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)
		rand2 := RandUniform(tc.prng, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)

		scalarBigint := bignum.NewInt(rand1)
		scalarBigint.Mul(scalarBigint, bignum.NewInt(rand2))

		tc.ringQ.MulScalar(polWant, rand1, polWant)
		tc.ringQ.MulScalar(polWant, rand2, polWant)
		tc.ringQ.MulScalarBigint(polTest, scalarBigint, polTest)

		require.True(t, tc.ringQ.Equal(polWant, polTest))
	})
}

func testExtendBasis(tc *testParams, t *testing.T) {

	N := tc.ringQ.N()

	t.Run(testString("ModUp/QToP", tc.ringQ), func(t *testing.T) {

		prng, _ := sampling.NewPRNG()

		basisextender := NewBasisExtender(tc.ringQ, tc.ringP)

		levelQ := tc.ringQ.Level() - 1
		levelP := tc.ringP.Level() - 1

		ringQ := tc.ringQ.AtLevel(levelQ)
		ringP := tc.ringP.AtLevel(levelP)

		Q := ringQ.Modulus()

		QHalf := new(big.Int).Set(Q)
		QHalf.Rsh(QHalf, 1)

		coeffs := make([]*big.Int, N)
		for i := 0; i < N; i++ {
			coeffs[i] = bignum.RandInt(prng, Q)
			coeffs[i].Sub(coeffs[i], QHalf)
		}

		PolQHave := ringQ.NewPoly()
		PolPTest := ringP.NewPoly()
		PolPWant := ringP.NewPoly()

		ringQ.SetCoefficientsBigint(coeffs, PolQHave)
		ringP.SetCoefficientsBigint(coeffs, PolPWant)

		basisextender.ModUpQtoP(levelQ, levelP, PolQHave, PolPTest)
		ringP.Reduce(PolPTest, PolPTest)

		for i := 0; i < PolPWant.Level()+1; i++ {
			for j := 0; j < N; j++ {
				require.Equal(t, PolPWant.Coeffs[i][j], PolPTest.Coeffs[i][j])
			}
		}
	})

	t.Run(testString("ModUp/PToQ", tc.ringQ), func(t *testing.T) {

		prng, _ := sampling.NewPRNG()

		basisextender := NewBasisExtender(tc.ringQ, tc.ringP)

		levelQ := tc.ringQ.Level() - 1
		levelP := tc.ringP.Level() - 1

		ringQ := tc.ringQ.AtLevel(levelQ)
		ringP := tc.ringP.AtLevel(levelP)

		P := ringP.Modulus()

		PHalf := new(big.Int).Set(P)
		PHalf.Rsh(PHalf, 1)

		coeffs := make([]*big.Int, N)
		for i := 0; i < N; i++ {
			coeffs[i] = bignum.RandInt(prng, P)
			coeffs[i].Sub(coeffs[i], PHalf)
		}

		PolQWant := ringQ.NewPoly()
		PolQTest := ringP.NewPoly()
		PolPHave := ringP.NewPoly()

		ringQ.SetCoefficientsBigint(coeffs, PolQWant)
		ringP.SetCoefficientsBigint(coeffs, PolPHave)

		basisextender.ModUpPtoQ(levelQ, levelP, PolPHave, PolQTest)
		ringQ.Reduce(PolQTest, PolQTest)

		for i := 0; i < PolQWant.Level()+1; i++ {
			for j := 0; j < N; j++ {
				require.Equal(t, PolQWant.Coeffs[i][j], PolQTest.Coeffs[i][j])
			}
		}
	})

	t.Run(testString("ModDown/QPToQ", tc.ringQ), func(t *testing.T) {

		prng, _ := sampling.NewPRNG()

		basisextender := NewBasisExtender(tc.ringQ, tc.ringP)

		levelQ := tc.ringQ.Level() - 1
		levelP := tc.ringP.Level() - 1

		ringQ := tc.ringQ.AtLevel(levelQ)
		ringP := tc.ringP.AtLevel(levelP)

		Q := ringQ.Modulus()
		P := ringP.Modulus()

		QP := new(big.Int).Mul(Q, P)

		coeffs := make([]*big.Int, N)
		for i := 0; i < N; i++ {
			coeffs[i] = bignum.RandInt(prng, QP)
			coeffs[i].Quo(coeffs[i], bignum.NewInt(10))
		}

		coeffsWant := make([]*big.Int, N)
		for i := range coeffs {
			coeffsWant[i] = new(big.Int).Set(coeffs[i])
			bignum.DivRound(coeffsWant[i], P, coeffsWant[i])
		}

		PolQHave := ringQ.NewPoly()
		PolPHave := ringP.NewPoly()
		PolQWant := ringP.NewPoly()

		ringQ.SetCoefficientsBigint(coeffs, PolQHave)
		ringP.SetCoefficientsBigint(coeffs, PolPHave)
		ringQ.SetCoefficientsBigint(coeffsWant, PolQWant)

		basisextender.ModDownQPtoQ(levelQ, levelP, PolQHave, PolPHave, PolQHave)
		ringQ.Reduce(PolQHave, PolQHave)

		for i := 0; i < PolQHave.Level()+1; i++ {
			for j := 0; j < N; j++ {
				require.Equal(t, PolQHave.Coeffs[i][j], PolQWant.Coeffs[i][j])
			}
		}
	})

	t.Run(testString("ModDown/QPToP", tc.ringQ), func(t *testing.T) {

		prng, _ := sampling.NewPRNG()

		basisextender := NewBasisExtender(tc.ringQ, tc.ringP)

		levelQ := tc.ringQ.Level() - 1
		levelP := tc.ringP.Level() - 1

		ringQ := tc.ringQ.AtLevel(levelQ)
		ringP := tc.ringP.AtLevel(levelP)

		Q := ringQ.Modulus()
		P := ringP.Modulus()

		QP := new(big.Int).Mul(Q, P)

		coeffs := make([]*big.Int, N)
		for i := 0; i < N; i++ {
			coeffs[i] = bignum.RandInt(prng, QP)
			coeffs[i].Quo(coeffs[i], bignum.NewInt(10))
		}

		coeffsWant := make([]*big.Int, N)
		for i := range coeffs {
			coeffsWant[i] = new(big.Int).Set(coeffs[i])
			bignum.DivRound(coeffsWant[i], Q, coeffsWant[i])
		}

		PolQHave := ringQ.NewPoly()
		PolPHave := ringP.NewPoly()
		PolPWant := ringP.NewPoly()

		ringQ.SetCoefficientsBigint(coeffs, PolQHave)
		ringP.SetCoefficientsBigint(coeffs, PolPHave)
		ringP.SetCoefficientsBigint(coeffsWant, PolPWant)

		basisextender.ModDownQPtoP(levelQ, levelP, PolQHave, PolPHave, PolPHave)
		ringP.Reduce(PolPHave, PolPHave)

		for i := 0; i < PolQHave.Level()+1; i++ {
			for j := 0; j < N; j++ {
				require.Equal(t, PolPHave.Coeffs[i][j], PolPWant.Coeffs[i][j])
			}
		}
	})
}

func testMultByMonomial(tc *testParams, t *testing.T) {

	t.Run(testString("MultByMonomial", tc.ringQ), func(t *testing.T) {

		p1 := tc.uniformSamplerQ.ReadNew()

		p3Test := tc.ringQ.NewPoly()
		p3Want := tc.ringQ.NewPoly()

		tc.ringQ.MultByMonomial(p1, 1, p3Test)
		tc.ringQ.MultByMonomial(p3Test, 8, p3Test)

		tc.ringQ.MultByMonomial(p1, 9, p3Want)

		require.Equal(t, p3Want.Coeffs[0][:tc.ringQ.N()], p3Test.Coeffs[0][:tc.ringQ.N()])
	})
}

func testShift(t *testing.T) {

	r, _ := NewRing(16, []uint64{97})
	p1, p2 := r.NewPoly(), r.NewPoly()

	for i := range p1.Coeffs[0] {
		p1.Coeffs[0][i] = uint64(i)
	}

	r.Shift(p1, 3, p2)
	require.Equal(t, p2.Coeffs[0], []uint64{3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2})

}
