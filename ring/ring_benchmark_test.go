package ring

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

func BenchmarkRing(b *testing.B) {

	var err error

	for _, params := range testParameters[:] {

		var tc *testParams
		if tc, err = genTestParams(params); err != nil {
			b.Fatal(err)
		}

		benchGenRing(tc, b)
		benchMarshalling(tc, b)
		benchSampling(tc, b)
		benchMontgomery(tc, b)
		benchMulCoeffs(tc, b)
		benchAddCoeffs(tc, b)
		benchSubCoeffs(tc, b)
		benchNegCoeffs(tc, b)
		benchMulScalar(tc, b)
		benchExtendBasis(tc, b)
		benchDivByLastModulus(tc, b)
		benchMRed(tc, b)
		benchBRed(tc, b)
		benchBRedAdd(tc, b)
	}
}

func benchGenRing(tc *testParams, b *testing.B) {

	b.Run(testString("GenRing", tc.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			if _, err := NewRing(tc.ringQ.N(), tc.ringQ.ModuliChain()); err != nil {
				b.Error(err)
			}
		}
	})
}

func benchMarshalling(tc *testParams, b *testing.B) {

	var err error

	p := tc.uniformSamplerQ.ReadNew()

	b.Run(testString("Marshalling/MarshalPoly", tc.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			if _, err = p.MarshalBinary(); err != nil {
				b.Error(err)
			}
		}
	})

	var data []byte
	if data, err = p.MarshalBinary(); err != nil {
		b.Error(err)
	}

	b.Run(testString("Marshalling/UnmarshalPoly", tc.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			if err = p.UnmarshalBinary(data); err != nil {
				b.Error(err)
			}
		}
	})
}

func benchSampling(tc *testParams, b *testing.B) {

	pol := tc.ringQ.NewPoly()

	b.Run(testString("Sampling/Gaussian", tc.ringQ), func(b *testing.B) {

		sampler, err := NewSampler(tc.prng, tc.ringQ, DiscreteGaussian{Sigma: DefaultSigma, Bound: DefaultBound}, false)
		require.NoError(b, err)

		for i := 0; i < b.N; i++ {
			sampler.Read(pol)
		}
	})

	b.Run(testString("Sampling/Ternary/0.3", tc.ringQ), func(b *testing.B) {

		sampler, err := NewSampler(tc.prng, tc.ringQ, Ternary{P: 1.0 / 3}, true)
		require.NoError(b, err)

		for i := 0; i < b.N; i++ {
			sampler.Read(pol)
		}
	})

	b.Run(testString("Sampling/Ternary/0.5", tc.ringQ), func(b *testing.B) {

		sampler, err := NewSampler(tc.prng, tc.ringQ, Ternary{P: 0.5}, true)
		require.NoError(b, err)

		for i := 0; i < b.N; i++ {
			sampler.Read(pol)
		}
	})

	b.Run(testString("Sampling/Ternary/sparse128", tc.ringQ), func(b *testing.B) {

		sampler, err := NewSampler(tc.prng, tc.ringQ, Ternary{H: 128}, true)
		require.NoError(b, err)

		for i := 0; i < b.N; i++ {
			sampler.Read(pol)
		}
	})

	b.Run(testString("Sampling/Uniform", tc.ringQ), func(b *testing.B) {

		sampler, err := NewSampler(tc.prng, tc.ringQ, Uniform{}, true)
		require.NoError(b, err)

		for i := 0; i < b.N; i++ {
			sampler.Read(pol)
		}
	})
}

func benchMontgomery(tc *testParams, b *testing.B) {

	p := tc.uniformSamplerQ.ReadNew()

	b.Run(testString("Montgomery/MForm", tc.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			tc.ringQ.MForm(p, p)
		}
	})

	b.Run(testString("Montgomery/InvMForm", tc.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			tc.ringQ.IMForm(p, p)
		}
	})
}

func benchMulCoeffs(tc *testParams, b *testing.B) {

	p0 := tc.uniformSamplerQ.ReadNew()
	p1 := tc.uniformSamplerQ.ReadNew()

	b.Run(testString("MulCoeffs/Montgomery", tc.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			tc.ringQ.MulCoeffsMontgomery(p0, p1, p0)
		}
	})

	b.Run(testString("MulCoeffs/MontgomeryLazy", tc.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			tc.ringQ.MulCoeffsMontgomeryLazy(p0, p1, p0)
		}
	})

	b.Run(testString("MulCoeffs/Barrett", tc.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			tc.ringQ.MulCoeffsBarrett(p0, p1, p0)
		}
	})

	b.Run(testString("MulCoeffs/BarrettLazy", tc.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			tc.ringQ.MulCoeffsBarrettLazy(p0, p1, p0)
		}
	})
}

func benchAddCoeffs(tc *testParams, b *testing.B) {

	p0 := tc.uniformSamplerQ.ReadNew()
	p1 := tc.uniformSamplerQ.ReadNew()

	b.Run(testString("AddCoeffs/Add", tc.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			tc.ringQ.Add(p0, p1, p0)
		}
	})

	b.Run(testString("AddCoeffs/AddLazy", tc.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			tc.ringQ.AddLazy(p0, p1, p0)
		}
	})
}

func benchSubCoeffs(tc *testParams, b *testing.B) {

	p0 := tc.uniformSamplerQ.ReadNew()
	p1 := tc.uniformSamplerQ.ReadNew()

	b.Run(testString("SubCoeffs/Sub", tc.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			tc.ringQ.Sub(p0, p1, p0)
		}
	})

	b.Run(testString("SubCoeffs/SubLazy", tc.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			tc.ringQ.SubLazy(p0, p1, p0)
		}
	})
}

func benchNegCoeffs(tc *testParams, b *testing.B) {

	p0 := tc.uniformSamplerQ.ReadNew()

	b.Run(testString("NegCoeffs", tc.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			tc.ringQ.Neg(p0, p0)
		}
	})
}

func benchMulScalar(tc *testParams, b *testing.B) {

	p := tc.uniformSamplerQ.ReadNew()

	rand1 := RandUniform(tc.prng, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)
	rand2 := RandUniform(tc.prng, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)

	scalarBigint := bignum.NewInt(rand1)
	scalarBigint.Mul(scalarBigint, bignum.NewInt(rand2))

	b.Run(testString("MulScalar/uint64", tc.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			tc.ringQ.MulScalar(p, rand1, p)
		}
	})

	b.Run(testString("MulScalar/big.Int", tc.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			tc.ringQ.MulScalarBigint(p, scalarBigint, p)
		}
	})
}

func benchExtendBasis(tc *testParams, b *testing.B) {

	basisExtender := NewBasisExtender(tc.ringQ, tc.ringP)

	p0 := tc.uniformSamplerQ.ReadNew()
	p1 := tc.uniformSamplerP.ReadNew()

	levelQ := tc.ringQ.MaxLevel()
	levelP := tc.ringP.MaxLevel()

	b.Run(fmt.Sprintf("ExtendBasis/ModUp/N=%d/limbsQ=%d/limbsP=%d", tc.ringQ.N(), tc.ringQ.ModuliChainLength(), tc.ringP.ModuliChainLength()), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			basisExtender.ModUpQtoP(levelQ, levelP, p0, p1)
		}
	})

	b.Run(fmt.Sprintf("ExtendBasis/ModDown/N=%d/limbsQ=%d/limbsP=%d", tc.ringQ.N(), tc.ringQ.ModuliChainLength(), tc.ringP.ModuliChainLength()), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			basisExtender.ModDownQPtoQ(levelQ, levelP, p0, p1, p0)
		}
	})

	b.Run(fmt.Sprintf("ExtendBasis/ModDownNTT/N=%d/limbsQ=%d/limbsP=%d", tc.ringQ.N(), tc.ringQ.ModuliChainLength(), tc.ringP.ModuliChainLength()), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			basisExtender.ModDownQPtoQNTT(levelQ, levelP, p0, p1, p0)
		}
	})
}

func benchDivByLastModulus(tc *testParams, b *testing.B) {

	p0 := tc.uniformSamplerQ.ReadNew()
	p1 := tc.ringQ.AtLevel(p0.Level() - 1).NewPoly()

	buff := tc.ringQ.NewPoly()

	b.Run(testString("DivByLastModulus/Floor", tc.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			tc.ringQ.DivFloorByLastModulus(p0, p1)
		}
	})

	b.Run(testString("DivByLastModulus/FloorNTT", tc.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			tc.ringQ.DivFloorByLastModulusNTT(p0, buff, p1)
		}
	})

	b.Run(testString("DivByLastModulus/Round", tc.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			tc.ringQ.DivRoundByLastModulus(p0, p1)
		}
	})

	b.Run(testString("DivByLastModulus/RoundNTT", tc.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			tc.ringQ.DivRoundByLastModulusNTT(p0, buff, p1)
		}
	})
}

func benchBRed(tc *testParams, b *testing.B) {

	var q, x, y uint64 = 1033576114481528833, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF

	brc := GenBRedConstant(q)

	b.ResetTimer()

	b.Run("BRed", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			x = BRed(x, y, q, brc)
		}
	})
}

func benchMRed(tc *testParams, b *testing.B) {

	var q, x, y uint64 = 1033576114481528833, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF

	y = MForm(y, q, GenBRedConstant(q))

	mrc := GenMRedConstant(q)

	b.ResetTimer()

	b.Run("MRed", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			x = MRed(x, y, q, mrc)
		}
	})
}

func benchBRedAdd(tc *testParams, b *testing.B) {

	var q, x uint64 = 1033576114481528833, 0xFFFFFFFFFFFFFFFF

	brc := GenBRedConstant(q)

	b.ResetTimer()

	b.Run("BRedAdd", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			BRedAdd(x, q, brc)
		}
	})
}
