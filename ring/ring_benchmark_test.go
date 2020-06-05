package ring

import (
	"fmt"
	"github.com/ldsec/lattigo/utils"
	"math/bits"
	"math/rand"
	"testing"
)

func BenchmarkRing(b *testing.B) {
	b.Run("GenRingContext", benchGenRingContext)
	b.Run("Marshalling", benchMarshalling)
	b.Run("Sampling", benchSampling)
	b.Run("Montgomery", benchMontgomeryForm)
	b.Run("NTT", benchNTT)
	b.Run("MulCoeffs", benchMulCoeffs)
	b.Run("AddCoeffs", benchAddCoeffs)
	b.Run("SubCoeffs", benchSubCoeffs)
	b.Run("NegCoeffs", benchNegCoeffs)
	b.Run("MulScalar", benchMulScalar)
	b.Run("ExtendBasis", benchExtendBasis)
	b.Run("DivByLastModulus", benchDivByLastModulus)
	b.Run("MRed", benchMRed)
	b.Run("BRed", benchBRed)
	b.Run("BRedAdd", benchBRedAdd)

}

func benchGenRingContext(b *testing.B) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		b.Run(testString("", context), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				genPolyContext(parameters[0])
			}

		})
	}
}

func benchMarshalling(b *testing.B) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, context)

		p := uniformSampler.SampleNew()

		b.Run(testString("Marshal/Poly/", context), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				p.MarshalBinary()
			}

		})

		data, _ := p.MarshalBinary()

		b.Run(testString("Unmarshal/Poly/", context), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				p.UnmarshalBinary(data)
			}
		})
	}
}

func benchSampling(b *testing.B) {

	sigma := 3.19

	bound := uint64(sigma * 6)

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		pol := context.NewPoly()

		b.Run(testString("Gaussian/PRNG/", context), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			gaussianSampler := NewGaussianSampler(prng, context)

			for i := 0; i < b.N; i++ {
				gaussianSampler.SampleLvl(uint64(len(context.Modulus)-1), pol, sigma, bound)
			}
		})

		b.Run(testString("Ternary/0.3/", context), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			ternarySampler := NewTernarySampler(prng, context)

			for i := 0; i < b.N; i++ {
				ternarySampler.Sample(pol, 1.0/3)
			}
		})

		b.Run(testString("Ternary/0.5/", context), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			ternarySampler := NewTernarySampler(prng, context)

			for i := 0; i < b.N; i++ {
				ternarySampler.Sample(pol, 0.5)
			}
		})

		b.Run(testString("Ternary/sparse128/", context), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			ternarySampler := NewTernarySampler(prng, context)

			for i := 0; i < b.N; i++ {
				ternarySampler.SampleSparse(pol, 128)
			}
		})

		b.Run(testString("Uniform/PRNG/", context), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, context)

			for i := 0; i < b.N; i++ {
				uniformSampler.Sample(pol)
			}
		})

	}
}

func benchMontgomeryForm(b *testing.B) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, context)

		p := uniformSampler.SampleNew()

		b.Run(testString("MForm/", context), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				context.MForm(p, p)
			}
		})

		b.Run(testString("InvMForm/", context), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				context.InvMForm(p, p)
			}
		})
	}
}

func benchNTT(b *testing.B) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, context)

		p := uniformSampler.SampleNew()

		b.Run(testString("NTT/", context), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				context.NTT(p, p)
			}
		})

		b.Run(testString("InvNTT/", context), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				context.InvNTT(p, p)
			}
		})

		b.Run(testString("NTTBarrett/", context), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				context.NTTBarrett(p, p)
			}
		})

		b.Run(testString("InvNTTBarrett/", context), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				context.InvNTTBarrett(p, p)
			}
		})
	}
}

func benchMulCoeffs(b *testing.B) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, context)

		p0 := uniformSampler.SampleNew()
		p1 := uniformSampler.SampleNew()

		b.Run(testString("Barrett/", context), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				context.MulCoeffs(p0, p1, p0)
			}
		})

		b.Run(testString("BarrettConstant/", context), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				context.MulCoeffsConstant(p0, p1, p0)
			}
		})

		b.Run(testString("Montgomery/", context), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				context.MulCoeffsMontgomery(p0, p1, p0)
			}
		})

		b.Run(testString("MontgomeryConstant/", context), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				context.MulCoeffsMontgomeryConstant(p0, p1, p0)
			}
		})
	}
}

func benchAddCoeffs(b *testing.B) {
	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, context)

		p0 := uniformSampler.SampleNew()
		p1 := uniformSampler.SampleNew()

		b.Run(testString("Add/", context), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				context.Add(p0, p1, p0)
			}
		})

		b.Run(testString("AddConstant/", context), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				context.AddNoMod(p0, p1, p0)
			}
		})
	}
}

func benchSubCoeffs(b *testing.B) {
	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, context)

		p0 := uniformSampler.SampleNew()
		p1 := uniformSampler.SampleNew()

		b.Run(testString("Sub/", context), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				context.Sub(p0, p1, p0)
			}
		})

		b.Run(testString("SubConstant/", context), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				context.SubNoMod(p0, p1, p0)
			}
		})
	}
}

func benchNegCoeffs(b *testing.B) {
	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, context)

		p0 := uniformSampler.SampleNew()

		b.Run(testString("Neg", context), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				context.Neg(p0, p0)
			}
		})
	}
}

func benchMulScalar(b *testing.B) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, context)

		p := uniformSampler.SampleNew()

		rand1 := RandUniform(prng, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)
		rand2 := RandUniform(prng, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)

		scalarBigint := NewUint(rand1)
		scalarBigint.Mul(scalarBigint, NewUint(rand2))

		b.Run(testString("uint64/", context), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				context.MulScalar(p, rand1, p)
			}
		})

		b.Run(testString("big.Int", context), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				context.MulScalarBigint(p, scalarBigint, p)
			}
		})
	}
}

func benchExtendBasis(b *testing.B) {
	for _, parameters := range testParams.polyParams {

		contextQ := genPolyContext(parameters[0])
		contextP := genPolyContext(parameters[1])

		rescaleParams := make([]uint64, len(contextP.Modulus))

		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSamplerQ := NewUniformSampler(prng, contextQ)
		uniformSamplerP := NewUniformSampler(prng, contextP)

		for i, pi := range contextP.Modulus {
			rescaleParams[i] = RandUniform(prng, pi, (1<<uint64(bits.Len64(pi)-1) - 1))
		}

		basisExtender := NewFastBasisExtender(contextQ, contextP)

		p0 := uniformSamplerQ.SampleNew()
		p1 := uniformSamplerP.SampleNew()

		level := uint64(len(contextQ.Modulus) - 1)

		b.Run(fmt.Sprintf("ModUp/N=%d/limbsQ=%d/limbsP=%d", contextQ.N, len(contextQ.Modulus), len(contextP.Modulus)), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				basisExtender.ModUpSplitQP(level, p0, p1)
			}
		})

		b.Run(fmt.Sprintf("ModDown/N=%d/limbsQ=%d/limbsP=%d", contextQ.N, len(contextQ.Modulus), len(contextP.Modulus)), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				basisExtender.ModDownSplitedPQ(level, p0, p1, p0)
			}
		})

		b.Run(fmt.Sprintf("ModDownNTT/N=%d/limbsQ=%d/limbsP=%d", contextQ.N, len(contextQ.Modulus), len(contextP.Modulus)), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				basisExtender.ModDownSplitedNTTPQ(level, p0, p1, p0)
			}
		})
	}
}

func benchDivByLastModulus(b *testing.B) {
	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		var p0 *Poly

		b.Run(testString("Floor/", context), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, context)

			for i := 0; i < b.N; i++ {

				b.StopTimer()
				p0 = uniformSampler.SampleNew()
				b.StartTimer()

				context.DivFloorByLastModulus(p0)
			}
		})

		b.Run(testString("FloorNTT/", context), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, context)

			for i := 0; i < b.N; i++ {

				b.StopTimer()
				p0 = uniformSampler.SampleNew()
				b.StartTimer()

				context.DivFloorByLastModulusNTT(p0)
			}
		})

		b.Run(testString("Round/", context), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, context)

			for i := 0; i < b.N; i++ {

				b.StopTimer()
				p0 = uniformSampler.SampleNew()
				b.StartTimer()

				context.DivRoundByLastModulus(p0)
			}
		})

		b.Run(testString("RoundNTT/", context), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, context)

			for i := 0; i < b.N; i++ {

				b.StopTimer()
				p0 = uniformSampler.SampleNew()
				b.StartTimer()

				context.DivRoundByLastModulusNTT(p0)
			}
		})
	}
}

func benchMRed(b *testing.B) {

	q := uint64(1033576114481528833)
	u := BRedParams(q)

	x := rand.Uint64() % q
	y := rand.Uint64() % q

	b.ResetTimer()

	b.Run(fmt.Sprintf("BRed"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			x = BRed(x, y, q, u)
		}
	})
}

func benchBRed(b *testing.B) {

	q := uint64(1033576114481528833)

	x := rand.Uint64() % q
	y := rand.Uint64() % q

	u := BRedParams(q)

	y = MForm(y, q, u)

	m := MRedParams(q)

	b.ResetTimer()

	b.Run(fmt.Sprintf("MRed"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			x = MRed(x, y, q, m)
		}
	})
}

func benchBRedAdd(b *testing.B) {

	q := uint64(1033576114481528833)
	u := BRedParams(q)

	x := rand.Uint64()

	b.ResetTimer()

	b.Run(fmt.Sprintf("BRedAdd"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			BRedAdd(x, q, u)
		}
	})
}
