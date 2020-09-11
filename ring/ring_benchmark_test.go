package ring

import (
	"fmt"
	"math/big"
	"math/bits"
	"math/rand"
	"testing"
)

func BenchmarkRing(b *testing.B) {

	var defaultParams []*Parameters

	if testing.Short() {
		defaultParams = DefaultParams[:3]
	} else {
		defaultParams = DefaultParams
	}

	for _, defaultParam := range defaultParams {

		if err := GenTestParams(defaultParam); err != nil {
			panic(err)
		}

		benchGenRing(b)
		benchMarshalling(b)
		benchSampling(b)
		benchMontgomery(b)
		benchNTT(b)
		benchMulCoeffs(b)
		benchAddCoeffs(b)
		benchSubCoeffs(b)
		benchNegCoeffs(b)
		benchMulScalar(b)
		benchExtendBasis(b)
		benchDivByLastModulus(b)
		benchDivByRNSBasis(b)
		benchMRed(b)
		benchBRed(b)
		benchBRedAdd(b)
	}
}

func benchGenRing(b *testing.B) {

	b.Run(testString("GenRing/", params.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			NewRing(params.ringQ.N, params.ringQ.Modulus)
		}
	})
}

func benchMarshalling(b *testing.B) {

	p := params.uniformSamplerQ.ReadNew()

	b.Run(testString("Marshalling/MarshalPoly/", params.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.MarshalBinary()
		}
	})

	data, _ := p.MarshalBinary()

	b.Run(testString("Marshalling/UnmarshalPoly/", params.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.UnmarshalBinary(data)
		}
	})
}

func benchSampling(b *testing.B) {

	pol := params.ringQ.NewPoly()

	b.Run(testString("Sampling/Gaussian/", params.ringQ), func(b *testing.B) {

		gaussianSampler := NewGaussianSampler(params.prng, params.ringQ, DefaultSigma, DefaultBound)

		for i := 0; i < b.N; i++ {
			gaussianSampler.ReadLvl(uint64(len(params.ringQ.Modulus)-1), pol)
		}
	})

	b.Run(testString("Sampling/Ternary/0.3/", params.ringQ), func(b *testing.B) {

		ternarySampler := NewTernarySampler(params.prng, params.ringQ, 1.0/3, true)

		for i := 0; i < b.N; i++ {
			ternarySampler.Read(pol)
		}
	})

	b.Run(testString("Sampling/Ternary/0.5/", params.ringQ), func(b *testing.B) {

		ternarySampler := NewTernarySampler(params.prng, params.ringQ, 0.5, true)

		for i := 0; i < b.N; i++ {
			ternarySampler.Read(pol)
		}
	})

	b.Run(testString("Sampling/Ternary/sparse128/", params.ringQ), func(b *testing.B) {

		ternarySampler := NewTernarySamplerSparse(params.prng, params.ringQ, 128, true)

		for i := 0; i < b.N; i++ {
			ternarySampler.Read(pol)
		}
	})

	b.Run(testString("Sampling/Uniform/", params.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			params.uniformSamplerQ.Read(pol)
		}
	})
}

func benchMontgomery(b *testing.B) {

	p := params.uniformSamplerQ.ReadNew()

	b.Run(testString("Montgomery/MForm/", params.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			params.ringQ.MForm(p, p)
		}
	})

	b.Run(testString("Montgomery/InvMForm/", params.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			params.ringQ.InvMForm(p, p)
		}
	})
}

func benchNTT(b *testing.B) {

	p := params.uniformSamplerQ.ReadNew()

	b.Run(testString("NTT/NTT/Montgomery/", params.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			params.ringQ.NTT(p, p)
		}
	})

	b.Run(testString("NTT/InvNTT/Montgomery/", params.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			params.ringQ.InvNTT(p, p)
		}
	})

	b.Run(testString("NTT/NTT/Barrett/", params.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			params.ringQ.NTTBarrett(p, p)
		}
	})

	b.Run(testString("NTT/InvNTT/Barrett/", params.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			params.ringQ.InvNTTBarrett(p, p)
		}
	})
}

func benchMulCoeffs(b *testing.B) {

	p0 := params.uniformSamplerQ.ReadNew()
	p1 := params.uniformSamplerQ.ReadNew()

	b.Run(testString("MulCoeffs/Montgomery/", params.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			params.ringQ.MulCoeffsMontgomery(p0, p1, p0)
		}
	})

	b.Run(testString("MulCoeffs/MontgomeryConstant/", params.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			params.ringQ.MulCoeffsMontgomeryConstant(p0, p1, p0)
		}
	})

	b.Run(testString("MulCoeffs/Barrett/", params.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			params.ringQ.MulCoeffs(p0, p1, p0)
		}
	})

	b.Run(testString("MulCoeffs/BarrettConstant/", params.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			params.ringQ.MulCoeffsConstant(p0, p1, p0)
		}
	})
}

func benchAddCoeffs(b *testing.B) {

	p0 := params.uniformSamplerQ.ReadNew()
	p1 := params.uniformSamplerQ.ReadNew()

	b.Run(testString("AddCoeffs/Add/", params.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			params.ringQ.Add(p0, p1, p0)
		}
	})

	b.Run(testString("AddCoeffs/AddConstant/", params.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			params.ringQ.AddNoMod(p0, p1, p0)
		}
	})
}

func benchSubCoeffs(b *testing.B) {

	p0 := params.uniformSamplerQ.ReadNew()
	p1 := params.uniformSamplerQ.ReadNew()

	b.Run(testString("SubCoeffs/Sub/", params.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			params.ringQ.Sub(p0, p1, p0)
		}
	})

	b.Run(testString("SubCoeffs/SubConstant/", params.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			params.ringQ.SubNoMod(p0, p1, p0)
		}
	})
}

func benchNegCoeffs(b *testing.B) {

	p0 := params.uniformSamplerQ.ReadNew()

	b.Run(testString("NegCoeffs", params.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			params.ringQ.Neg(p0, p0)
		}
	})
}

func benchMulScalar(b *testing.B) {

	p := params.uniformSamplerQ.ReadNew()

	rand1 := RandUniform(params.prng, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)
	rand2 := RandUniform(params.prng, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)

	scalarBigint := NewUint(rand1)
	scalarBigint.Mul(scalarBigint, NewUint(rand2))

	b.Run(testString("MulScalar/uint64/", params.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			params.ringQ.MulScalar(p, rand1, p)
		}
	})

	b.Run(testString("MulScalar/big.Int/", params.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			params.ringQ.MulScalarBigint(p, scalarBigint, p)
		}
	})
}

func benchExtendBasis(b *testing.B) {

	rescaleParams := make([]uint64, len(params.ringP.Modulus))

	for i, pi := range params.ringP.Modulus {
		rescaleParams[i] = RandUniform(params.prng, pi, (1<<uint64(bits.Len64(pi)-1) - 1))
	}

	basisExtender := NewFastBasisExtender(params.ringQ, params.ringP)

	p0 := params.uniformSamplerQ.ReadNew()
	p1 := params.uniformSamplerP.ReadNew()

	level := uint64(len(params.ringQ.Modulus) - 1)

	b.Run(fmt.Sprintf("ExtendBasis/ModUp/N=%d/limbsQ=%d/limbsP=%d", params.ringQ.N, len(params.ringQ.Modulus), len(params.ringP.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			basisExtender.ModUpSplitQP(level, p0, p1)
		}
	})

	b.Run(fmt.Sprintf("ExtendBasis/ModDown/N=%d/limbsQ=%d/limbsP=%d", params.ringQ.N, len(params.ringQ.Modulus), len(params.ringP.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			basisExtender.ModDownSplitPQ(level, p0, p1, p0)
		}
	})

	b.Run(fmt.Sprintf("ExtendBasis/ModDownNTT/N=%d/limbsQ=%d/limbsP=%d", params.ringQ.N, len(params.ringQ.Modulus), len(params.ringP.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			basisExtender.ModDownSplitNTTPQ(level, p0, p1, p0)
		}
	})
}

func benchDivByLastModulus(b *testing.B) {

	var p0 *Poly

	b.Run(testString("DivByLastModulus/Floor/", params.ringQ), func(b *testing.B) {

		for i := 0; i < b.N; i++ {

			b.StopTimer()
			p0 = params.uniformSamplerQ.ReadNew()
			b.StartTimer()

			params.ringQ.DivFloorByLastModulus(p0)
		}
	})

	b.Run(testString("DivByLastModulus/FloorNTT/", params.ringQ), func(b *testing.B) {

		for i := 0; i < b.N; i++ {

			b.StopTimer()
			p0 = params.uniformSamplerQ.ReadNew()
			b.StartTimer()

			params.ringQ.DivFloorByLastModulusNTT(p0)
		}
	})

	b.Run(testString("DivByLastModulus/Round/", params.ringQ), func(b *testing.B) {

		for i := 0; i < b.N; i++ {

			b.StopTimer()
			p0 = params.uniformSamplerQ.ReadNew()
			b.StartTimer()

			params.ringQ.DivRoundByLastModulus(p0)
		}
	})

	b.Run(testString("DivByLastModulus/RoundNTT/", params.ringQ), func(b *testing.B) {

		for i := 0; i < b.N; i++ {

			b.StopTimer()
			p0 = params.uniformSamplerQ.ReadNew()
			b.StartTimer()

			params.ringQ.DivRoundByLastModulusNTT(p0)
		}
	})
}

func benchDivByRNSBasis(b *testing.B) {

	b.Run(testString("DivByRNSBasis/Simple/DivByQOverTRounded/reconstructAndScale/", params.ringQ), func(b *testing.B) {

		rescaler := NewSimpleScaler(T, params.ringQ)

		coeffs := make([]*big.Int, params.ringQ.N)
		for i := uint64(0); i < params.ringQ.N; i++ {
			coeffs[i] = RandInt(params.ringQ.ModulusBigint)
		}

		tmp0 := params.ringQ.NewPoly()
		tmp1 := params.ringQ.NewPoly()

		params.ringQ.SetCoefficientsBigint(coeffs, tmp0)

		for i := 0; i < b.N; i++ {
			rescaler.reconstructAndScale(tmp0, tmp1)
		}
	})

	b.Run(testString("DivByRNSBasis/Simple/DivByQOverTRounded/reconstructThenScale/", params.ringQ), func(b *testing.B) {

		rescaler := NewSimpleScaler(T, params.ringQ)

		coeffs := make([]*big.Int, params.ringQ.N)
		for i := uint64(0); i < params.ringQ.N; i++ {
			coeffs[i] = RandInt(params.ringQ.ModulusBigint)
		}

		tmp0 := params.ringQ.NewPoly()
		tmp1 := params.ringQ.NewPoly()

		params.ringQ.SetCoefficientsBigint(coeffs, tmp0)

		for i := 0; i < b.N; i++ {
			rescaler.reconstructThenScale(tmp0, tmp1)
		}
	})

	b.Run(testString("DivByRNSBasis/RNS/DivByQOverTRounded/", params.ringQ), func(b *testing.B) {

		coeffs := make([]*big.Int, params.ringQ.N)
		for i := uint64(0); i < params.ringQ.N; i++ {
			coeffs[i] = RandInt(params.ringQ.ModulusBigint)
		}

		scaler := NewRNSScaler(T, params.ringQ)
		polyQ := params.ringQ.NewPoly()
		polyT := NewPoly(params.ringQ.N, 1)
		params.ringQ.SetCoefficientsBigint(coeffs, polyQ)

		for i := 0; i < b.N; i++ {
			scaler.DivByQOverTRounded(polyQ, polyT)
		}
	})
}

func benchBRed(b *testing.B) {

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

func benchMRed(b *testing.B) {

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
