package ring

import (
	"fmt"
	"math/big"
	"math/bits"
	"math/rand"
	"testing"

	"github.com/ldsec/lattigo/utils"
)

func BenchmarkRing(b *testing.B) {
	b.Run("GenRing", benchGenRing)
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
	b.Run("DivByRNSBasis", benchDivByRNSBasis)
	b.Run("MRed", benchMRed)
	b.Run("BRed", benchBRed)
	b.Run("BRedAdd", benchBRedAdd)

}

func benchGenRing(b *testing.B) {

	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])

		b.Run(testString("", ringQ), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				getRing(parameters[0])
			}

		})
	}
}

func benchMarshalling(b *testing.B) {

	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, ringQ)

		p := uniformSampler.ReadNew()

		b.Run(testString("Marshal/Poly/", ringQ), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				p.MarshalBinary()
			}

		})

		data, _ := p.MarshalBinary()

		b.Run(testString("Unmarshal/Poly/", ringQ), func(b *testing.B) {
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

		ringQ := getRing(parameters[0])

		pol := ringQ.NewPoly()

		b.Run(testString("Gaussian/PRNG/", ringQ), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			gaussianSampler := NewGaussianSampler(prng, ringQ, sigma, bound)

			for i := 0; i < b.N; i++ {
				gaussianSampler.ReadLvl(uint64(len(ringQ.Modulus)-1), pol)
			}
		})

		b.Run(testString("Ternary/0.3/", ringQ), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			ternarySampler := NewTernarySampler(prng, ringQ, 1.0/3, true)

			for i := 0; i < b.N; i++ {
				ternarySampler.Read(pol)
			}
		})

		b.Run(testString("Ternary/0.5/", ringQ), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			ternarySampler := NewTernarySampler(prng, ringQ, 0.5, true)

			for i := 0; i < b.N; i++ {
				ternarySampler.Read(pol)
			}
		})

		b.Run(testString("Ternary/sparse128/", ringQ), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			ternarySampler := NewTernarySamplerSparse(prng, ringQ, 128, true)

			for i := 0; i < b.N; i++ {
				ternarySampler.Read(pol)
			}
		})

		b.Run(testString("Uniform/PRNG/", ringQ), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, ringQ)

			for i := 0; i < b.N; i++ {
				uniformSampler.Read(pol)
			}
		})

	}
}

func benchMontgomeryForm(b *testing.B) {

	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, ringQ)

		p := uniformSampler.ReadNew()

		b.Run(testString("MForm/", ringQ), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				ringQ.MForm(p, p)
			}
		})

		b.Run(testString("InvMForm/", ringQ), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				ringQ.InvMForm(p, p)
			}
		})
	}
}

func benchNTT(b *testing.B) {

	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])

		var NTT func(*Poly, *Poly)
		if ringQ.N == 16384 {
			NTT = ringQ.NTT
		} else {
			NTT = ringQ.NTT
		}

		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, ringQ)

		p := uniformSampler.ReadNew()

		b.Run(testString("NTT/", ringQ), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				NTT(p, p)
			}
		})

		b.Run(testString("InvNTT/", ringQ), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				ringQ.InvNTT(p, p)
			}
		})

		b.Run(testString("NTTBarrett/", ringQ), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				ringQ.NTTBarrett(p, p)
			}
		})

		b.Run(testString("InvNTTBarrett/", ringQ), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				ringQ.InvNTTBarrett(p, p)
			}
		})
	}
}

func benchMulCoeffs(b *testing.B) {

	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, ringQ)

		p0 := uniformSampler.ReadNew()
		p1 := uniformSampler.ReadNew()

		b.Run(testString("Barrett/", ringQ), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				ringQ.MulCoeffs(p0, p1, p0)
			}
		})

		b.Run(testString("BarrettConstant/", ringQ), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				ringQ.MulCoeffsConstant(p0, p1, p0)
			}
		})

		b.Run(testString("Montgomery/", ringQ), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				ringQ.MulCoeffsMontgomery(p0, p1, p0)
			}
		})

		b.Run(testString("MontgomeryConstant/", ringQ), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				ringQ.MulCoeffsMontgomeryConstant(p0, p1, p0)
			}
		})
	}
}

func benchAddCoeffs(b *testing.B) {
	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, ringQ)

		p0 := uniformSampler.ReadNew()
		p1 := uniformSampler.ReadNew()

		b.Run(testString("Add/", ringQ), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				ringQ.Add(p0, p1, p0)
			}
		})

		b.Run(testString("AddConstant/", ringQ), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				ringQ.AddNoMod(p0, p1, p0)
			}
		})
	}
}

func benchSubCoeffs(b *testing.B) {
	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, ringQ)

		p0 := uniformSampler.ReadNew()
		p1 := uniformSampler.ReadNew()

		b.Run(testString("Sub/", ringQ), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				ringQ.Sub(p0, p1, p0)
			}
		})

		b.Run(testString("SubConstant/", ringQ), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				ringQ.SubNoMod(p0, p1, p0)
			}
		})
	}
}

func benchNegCoeffs(b *testing.B) {
	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, ringQ)

		p0 := uniformSampler.ReadNew()

		b.Run(testString("Neg", ringQ), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				ringQ.Neg(p0, p0)
			}
		})
	}
}

func benchMulScalar(b *testing.B) {

	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, ringQ)

		p := uniformSampler.ReadNew()

		rand1 := RandUniform(prng, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)
		rand2 := RandUniform(prng, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)

		scalarBigint := NewUint(rand1)
		scalarBigint.Mul(scalarBigint, NewUint(rand2))

		b.Run(testString("uint64/", ringQ), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				ringQ.MulScalar(p, rand1, p)
			}
		})

		b.Run(testString("big.Int", ringQ), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				ringQ.MulScalarBigint(p, scalarBigint, p)
			}
		})
	}
}

func benchExtendBasis(b *testing.B) {
	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])
		ringP := getRing(parameters[1])

		rescaleParams := make([]uint64, len(ringP.Modulus))

		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSamplerQ := NewUniformSampler(prng, ringQ)
		uniformSamplerP := NewUniformSampler(prng, ringP)

		for i, pi := range ringP.Modulus {
			rescaleParams[i] = RandUniform(prng, pi, (1<<uint64(bits.Len64(pi)-1) - 1))
		}

		basisExtender := NewFastBasisExtender(ringQ, ringP)

		p0 := uniformSamplerQ.ReadNew()
		p1 := uniformSamplerP.ReadNew()

		level := uint64(len(ringQ.Modulus) - 1)

		b.Run(fmt.Sprintf("ModUp/N=%d/limbsQ=%d/limbsP=%d", ringQ.N, len(ringQ.Modulus), len(ringP.Modulus)), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				basisExtender.ModUpSplitQP(level, p0, p1)
			}
		})

		b.Run(fmt.Sprintf("ModDown/N=%d/limbsQ=%d/limbsP=%d", ringQ.N, len(ringQ.Modulus), len(ringP.Modulus)), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				basisExtender.ModDownSplitPQ(level, p0, p1, p0)
			}
		})

		b.Run(fmt.Sprintf("ModDownNTT/N=%d/limbsQ=%d/limbsP=%d", ringQ.N, len(ringQ.Modulus), len(ringP.Modulus)), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				basisExtender.ModDownSplitNTTPQ(level, p0, p1, p0)
			}
		})
	}
}

func benchDivByLastModulus(b *testing.B) {
	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])

		var p0 *Poly

		b.Run(testString("Floor/", ringQ), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, ringQ)

			for i := 0; i < b.N; i++ {

				b.StopTimer()
				p0 = uniformSampler.ReadNew()
				b.StartTimer()

				ringQ.DivFloorByLastModulus(p0)
			}
		})

		b.Run(testString("FloorNTT/", ringQ), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, ringQ)

			for i := 0; i < b.N; i++ {

				b.StopTimer()
				p0 = uniformSampler.ReadNew()
				b.StartTimer()

				ringQ.DivFloorByLastModulusNTT(p0)
			}
		})

		b.Run(testString("Round/", ringQ), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, ringQ)

			for i := 0; i < b.N; i++ {

				b.StopTimer()
				p0 = uniformSampler.ReadNew()
				b.StartTimer()

				ringQ.DivRoundByLastModulus(p0)
			}
		})

		b.Run(testString("RoundNTT/", ringQ), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, ringQ)

			for i := 0; i < b.N; i++ {

				b.StopTimer()
				p0 = uniformSampler.ReadNew()
				b.StartTimer()

				ringQ.DivRoundByLastModulusNTT(p0)
			}
		})
	}
}

func benchDivByRNSBasis(b *testing.B) {

	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])

		b.Run(testString("SimpleScaler/DivByQOverTRounded/reconstructAndScale/", ringQ), func(b *testing.B) {

			rescaler := NewSimpleScaler(testParams.T, ringQ)

			coeffs := make([]*big.Int, ringQ.N)
			for i := uint64(0); i < ringQ.N; i++ {
				coeffs[i] = RandInt(ringQ.ModulusBigint)
			}

			tmp0 := ringQ.NewPoly()
			tmp1 := ringQ.NewPoly()

			ringQ.SetCoefficientsBigint(coeffs, tmp0)

			for i := 0; i < b.N; i++ {
				rescaler.reconstructAndScale(tmp0, tmp1)
			}
		})

		b.Run(testString("SimpleScaler/DivByQOverTRounded/reconstructThenScale/", ringQ), func(b *testing.B) {

			rescaler := NewSimpleScaler(testParams.T, ringQ)

			coeffs := make([]*big.Int, ringQ.N)
			for i := uint64(0); i < ringQ.N; i++ {
				coeffs[i] = RandInt(ringQ.ModulusBigint)
			}

			tmp0 := ringQ.NewPoly()
			tmp1 := ringQ.NewPoly()

			ringQ.SetCoefficientsBigint(coeffs, tmp0)

			for i := 0; i < b.N; i++ {
				rescaler.reconstructThenScale(tmp0, tmp1)
			}
		})

		b.Run(testString("RNSScaler/DivByQOverTRounded/", ringQ), func(b *testing.B) {

			coeffs := make([]*big.Int, ringQ.N)
			for i := uint64(0); i < ringQ.N; i++ {
				coeffs[i] = RandInt(ringQ.ModulusBigint)
			}

			scaler := NewRNSScaler(testParams.T, ringQ)
			polyQ := ringQ.NewPoly()
			polyT := NewPoly(ringQ.N, 1)
			ringQ.SetCoefficientsBigint(coeffs, polyQ)

			for i := 0; i < b.N; i++ {
				scaler.DivByQOverTRounded(polyQ, polyT)
			}
		})
	}
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
