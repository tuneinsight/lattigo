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
	b.Run("DivByRNSBasis", benchDivByRNSBasis)
	b.Run("MRed", benchMRed)
	b.Run("BRed", benchBRed)
	b.Run("BRedAdd", benchBRedAdd)

}

func benchGenRingContext(b *testing.B) {

	for _, parameters := range testParams.polyParams {

		r := getRing(parameters[0])

		b.Run(testString("", r), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				getRing(parameters[0])
			}

		})
	}
}

func benchMarshalling(b *testing.B) {

	for _, parameters := range testParams.polyParams {

		r := getRing(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, r)

		p := uniformSampler.ReadNew()

		b.Run(testString("Marshal/Poly/", r), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				p.MarshalBinary()
			}

		})

		data, _ := p.MarshalBinary()

		b.Run(testString("Unmarshal/Poly/", r), func(b *testing.B) {
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

		r := getRing(parameters[0])

		pol := r.NewPoly()

		b.Run(testString("Gaussian/PRNG/", r), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			gaussianSampler := NewGaussianSampler(prng, r, sigma, bound)

			for i := 0; i < b.N; i++ {
				gaussianSampler.ReadLvl(uint64(len(r.Modulus)-1), pol)
			}
		})

		b.Run(testString("Ternary/0.3/", r), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			ternarySampler := NewTernarySampler(prng, r, 1.0/3, true)

			for i := 0; i < b.N; i++ {
				ternarySampler.Read(pol)
			}
		})

		b.Run(testString("Ternary/0.5/", r), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			ternarySampler := NewTernarySampler(prng, r, 0.5, true)

			for i := 0; i < b.N; i++ {
				ternarySampler.Read(pol)
			}
		})

		b.Run(testString("Ternary/sparse128/", r), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			ternarySampler := NewTernarySamplerSparse(prng, r, 128, true)

			for i := 0; i < b.N; i++ {
				ternarySampler.Read(pol)
			}
		})

		b.Run(testString("Uniform/PRNG/", r), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, r)

			for i := 0; i < b.N; i++ {
				uniformSampler.Read(pol)
			}
		})

	}
}

func benchMontgomeryForm(b *testing.B) {

	for _, parameters := range testParams.polyParams {

		r := getRing(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, r)

		p := uniformSampler.ReadNew()

		b.Run(testString("MForm/", r), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				r.MForm(p, p)
			}
		})

		b.Run(testString("InvMForm/", r), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				r.InvMForm(p, p)
			}
		})
	}
}

func benchNTT(b *testing.B) {

	for _, parameters := range testParams.polyParams {

		r := getRing(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, r)

		p := uniformSampler.ReadNew()

		b.Run(testString("NTT/", r), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				r.NTT(p, p)
			}
		})

		b.Run(testString("InvNTT/", r), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				r.InvNTT(p, p)
			}
		})

		b.Run(testString("NTTBarrett/", r), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				r.NTTBarrett(p, p)
			}
		})

		b.Run(testString("InvNTTBarrett/", r), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				r.InvNTTBarrett(p, p)
			}
		})
	}
}

func benchMulCoeffs(b *testing.B) {

	for _, parameters := range testParams.polyParams {

		r := getRing(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, r)

		p0 := uniformSampler.ReadNew()
		p1 := uniformSampler.ReadNew()

		b.Run(testString("Barrett/", r), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				r.MulCoeffs(p0, p1, p0)
			}
		})

		b.Run(testString("BarrettConstant/", r), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				r.MulCoeffsConstant(p0, p1, p0)
			}
		})

		b.Run(testString("Montgomery/", r), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				r.MulCoeffsMontgomery(p0, p1, p0)
			}
		})

		b.Run(testString("MontgomeryConstant/", r), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				r.MulCoeffsMontgomeryConstant(p0, p1, p0)
			}
		})
	}
}

func benchAddCoeffs(b *testing.B) {
	for _, parameters := range testParams.polyParams {

		r := getRing(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, r)

		p0 := uniformSampler.ReadNew()
		p1 := uniformSampler.ReadNew()

		b.Run(testString("Add/", r), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				r.Add(p0, p1, p0)
			}
		})

		b.Run(testString("AddConstant/", r), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				r.AddNoMod(p0, p1, p0)
			}
		})
	}
}

func benchSubCoeffs(b *testing.B) {
	for _, parameters := range testParams.polyParams {

		r := getRing(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, r)

		p0 := uniformSampler.ReadNew()
		p1 := uniformSampler.ReadNew()

		b.Run(testString("Sub/", r), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				r.Sub(p0, p1, p0)
			}
		})

		b.Run(testString("SubConstant/", r), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				r.SubNoMod(p0, p1, p0)
			}
		})
	}
}

func benchNegCoeffs(b *testing.B) {
	for _, parameters := range testParams.polyParams {

		r := getRing(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, r)

		p0 := uniformSampler.ReadNew()

		b.Run(testString("Neg", r), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				r.Neg(p0, p0)
			}
		})
	}
}

func benchMulScalar(b *testing.B) {

	for _, parameters := range testParams.polyParams {

		r := getRing(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, r)

		p := uniformSampler.ReadNew()

		rand1 := RandUniform(prng, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)
		rand2 := RandUniform(prng, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)

		scalarBigint := NewUint(rand1)
		scalarBigint.Mul(scalarBigint, NewUint(rand2))

		b.Run(testString("uint64/", r), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				r.MulScalar(p, rand1, p)
			}
		})

		b.Run(testString("big.Int", r), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				r.MulScalarBigint(p, scalarBigint, p)
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
				basisExtender.ModDownSplitedPQ(level, p0, p1, p0)
			}
		})

		b.Run(fmt.Sprintf("ModDownNTT/N=%d/limbsQ=%d/limbsP=%d", ringQ.N, len(ringQ.Modulus), len(ringP.Modulus)), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				basisExtender.ModDownSplitedNTTPQ(level, p0, p1, p0)
			}
		})
	}
}

func benchDivByLastModulus(b *testing.B) {
	for _, parameters := range testParams.polyParams {

		r := getRing(parameters[0])

		var p0 *Poly

		b.Run(testString("Floor/", r), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, r)

			for i := 0; i < b.N; i++ {

				b.StopTimer()
				p0 = uniformSampler.ReadNew()
				b.StartTimer()

				r.DivFloorByLastModulus(p0)
			}
		})

		b.Run(testString("FloorNTT/", r), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, r)

			for i := 0; i < b.N; i++ {

				b.StopTimer()
				p0 = uniformSampler.ReadNew()
				b.StartTimer()

				r.DivFloorByLastModulusNTT(p0)
			}
		})

		b.Run(testString("Round/", r), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, r)

			for i := 0; i < b.N; i++ {

				b.StopTimer()
				p0 = uniformSampler.ReadNew()
				b.StartTimer()

				r.DivRoundByLastModulus(p0)
			}
		})

		b.Run(testString("RoundNTT/", r), func(b *testing.B) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, r)

			for i := 0; i < b.N; i++ {

				b.StopTimer()
				p0 = uniformSampler.ReadNew()
				b.StartTimer()

				r.DivRoundByLastModulusNTT(p0)
			}
		})
	}
}

func benchDivByRNSBasis(b *testing.B) {

	for _, parameters := range testParams.polyParams {

		r := getRing(parameters[0])

		b.Run(testString("SimpleScaler/DivByQOverTRounded/reconstructAndScale/", r), func(b *testing.B) {

			rescaler := NewSimpleScaler(testParams.T, r)

			coeffs := make([]*big.Int, r.N)
			for i := uint64(0); i < r.N; i++ {
				coeffs[i] = RandInt(r.ModulusBigint)
			}

			tmp0 := r.NewPoly()
			tmp1 := r.NewPoly()

			r.SetCoefficientsBigint(coeffs, tmp0)

			for i := 0; i < b.N; i++ {
				rescaler.reconstructAndScale(tmp0, tmp1)
			}
		})

		b.Run(testString("SimpleScaler/DivByQOverTRounded/reconstructThenScale/", r), func(b *testing.B) {

			rescaler := NewSimpleScaler(testParams.T, r)

			coeffs := make([]*big.Int, r.N)
			for i := uint64(0); i < r.N; i++ {
				coeffs[i] = RandInt(r.ModulusBigint)
			}

			tmp0 := r.NewPoly()
			tmp1 := r.NewPoly()

			r.SetCoefficientsBigint(coeffs, tmp0)

			for i := 0; i < b.N; i++ {
				rescaler.reconstructThenScale(tmp0, tmp1)
			}
		})

		b.Run(testString("RNSScaler/DivByQOverTRounded/", r), func(b *testing.B) {

			coeffs := make([]*big.Int, r.N)
			for i := uint64(0); i < r.N; i++ {
				coeffs[i] = RandInt(r.ModulusBigint)
			}

			scaler := NewRNSScaler(testParams.T, r)
			polyQ := r.NewPoly()
			polyT := NewPoly(r.N, 1)
			r.SetCoefficientsBigint(coeffs, polyQ)

			for i := 0; i < b.N; i++ {
				scaler.DivByQOverTRounded(polyQ, polyT)
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
