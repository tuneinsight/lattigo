package ring

import (
	"fmt"
	"math/rand"
	"testing"
)

func BenchmarkPolynomial(b *testing.B) {

	for i := uint64(0); i < 1; i++ {

		N := uint64(1 << (12 + i))
		T := uint64(65537)

		Qi := Qi60[uint64(len(Qi60))-2<<i:]

		Pi := Pi60[uint64(len(Pi60))-((2<<i)+1):]

		sigma := 3.19

		contextT := NewContext()
		contextT.SetParameters(N, []uint64{T})
		contextT.GenNTTParams()

		contextQ := NewContext()
		contextQ.SetParameters(N, Qi)
		contextQ.GenNTTParams()

		contextP := NewContext()
		contextP.SetParameters(N, Pi)
		contextP.GenNTTParams()

		contextQP := NewContext()
		contextQP.Merge(contextQ, contextP)

		benchmarkContext(N, Qi, b)

		benchmarkKYSGaussPoly(sigma, contextQ, b)

		benchmarkTernaryPoly(contextQ, b)

		benchmarkUniformPoly(contextQ, b)

		benchmarkMForm(contextQ, b)

		benchmarkNTT(contextQ, b)

		benchmarkInvNTT(contextQ, b)

		benchmarkMulScalar(contextQ, b)

		benchmarkNeg(contextQ, b)

		benchmarkSub(contextQ, b)

		benchmarkAdd(contextQ, b)

		benchmarkMulCoeffs(contextQ, b)

		benchmarkMulCoeffsMontgomery(contextQ, b)

		benchmarkMulPoly(contextQ, b)

		benchmarkMulPolyMontgomery(contextQ, b)

		//benchmarkMulPolyNaiveMontgomery(contextQ, b)

		benchmarkExtendBasis(contextQ, contextP, contextQP, b)

		benchmarkSimpleScalerScale(T, contextQ, b)

		benchmarkComplexScalerScale(T, contextQ, contextP, contextQP, b)

		benchmarkMarshaler(contextQ, b)

		benchmarkUnMarshaler(contextQ, b)

		benchmarkBRed(b)

		benchmarkBRedAdd(b)

		benchmarkMRed(b)

	}
}

func benchmarkContext(N uint64, Qi []uint64, b *testing.B) {
	var context *Context
	b.Run(fmt.Sprintf("N=%d/limbs=%d/Context", N, len(Qi)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context = NewContext()
			context.SetParameters(N, Qi)
			context.GenNTTParams()
		}
	})
}

func benchmarkUnMarshaler(context *Context, b *testing.B) {

	p := context.NewUniformPoly()
	b.Run(fmt.Sprintf("N=%d/limbs=%d/MarshalBinary", context.N, len(context.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.MarshalBinary()
		}
	})
}

func benchmarkMarshaler(context *Context, b *testing.B) {

	p := context.NewUniformPoly()

	data, _ := p.MarshalBinary()
	b.Run(fmt.Sprintf("N=%d/limbs=%d/UnMarshalBinary", context.N, len(context.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.UnMarshalBinary(data)
		}
	})
}

func benchmarkKYSGaussPoly(sigma float64, context *Context, b *testing.B) {

	bound := int(sigma * 6)

	KYS := context.NewKYSampler(sigma, bound)

	pol := context.NewPoly()

	b.Run(fmt.Sprintf("N=%d/limbs=%d/KYS.Sample", context.N, len(context.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			KYS.Sample(pol)
		}
	})

	b.Run(fmt.Sprintf("N=%d/limbs=%d/KYS.SampleNTT", context.N, len(context.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			KYS.SampleNTT(pol)
		}
	})
}

func benchmarkTernaryPoly(context *Context, b *testing.B) {

	ternarySampler := context.NewTernarySampler()

	pol := context.NewPoly()

	b.Run(fmt.Sprintf("N=%d/limbs=%d/SampleTernary(0.5)", context.N, len(context.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ternarySampler.Sample(0.5, pol)
		}
	})

	b.Run(fmt.Sprintf("N=%d/limbs=%d/SampleTernary(0.5)MontgomeryNTT", context.N, len(context.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ternarySampler.SampleMontgomeryNTT(0.5, pol)
		}
	})

	b.Run(fmt.Sprintf("N=%d/limbs=%d/SampleTernary(1/3)", context.N, len(context.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ternarySampler.Sample(1.0/3, pol)
		}
	})

	b.Run(fmt.Sprintf("N=%d/limbs=%d/SampleTernary(1/3)MontgomeryNTT", context.N, len(context.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ternarySampler.SampleMontgomeryNTT(1.0/3, pol)

		}
	})
}

func benchmarkUniformPoly(context *Context, b *testing.B) {

	b.Run(fmt.Sprintf("N=%d/limbs=%d/NewUniformPoly", context.N, len(context.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context.NewUniformPoly()
		}
	})
}

func benchmarkMForm(context *Context, b *testing.B) {

	p := context.NewUniformPoly()

	b.Run(fmt.Sprintf("N=%d/limbs=%d/MForm", context.N, len(context.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context.MForm(p, p)
		}
	})
}

func benchmarkNTT(context *Context, b *testing.B) {

	p := context.NewUniformPoly()
	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/limbs=%d/NTT", context.N, len(context.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context.NTT(p, p)
		}
	})
}

func benchmarkInvNTT(context *Context, b *testing.B) {

	p := context.NewUniformPoly()

	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/limbs=%d/InvNTT", context.N, len(context.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context.InvNTT(p, p)
		}
	})
}

func benchmarkMulCoeffs(context *Context, b *testing.B) {

	p := context.NewUniformPoly()

	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/limbs=%d/MulCoeffs", context.N, len(context.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context.MulCoeffs(p, p, p)
		}
	})
}

func benchmarkMulCoeffsMontgomery(context *Context, b *testing.B) {

	p := context.NewUniformPoly()

	context.MForm(p, p)

	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/limbs=%d/MulCoeffs_Montgomery", context.N, len(context.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context.MulCoeffsMontgomery(p, p, p)
		}
	})
}

func benchmarkMulPoly(context *Context, b *testing.B) {

	p := context.NewUniformPoly()

	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/limbs=%d/MulPoly", context.N, len(context.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context.MulPoly(p, p, p)
		}
	})
}

func benchmarkMulPolyMontgomery(context *Context, b *testing.B) {

	p := context.NewUniformPoly()

	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/limbs=%d/MulPoly_Montgomery", context.N, len(context.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context.MulPolyMontgomery(p, p, p)
		}
	})
}

func benchmarkMulPolyNaiveMontgomery(context *Context, b *testing.B) {

	p := context.NewUniformPoly()

	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/limbs=%d/MulPoly_Naive_Montgomery", context.N, len(context.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context.MulPolyNaiveMontgomery(p, p, p)
		}
	})
}

func benchmarkAdd(context *Context, b *testing.B) {

	p := context.NewUniformPoly()

	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/limbs=%d/Add", context.N, len(context.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context.Add(p, p, p)
		}
	})
}

func benchmarkSub(context *Context, b *testing.B) {

	p := context.NewUniformPoly()

	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/limbs=%d/Sub", context.N, len(context.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context.Sub(p, p, p)
		}
	})
}

func benchmarkNeg(context *Context, b *testing.B) {

	p := context.NewUniformPoly()

	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/limbs=%d/Neg", context.N, len(context.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context.Neg(p, p)
		}
	})
}

func benchmarkMulScalar(context *Context, b *testing.B) {

	p := context.NewUniformPoly()

	scalar := uint64(12345678987654321)

	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/limbs=%d/MulScalar", context.N, len(context.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context.MulScalar(p, scalar, p)
		}
	})
}

func benchmarkExtendBasis(contextQ, contextP, contextQP *Context, b *testing.B) {

	BasisExtenderQP := NewBasisExtender(contextQ, contextP)

	p0 := contextQ.NewUniformPoly()
	p1 := contextQP.NewPoly()

	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/limbs=%d+%d/ExtendBasis", contextQ.N, len(contextQ.Modulus), len(contextP.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			BasisExtenderQP.ExtendBasis(p0, p1)
		}
	})

	b.Run(fmt.Sprintf("N=%d/limbs=%d+%d/ExtendBasis_Approximate", contextQ.N, len(contextQ.Modulus), len(contextP.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			BasisExtenderQP.ExtendBasisApproximate(p0, p1)
		}
	})
}

func benchmarkSimpleScalerScale(T uint64, context *Context, b *testing.B) {

	SimpleScaler := NewSimpleScaler(T, context)

	p0 := context.NewUniformPoly()
	p1 := context.NewPoly()

	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/limbs=%d/SimpleScaling", context.N, len(context.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			SimpleScaler.Scale(p0, p1)
		}
	})
}

func benchmarkComplexScalerScale(T uint64, contextQ, contextP, contextQP *Context, b *testing.B) {

	ComplexScalerQP := NewComplexScaler(T, contextQ, contextP)

	p0 := contextQP.NewUniformPoly()
	p1 := contextQ.NewUniformPoly()

	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/limbs=%d+%d/ComplexScaling", contextQ.N, len(contextP.Modulus), len(contextQ.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ComplexScalerQP.Scale(p0, p1)

		}
	})
}

func benchmarkBRed(b *testing.B) {

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

func benchmarkBRedAdd(b *testing.B) {

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

func benchmarkMRed(b *testing.B) {

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
