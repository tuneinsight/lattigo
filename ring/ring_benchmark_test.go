package ring

import (
	"fmt"
	"math/rand"
	"testing"
)

func Benchmark_POLYNOMIAL(b *testing.B) {

	for i := uint64(0); i < 1; i++ {

		N := uint64(2 << (12 + i))
		T := uint64(65537)

		N = 1 << 16

		Qi := Qi60[uint64(len(Qi60))-2<<i:]

		Qi = Qi60[uint64(len(Qi60))-32:]

		Pi := Pi60[uint64(len(Pi60))-((2<<i)+1):]

		sigma := 3.19

		contextT := NewContext()
		contextT.SetParameters(N, []uint64{T})
		contextT.ValidateParameters()

		contextQ := NewContext()
		contextQ.SetParameters(N, Qi)
		contextQ.ValidateParameters()

		contextP := NewContext()
		contextP.SetParameters(N, Pi)
		contextP.ValidateParameters()

		contextQP := NewContext()
		contextQP.Merge(contextQ, contextP)

		benchmark_Context(N, Qi, b)

		benchmark_KYSGaussPoly(sigma, contextQ, b)

		benchmark_TernaryPoly(contextQ, b)

		benchmark_UniformPoly(contextQ, b)

		benchmark_MForm(contextQ, b)

		benchmark_NTT(contextQ, b)

		benchmark_InvNTT(contextQ, b)

		benchmark_MulScalar(contextQ, b)

		benchmark_Neg(contextQ, b)

		benchmark_Sub(contextQ, b)

		benchmark_Add(contextQ, b)

		benchmark_MulCoeffs(contextQ, b)

		benchmark_MulCoeffsMontgomery(contextQ, b)

		benchmark_MulPoly(contextQ, b)

		benchmark_MulPolyMontgomery(contextQ, b)

		benchmark_MulPolyNaiveMontgomery(contextQ, b)

		benchmark_ExtendBasis(contextQ, contextP, contextQP, b)

		benchmark_SimpleScaling(T, contextQ, b)

		benchmark_ComplexScaling(T, contextQ, contextP, contextQP, b)

		benchmark_Marshaler(contextQ, b)

		benchmark_UnMarshaler(contextQ, b)

		benchmark_BRed(b)

		benchmark_BRedAdd(b)

		benchmark_MRed(b)

	}
}

func benchmark_Context(N uint64, Qi []uint64, b *testing.B) {
	var context *Context
	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/Context", N, len(Qi), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context = NewContext()
			context.SetParameters(N, Qi)
			context.ValidateParameters()
		}
	})
}

func benchmark_UnMarshaler(context *Context, b *testing.B) {

	p := context.NewUniformPoly()
	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/MarshalBinary", context.N, len(context.Modulus), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.MarshalBinary()
		}
	})
}

func benchmark_Marshaler(context *Context, b *testing.B) {

	p := context.NewUniformPoly()

	data, _ := p.MarshalBinary()
	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/UnMarshalBinary", context.N, len(context.Modulus), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.UnMarshalBinary(data)
		}
	})
}

func benchmark_KYSGaussPoly(sigma float64, context *Context, b *testing.B) {

	bound := int(sigma * 6)

	KYS := context.NewKYSampler(sigma, bound)

	pol := context.NewPoly()

	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/KYS.Sample", context.N, len(context.Modulus), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			KYS.Sample(pol)
		}
	})

	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/KYS.SampleNTT", context.N, len(context.Modulus), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			KYS.SampleNTT(pol)
		}
	})
}

func benchmark_TernaryPoly(context *Context, b *testing.B) {

	ternarySampler := context.NewTernarySampler()

	pol := context.NewPoly()

	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/SampleTernary(0.5)", context.N, len(context.Modulus), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ternarySampler.Sample(0.5, pol)
		}
	})

	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/SampleTernary(0.5)MontgomeryNTT", context.N, len(context.Modulus), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ternarySampler.SampleMontgomeryNTT(0.5, pol)
		}
	})

	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/SampleTernary(1/3)", context.N, len(context.Modulus), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ternarySampler.Sample(1.0/3, pol)
		}
	})

	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/SampleTernary(1/3)MontgomeryNTT", context.N, len(context.Modulus), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ternarySampler.SampleMontgomeryNTT(1.0/3, pol)

		}
	})

	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/SampleTernary(0.66)", context.N, len(context.Modulus), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			SampleTernary(context.N, 0.5)
		}
	})

}

func benchmark_UniformPoly(context *Context, b *testing.B) {

	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/NewUniformPoly", context.N, len(context.Modulus), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context.NewUniformPoly()
		}
	})
}

func benchmark_MForm(context *Context, b *testing.B) {

	p := context.NewUniformPoly()

	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/MForm", context.N, len(context.Modulus), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context.MForm(p, p)
		}
	})
}

func benchmark_NTT(context *Context, b *testing.B) {

	p := context.NewUniformPoly()
	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/NTT", context.N, len(context.Modulus), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context.NTT(p, p)
		}
	})
}

func benchmark_InvNTT(context *Context, b *testing.B) {

	p := context.NewUniformPoly()

	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/InvNTT", context.N, len(context.Modulus), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context.InvNTT(p, p)
		}
	})
}

func benchmark_MulCoeffs(context *Context, b *testing.B) {

	p := context.NewUniformPoly()

	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/MulCoeffs", context.N, len(context.Modulus), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context.MulCoeffs(p, p, p)
		}
	})
}

func benchmark_MulCoeffsMontgomery(context *Context, b *testing.B) {

	p := context.NewUniformPoly()

	context.MForm(p, p)

	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/MulCoeffs_Montgomery", context.N, len(context.Modulus), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context.MulCoeffsMontgomery(p, p, p)
		}
	})
}

func benchmark_MulPoly(context *Context, b *testing.B) {

	p := context.NewUniformPoly()

	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/MulPoly", context.N, len(context.Modulus), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context.MulPoly(p, p, p)
		}
	})
}

func benchmark_MulPolyMontgomery(context *Context, b *testing.B) {

	p := context.NewUniformPoly()

	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/MulPoly_Montgomery", context.N, len(context.Modulus), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context.MulPolyMontgomery(p, p, p)
		}
	})
}

func benchmark_MulPolyNaiveMontgomery(context *Context, b *testing.B) {

	p := context.NewUniformPoly()

	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/MulPoly_Naive_Montgomery", context.N, len(context.Modulus), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context.MulPolyNaiveMontgomery(p, p, p)
		}
	})
}

func benchmark_Add(context *Context, b *testing.B) {

	p := context.NewUniformPoly()

	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/Add", context.N, len(context.Modulus), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context.Add(p, p, p)
		}
	})
}

func benchmark_Sub(context *Context, b *testing.B) {

	p := context.NewUniformPoly()

	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/Sub", context.N, len(context.Modulus), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context.Sub(p, p, p)
		}
	})
}

func benchmark_Neg(context *Context, b *testing.B) {

	p := context.NewUniformPoly()

	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/Neg", context.N, len(context.Modulus), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context.Neg(p, p)
		}
	})
}

func benchmark_MulScalar(context *Context, b *testing.B) {

	p := context.NewUniformPoly()

	scalar := uint64(12345678987654321)

	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/MulScalar", context.N, len(context.Modulus), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			context.MulScalar(p, scalar, p)
		}
	})
}

func benchmark_ExtendBasis(contextQ, contextP, contextQP *Context, b *testing.B) {

	BasisExtenderQP, _ := NewBasisExtender(contextQ, contextP)

	p0 := contextQ.NewUniformPoly()
	p1 := contextQP.NewPoly()

	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/Pi=%dx%dbit/ExtendBasis", contextQ.N, len(contextQ.Modulus), 60, len(contextP.Modulus), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			BasisExtenderQP.ExtendBasis(p0, p1)
		}
	})

	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/Pi=%dx%dbit/ExtendBasis_Approximate", contextQ.N, len(contextQ.Modulus), 60, len(contextP.Modulus), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			BasisExtenderQP.ExtendBasisApproximate(p0, p1)
		}
	})
}

func benchmark_SimpleScaling(T uint64, context *Context, b *testing.B) {

	SimpleScaler, _ := NewSimpleScaler(T, context)

	p0 := context.NewUniformPoly()
	p1 := context.NewPoly()

	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/SimpleScaling", context.N, len(context.Modulus), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			SimpleScaler.Scale(p0, p1)
		}
	})
}

func benchmark_ComplexScaling(T uint64, contextQ, contextP, contextQP *Context, b *testing.B) {

	ComplexScalerQP, _ := NewComplexScaler(T, contextQ, contextP)

	p0 := contextQP.NewUniformPoly()
	p1 := contextQ.NewUniformPoly()

	b.ResetTimer()

	b.Run(fmt.Sprintf("N=%d/Qi=%dx%dbit/Pi=%dx%dbit/ComplexScaling", contextQ.N, len(contextP.Modulus), 60, len(contextQ.Modulus), 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ComplexScalerQP.Scale(p0, p1)

		}
	})
}

func benchmark_BRed(b *testing.B) {

	q := uint64(1033576114481528833)
	u := BRedParams(q)

	x := rand.Uint64() % q
	y := rand.Uint64() % q

	b.ResetTimer()

	b.Run(fmt.Sprintf("Qi=%d/BRed", 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			x = BRed(x, y, q, u)
		}
	})
}

func benchmark_BRedAdd(b *testing.B) {

	q := uint64(1033576114481528833)
	u := BRedParams(q)

	x := rand.Uint64()

	b.ResetTimer()

	b.Run(fmt.Sprintf("Qi=%dbit/BRedAdd", 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			BRedAdd(x, q, u)
		}
	})
}

func benchmark_MRed(b *testing.B) {

	q := uint64(1033576114481528833)

	x := rand.Uint64() % q
	y := rand.Uint64() % q

	u := BRedParams(q)

	y = MForm(y, q, u)

	m := MRedParams(q)

	b.ResetTimer()

	b.Run(fmt.Sprintf("Qi=%dbit/MRed", 60), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			x = MRed(x, y, q, m)
		}
	})
}
