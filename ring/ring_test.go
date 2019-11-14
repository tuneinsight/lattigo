package ring

import (
	"fmt"
	"log"
	"math/bits"
	"math/rand"
	"testing"
	"time"
)

func Test_Polynomial(t *testing.T) {

	rand.Seed(time.Now().UnixNano())

	for i := uint64(0); i < 1; i++ {

		N := uint64(2 << (12 + i))
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

		test_PRNG(contextQ, t)

		// ok!
		test_GenerateNTTPrimes(N, Qi[0], t)

		// ok!
		test_ImportExportPolyString(contextQ, t)

		// ok!
		test_Marshaler(contextQ, t)

		// TODO : check that the coefficients are within the bound
		test_GaussianPoly(sigma, contextQ, t)

		// ok!
		test_BRed(contextQ, t)

		// ok!
		test_MRed(contextQ, t)

		// ok!
		test_Rescale(contextQ, t)

		// ok!
		test_MulScalarBigint(contextQ, t)

		// ok!
		test_Shift(contextQ, t)

		// ok!
		test_GaloisShift(contextQ, t)

		// ok!
		test_MForm(contextQ, t)

		// ok!
		test_MulPoly(contextQ, t)

		// ok!
		test_MulPoly_Montgomery(contextQ, t)

		// ok!
		test_ExtendBasis(contextQ, contextP, contextQP, t)

		// ok!
		test_SimpleScaling(T, contextQ, contextP, t)

		test_MultByMonomial(contextQ, t)
	}
}

func test_PRNG(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("PRNG"), func(t *testing.T) {

		Ha, _ := NewPRNG(nil)
		Hb, _ := NewPRNG(nil)

		// Random 32 byte seed
		seed1 := []byte{0x48, 0xc3, 0x31, 0x12, 0x74, 0x98, 0xd3, 0xf2,
			0x7b, 0x15, 0x15, 0x9b, 0x50, 0xc4, 0x9c, 0x00,
			0x7d, 0xa5, 0xea, 0x68, 0x1f, 0xed, 0x4f, 0x99,
			0x54, 0xc0, 0x52, 0xc0, 0x75, 0xff, 0xf7, 0x5c}

		// New reseed of the PRNG after one clock cycle with the seed1
		seed2 := []byte{250, 228, 6, 63, 97, 110, 68, 153,
			147, 236, 236, 37, 152, 89, 129, 32,
			185, 5, 221, 180, 160, 217, 247, 201,
			211, 188, 160, 163, 176, 83, 83, 138}

		Ha.Seed(seed1)
		Hb.Seed(append(seed1, seed2...)) //Append works since blake2b hashes blocks of 512 bytes

		Ha.SetClock(256)
		Hb.SetClock(255)

		a := Ha.Clock()
		b := Hb.Clock()

		for i := 0; i < 32; i++ {
			if a[i] != b[i] {
				t.Errorf("error : error prng")
				break
			}
		}

		crs_generator_1, _ := NewCRPGenerator(nil, context)
		crs_generator_2, _ := NewCRPGenerator(nil, context)

		crs_generator_1.Seed(seed1)
		crs_generator_2.Seed(append(seed1, seed2...)) //Append works since blake2b hashes blocks of 512 bytes

		crs_generator_1.SetClock(256)
		crs_generator_2.SetClock(255)

		p0 := crs_generator_1.Clock()
		p1 := crs_generator_2.Clock()

		if context.Equal(p0, p1) != true {
			t.Errorf("error : crs prng generator")
		}
	})

}

func test_GenerateNTTPrimes(N, Qi uint64, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/GenerateNTTPrimes", N), func(t *testing.T) {

		primes, err := GenerateNTTPrimes(N, Qi, 100, 60, true)

		if err != nil {
			t.Errorf("error : generateNTTPrimes")
		}

		for _, q := range primes {
			if bits.Len64(q) != 60 || q&((N<<1)-1) != 1 || IsPrime(q) != true {
				t.Errorf("error : GenerateNTTPrimes for q = %v", q)
				break
			}
		}
	})
}

func test_ImportExportPolyString(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d/ImportExportPolyString", context.N, len(context.Modulus)), func(t *testing.T) {

		p0 := context.NewUniformPoly()
		p1 := context.NewPoly()

		context.SetCoefficientsString(context.PolyToString(p0), p1)

		if context.Equal(p0, p1) != true {
			t.Errorf("error : import/export ring from/to string")
		}
	})
}

func test_Rescale(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d/DivFloorByLastModulusMany", context.N, len(context.Modulus)), func(t *testing.T) {

		coeffs := make([]*Int, context.N)
		for i := uint64(0); i < context.N; i++ {
			coeffs[i] = RandInt(context.ModulusBigint)
			coeffs[i].Div(coeffs[i], NewUint(10))
		}

		nbRescals := len(context.Modulus) - 1

		coeffsWant := make([]*Int, context.N)
		for i := range coeffs {
			coeffsWant[i] = coeffs[i].Copy()
			for j := 0; j < nbRescals; j++ {
				coeffsWant[i].Div(coeffsWant[i], NewUint(context.Modulus[len(context.Modulus)-1-j]))
			}
		}

		polTest := context.NewPoly()
		polWant := context.NewPoly()

		context.SetCoefficientsBigint(coeffs, polTest)
		context.SetCoefficientsBigint(coeffsWant, polWant)

		context.DivFloorByLastModulusMany(polTest, uint64(nbRescals))
		state := true
		for i := uint64(0); i < context.N && state; i++ {
			for j := 0; j < len(context.Modulus)-nbRescals && state; j++ {
				if polWant.Coeffs[j][i] != polTest.Coeffs[j][i] {
					t.Errorf("error : coeff %v Qi%v = %s, want %v have %v", i, j, coeffs[i].String(), polWant.Coeffs[j][i], polTest.Coeffs[j][i])
					state = false
				}
			}
		}
	})

	t.Run(fmt.Sprintf("N=%d/limbs=%d/DivRoundByLastModulusMany", context.N, len(context.Modulus)), func(t *testing.T) {

		coeffs := make([]*Int, context.N)
		for i := uint64(0); i < context.N; i++ {
			coeffs[i] = RandInt(context.ModulusBigint)
			coeffs[i].Div(coeffs[i], NewUint(10))
		}

		nbRescals := len(context.Modulus) - 1

		coeffsWant := make([]*Int, context.N)
		for i := range coeffs {
			coeffsWant[i] = coeffs[i].Copy()
			for j := 0; j < nbRescals; j++ {
				coeffsWant[i].DivRound(coeffsWant[i], NewUint(context.Modulus[len(context.Modulus)-1-j]))
			}
		}

		polTest := context.NewPoly()
		polWant := context.NewPoly()

		context.SetCoefficientsBigint(coeffs, polTest)
		context.SetCoefficientsBigint(coeffsWant, polWant)

		context.DivRoundByLastModulusMany(polTest, uint64(nbRescals))
		state := true
		for i := uint64(0); i < context.N && state; i++ {
			for j := 0; j < len(context.Modulus)-nbRescals && state; j++ {
				if polWant.Coeffs[j][i] != polTest.Coeffs[j][i] {
					t.Errorf("error : coeff %v Qi%v = %s, want %v have %v", i, j, coeffs[i].String(), polWant.Coeffs[j][i], polTest.Coeffs[j][i])
					state = false
				}
			}
		}
	})
}

func test_Marshaler(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d/MarshalContext", context.N, len(context.Modulus)), func(t *testing.T) {

		data, _ := context.MarshalBinary()

		contextTest := NewContext()
		contextTest.UnMarshalBinary(data)

		if contextTest.N != context.N {
			t.Errorf("ERROR encoding/decoding N")
		}

		for i := range context.Modulus {
			if contextTest.Modulus[i] != context.Modulus[i] {
				t.Errorf("ERROR encoding/decoding Modulus")
			}
		}
	})

	t.Run(fmt.Sprintf("N=%d/limbs=%d/MarshalPoly", context.N, len(context.Modulus)), func(t *testing.T) {

		p := context.NewUniformPoly()
		pTest := context.NewPoly()

		data, _ := p.MarshalBinary()

		pTest, _ = pTest.UnMarshalBinary(data)

		for i := range context.Modulus {
			for j := uint64(0); j < context.N; j++ {
				if p.Coeffs[i][j] != pTest.Coeffs[i][j] {
					t.Errorf("PolyBytes Import Error : want %v, have %v", p.Coeffs[i][j], pTest.Coeffs[i][j])
				}
			}
		}
	})
}

func test_GaussianPoly(sigma float64, context *Context, t *testing.T) {

	bound := int(sigma * 6)
	KYS := context.NewKYSampler(sigma, bound)

	pol := context.NewPoly()

	t.Run(fmt.Sprintf("N=%d/limbs=%d/NewGaussPoly", context.N, len(context.Modulus)), func(t *testing.T) {
		KYS.Sample(pol)
	})

	countOne := 0
	countZer := 0
	countMOn := 0
	t.Run(fmt.Sprintf("N=%d/limbs=%d/NewTernaryPoly", context.N, len(context.Modulus)), func(t *testing.T) {
		if err := context.SampleTernary(pol, 1.0/3); err != nil {
			log.Fatal(err)
		}

		//fmt.Println(pol.Coeffs[0])

		for i := range pol.Coeffs[0] {
			if pol.Coeffs[0][i] == context.Modulus[0]-1 {
				countMOn += 1
			}

			if pol.Coeffs[0][i] == 0 {
				countZer += 1
			}

			if pol.Coeffs[0][i] == 1 {
				countOne += 1
			}
		}

		//fmt.Println("-1 :", countMOn)
		//fmt.Println(" 0 :", countZer)
		//fmt.Println(" 1 :", countOne)
	})
}

func test_BRed(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d/BRed", context.N, len(context.Modulus)), func(t *testing.T) {
		for j, q := range context.Modulus {

			bigQ := NewUint(q)

			for i := 0; i < 65536; i++ {
				x := rand.Uint64() % q
				y := rand.Uint64() % q

				result := NewUint(x)
				result.Mul(result, NewUint(y))
				result.Mod(result, bigQ)

				test := BRed(x, y, q, context.bredParams[j])
				want := result.Uint64()

				if test != want {
					t.Errorf("error : 128bit barrett multiplication, x = %v, y=%v, have = %v, want =%v", x, y, test, want)
					break
				}
			}
		}
	})
}

func test_MRed(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d/MRed", context.N, len(context.Modulus)), func(t *testing.T) {
		for j := range context.Modulus {

			q := context.Modulus[j]

			bigQ := NewUint(q)

			for i := 0; i < 65536; i++ {

				x := rand.Uint64() % q
				y := rand.Uint64() % q

				result := NewUint(x)
				result.Mul(result, NewUint(y))
				result.Mod(result, bigQ)

				test := MRed(x, MForm(y, q, context.bredParams[j]), q, context.mredParams[j])
				want := result.Uint64()

				if test != want {
					t.Errorf("error : 128bit montgomery multiplication, x = %v, y=%v, have = %v, want =%v", x, y, test, want)
					break
				}

			}
		}
	})
}

func test_Shift(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d/Shift", context.N, len(context.Modulus)), func(t *testing.T) {
		pWant := context.NewUniformPoly()
		pTest := context.NewPoly()

		context.Shift(pTest, 1, pWant)
		for i := range context.Modulus {
			for j := uint64(0); j < context.N; j++ {
				if pTest.Coeffs[i][(j+1)%context.N] != pWant.Coeffs[i][j] {
					t.Errorf("error : ShiftR")
					break
				}
			}
		}
	})
}

func test_GaloisShift(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d/GaloisShift", context.N, len(context.Modulus)), func(t *testing.T) {

		pWant := context.NewUniformPoly()
		pTest := pWant.CopyNew()

		context.BitReverse(pTest, pTest)
		context.InvNTT(pTest, pTest)
		context.Rotate(pTest, 1, pTest)
		context.NTT(pTest, pTest)
		context.BitReverse(pTest, pTest)
		context.Reduce(pTest, pTest)

		context.Shift(pWant, 1, pWant)

		for i := range context.Modulus {
			for j := uint64(0); j < context.N; j++ {
				if pTest.Coeffs[i][j] != pWant.Coeffs[i][j] {
					t.Errorf("error : GaloisShiftR")
					break
				}
			}
		}
	})
}

func test_MForm(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d/MForm", context.N, len(context.Modulus)), func(t *testing.T) {

		polWant := context.NewUniformPoly()
		polTest := context.NewPoly()

		context.MForm(polWant, polTest)
		context.InvMForm(polTest, polTest)

		if context.Equal(polWant, polTest) != true {
			t.Errorf("error : 128bit MForm/InvMForm")
		}
	})

}

func test_MulScalarBigint(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d/MulScalarBigint", context.N, len(context.Modulus)), func(t *testing.T) {

		polWant := context.NewUniformPoly()
		polTest := polWant.CopyNew()

		rand1 := RandUniform(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)
		rand2 := RandUniform(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)

		scalarBigint := NewUint(rand1)
		scalarBigint.Mul(scalarBigint, NewUint(rand2))

		context.MulScalar(polWant, rand1, polWant)
		context.MulScalar(polWant, rand2, polWant)
		context.MulScalarBigint(polTest, scalarBigint, polTest)

		if context.Equal(polWant, polTest) != true {
			t.Errorf("error : mulScalarBigint")
		}
	})

}

func test_MulPoly(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d/MulPoly", context.N, len(context.Modulus)), func(t *testing.T) {

		p1 := context.NewUniformPoly()
		p2 := context.NewUniformPoly()
		p3Test := context.NewPoly()
		p3Want := context.NewPoly()

		context.Reduce(p1, p1)
		context.Reduce(p2, p2)

		context.MulPolyNaive(p1, p2, p3Want)

		context.MulPoly(p1, p2, p3Test)

		for i := uint64(0); i < context.N; i++ {
			if p3Want.Coeffs[0][i] != p3Test.Coeffs[0][i] {
				t.Errorf("ERROR MUL COEFF %v, want %v - has %v", i, p3Want.Coeffs[0][i], p3Test.Coeffs[0][i])
				break
			}
		}
	})
}

func test_MulPoly_Montgomery(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d/MulPoly_Montgomery", context.N, len(context.Modulus)), func(t *testing.T) {
		p1 := context.NewUniformPoly()
		p2 := context.NewUniformPoly()
		p3Test := context.NewPoly()
		p3Want := context.NewPoly()

		context.MForm(p1, p1)
		context.MForm(p2, p2)

		context.MulPolyNaiveMontgomery(p1, p2, p3Want)
		context.MulPolyMontgomery(p1, p2, p3Test)

		context.InvMForm(p3Test, p3Test)
		context.InvMForm(p3Want, p3Want)

		for i := uint64(0); i < context.N; i++ {
			if p3Want.Coeffs[0][i] != p3Test.Coeffs[0][i] {
				t.Errorf("ERROR MUL COEFF %v, want %v - has %v", i, p3Want.Coeffs[0][i], p3Test.Coeffs[0][i])
				break
			}
		}
	})
}

func test_ExtendBasis(contextQ, contextP, contextQP *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d+%d/ExtendBasis", contextQ.N, len(contextQ.Modulus), len(contextP.Modulus)), func(t *testing.T) {

		basisextender := NewBasisExtender(contextQ, contextP)

		coeffs := make([]*Int, contextQ.N)
		for i := uint64(0); i < contextQ.N; i++ {
			coeffs[i] = RandInt(contextQ.ModulusBigint)
		}

		PolTest := contextQ.NewPoly()
		PolWant := contextQP.NewPoly()

		contextQ.SetCoefficientsBigint(coeffs, PolTest)
		contextQP.SetCoefficientsBigint(coeffs, PolWant)

		basisextender.ExtendBasis(PolTest, PolTest)

		for i := range contextQP.Modulus {
			for j := uint64(0); j < contextQ.N; j++ {
				if PolTest.Coeffs[i][j] != PolWant.Coeffs[i][j] {
					t.Errorf("error extendBasis, want %v - has %v", PolTest.Coeffs[i][j], PolWant.Coeffs[i][j])
					break
				}
			}
		}
	})
}

func test_SimpleScaling(T uint64, contextT, contextQ *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d/T=%v/SimpleScaling", contextQ.N, len(contextQ.Modulus), T), func(t *testing.T) {

		rescaler := NewSimpleScaler(T, contextQ)

		coeffs := make([]*Int, contextQ.N)
		for i := uint64(0); i < contextQ.N; i++ {
			coeffs[i] = RandInt(contextQ.ModulusBigint)
		}

		coeffsWant := make([]*Int, contextQ.N)
		for i := range coeffs {
			coeffsWant[i] = Copy(coeffs[i])
			coeffsWant[i].Mul(coeffsWant[i], NewUint(T))
			coeffsWant[i].DivRound(coeffsWant[i], contextQ.ModulusBigint)
			coeffsWant[i].Mod(coeffsWant[i], NewUint(T))
		}

		PolTest := contextQ.NewPoly()
		PolWant := contextT.NewPoly()

		contextQ.SetCoefficientsBigint(coeffs, PolTest)

		rescaler.Scale(PolTest, PolTest)

		contextT.SetCoefficientsBigint(coeffsWant, PolWant)

		for i := uint64(0); i < contextT.N; i++ {
			if PolWant.Coeffs[0][i] != PolTest.Coeffs[0][i] {
				t.Errorf("error : simple scaling, want %v have %v", PolWant.Coeffs[0][i], PolTest.Coeffs[0][i])
				break
			}
		}
	})
}

func test_MultByMonomial(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/limbs=%d/MultByMonomial", context.N, len(context.Modulus)), func(t *testing.T) {

		p1 := context.NewUniformPoly()

		p3Test := context.NewPoly()
		p3Want := context.NewPoly()

		context.MultByMonomial(p1, 1, p3Test)
		context.MultByMonomial(p3Test, 8, p3Test)

		context.MultByMonomial(p1, 9, p3Want)

		for i := uint64(0); i < context.N; i++ {
			if p3Want.Coeffs[0][i] != p3Test.Coeffs[0][i] {
				t.Errorf("ERROR MUL BY MONOMIAL %v, want %v - has %v", i, p3Want.Coeffs[0][i], p3Test.Coeffs[0][i])
				break
			}
		}
	})
}
