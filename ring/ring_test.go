package ring

import (
	"fmt"
	"math"
	"math/big"
	"math/bits"
	"math/rand"
	"testing"
	"time"
)

type PolynomialTestParams struct {
	T          uint64
	polyParams [][2]*Parameters
	sigma      float64
}

func testString(opname string, context *Context) string {
	return fmt.Sprintf("%sN=%d/limbs=%d", opname, context.N, len(context.Modulus))
}

var testParams = new(PolynomialTestParams)

func init() {
	rand.Seed(time.Now().UnixNano())

	testParams.T = 0x3ee0001

	testParams.polyParams = [][2]*Parameters{
		{DefaultParamsQi[12], DefaultParamsPi[12]},
		{DefaultParamsQi[13], DefaultParamsPi[13]},
		{DefaultParamsQi[14], DefaultParamsPi[14]},
		//{DefaultParamsQi[15], DefaultParamsPi[15]},
		//{DefaultParamsQi[16], DefaultParamsPi[16]},
	}

	testParams.sigma = 3.19
}

func testRing(t *testing.T) {
	t.Run("PRNG", testPRNG)
	t.Run("GenerateNTTPrimes", testGenerateNTTPrimes)
	t.Run("ImportExportPolyString", testImportExportPolyString)
	t.Run("DivFloorByLastModulusMany", testDivFloorByLastModulusMany)
	t.Run("DivRoundByLastModulusMany", testDivRoundByLastModulusMany)
	t.Run("MarshalBinary", testMarshalBinary)
	t.Run("GaussianSampler", testGaussianSampler)
	t.Run("TernarySampler", testTernarySampler)
	t.Run("GaloisShift", testGaloisShift)
	t.Run("BRed", testBRed)
	t.Run("MRed", testMRed)
	t.Run("MulScalarBigint", testMulScalarBigint)
	t.Run("MulPoly", testMulPoly)
	t.Run("ExtendBasis", testExtendBasis)
	t.Run("SimpleScaling", testSimpleScaling)
	t.Run("MultByMonomial", testMultByMonomial)
}

func genPolyContext(params *Parameters) (context *Context) {
	context = NewContext()
	context.SetParameters(params.N, params.Moduli)
	context.GenNTTParams()
	return context
}

func testPRNG(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		t.Run(testString("", context), func(t *testing.T) {

			crsGenerator1 := NewCRPGenerator(nil, context)
			crsGenerator2 := NewCRPGenerator(nil, context)

			crsGenerator1.Seed(nil)
			crsGenerator2.Seed(nil)

			crsGenerator1.SetClock(256)
			crsGenerator2.SetClock(256)

			p0 := crsGenerator1.ClockNew()
			p1 := crsGenerator2.ClockNew()

			if context.Equal(p0, p1) != true {
				t.Errorf("crs prng generator")
			}
		})
	}
}

func testGenerateNTTPrimes(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		t.Run(testString("", context), func(t *testing.T) {

			primes, err := GenerateNTTPrimes(context.N, context.Modulus[0], uint64(len(context.Modulus)), 60, true)

			if err != nil {
				t.Errorf("generateNTTPrimes")
			}

			for _, q := range primes {
				if bits.Len64(q) != 60 || q&((context.N<<1)-1) != 1 || IsPrime(q) != true {
					t.Errorf("GenerateNTTPrimes for q = %v", q)
					break
				}
			}
		})
	}
}

func testImportExportPolyString(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		t.Run(testString("", context), func(t *testing.T) {

			p0 := context.NewUniformPoly()
			p1 := context.NewPoly()

			context.SetCoefficientsString(context.PolyToString(p0), p1)

			if context.Equal(p0, p1) != true {
				t.Errorf("import/export ring from/to string")
			}
		})
	}
}

func testDivFloorByLastModulusMany(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		t.Run(testString("", context), func(t *testing.T) {

			coeffs := make([]*big.Int, context.N)
			for i := uint64(0); i < context.N; i++ {
				coeffs[i] = RandInt(context.ModulusBigint)
				coeffs[i].Quo(coeffs[i], NewUint(10))
			}

			nbRescals := len(context.Modulus) - 1

			coeffsWant := make([]*big.Int, context.N)
			for i := range coeffs {
				coeffsWant[i] = new(big.Int).Set(coeffs[i])
				for j := 0; j < nbRescals; j++ {
					coeffsWant[i].Quo(coeffsWant[i], NewUint(context.Modulus[len(context.Modulus)-1-j]))
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
	}
}

func testDivRoundByLastModulusMany(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		t.Run(testString("", context), func(t *testing.T) {

			coeffs := make([]*big.Int, context.N)
			for i := uint64(0); i < context.N; i++ {
				coeffs[i] = RandInt(context.ModulusBigint)
				coeffs[i].Quo(coeffs[i], NewUint(10))
			}

			nbRescals := len(context.Modulus) - 1

			coeffsWant := make([]*big.Int, context.N)
			for i := range coeffs {
				coeffsWant[i] = new(big.Int).Set(coeffs[i])
				for j := 0; j < nbRescals; j++ {
					DivRound(coeffsWant[i], NewUint(context.Modulus[len(context.Modulus)-1-j]), coeffsWant[i])
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
}

func testMarshalBinary(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		t.Run(testString("Context/", context), func(t *testing.T) {

			data, _ := context.MarshalBinary()

			contextTest := NewContext()
			contextTest.UnmarshalBinary(data)

			if contextTest.N != context.N {
				t.Errorf("ERROR encoding/decoding N")
			}

			for i := range context.Modulus {
				if contextTest.Modulus[i] != context.Modulus[i] {
					t.Errorf("ERROR encoding/decoding Modulus")
				}
			}
		})

		t.Run(testString("Poly/", context), func(t *testing.T) {

			p := context.NewUniformPoly()
			pTest := context.NewPoly()

			data, _ := p.MarshalBinary()

			_ = pTest.UnmarshalBinary(data)

			for i := range context.Modulus {
				for j := uint64(0); j < context.N; j++ {
					if p.Coeffs[i][j] != pTest.Coeffs[i][j] {
						t.Errorf("PolyBytes Import Error : want %v, have %v", p.Coeffs[i][j], pTest.Coeffs[i][j])
					}
				}
			}
		})
	}
}

func testGaussianSampler(t *testing.T) {

	sigma := testParams.sigma
	bound := int(sigma * 6)

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		t.Run(testString("", context), func(t *testing.T) {

			KYS := context.NewKYSampler(sigma, bound)

			pol := context.NewPoly()

			KYS.Sample(pol)

			state := true

			for i := uint64(0); i < context.N; i++ {
				for j, qi := range context.Modulus {

					if (uint64(bound) < pol.Coeffs[j][i]) && (pol.Coeffs[j][i] < (qi - uint64(bound))) {
						state = false
						break
					}
				}

				if !state {
					break
				}
			}

			if !state {
				t.Errorf("GaussianSampler")
			}
		})
	}
}

func testTernarySampler(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		t.Run(testString("", context), func(t *testing.T) {

			countOne := uint64(0)
			countZer := uint64(0)
			countMOn := uint64(0)

			pol := context.NewPoly()

			rho := 1.0 / 3

			context.SampleTernary(pol, rho)

			for i := range pol.Coeffs[0] {
				if pol.Coeffs[0][i] == context.Modulus[0]-1 {
					countMOn++
				}

				if pol.Coeffs[0][i] == 0 {
					countZer++
				}

				if pol.Coeffs[0][i] == 1 {
					countOne++
				}
			}

			ratio := math.Round(float64(countOne+countMOn)/float64(countZer)*100.0) / 100.0

			if ((1-rho)/rho)*0.90 > ratio || ((1-rho)/rho)*1.10 < ratio {
				fmt.Println(float64(countOne+countMOn) / float64(countZer))
				t.Errorf("TernarySampler")
			}
		})
	}
}

func testBRed(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		t.Run(testString("", context), func(t *testing.T) {

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
						t.Errorf("128bit barrett multiplication, x = %v, y=%v, have = %v, want =%v", x, y, test, want)
						break
					}
				}
			}
		})
	}
}

func testMRed(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		t.Run(testString("", context), func(t *testing.T) {

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
						t.Errorf("128bit montgomery multiplication, x = %v, y=%v, have = %v, want =%v", x, y, test, want)
						break
					}

				}
			}
		})
	}
}

func testGaloisShift(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		t.Run(testString("", context), func(t *testing.T) {

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
						t.Errorf("GaloisShiftR")
						break
					}
				}
			}
		})
	}
}

func testMForm(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		t.Run(testString("", context), func(t *testing.T) {

			polWant := context.NewUniformPoly()
			polTest := context.NewPoly()

			context.MForm(polWant, polTest)
			context.InvMForm(polTest, polTest)

			if context.Equal(polWant, polTest) != true {
				t.Errorf("128bit MForm/InvMForm")
			}
		})
	}
}

func testMulScalarBigint(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		t.Run(testString("", context), func(t *testing.T) {

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
}

func testMulPoly(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		p1 := context.NewUniformPoly()
		p2 := context.NewUniformPoly()
		p3Test := context.NewPoly()
		p3Want := context.NewPoly()

		context.Reduce(p1, p1)
		context.Reduce(p2, p2)

		context.MulPolyNaive(p1, p2, p3Want)

		t.Run(testString("MulPolyBarrett/", context), func(t *testing.T) {

			context.MulPoly(p1, p2, p3Test)

			for i := uint64(0); i < context.N; i++ {
				if p3Want.Coeffs[0][i] != p3Test.Coeffs[0][i] {
					t.Errorf("MULPOLY COEFF %v, want %v - has %v", i, p3Want.Coeffs[0][i], p3Test.Coeffs[0][i])
					break
				}
			}
		})

		t.Run(testString("MulPolyMontgomery/", context), func(t *testing.T) {

			context.MForm(p1, p1)
			context.MForm(p2, p2)

			context.MulPolyMontgomery(p1, p2, p3Test)

			context.InvMForm(p3Test, p3Test)

			for i := uint64(0); i < context.N; i++ {
				if p3Want.Coeffs[0][i] != p3Test.Coeffs[0][i] {
					t.Errorf("MULPOLY COEFF %v, want %v - has %v", i, p3Want.Coeffs[0][i], p3Test.Coeffs[0][i])
					break
				}
			}
		})
	}
}

func testExtendBasis(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		contextQ := genPolyContext(parameters[0])
		contextP := genPolyContext(parameters[0])

		t.Run(testString("", contextQ), func(t *testing.T) {

			basisextender := NewBasisExtender(contextQ, contextP)

			coeffs := make([]*big.Int, contextQ.N)
			for i := uint64(0); i < contextQ.N; i++ {
				coeffs[i] = RandInt(contextQ.ModulusBigint)
			}

			Pol := contextQ.NewPoly()
			PolTest := contextP.NewPoly()
			PolWant := contextP.NewPoly()

			contextQ.SetCoefficientsBigint(coeffs, Pol)
			contextP.SetCoefficientsBigint(coeffs, PolWant)

			basisextender.ExtendBasisSplit(Pol, PolTest)

			for i := range contextP.Modulus {
				for j := uint64(0); j < contextQ.N; j++ {
					if PolTest.Coeffs[i][j] != PolWant.Coeffs[i][j] {
						t.Errorf("error extendBasis, have %v - want %v", PolTest.Coeffs[i][j], PolWant.Coeffs[i][j])
						break
					}
				}
			}
		})
	}
}

func testSimpleScaling(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		t.Run(testString("", context), func(t *testing.T) {

			rescaler := NewSimpleScaler(testParams.T, context)

			coeffs := make([]*big.Int, context.N)
			for i := uint64(0); i < context.N; i++ {
				coeffs[i] = RandInt(context.ModulusBigint)
			}

			coeffsWant := make([]*big.Int, context.N)
			for i := range coeffs {
				coeffsWant[i] = new(big.Int).Set(coeffs[i])
				coeffsWant[i].Mul(coeffsWant[i], NewUint(testParams.T))
				DivRound(coeffsWant[i], context.ModulusBigint, coeffsWant[i])
				coeffsWant[i].Mod(coeffsWant[i], NewUint(testParams.T))
			}

			PolTest := context.NewPoly()

			context.SetCoefficientsBigint(coeffs, PolTest)

			rescaler.Scale(PolTest, PolTest)

			for i := uint64(0); i < context.N; i++ {
				if PolTest.Coeffs[0][i] != coeffsWant[i].Uint64() {
					t.Errorf("simple scaling, coeffs %v want %v have %v", i, coeffsWant[i].Uint64(), PolTest.Coeffs[0][i])
					break
				}
			}
		})
	}
}

func testMultByMonomial(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		t.Run(testString("", context), func(t *testing.T) {

			p1 := context.NewUniformPoly()

			p3Test := context.NewPoly()
			p3Want := context.NewPoly()

			context.MultByMonomial(p1, 1, p3Test)
			context.MultByMonomial(p3Test, 8, p3Test)

			context.MultByMonomial(p1, 9, p3Want)

			for i := uint64(0); i < context.N; i++ {
				if p3Want.Coeffs[0][i] != p3Test.Coeffs[0][i] {
					t.Errorf("MULT BY MONOMIAL %v, want %v - has %v", i, p3Want.Coeffs[0][i], p3Test.Coeffs[0][i])
					break
				}
			}
		})
	}
}
