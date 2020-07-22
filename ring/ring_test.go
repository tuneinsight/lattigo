package ring

import (
	"fmt"
	"math"
	"math/big"
	"math/bits"
	"math/rand"
	"testing"
	"time"

	"github.com/ldsec/lattigo/utils"

	"github.com/stretchr/testify/require"
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

func TestRing(t *testing.T) {
	t.Run("PRNG", testPRNG)
	t.Run("GenerateNTTPrimes", testGenerateNTTPrimes)
	t.Run("ImportExportPolyString", testImportExportPolyString)
	t.Run("DivFloorByLastModulusMany", testDivFloorByLastModulusMany)
	t.Run("DivRoundByLastModulusMany", testDivRoundByLastModulusMany)
	t.Run("MarshalBinary", testMarshalBinary)
	t.Run("UniformSampler", testUniformSampler)
	t.Run("GaussianSampler", testGaussianSampler)
	t.Run("TernarySampler", testTernarySampler)
	t.Run("GaloisShift", testGaloisShift)
	t.Run("BRed", testBRed)
	t.Run("MRed", testMRed)
	t.Run("MulScalarBigint", testMulScalarBigint)
	t.Run("MulPoly", testMulPoly)
	t.Run("ExtendBasis", testExtendBasis)
	t.Run("Scaling", testScaling)
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
		sum := make([]byte, context.N)

		t.Run(testString("", context), func(t *testing.T) {
			prng1, err := utils.NewKeyedPRNG(nil)
			if err != nil {
				panic(err)
			}
			prng2, err := utils.NewKeyedPRNG(nil)
			if err != nil {
				panic(err)
			}

			prng1.SetClock(sum, 256)
			prng2.SetClock(sum, 256)

			crsGenerator1 := NewUniformSampler(prng1, context)
			crsGenerator2 := NewUniformSampler(prng2, context)

			p0 := crsGenerator1.ReadNew()
			p1 := crsGenerator2.ReadNew()

			require.True(t, context.Equal(p0, p1))
		})
	}
}

func testGenerateNTTPrimes(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		t.Run(testString("", context), func(t *testing.T) {

			primes := GenerateNTTPrimes(55, uint64(bits.Len64(context.N)-1), uint64(len(context.Modulus)))

			for _, q := range primes {
				require.Equal(t, q&((context.N<<1)-1), uint64(1))
				require.True(t, IsPrime(q), q)
			}
		})
	}
}

func testImportExportPolyString(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, context)

		t.Run(testString("", context), func(t *testing.T) {

			p0 := uniformSampler.ReadNew()
			p1 := context.NewPoly()

			context.SetCoefficientsString(context.PolyToString(p0), p1)

			require.True(t, context.Equal(p0, p1))
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
			for i := uint64(0); i < context.N; i++ {
				for j := 0; j < len(context.Modulus)-nbRescals; j++ {
					require.Equalf(t, polWant.Coeffs[j][i], polTest.Coeffs[j][i], "coeff %v Qi%v = %s", i, j, coeffs[i].String())
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
			for i := uint64(0); i < context.N; i++ {
				for j := 0; j < len(context.Modulus)-nbRescals; j++ {
					require.Equalf(t, polWant.Coeffs[j][i], polTest.Coeffs[j][i], "coeff %v Qi%v = %s", i, j, coeffs[i].String())
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

			require.Equal(t, contextTest.N, context.N)
			require.Equal(t, contextTest.Modulus, context.Modulus)
		})

		t.Run(testString("Poly/", context), func(t *testing.T) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, context)

			p := uniformSampler.ReadNew()
			pTest := context.NewPoly()

			data, _ := p.MarshalBinary()

			_ = pTest.UnmarshalBinary(data)

			for i := range context.Modulus {
				require.Equal(t, p.Coeffs[i][:context.N], pTest.Coeffs[i][:context.N])
			}
		})
	}
}

func testUniformSampler(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, context)

		t.Run(testString("Read", context), func(t *testing.T) {
			pol := context.NewPoly()
			uniformSampler.Read(pol)
			for i := uint64(0); i < context.N; i++ {
				for j, qi := range context.Modulus {
					require.False(t, pol.Coeffs[j][i] > qi)
				}
			}
		})

		t.Run(testString("ReadNew", context), func(t *testing.T) {
			pol := uniformSampler.ReadNew()
			for i := uint64(0); i < context.N; i++ {
				for j, qi := range context.Modulus {
					require.False(t, pol.Coeffs[j][i] > qi)
				}
			}
		})

	}
}

func testGaussianSampler(t *testing.T) {

	sigma := testParams.sigma
	bound := uint64(sigma * 6)

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		t.Run(testString("", context), func(t *testing.T) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			gaussianSampler := NewGaussianSampler(prng, context, sigma, bound)
			pol := gaussianSampler.ReadNew()

			for i := uint64(0); i < context.N; i++ {
				for j, qi := range context.Modulus {
					require.False(t, uint64(bound) < pol.Coeffs[j][i] && pol.Coeffs[j][i] < (qi-uint64(bound)))
				}
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

			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			ternarySampler := NewTernarySampler(prng, context, rho, false)

			ternarySampler.Read(pol)

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

			threshold := 0.075

			ratio := math.Round(float64(countOne+countMOn)/float64(countZer)*100.0) / 100.0

			min := ((1 - rho) / rho) * (1.0 - threshold)
			max := ((1 - rho) / rho) * (1.0 + threshold)

			require.Greater(t, ratio, min)
			require.Less(t, ratio, max)
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

					require.Equalf(t, BRed(x, y, q, context.bredParams[j]), result.Uint64(), "x = %v, y=%v", x, y)
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

					require.Equalf(t, MRed(x, MForm(y, q, context.bredParams[j]), q, context.mredParams[j]), result.Uint64(), "x = %v, y=%v", x, y)
				}
			}
		})
	}
}

func testGaloisShift(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		t.Run(testString("", context), func(t *testing.T) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, context)

			pWant := uniformSampler.ReadNew()
			pTest := pWant.CopyNew()

			context.BitReverse(pTest, pTest)
			context.InvNTT(pTest, pTest)
			context.Rotate(pTest, 1, pTest)
			context.NTT(pTest, pTest)
			context.BitReverse(pTest, pTest)
			context.Reduce(pTest, pTest)

			context.Shift(pWant, 1, pWant)

			for i := range context.Modulus {
				require.Equal(t, pTest.Coeffs[i][:context.N], pWant.Coeffs[i][:context.N])
			}
		})
	}
}

func testMForm(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		t.Run(testString("", context), func(t *testing.T) {

			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, context)

			polWant := uniformSampler.ReadNew()
			polTest := context.NewPoly()

			context.MForm(polWant, polTest)
			context.InvMForm(polTest, polTest)

			require.True(t, context.Equal(polWant, polTest))
		})
	}
}

func testMulScalarBigint(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		t.Run(testString("", context), func(t *testing.T) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, context)

			polWant := uniformSampler.ReadNew()
			polTest := polWant.CopyNew()

			rand1 := RandUniform(prng, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)
			rand2 := RandUniform(prng, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)

			scalarBigint := NewUint(rand1)
			scalarBigint.Mul(scalarBigint, NewUint(rand2))

			context.MulScalar(polWant, rand1, polWant)
			context.MulScalar(polWant, rand2, polWant)
			context.MulScalarBigint(polTest, scalarBigint, polTest)

			require.True(t, context.Equal(polWant, polTest))
		})
	}
}

func testMulPoly(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, context)

		p1 := uniformSampler.ReadNew()
		p2 := uniformSampler.ReadNew()
		p3Test := context.NewPoly()
		p3Want := context.NewPoly()

		context.Reduce(p1, p1)
		context.Reduce(p2, p2)

		context.MulPolyNaive(p1, p2, p3Want)

		t.Run(testString("MulPolyBarrett/", context), func(t *testing.T) {

			context.MulPoly(p1, p2, p3Test)

			require.Equal(t, p3Want.Coeffs[0][:context.N], p3Test.Coeffs[0][:context.N])
		})

		t.Run(testString("MulPolyMontgomery/", context), func(t *testing.T) {

			context.MForm(p1, p1)
			context.MForm(p2, p2)

			context.MulPolyMontgomery(p1, p2, p3Test)

			context.InvMForm(p3Test, p3Test)

			require.Equal(t, p3Want.Coeffs[0][:context.N], p3Test.Coeffs[0][:context.N])
		})
	}
}

func testExtendBasis(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		contextQ := genPolyContext(parameters[0])
		contextP := genPolyContext(parameters[0])

		t.Run(testString("", contextQ), func(t *testing.T) {

			basisextender := NewFastBasisExtender(contextQ, contextP)

			coeffs := make([]*big.Int, contextQ.N)
			for i := uint64(0); i < contextQ.N; i++ {
				coeffs[i] = RandInt(contextQ.ModulusBigint)
			}

			Pol := contextQ.NewPoly()
			PolTest := contextP.NewPoly()
			PolWant := contextP.NewPoly()

			contextQ.SetCoefficientsBigint(coeffs, Pol)
			contextP.SetCoefficientsBigint(coeffs, PolWant)

			basisextender.ModUpSplitQP(uint64(len(contextQ.Modulus)-1), Pol, PolTest)

			for i := range contextP.Modulus {
				require.Equal(t, PolTest.Coeffs[i][:contextQ.N], PolWant.Coeffs[i][:contextQ.N])
			}
		})
	}
}

func testScaling(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		t.Run(testString("SimpleScaling", context), func(t *testing.T) {

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

			rescaler.DivByQOverTRounded(PolTest, PolTest)

			for i := uint64(0); i < context.N; i++ {
				require.Equal(t, PolTest.Coeffs[0][i], coeffsWant[i].Uint64())
			}
		})

		t.Run(testString("RNSScaling", context), func(t *testing.T) {

			scaler := NewRNSScaler(testParams.T, context)

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

			polyQ := context.NewPoly()
			polyT := NewPoly(context.N, 1)
			context.SetCoefficientsBigint(coeffs, polyQ)

			scaler.DivByQOverTRounded(polyQ, polyT)

			for i := uint64(0); i < context.N; i++ {
				require.Equal(t, polyT.Coeffs[0][i], coeffsWant[i].Uint64())
			}
		})
	}
}

func testMultByMonomial(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		context := genPolyContext(parameters[0])

		t.Run(testString("", context), func(t *testing.T) {

			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, context)

			p1 := uniformSampler.ReadNew()

			p3Test := context.NewPoly()
			p3Want := context.NewPoly()

			context.MultByMonomial(p1, 1, p3Test)
			context.MultByMonomial(p3Test, 8, p3Test)

			context.MultByMonomial(p1, 9, p3Want)

			require.Equal(t, p3Want.Coeffs[0][:context.N], p3Test.Coeffs[0][:context.N])
		})
	}
}
