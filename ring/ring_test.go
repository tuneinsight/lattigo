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

func testString(opname string, ringQ *Ring) string {
	return fmt.Sprintf("%sN=%d/limbs=%d", opname, ringQ.N, len(ringQ.Modulus))
}

var testParams = new(PolynomialTestParams)

func init() {
	rand.Seed(time.Now().UnixNano())

	testParams.T = 0x3ee0001

	testParams.polyParams = [][2]*Parameters{
		//{DefaultParamsQi[12], DefaultParamsPi[12]},
		//{DefaultParamsQi[13], DefaultParamsPi[13]},
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

func getRing(params *Parameters) *Ring {
	r, err := NewRing(params.N, params.Moduli)
	if err != nil {
		panic(err)
	}
	return r
}

func testPRNG(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])
		sum := make([]byte, ringQ.N)

		t.Run(testString("", ringQ), func(t *testing.T) {
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

			crsGenerator1 := NewUniformSampler(prng1, ringQ)
			crsGenerator2 := NewUniformSampler(prng2, ringQ)

			p0 := crsGenerator1.ReadNew()
			p1 := crsGenerator2.ReadNew()

			require.True(t, ringQ.Equal(p0, p1))
		})
	}
}

func testGenerateNTTPrimes(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])

		t.Run(testString("", ringQ), func(t *testing.T) {

			primes := GenerateNTTPrimes(55, uint64(bits.Len64(ringQ.N)-1), uint64(len(ringQ.Modulus)))

			for _, q := range primes {
				require.Equal(t, q&((ringQ.N<<1)-1), uint64(1))
				require.True(t, IsPrime(q), q)
			}
		})
	}
}

func testImportExportPolyString(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, ringQ)

		t.Run(testString("", ringQ), func(t *testing.T) {

			p0 := uniformSampler.ReadNew()
			p1 := ringQ.NewPoly()

			ringQ.SetCoefficientsString(ringQ.PolyToString(p0), p1)

			require.True(t, ringQ.Equal(p0, p1))
		})
	}
}

func testDivFloorByLastModulusMany(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])

		t.Run(testString("", ringQ), func(t *testing.T) {

			coeffs := make([]*big.Int, ringQ.N)
			for i := uint64(0); i < ringQ.N; i++ {
				coeffs[i] = RandInt(ringQ.ModulusBigint)
				coeffs[i].Quo(coeffs[i], NewUint(10))
			}

			nbRescals := len(ringQ.Modulus) - 1

			coeffsWant := make([]*big.Int, ringQ.N)
			for i := range coeffs {
				coeffsWant[i] = new(big.Int).Set(coeffs[i])
				for j := 0; j < nbRescals; j++ {
					coeffsWant[i].Quo(coeffsWant[i], NewUint(ringQ.Modulus[len(ringQ.Modulus)-1-j]))
				}
			}

			polTest := ringQ.NewPoly()
			polWant := ringQ.NewPoly()

			ringQ.SetCoefficientsBigint(coeffs, polTest)
			ringQ.SetCoefficientsBigint(coeffsWant, polWant)

			ringQ.DivFloorByLastModulusMany(polTest, uint64(nbRescals))
			for i := uint64(0); i < ringQ.N; i++ {
				for j := 0; j < len(ringQ.Modulus)-nbRescals; j++ {
					require.Equalf(t, polWant.Coeffs[j][i], polTest.Coeffs[j][i], "coeff %v Qi%v = %s", i, j, coeffs[i].String())
				}
			}
		})
	}
}

func testDivRoundByLastModulusMany(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])

		t.Run(testString("", ringQ), func(t *testing.T) {

			coeffs := make([]*big.Int, ringQ.N)
			for i := uint64(0); i < ringQ.N; i++ {
				coeffs[i] = RandInt(ringQ.ModulusBigint)
				coeffs[i].Quo(coeffs[i], NewUint(10))
			}

			nbRescals := len(ringQ.Modulus) - 1

			coeffsWant := make([]*big.Int, ringQ.N)
			for i := range coeffs {
				coeffsWant[i] = new(big.Int).Set(coeffs[i])
				for j := 0; j < nbRescals; j++ {
					DivRound(coeffsWant[i], NewUint(ringQ.Modulus[len(ringQ.Modulus)-1-j]), coeffsWant[i])
				}
			}

			polTest := ringQ.NewPoly()
			polWant := ringQ.NewPoly()

			ringQ.SetCoefficientsBigint(coeffs, polTest)
			ringQ.SetCoefficientsBigint(coeffsWant, polWant)

			ringQ.DivRoundByLastModulusMany(polTest, uint64(nbRescals))
			for i := uint64(0); i < ringQ.N; i++ {
				for j := 0; j < len(ringQ.Modulus)-nbRescals; j++ {
					require.Equalf(t, polWant.Coeffs[j][i], polTest.Coeffs[j][i], "coeff %v Qi%v = %s", i, j, coeffs[i].String())
				}
			}
		})
	}
}

func testMarshalBinary(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])

		t.Run(testString("Ring/", ringQ), func(t *testing.T) {

			data, _ := ringQ.MarshalBinary()

			ringQTest := new(Ring)
			ringQTest.UnmarshalBinary(data)

			require.Equal(t, ringQTest.N, ringQ.N)
			require.Equal(t, ringQTest.Modulus, ringQ.Modulus)
		})

		t.Run(testString("Poly/", ringQ), func(t *testing.T) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, ringQ)

			p := uniformSampler.ReadNew()
			pTest := ringQ.NewPoly()

			data, _ := p.MarshalBinary()

			_ = pTest.UnmarshalBinary(data)

			for i := range ringQ.Modulus {
				require.Equal(t, p.Coeffs[i][:ringQ.N], pTest.Coeffs[i][:ringQ.N])
			}
		})
	}
}

func testUniformSampler(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])

		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, ringQ)

		t.Run(testString("Read", ringQ), func(t *testing.T) {
			pol := ringQ.NewPoly()
			uniformSampler.Read(pol)
			for i := uint64(0); i < ringQ.N; i++ {
				for j, qi := range ringQ.Modulus {
					require.False(t, pol.Coeffs[j][i] > qi)
				}
			}
		})

		t.Run(testString("ReadNew", ringQ), func(t *testing.T) {
			pol := uniformSampler.ReadNew()
			for i := uint64(0); i < ringQ.N; i++ {
				for j, qi := range ringQ.Modulus {
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

		ringQ := getRing(parameters[0])

		t.Run(testString("", ringQ), func(t *testing.T) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			gaussianSampler := NewGaussianSampler(prng, ringQ, sigma, bound)
			pol := gaussianSampler.ReadNew()

			for i := uint64(0); i < ringQ.N; i++ {
				for j, qi := range ringQ.Modulus {
					require.False(t, uint64(bound) < pol.Coeffs[j][i] && pol.Coeffs[j][i] < (qi-uint64(bound)))
				}
			}
		})
	}
}

func testTernarySampler(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])

		t.Run(testString("", ringQ), func(t *testing.T) {

			countOne := uint64(0)
			countZer := uint64(0)
			countMOn := uint64(0)

			pol := ringQ.NewPoly()

			rho := 1.0 / 3

			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			ternarySampler := NewTernarySampler(prng, ringQ, rho, false)

			ternarySampler.Read(pol)

			for i := range pol.Coeffs[0] {
				if pol.Coeffs[0][i] == ringQ.Modulus[0]-1 {
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

		ringQ := getRing(parameters[0])

		t.Run(testString("", ringQ), func(t *testing.T) {

			for j, q := range ringQ.Modulus {

				bigQ := NewUint(q)

				for i := 0; i < 65536; i++ {
					x := rand.Uint64() % q
					y := rand.Uint64() % q

					result := NewUint(x)
					result.Mul(result, NewUint(y))
					result.Mod(result, bigQ)

					require.Equalf(t, BRed(x, y, q, ringQ.BredParams[j]), result.Uint64(), "x = %v, y=%v", x, y)
				}
			}
		})
	}
}

func testMRed(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])

		t.Run(testString("", ringQ), func(t *testing.T) {

			for j := range ringQ.Modulus {

				q := ringQ.Modulus[j]

				bigQ := NewUint(q)

				for i := 0; i < 65536; i++ {

					x := rand.Uint64() % q
					y := rand.Uint64() % q

					result := NewUint(x)
					result.Mul(result, NewUint(y))
					result.Mod(result, bigQ)

					require.Equalf(t, MRed(x, MForm(y, q, ringQ.BredParams[j]), q, ringQ.MredParams[j]), result.Uint64(), "x = %v, y=%v", x, y)
				}
			}
		})
	}
}

func testGaloisShift(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])

		t.Run(testString("", ringQ), func(t *testing.T) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, ringQ)

			pWant := uniformSampler.ReadNew()
			pTest := pWant.CopyNew()

			ringQ.BitReverse(pTest, pTest)
			ringQ.InvNTT(pTest, pTest)
			ringQ.Rotate(pTest, 1, pTest)
			ringQ.NTT(pTest, pTest)
			ringQ.BitReverse(pTest, pTest)
			ringQ.Reduce(pTest, pTest)

			ringQ.Shift(pWant, 1, pWant)

			for i := range ringQ.Modulus {
				require.Equal(t, pTest.Coeffs[i][:ringQ.N], pWant.Coeffs[i][:ringQ.N])
			}
		})
	}
}

func testMForm(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])

		t.Run(testString("", ringQ), func(t *testing.T) {

			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, ringQ)

			polWant := uniformSampler.ReadNew()
			polTest := ringQ.NewPoly()

			ringQ.MForm(polWant, polTest)
			ringQ.InvMForm(polTest, polTest)

			require.True(t, ringQ.Equal(polWant, polTest))
		})
	}
}

func testMulScalarBigint(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])

		t.Run(testString("", ringQ), func(t *testing.T) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, ringQ)

			polWant := uniformSampler.ReadNew()
			polTest := polWant.CopyNew()

			rand1 := RandUniform(prng, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)
			rand2 := RandUniform(prng, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)

			scalarBigint := NewUint(rand1)
			scalarBigint.Mul(scalarBigint, NewUint(rand2))

			ringQ.MulScalar(polWant, rand1, polWant)
			ringQ.MulScalar(polWant, rand2, polWant)
			ringQ.MulScalarBigint(polTest, scalarBigint, polTest)

			require.True(t, ringQ.Equal(polWant, polTest))
		})
	}
}

func testMulPoly(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, ringQ)

		p1 := uniformSampler.ReadNew()
		p2 := uniformSampler.ReadNew()
		p3Test := ringQ.NewPoly()
		p3Want := ringQ.NewPoly()

		ringQ.Reduce(p1, p1)
		ringQ.Reduce(p2, p2)

		ringQ.MulPolyNaive(p1, p2, p3Want)

		t.Run(testString("MulPolyBarrett/", ringQ), func(t *testing.T) {

			ringQ.MulPoly(p1, p2, p3Test)

			require.Equal(t, p3Want.Coeffs[0][:ringQ.N], p3Test.Coeffs[0][:ringQ.N])
		})

		t.Run(testString("MulPolyMontgomery/", ringQ), func(t *testing.T) {

			ringQ.MForm(p1, p1)
			ringQ.MForm(p2, p2)

			ringQ.MulPolyMontgomery(p1, p2, p3Test)

			ringQ.InvMForm(p3Test, p3Test)

			require.Equal(t, p3Want.Coeffs[0][:ringQ.N], p3Test.Coeffs[0][:ringQ.N])
		})
	}
}

func testExtendBasis(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])
		ringP := getRing(parameters[0])

		t.Run(testString("", ringQ), func(t *testing.T) {

			basisextender := NewFastBasisExtender(ringQ, ringP)

			coeffs := make([]*big.Int, ringQ.N)
			for i := uint64(0); i < ringQ.N; i++ {
				coeffs[i] = RandInt(ringQ.ModulusBigint)
			}

			Pol := ringQ.NewPoly()
			PolTest := ringP.NewPoly()
			PolWant := ringP.NewPoly()

			ringQ.SetCoefficientsBigint(coeffs, Pol)
			ringP.SetCoefficientsBigint(coeffs, PolWant)

			basisextender.ModUpSplitQP(uint64(len(ringQ.Modulus)-1), Pol, PolTest)

			for i := range ringP.Modulus {
				require.Equal(t, PolTest.Coeffs[i][:ringQ.N], PolWant.Coeffs[i][:ringQ.N])
			}
		})
	}
}

func testScaling(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])

		t.Run(testString("SimpleScaling", ringQ), func(t *testing.T) {

			rescaler := NewSimpleScaler(testParams.T, ringQ)

			coeffs := make([]*big.Int, ringQ.N)
			for i := uint64(0); i < ringQ.N; i++ {
				coeffs[i] = RandInt(ringQ.ModulusBigint)
			}

			coeffsWant := make([]*big.Int, ringQ.N)
			for i := range coeffs {
				coeffsWant[i] = new(big.Int).Set(coeffs[i])
				coeffsWant[i].Mul(coeffsWant[i], NewUint(testParams.T))
				DivRound(coeffsWant[i], ringQ.ModulusBigint, coeffsWant[i])
				coeffsWant[i].Mod(coeffsWant[i], NewUint(testParams.T))
			}

			PolTest := ringQ.NewPoly()

			ringQ.SetCoefficientsBigint(coeffs, PolTest)

			rescaler.DivByQOverTRounded(PolTest, PolTest)

			for i := uint64(0); i < ringQ.N; i++ {
				require.Equal(t, PolTest.Coeffs[0][i], coeffsWant[i].Uint64())
			}
		})

		t.Run(testString("RNSScaling", ringQ), func(t *testing.T) {

			scaler := NewRNSScaler(testParams.T, ringQ)

			coeffs := make([]*big.Int, ringQ.N)
			for i := uint64(0); i < ringQ.N; i++ {
				coeffs[i] = RandInt(ringQ.ModulusBigint)
			}

			coeffsWant := make([]*big.Int, ringQ.N)
			for i := range coeffs {
				coeffsWant[i] = new(big.Int).Set(coeffs[i])
				coeffsWant[i].Mul(coeffsWant[i], NewUint(testParams.T))
				DivRound(coeffsWant[i], ringQ.ModulusBigint, coeffsWant[i])
				coeffsWant[i].Mod(coeffsWant[i], NewUint(testParams.T))
			}

			polyQ := ringQ.NewPoly()
			polyT := NewPoly(ringQ.N, 1)
			ringQ.SetCoefficientsBigint(coeffs, polyQ)

			scaler.DivByQOverTRounded(polyQ, polyT)

			for i := uint64(0); i < ringQ.N; i++ {
				require.Equal(t, polyT.Coeffs[0][i], coeffsWant[i].Uint64())
			}
		})
	}
}

func testMultByMonomial(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		ringQ := getRing(parameters[0])

		t.Run(testString("", ringQ), func(t *testing.T) {

			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, ringQ)

			p1 := uniformSampler.ReadNew()

			p3Test := ringQ.NewPoly()
			p3Want := ringQ.NewPoly()

			ringQ.MultByMonomial(p1, 1, p3Test)
			ringQ.MultByMonomial(p3Test, 8, p3Test)

			ringQ.MultByMonomial(p1, 9, p3Want)

			require.Equal(t, p3Want.Coeffs[0][:ringQ.N], p3Test.Coeffs[0][:ringQ.N])
		})
	}
}
