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

func testString(opname string, r *Ring) string {
	return fmt.Sprintf("%sN=%d/limbs=%d", opname, r.N, len(r.Modulus))
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

func getRing(params *Parameters) (r *Ring) {
	var err error
	if r, err = NewRing(params.N, params.Moduli); err != nil {
		panic(err)
	}
	return r
}

func testPRNG(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		r := getRing(parameters[0])
		sum := make([]byte, r.N)

		t.Run(testString("", r), func(t *testing.T) {
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

			crsGenerator1 := NewUniformSampler(prng1, r)
			crsGenerator2 := NewUniformSampler(prng2, r)

			p0 := crsGenerator1.ReadNew()
			p1 := crsGenerator2.ReadNew()

			require.True(t, r.Equal(p0, p1))
		})
	}
}

func testGenerateNTTPrimes(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		r := getRing(parameters[0])

		t.Run(testString("", r), func(t *testing.T) {

			primes := GenerateNTTPrimes(55, uint64(bits.Len64(r.N)-1), uint64(len(r.Modulus)))

			for _, q := range primes {
				require.Equal(t, q&((r.N<<1)-1), uint64(1))
				require.True(t, IsPrime(q), q)
			}
		})
	}
}

func testImportExportPolyString(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		r := getRing(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, r)

		t.Run(testString("", r), func(t *testing.T) {

			p0 := uniformSampler.ReadNew()
			p1 := r.NewPoly()

			r.SetCoefficientsString(r.PolyToString(p0), p1)

			require.True(t, r.Equal(p0, p1))
		})
	}
}

func testDivFloorByLastModulusMany(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		r := getRing(parameters[0])

		t.Run(testString("", r), func(t *testing.T) {

			coeffs := make([]*big.Int, r.N)
			for i := uint64(0); i < r.N; i++ {
				coeffs[i] = RandInt(r.ModulusBigint)
				coeffs[i].Quo(coeffs[i], NewUint(10))
			}

			nbRescals := len(r.Modulus) - 1

			coeffsWant := make([]*big.Int, r.N)
			for i := range coeffs {
				coeffsWant[i] = new(big.Int).Set(coeffs[i])
				for j := 0; j < nbRescals; j++ {
					coeffsWant[i].Quo(coeffsWant[i], NewUint(r.Modulus[len(r.Modulus)-1-j]))
				}
			}

			polTest := r.NewPoly()
			polWant := r.NewPoly()

			r.SetCoefficientsBigint(coeffs, polTest)
			r.SetCoefficientsBigint(coeffsWant, polWant)

			r.DivFloorByLastModulusMany(polTest, uint64(nbRescals))
			for i := uint64(0); i < r.N; i++ {
				for j := 0; j < len(r.Modulus)-nbRescals; j++ {
					require.Equalf(t, polWant.Coeffs[j][i], polTest.Coeffs[j][i], "coeff %v Qi%v = %s", i, j, coeffs[i].String())
				}
			}
		})
	}
}

func testDivRoundByLastModulusMany(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		r := getRing(parameters[0])

		t.Run(testString("", r), func(t *testing.T) {

			coeffs := make([]*big.Int, r.N)
			for i := uint64(0); i < r.N; i++ {
				coeffs[i] = RandInt(r.ModulusBigint)
				coeffs[i].Quo(coeffs[i], NewUint(10))
			}

			nbRescals := len(r.Modulus) - 1

			coeffsWant := make([]*big.Int, r.N)
			for i := range coeffs {
				coeffsWant[i] = new(big.Int).Set(coeffs[i])
				for j := 0; j < nbRescals; j++ {
					DivRound(coeffsWant[i], NewUint(r.Modulus[len(r.Modulus)-1-j]), coeffsWant[i])
				}
			}

			polTest := r.NewPoly()
			polWant := r.NewPoly()

			r.SetCoefficientsBigint(coeffs, polTest)
			r.SetCoefficientsBigint(coeffsWant, polWant)

			r.DivRoundByLastModulusMany(polTest, uint64(nbRescals))
			for i := uint64(0); i < r.N; i++ {
				for j := 0; j < len(r.Modulus)-nbRescals; j++ {
					require.Equalf(t, polWant.Coeffs[j][i], polTest.Coeffs[j][i], "coeff %v Qi%v = %s", i, j, coeffs[i].String())
				}
			}
		})
	}
}

func testMarshalBinary(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		r := getRing(parameters[0])

		t.Run(testString("r/", r), func(t *testing.T) {

			data, _ := r.MarshalBinary()

			contextTest := new(Ring)
			contextTest.UnmarshalBinary(data)

			require.Equal(t, contextTest.N, r.N)
			require.Equal(t, contextTest.Modulus, r.Modulus)
		})

		t.Run(testString("Poly/", r), func(t *testing.T) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, r)

			p := uniformSampler.ReadNew()
			pTest := r.NewPoly()

			data, _ := p.MarshalBinary()

			_ = pTest.UnmarshalBinary(data)

			for i := range r.Modulus {
				require.Equal(t, p.Coeffs[i][:r.N], pTest.Coeffs[i][:r.N])
			}
		})
	}
}

func testUniformSampler(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		r := getRing(parameters[0])

		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, r)

		t.Run(testString("Read", r), func(t *testing.T) {
			pol := r.NewPoly()
			uniformSampler.Read(pol)
			for i := uint64(0); i < r.N; i++ {
				for j, qi := range r.Modulus {
					require.False(t, pol.Coeffs[j][i] > qi)
				}
			}
		})

		t.Run(testString("ReadNew", r), func(t *testing.T) {
			pol := uniformSampler.ReadNew()
			for i := uint64(0); i < r.N; i++ {
				for j, qi := range r.Modulus {
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

		r := getRing(parameters[0])

		t.Run(testString("", r), func(t *testing.T) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			gaussianSampler := NewGaussianSampler(prng, r, sigma, bound)
			pol := gaussianSampler.ReadNew()

			for i := uint64(0); i < r.N; i++ {
				for j, qi := range r.Modulus {
					require.False(t, uint64(bound) < pol.Coeffs[j][i] && pol.Coeffs[j][i] < (qi-uint64(bound)))
				}
			}
		})
	}
}

func testTernarySampler(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		r := getRing(parameters[0])

		t.Run(testString("", r), func(t *testing.T) {

			countOne := uint64(0)
			countZer := uint64(0)
			countMOn := uint64(0)

			pol := r.NewPoly()

			rho := 1.0 / 3

			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			ternarySampler := NewTernarySampler(prng, r, rho, false)

			ternarySampler.Read(pol)

			for i := range pol.Coeffs[0] {
				if pol.Coeffs[0][i] == r.Modulus[0]-1 {
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

		r := getRing(parameters[0])

		t.Run(testString("", r), func(t *testing.T) {

			for j, q := range r.Modulus {

				bigQ := NewUint(q)

				for i := 0; i < 65536; i++ {
					x := rand.Uint64() % q
					y := rand.Uint64() % q

					result := NewUint(x)
					result.Mul(result, NewUint(y))
					result.Mod(result, bigQ)

					require.Equalf(t, BRed(x, y, q, r.BredParams[j]), result.Uint64(), "x = %v, y=%v", x, y)
				}
			}
		})
	}
}

func testMRed(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		r := getRing(parameters[0])

		t.Run(testString("", r), func(t *testing.T) {

			for j := range r.Modulus {

				q := r.Modulus[j]

				bigQ := NewUint(q)

				for i := 0; i < 65536; i++ {

					x := rand.Uint64() % q
					y := rand.Uint64() % q

					result := NewUint(x)
					result.Mul(result, NewUint(y))
					result.Mod(result, bigQ)

					require.Equalf(t, MRed(x, MForm(y, q, r.BredParams[j]), q, r.MredParams[j]), result.Uint64(), "x = %v, y=%v", x, y)
				}
			}
		})
	}
}

func testGaloisShift(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		r := getRing(parameters[0])

		t.Run(testString("", r), func(t *testing.T) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, r)

			pWant := uniformSampler.ReadNew()
			pTest := pWant.CopyNew()

			r.BitReverse(pTest, pTest)
			r.InvNTT(pTest, pTest)
			r.Rotate(pTest, 1, pTest)
			r.NTT(pTest, pTest)
			r.BitReverse(pTest, pTest)
			r.Reduce(pTest, pTest)

			r.Shift(pWant, 1, pWant)

			for i := range r.Modulus {
				require.Equal(t, pTest.Coeffs[i][:r.N], pWant.Coeffs[i][:r.N])
			}
		})
	}
}

func testMForm(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		r := getRing(parameters[0])

		t.Run(testString("", r), func(t *testing.T) {

			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, r)

			polWant := uniformSampler.ReadNew()
			polTest := r.NewPoly()

			r.MForm(polWant, polTest)
			r.InvMForm(polTest, polTest)

			require.True(t, r.Equal(polWant, polTest))
		})
	}
}

func testMulScalarBigint(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		r := getRing(parameters[0])

		t.Run(testString("", r), func(t *testing.T) {
			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, r)

			polWant := uniformSampler.ReadNew()
			polTest := polWant.CopyNew()

			rand1 := RandUniform(prng, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)
			rand2 := RandUniform(prng, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)

			scalarBigint := NewUint(rand1)
			scalarBigint.Mul(scalarBigint, NewUint(rand2))

			r.MulScalar(polWant, rand1, polWant)
			r.MulScalar(polWant, rand2, polWant)
			r.MulScalarBigint(polTest, scalarBigint, polTest)

			require.True(t, r.Equal(polWant, polTest))
		})
	}
}

func testMulPoly(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		r := getRing(parameters[0])
		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		uniformSampler := NewUniformSampler(prng, r)

		p1 := uniformSampler.ReadNew()
		p2 := uniformSampler.ReadNew()
		p3Test := r.NewPoly()
		p3Want := r.NewPoly()

		r.Reduce(p1, p1)
		r.Reduce(p2, p2)

		r.MulPolyNaive(p1, p2, p3Want)

		t.Run(testString("MulPolyBarrett/", r), func(t *testing.T) {

			r.MulPoly(p1, p2, p3Test)

			require.Equal(t, p3Want.Coeffs[0][:r.N], p3Test.Coeffs[0][:r.N])
		})

		t.Run(testString("MulPolyMontgomery/", r), func(t *testing.T) {

			r.MForm(p1, p1)
			r.MForm(p2, p2)

			r.MulPolyMontgomery(p1, p2, p3Test)

			r.InvMForm(p3Test, p3Test)

			require.Equal(t, p3Want.Coeffs[0][:r.N], p3Test.Coeffs[0][:r.N])
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

		r := getRing(parameters[0])

		t.Run(testString("SimpleScaling", r), func(t *testing.T) {

			rescaler := NewSimpleScaler(testParams.T, r)

			coeffs := make([]*big.Int, r.N)
			for i := uint64(0); i < r.N; i++ {
				coeffs[i] = RandInt(r.ModulusBigint)
			}

			coeffsWant := make([]*big.Int, r.N)
			for i := range coeffs {
				coeffsWant[i] = new(big.Int).Set(coeffs[i])
				coeffsWant[i].Mul(coeffsWant[i], NewUint(testParams.T))
				DivRound(coeffsWant[i], r.ModulusBigint, coeffsWant[i])
				coeffsWant[i].Mod(coeffsWant[i], NewUint(testParams.T))
			}

			PolTest := r.NewPoly()

			r.SetCoefficientsBigint(coeffs, PolTest)

			rescaler.DivByQOverTRounded(PolTest, PolTest)

			for i := uint64(0); i < r.N; i++ {
				require.Equal(t, PolTest.Coeffs[0][i], coeffsWant[i].Uint64())
			}
		})

		t.Run(testString("RNSScaling", r), func(t *testing.T) {

			scaler := NewRNSScaler(testParams.T, r)

			coeffs := make([]*big.Int, r.N)
			for i := uint64(0); i < r.N; i++ {
				coeffs[i] = RandInt(r.ModulusBigint)
			}

			coeffsWant := make([]*big.Int, r.N)
			for i := range coeffs {
				coeffsWant[i] = new(big.Int).Set(coeffs[i])
				coeffsWant[i].Mul(coeffsWant[i], NewUint(testParams.T))
				DivRound(coeffsWant[i], r.ModulusBigint, coeffsWant[i])
				coeffsWant[i].Mod(coeffsWant[i], NewUint(testParams.T))
			}

			polyQ := r.NewPoly()
			polyT := NewPoly(r.N, 1)
			r.SetCoefficientsBigint(coeffs, polyQ)

			scaler.DivByQOverTRounded(polyQ, polyT)

			for i := uint64(0); i < r.N; i++ {
				require.Equal(t, polyT.Coeffs[0][i], coeffsWant[i].Uint64())
			}
		})
	}
}

func testMultByMonomial(t *testing.T) {

	for _, parameters := range testParams.polyParams {

		r := getRing(parameters[0])

		t.Run(testString("", r), func(t *testing.T) {

			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			uniformSampler := NewUniformSampler(prng, r)

			p1 := uniformSampler.ReadNew()

			p3Test := r.NewPoly()
			p3Want := r.NewPoly()

			r.MultByMonomial(p1, 1, p3Test)
			r.MultByMonomial(p3Test, 8, p3Test)

			r.MultByMonomial(p1, 9, p3Want)

			require.Equal(t, p3Want.Coeffs[0][:r.N], p3Test.Coeffs[0][:r.N])
		})
	}
}
