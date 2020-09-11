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

var T = uint64(0x3ee0001)
var DefaultSigma = float64(3.2)
var DefaultBound = uint64(6 * DefaultSigma)

func testString(opname string, ringQ *Ring) string {
	return fmt.Sprintf("%sN=%d/limbs=%d", opname, ringQ.N, len(ringQ.Modulus))
}

var params = new(testParams)

type testParams struct {
	ringQ           *Ring
	ringP           *Ring
	prng            utils.PRNG
	uniformSamplerQ *UniformSampler
	uniformSamplerP *UniformSampler
}

func GenTestParams(defaultParams *Parameters) (err error) {

	if params.ringQ, err = NewRing(1<<defaultParams.logN, defaultParams.qi); err != nil {
		return err
	}
	if params.ringP, err = NewRing(1<<defaultParams.logN, defaultParams.pi); err != nil {
		return err
	}
	if params.prng, err = utils.NewPRNG(); err != nil {
		return err
	}
	params.uniformSamplerQ = NewUniformSampler(params.prng, params.ringQ)
	params.uniformSamplerP = NewUniformSampler(params.prng, params.ringP)
	return nil
}

func TestRing(t *testing.T) {

	rand.Seed(time.Now().UnixNano())

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

		testPRNG(t)
		testGenerateNTTPrimes(t)
		testImportExportPolyString(t)
		testDivFloorByLastModulusMany(t)
		testDivRoundByLastModulusMany(t)
		testMarshalBinary(t)
		testUniformSampler(t)
		testGaussianSampler(t)
		testTernarySampler(t)
		testGaloisShift(t)
		testModularReduction(t)
		testMulScalarBigint(t)
		testMulPoly(t)
		testExtendBasis(t)
		testScaling(t)
		testMultByMonomial(t)
	}
}

func testPRNG(t *testing.T) {

	sum := make([]byte, params.ringQ.N)
	t.Run(testString("PRNG/", params.ringQ), func(t *testing.T) {
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

		crsGenerator1 := NewUniformSampler(prng1, params.ringQ)
		crsGenerator2 := NewUniformSampler(prng2, params.ringQ)

		p0 := crsGenerator1.ReadNew()
		p1 := crsGenerator2.ReadNew()

		require.True(t, params.ringQ.Equal(p0, p1))
	})

}

func testGenerateNTTPrimes(t *testing.T) {

	t.Run(testString("GenerateNTTPrimes/", params.ringQ), func(t *testing.T) {

		primes := GenerateNTTPrimes(55, uint64(bits.Len64(params.ringQ.N)-1), uint64(len(params.ringQ.Modulus)))

		for _, q := range primes {
			require.Equal(t, q&((params.ringQ.N<<1)-1), uint64(1))
			require.True(t, IsPrime(q), q)
		}
	})
}

func testImportExportPolyString(t *testing.T) {

	t.Run(testString("ImportExportPolyString/", params.ringQ), func(t *testing.T) {

		p0 := params.uniformSamplerQ.ReadNew()
		p1 := params.ringQ.NewPoly()

		params.ringQ.SetCoefficientsString(params.ringQ.PolyToString(p0), p1)

		require.True(t, params.ringQ.Equal(p0, p1))
	})
}

func testDivFloorByLastModulusMany(t *testing.T) {

	t.Run(testString("DivFloorByLastModulusMany/", params.ringQ), func(t *testing.T) {

		coeffs := make([]*big.Int, params.ringQ.N)
		for i := uint64(0); i < params.ringQ.N; i++ {
			coeffs[i] = RandInt(params.ringQ.ModulusBigint)
			coeffs[i].Quo(coeffs[i], NewUint(10))
		}

		nbRescals := len(params.ringQ.Modulus) - 1

		coeffsWant := make([]*big.Int, params.ringQ.N)
		for i := range coeffs {
			coeffsWant[i] = new(big.Int).Set(coeffs[i])
			for j := 0; j < nbRescals; j++ {
				coeffsWant[i].Quo(coeffsWant[i], NewUint(params.ringQ.Modulus[len(params.ringQ.Modulus)-1-j]))
			}
		}

		polTest := params.ringQ.NewPoly()
		polWant := params.ringQ.NewPoly()

		params.ringQ.SetCoefficientsBigint(coeffs, polTest)
		params.ringQ.SetCoefficientsBigint(coeffsWant, polWant)

		params.ringQ.DivFloorByLastModulusMany(polTest, uint64(nbRescals))
		for i := uint64(0); i < params.ringQ.N; i++ {
			for j := 0; j < len(params.ringQ.Modulus)-nbRescals; j++ {
				require.Equalf(t, polWant.Coeffs[j][i], polTest.Coeffs[j][i], "coeff %v Qi%v = %s", i, j, coeffs[i].String())
			}
		}
	})
}

func testDivRoundByLastModulusMany(t *testing.T) {

	t.Run(testString("DivRoundByLastModulusMany/", params.ringQ), func(t *testing.T) {

		coeffs := make([]*big.Int, params.ringQ.N)
		for i := uint64(0); i < params.ringQ.N; i++ {
			coeffs[i] = RandInt(params.ringQ.ModulusBigint)
			coeffs[i].Quo(coeffs[i], NewUint(10))
		}

		nbRescals := len(params.ringQ.Modulus) - 1

		coeffsWant := make([]*big.Int, params.ringQ.N)
		for i := range coeffs {
			coeffsWant[i] = new(big.Int).Set(coeffs[i])
			for j := 0; j < nbRescals; j++ {
				DivRound(coeffsWant[i], NewUint(params.ringQ.Modulus[len(params.ringQ.Modulus)-1-j]), coeffsWant[i])
			}
		}

		polTest := params.ringQ.NewPoly()
		polWant := params.ringQ.NewPoly()

		params.ringQ.SetCoefficientsBigint(coeffs, polTest)
		params.ringQ.SetCoefficientsBigint(coeffsWant, polWant)

		params.ringQ.DivRoundByLastModulusMany(polTest, uint64(nbRescals))
		for i := uint64(0); i < params.ringQ.N; i++ {
			for j := 0; j < len(params.ringQ.Modulus)-nbRescals; j++ {
				require.Equalf(t, polWant.Coeffs[j][i], polTest.Coeffs[j][i], "coeff %v Qi%v = %s", i, j, coeffs[i].String())
			}
		}
	})
}

func testMarshalBinary(t *testing.T) {

	t.Run(testString("MarshalBinary/Ring/", params.ringQ), func(t *testing.T) {

		data, _ := params.ringQ.MarshalBinary()

		ringQTest := new(Ring)
		ringQTest.UnmarshalBinary(data)

		require.Equal(t, ringQTest.N, params.ringQ.N)
		require.Equal(t, ringQTest.Modulus, params.ringQ.Modulus)
	})

	t.Run(testString("MarshalBinary/Poly/", params.ringQ), func(t *testing.T) {

		p := params.uniformSamplerQ.ReadNew()
		pTest := params.ringQ.NewPoly()

		data, _ := p.MarshalBinary()

		_ = pTest.UnmarshalBinary(data)

		for i := range params.ringQ.Modulus {
			require.Equal(t, p.Coeffs[i][:params.ringQ.N], pTest.Coeffs[i][:params.ringQ.N])
		}
	})
}

func testUniformSampler(t *testing.T) {

	t.Run(testString("UniformSampler/Read/", params.ringQ), func(t *testing.T) {
		pol := params.ringQ.NewPoly()
		params.uniformSamplerQ.Read(pol)
		for i := uint64(0); i < params.ringQ.N; i++ {
			for j, qi := range params.ringQ.Modulus {
				require.False(t, pol.Coeffs[j][i] > qi)
			}
		}
	})

	t.Run(testString("UniformSampler/ReadNew/", params.ringQ), func(t *testing.T) {
		pol := params.uniformSamplerQ.ReadNew()
		for i := uint64(0); i < params.ringQ.N; i++ {
			for j, qi := range params.ringQ.Modulus {
				require.False(t, pol.Coeffs[j][i] > qi)
			}
		}
	})
}

func testGaussianSampler(t *testing.T) {

	t.Run(testString("GaussianSampler/", params.ringQ), func(t *testing.T) {
		gaussianSampler := NewGaussianSampler(params.prng, params.ringQ, DefaultSigma, DefaultBound)
		pol := gaussianSampler.ReadNew()

		for i := uint64(0); i < params.ringQ.N; i++ {
			for j, qi := range params.ringQ.Modulus {
				require.False(t, uint64(DefaultBound) < pol.Coeffs[j][i] && pol.Coeffs[j][i] < (qi-uint64(DefaultBound)))
			}
		}
	})
}

func testTernarySampler(t *testing.T) {

	t.Run(testString("TernarySampler/", params.ringQ), func(t *testing.T) {

		countOne := uint64(0)
		countZer := uint64(0)
		countMOn := uint64(0)

		pol := params.ringQ.NewPoly()

		rho := 1.0 / 3

		prng, err := utils.NewPRNG()
		if err != nil {
			panic(err)
		}
		ternarySampler := NewTernarySampler(prng, params.ringQ, rho, false)

		ternarySampler.Read(pol)

		for i := range pol.Coeffs[0] {
			if pol.Coeffs[0][i] == params.ringQ.Modulus[0]-1 {
				countMOn++
			}

			if pol.Coeffs[0][i] == 0 {
				countZer++
			}

			if pol.Coeffs[0][i] == 1 {
				countOne++
			}
		}

		threshold := 0.066

		ratio := math.Round(float64(countOne+countMOn)/float64(countZer)*100.0) / 100.0

		min := ((1 - rho) / rho) * (1.0 - threshold)
		max := ((1 - rho) / rho) * (1.0 + threshold)

		require.Greater(t, ratio, min)
		require.Less(t, ratio, max)
	})
}

func testModularReduction(t *testing.T) {

	t.Run(testString("ModularReduction/BRed/", params.ringQ), func(t *testing.T) {

		for j, q := range params.ringQ.Modulus {

			bigQ := NewUint(q)

			for i := 0; i < 65536; i++ {
				x := rand.Uint64() % q
				y := rand.Uint64() % q

				result := NewUint(x)
				result.Mul(result, NewUint(y))
				result.Mod(result, bigQ)

				require.Equalf(t, BRed(x, y, q, params.ringQ.BredParams[j]), result.Uint64(), "x = %v, y=%v", x, y)
			}
		}
	})

	t.Run(testString("ModularReduction/MRed/", params.ringQ), func(t *testing.T) {

		for j := range params.ringQ.Modulus {

			q := params.ringQ.Modulus[j]

			bigQ := NewUint(q)

			for i := 0; i < 65536; i++ {

				x := rand.Uint64() % q
				y := rand.Uint64() % q

				result := NewUint(x)
				result.Mul(result, NewUint(y))
				result.Mod(result, bigQ)

				require.Equalf(t, MRed(x, MForm(y, q, params.ringQ.BredParams[j]), q, params.ringQ.MredParams[j]), result.Uint64(), "x = %v, y=%v", x, y)
			}
		}
	})
}

func testGaloisShift(t *testing.T) {

	t.Run(testString("GaloisShift/", params.ringQ), func(t *testing.T) {

		pWant := params.uniformSamplerQ.ReadNew()
		pTest := pWant.CopyNew()

		params.ringQ.BitReverse(pTest, pTest)
		params.ringQ.InvNTT(pTest, pTest)
		params.ringQ.Rotate(pTest, 1, pTest)
		params.ringQ.NTT(pTest, pTest)
		params.ringQ.BitReverse(pTest, pTest)
		params.ringQ.Reduce(pTest, pTest)

		params.ringQ.Shift(pWant, 1, pWant)

		for i := range params.ringQ.Modulus {
			require.Equal(t, pTest.Coeffs[i][:params.ringQ.N], pWant.Coeffs[i][:params.ringQ.N])
		}
	})
}

func testMForm(t *testing.T) {

	t.Run(testString("MForm/", params.ringQ), func(t *testing.T) {

		polWant := params.uniformSamplerQ.ReadNew()
		polTest := params.ringQ.NewPoly()

		params.ringQ.MForm(polWant, polTest)
		params.ringQ.InvMForm(polTest, polTest)

		require.True(t, params.ringQ.Equal(polWant, polTest))
	})
}

func testMulScalarBigint(t *testing.T) {

	t.Run(testString("MulScalarBigint/", params.ringQ), func(t *testing.T) {

		polWant := params.uniformSamplerQ.ReadNew()
		polTest := polWant.CopyNew()

		rand1 := RandUniform(params.prng, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)
		rand2 := RandUniform(params.prng, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)

		scalarBigint := NewUint(rand1)
		scalarBigint.Mul(scalarBigint, NewUint(rand2))

		params.ringQ.MulScalar(polWant, rand1, polWant)
		params.ringQ.MulScalar(polWant, rand2, polWant)
		params.ringQ.MulScalarBigint(polTest, scalarBigint, polTest)

		require.True(t, params.ringQ.Equal(polWant, polTest))
	})
}

func testMulPoly(t *testing.T) {

	p1 := params.uniformSamplerQ.ReadNew()
	p2 := params.uniformSamplerQ.ReadNew()
	p3Test := params.ringQ.NewPoly()
	p3Want := params.ringQ.NewPoly()

	params.ringQ.Reduce(p1, p1)
	params.ringQ.Reduce(p2, p2)

	params.ringQ.MulPolyNaive(p1, p2, p3Want)

	t.Run(testString("MulPoly/Barrett/", params.ringQ), func(t *testing.T) {

		params.ringQ.MulPoly(p1, p2, p3Test)

		require.Equal(t, p3Want.Coeffs[0][:params.ringQ.N], p3Test.Coeffs[0][:params.ringQ.N])
	})

	t.Run(testString("MulPoly/Montgomery/", params.ringQ), func(t *testing.T) {

		params.ringQ.MForm(p1, p1)
		params.ringQ.MForm(p2, p2)

		params.ringQ.MulPolyMontgomery(p1, p2, p3Test)

		params.ringQ.InvMForm(p3Test, p3Test)

		require.Equal(t, p3Want.Coeffs[0][:params.ringQ.N], p3Test.Coeffs[0][:params.ringQ.N])
	})
}

func testExtendBasis(t *testing.T) {

	t.Run(testString("ExtendBasis/", params.ringQ), func(t *testing.T) {

		basisextender := NewFastBasisExtender(params.ringQ, params.ringP)

		coeffs := make([]*big.Int, params.ringQ.N)
		for i := uint64(0); i < params.ringQ.N; i++ {
			coeffs[i] = RandInt(params.ringQ.ModulusBigint)
		}

		Pol := params.ringQ.NewPoly()
		PolTest := params.ringP.NewPoly()
		PolWant := params.ringP.NewPoly()

		params.ringQ.SetCoefficientsBigint(coeffs, Pol)
		params.ringP.SetCoefficientsBigint(coeffs, PolWant)

		basisextender.ModUpSplitQP(uint64(len(params.ringQ.Modulus)-1), Pol, PolTest)

		for i := range params.ringP.Modulus {
			require.Equal(t, PolTest.Coeffs[i][:params.ringQ.N], PolWant.Coeffs[i][:params.ringQ.N])
		}
	})
}

func testScaling(t *testing.T) {

	t.Run(testString("Scaling/Simple/", params.ringQ), func(t *testing.T) {

		rescaler := NewSimpleScaler(T, params.ringQ)

		coeffs := make([]*big.Int, params.ringQ.N)
		for i := uint64(0); i < params.ringQ.N; i++ {
			coeffs[i] = RandInt(params.ringQ.ModulusBigint)
		}

		coeffsWant := make([]*big.Int, params.ringQ.N)
		for i := range coeffs {
			coeffsWant[i] = new(big.Int).Set(coeffs[i])
			coeffsWant[i].Mul(coeffsWant[i], NewUint(T))
			DivRound(coeffsWant[i], params.ringQ.ModulusBigint, coeffsWant[i])
			coeffsWant[i].Mod(coeffsWant[i], NewUint(T))
		}

		PolTest := params.ringQ.NewPoly()

		params.ringQ.SetCoefficientsBigint(coeffs, PolTest)

		rescaler.DivByQOverTRounded(PolTest, PolTest)

		for i := uint64(0); i < params.ringQ.N; i++ {
			require.Equal(t, PolTest.Coeffs[0][i], coeffsWant[i].Uint64())
		}
	})

	t.Run(testString("Scaling/RNS", params.ringQ), func(t *testing.T) {

		scaler := NewRNSScaler(T, params.ringQ)

		coeffs := make([]*big.Int, params.ringQ.N)
		for i := uint64(0); i < params.ringQ.N; i++ {
			coeffs[i] = RandInt(params.ringQ.ModulusBigint)
		}

		coeffsWant := make([]*big.Int, params.ringQ.N)
		for i := range coeffs {
			coeffsWant[i] = new(big.Int).Set(coeffs[i])
			coeffsWant[i].Mul(coeffsWant[i], NewUint(T))
			DivRound(coeffsWant[i], params.ringQ.ModulusBigint, coeffsWant[i])
			coeffsWant[i].Mod(coeffsWant[i], NewUint(T))
		}

		polyQ := params.ringQ.NewPoly()
		polyT := NewPoly(params.ringQ.N, 1)
		params.ringQ.SetCoefficientsBigint(coeffs, polyQ)

		scaler.DivByQOverTRounded(polyQ, polyT)

		for i := uint64(0); i < params.ringQ.N; i++ {
			require.Equal(t, polyT.Coeffs[0][i], coeffsWant[i].Uint64())
		}
	})
}

func testMultByMonomial(t *testing.T) {

	t.Run(testString("MultByMonomial/", params.ringQ), func(t *testing.T) {

		p1 := params.uniformSamplerQ.ReadNew()

		p3Test := params.ringQ.NewPoly()
		p3Want := params.ringQ.NewPoly()

		params.ringQ.MultByMonomial(p1, 1, p3Test)
		params.ringQ.MultByMonomial(p3Test, 8, p3Test)

		params.ringQ.MultByMonomial(p1, 9, p3Want)

		require.Equal(t, p3Want.Coeffs[0][:params.ringQ.N], p3Test.Coeffs[0][:params.ringQ.N])
	})
}
