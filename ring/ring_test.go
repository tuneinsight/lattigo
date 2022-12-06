package ring

import (
	"flag"
	"fmt"
	"math/big"
	"testing"

	"github.com/tuneinsight/lattigo/v4/utils"

	"github.com/stretchr/testify/require"
)

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters). Overrides -short and requires -timeout=0.")

var T = uint64(0x3ee0001)
var DefaultSigma = float64(3.2)
var DefaultBound = int(6 * DefaultSigma)

func testString(opname string, ringQ *Ring) string {
	return fmt.Sprintf("%s/N=%d/limbs=%d", opname, ringQ.N(), ringQ.NbModuli())
}

type testParams struct {
	ringQ           *Ring
	ringP           *Ring
	prng            utils.PRNG
	uniformSamplerQ *UniformSampler
	uniformSamplerP *UniformSampler
}

func genTestParams(defaultParams Parameters) (tc *testParams, err error) {

	tc = new(testParams)

	if tc.ringQ, err = NewRing(1<<defaultParams.logN, defaultParams.qi); err != nil {
		return nil, err
	}
	if tc.ringP, err = NewRing(1<<defaultParams.logN, defaultParams.pi); err != nil {
		return nil, err
	}
	if tc.prng, err = utils.NewPRNG(); err != nil {
		return nil, err
	}
	tc.uniformSamplerQ = NewUniformSampler(tc.prng, tc.ringQ)
	tc.uniformSamplerP = NewUniformSampler(tc.prng, tc.ringP)
	return
}

func TestRing(t *testing.T) {

	var err error

	var defaultParams = DefaultParams[0:4] // the default test
	if testing.Short() {
		defaultParams = DefaultParams[0:2] // the short test suite
	}
	if *flagLongTest {
		defaultParams = DefaultParams // the long test suite
	}

	testNewRing(t)
	testShift(t)
	for _, defaultParam := range defaultParams[:] {

		var tc *testParams
		if tc, err = genTestParams(defaultParam); err != nil {
			t.Fatal(err)
		}

		testNTTConjugateInvariant(tc, t)
		testPRNG(tc, t)
		testGenerateNTTPrimes(tc, t)
		testDivFloorByLastModulusMany(tc, t)
		testDivRoundByLastModulusMany(tc, t)
		testMarshalBinary(tc, t)
		testUniformSampler(tc, t)
		testGaussianSampler(tc, t)
		testTernarySampler(tc, t)
		testModularReduction(tc, t)
		testMForm(tc, t)
		testMulScalarBigint(tc, t)
		testExtendBasis(tc, t)
		testMultByMonomial(tc, t)

	}
}

func testNTTConjugateInvariant(tc *testParams, t *testing.T) {

	t.Run(testString("NTTConjugateInvariant", tc.ringQ), func(t *testing.T) {

		ringQ := tc.ringQ
		Q := ringQ.Moduli()
		N := ringQ.N()
		ringQ2N, _ := NewRing(N<<1, Q)
		ringQConjugateInvariant, _ := NewRingFromType(N, Q, ConjugateInvariant)

		sampler := NewUniformSampler(tc.prng, ringQ)
		p1 := sampler.ReadNew()
		p2 := ringQ2N.NewPoly()

		for i, qi := range Q {
			copy(p2.Coeffs[i], p1.Coeffs[i])
			for j := 1; j < N; j++ {
				p2.Coeffs[i][N*2-j] = qi - p2.Coeffs[i][j]
			}
		}

		ringQ2N.NTT(p2, p2)
		ringQ2N.MForm(p2, p2)
		ringQ2N.MulCoeffsMontgomery(p2, p2, p2)
		ringQ2N.InvMForm(p2, p2)
		ringQ2N.InvNTT(p2, p2)

		p1tmp := ringQ2N.NewPoly()

		ringQConjugateInvariant.NTT(p1, p1tmp)
		ringQConjugateInvariant.MForm(p1tmp, p1tmp)
		ringQConjugateInvariant.MulCoeffsMontgomery(p1tmp, p1tmp, p1tmp)
		ringQConjugateInvariant.InvMForm(p1tmp, p1tmp)
		ringQConjugateInvariant.InvNTT(p1tmp, p1)

		for j := range Q {
			for i := 0; i < N; i++ {
				require.Equal(t, p1.Coeffs[j][i], p2.Coeffs[j][i])
			}
		}
	})
}

func testNewRing(t *testing.T) {
	t.Run("NewRing", func(t *testing.T) {
		r, err := NewRing(0, nil)
		require.Nil(t, r)
		require.Error(t, err)

		r, err = NewRing(0, []uint64{})
		require.Nil(t, r)
		require.Error(t, err)

		r, err = NewRing(4, []uint64{})
		require.Nil(t, r)
		require.Error(t, err)

		r, err = NewRing(8, []uint64{})
		require.Nil(t, r)
		require.Error(t, err)

		r, err = NewRing(16, []uint64{7}) // Passing non NTT-enabling coeff modulus
		require.NotNil(t, r)              // Should still return a Ring instance
		require.Error(t, err)             // Should also return an error due to non NTT

		r, err = NewRing(16, []uint64{4}) // Passing non prime moduli
		require.NotNil(t, r)              // Should still return a Ring instance
		require.Error(t, err)             // Should also return an error due to non NTT

		r, err = NewRing(16, []uint64{97, 7}) // Passing a NTT-enabling and a non NTT-enabling coeff modulus
		require.NotNil(t, r)                  // Should still return a Ring instance
		require.Error(t, err)                 // Should also return an error due to non NTT

		r, err = NewRing(16, []uint64{97, 97}) // Passing non CRT-enabling coeff modulus
		require.Nil(t, r)                      // Should not return a Ring instance
		require.Error(t, err)

		r, err = NewRing(16, []uint64{97}) // Passing NTT-enabling coeff modulus
		require.NotNil(t, r)
		require.NoError(t, err)

	})
}

func testPRNG(tc *testParams, t *testing.T) {

	t.Run(testString("PRNG", tc.ringQ), func(t *testing.T) {

		var err error

		var prng1, prng2 utils.PRNG

		if prng1, err = utils.NewKeyedPRNG(nil); err != nil {
			t.Error(err)
		}

		if prng2, err = utils.NewKeyedPRNG(nil); err != nil {
			t.Error(err)
		}

		crsGenerator1 := NewUniformSampler(prng1, tc.ringQ)
		crsGenerator2 := NewUniformSampler(prng2, tc.ringQ)

		p0 := crsGenerator1.ReadNew()
		p1 := crsGenerator2.ReadNew()

		require.True(t, tc.ringQ.Equal(p0, p1))
	})

}

func testGenerateNTTPrimes(tc *testParams, t *testing.T) {

	t.Run(testString("GenerateNTTPrimes", tc.ringQ), func(t *testing.T) {

		NthRoot := tc.ringQ.N() << 1

		primes := GenerateNTTPrimes(55, NthRoot, tc.ringQ.NbModuli())

		for _, q := range primes {
			require.Equal(t, q&uint64(NthRoot-1), uint64(1))
			require.True(t, IsPrime(q), q)
		}
	})
}

func testDivFloorByLastModulusMany(tc *testParams, t *testing.T) {

	t.Run(testString("DivFloorByLastModulusMany", tc.ringQ), func(t *testing.T) {

		N := tc.ringQ.N()

		level := tc.ringQ.MaxLevel()
		coeffs := make([]*big.Int, N)
		for i := 0; i < N; i++ {
			coeffs[i] = RandInt(tc.ringQ.ModulusAtLevel[level])
			coeffs[i].Quo(coeffs[i], NewUint(10))
		}

		nbRescales := level

		coeffsWant := make([]*big.Int, N)
		for i := range coeffs {
			coeffsWant[i] = new(big.Int).Set(coeffs[i])
			for j := 0; j < nbRescales; j++ {
				coeffsWant[i].Quo(coeffsWant[i], NewUint(tc.ringQ.Tables[level-j].Modulus))
			}
		}

		polTest0 := tc.ringQ.NewPoly()
		polTest1 := tc.ringQ.NewPoly()
		polWant := tc.ringQ.NewPoly()
		buff := tc.ringQ.NewPoly()

		tc.ringQ.SetCoefficientsBigintLvl(level, coeffs, polTest0)
		tc.ringQ.SetCoefficientsBigintLvl(level, coeffsWant, polWant)

		tc.ringQ.DivFloorByLastModulusManyLvl(polTest0.Level(), nbRescales, polTest0, buff, polTest1)
		for i := 0; i < N; i++ {
			for j := 0; j < polTest0.Level()-nbRescales+1; j++ {
				require.Equalf(t, polWant.Coeffs[j][i], polTest1.Coeffs[j][i], "coeff %v Qi%v = %s", i, j, coeffs[i].String())
			}
		}
	})
}

func testDivRoundByLastModulusMany(tc *testParams, t *testing.T) {

	t.Run(testString("DivRoundByLastModulusMany", tc.ringQ), func(t *testing.T) {

		N := tc.ringQ.N()

		level := tc.ringQ.MaxLevel()
		coeffs := make([]*big.Int, N)
		for i := 0; i < N; i++ {
			coeffs[i] = RandInt(tc.ringQ.ModulusAtLevel[level])
			coeffs[i].Quo(coeffs[i], NewUint(10))
		}

		nbRescals := level

		coeffsWant := make([]*big.Int, N)
		for i := range coeffs {
			coeffsWant[i] = new(big.Int).Set(coeffs[i])
			for j := 0; j < nbRescals; j++ {
				DivRound(coeffsWant[i], NewUint(tc.ringQ.Tables[level-j].Modulus), coeffsWant[i])
			}
		}

		polTest0 := tc.ringQ.NewPoly()
		polTest1 := tc.ringQ.NewPoly()
		polWant := tc.ringQ.NewPoly()
		buff := tc.ringQ.NewPoly()

		tc.ringQ.SetCoefficientsBigintLvl(level, coeffs, polTest0)
		tc.ringQ.SetCoefficientsBigintLvl(level, coeffsWant, polWant)

		tc.ringQ.DivRoundByLastModulusManyLvl(polTest0.Level(), nbRescals, polTest0, buff, polTest1)
		for i := 0; i < N; i++ {
			for j := 0; j < polTest0.Level()-nbRescals+1; j++ {
				require.Equalf(t, polWant.Coeffs[j][i], polTest1.Coeffs[j][i], "coeff %v Qi%v = %s", i, j, coeffs[i].String())
			}
		}
	})
}

func testMarshalBinary(tc *testParams, t *testing.T) {

	t.Run(testString("MarshalBinary/Ring", tc.ringQ), func(t *testing.T) {

		var err error

		var data []byte
		if data, err = tc.ringQ.MarshalBinary(); err != nil {
			t.Error(err)
		}

		ringQTest := new(Ring)
		if err = ringQTest.UnmarshalBinary(data); err != nil {
			t.Error(err)
		}

		require.Equal(t, ringQTest, tc.ringQ)
	})

	t.Run(testString("MarshalBinary/Poly", tc.ringQ), func(t *testing.T) {

		var err error

		p := tc.uniformSamplerQ.ReadNew()

		var data []byte
		if data, err = p.MarshalBinary(); err != nil {
			t.Error(err)
		}

		pTest := new(Poly)
		if err = pTest.UnmarshalBinary(data); err != nil {
			t.Error(err)
		}

		for i := range tc.ringQ.Tables {
			require.Equal(t, p.Coeffs[i][:tc.ringQ.N()], pTest.Coeffs[i][:tc.ringQ.N()])
		}
	})
}

func testUniformSampler(tc *testParams, t *testing.T) {

	N := tc.ringQ.N()

	t.Run(testString("UniformSampler/Read", tc.ringQ), func(t *testing.T) {
		pol := tc.ringQ.NewPoly()
		tc.uniformSamplerQ.Read(pol)

		for i, qi := range tc.ringQ.Moduli() {
			coeffs := pol.Coeffs[i]
			for j := 0; j < N; j++ {
				require.False(t, coeffs[j] > qi)
			}
		}
	})

	t.Run(testString("UniformSampler/ReadNew", tc.ringQ), func(t *testing.T) {
		pol := tc.uniformSamplerQ.ReadNew()

		for i, qi := range tc.ringQ.Moduli() {
			coeffs := pol.Coeffs[i]
			for j := 0; j < N; j++ {
				require.False(t, coeffs[j] > qi)
			}
		}
	})
}

func testGaussianSampler(tc *testParams, t *testing.T) {

	N := tc.ringQ.N()

	t.Run(testString("GaussianSampler", tc.ringQ), func(t *testing.T) {
		gaussianSampler := NewGaussianSampler(tc.prng, tc.ringQ, DefaultSigma, DefaultBound)
		pol := gaussianSampler.ReadNew()

		bound := uint64(DefaultBound)
		for i, qi := range tc.ringQ.Moduli() {
			coeffs := pol.Coeffs[i]
			negbound := qi - uint64(DefaultBound)
			for j := 0; j < N; j++ {
				require.False(t, bound < coeffs[j] && coeffs[j] < negbound)
			}
		}
	})
}

func testTernarySampler(tc *testParams, t *testing.T) {

	for _, p := range []float64{.5, 1. / 3., 128. / 65536.} {
		t.Run(testString(fmt.Sprintf("TernarySampler/p=%1.2f", p), tc.ringQ), func(t *testing.T) {

			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			ternarySampler := NewTernarySampler(prng, tc.ringQ, p, false)

			pol := ternarySampler.ReadNew()
			for i, qi := range tc.ringQ.Moduli() {
				minOne := qi - 1
				for _, c := range pol.Coeffs[i] {
					require.True(t, c == 0 || c == minOne || c == 1)
				}
			}
		})
	}

	for _, p := range []int{0, 64, 96, 128, 256} {
		t.Run(testString(fmt.Sprintf("TernarySampler/hw=%d", p), tc.ringQ), func(t *testing.T) {

			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}

			ternarySampler := NewTernarySamplerWithHammingWeight(prng, tc.ringQ, p, false)

			checkPoly := func(pol *Poly) {
				for i := range tc.ringQ.Tables {
					hw := 0
					for _, c := range pol.Coeffs[i] {
						if c != 0 {
							hw++
						}
					}
					require.True(t, hw == p)
				}
			}

			pol := ternarySampler.ReadNew()

			checkPoly(pol)

			ternarySampler.Read(pol)

			checkPoly(pol)
		})
	}
}

func testModularReduction(tc *testParams, t *testing.T) {

	t.Run(testString("ModularReduction/BRed", tc.ringQ), func(t *testing.T) {

		var x, y uint64
		var bigQ, result *big.Int

		for j, q := range tc.ringQ.Moduli() {

			bigQ = NewUint(q)

			bredParams := tc.ringQ.Tables[j].BRedParams

			x = 1
			y = 1

			result = NewUint(x)
			result.Mul(result, NewUint(y))
			result.Mod(result, bigQ)

			require.Equalf(t, BRed(x, y, q, bredParams), result.Uint64(), "x = %v, y=%v", x, y)

			x = 1
			y = q - 1

			result = NewUint(x)
			result.Mul(result, NewUint(y))
			result.Mod(result, bigQ)

			require.Equalf(t, BRed(x, y, q, bredParams), result.Uint64(), "x = %v, y=%v", x, y)

			x = 1
			y = 0xFFFFFFFFFFFFFFFF

			result = NewUint(x)
			result.Mul(result, NewUint(y))
			result.Mod(result, bigQ)

			require.Equalf(t, BRed(x, y, q, bredParams), result.Uint64(), "x = %v, y=%v", x, y)

			x = q - 1
			y = q - 1

			result = NewUint(x)
			result.Mul(result, NewUint(y))
			result.Mod(result, bigQ)

			require.Equalf(t, BRed(x, y, q, bredParams), result.Uint64(), "x = %v, y=%v", x, y)

			x = q - 1
			y = 0xFFFFFFFFFFFFFFFF

			result = NewUint(x)
			result.Mul(result, NewUint(y))
			result.Mod(result, bigQ)

			require.Equalf(t, BRed(x, y, q, bredParams), result.Uint64(), "x = %v, y=%v", x, y)

			x = 0xFFFFFFFFFFFFFFFF
			y = 0xFFFFFFFFFFFFFFFF

			result = NewUint(x)
			result.Mul(result, NewUint(y))
			result.Mod(result, bigQ)

			require.Equalf(t, BRed(x, y, q, bredParams), result.Uint64(), "x = %v, y=%v", x, y)
		}
	})

	t.Run(testString("ModularReduction/MRed", tc.ringQ), func(t *testing.T) {

		var x, y uint64
		var bigQ, result *big.Int

		for j, q := range tc.ringQ.Moduli() {

			bigQ = NewUint(q)

			bredParams := tc.ringQ.Tables[j].BRedParams
			mredparams := tc.ringQ.Tables[j].MRedParams

			x = 1
			y = 1

			result = NewUint(x)
			result.Mul(result, NewUint(y))
			result.Mod(result, bigQ)

			require.Equalf(t, MRed(x, MForm(y, q, bredParams), q, mredparams), result.Uint64(), "x = %v, y=%v", x, y)

			x = 1
			y = q - 1

			result = NewUint(x)
			result.Mul(result, NewUint(y))
			result.Mod(result, bigQ)

			require.Equalf(t, MRed(x, MForm(y, q, bredParams), q, mredparams), result.Uint64(), "x = %v, y=%v", x, y)

			x = 1
			y = 0xFFFFFFFFFFFFFFFF

			result = NewUint(x)
			result.Mul(result, NewUint(y))
			result.Mod(result, bigQ)

			require.Equalf(t, MRed(x, MForm(y, q, bredParams), q, mredparams), result.Uint64(), "x = %v, y=%v", x, y)

			x = q - 1
			y = q - 1

			result = NewUint(x)
			result.Mul(result, NewUint(y))
			result.Mod(result, bigQ)

			require.Equalf(t, MRed(x, MForm(y, q, bredParams), q, mredparams), result.Uint64(), "x = %v, y=%v", x, y)

			x = q - 1
			y = 0xFFFFFFFFFFFFFFFF

			result = NewUint(x)
			result.Mul(result, NewUint(y))
			result.Mod(result, bigQ)

			require.Equalf(t, MRed(x, MForm(y, q, bredParams), q, mredparams), result.Uint64(), "x = %v, y=%v", x, y)

			x = 0xFFFFFFFFFFFFFFFF
			y = 0xFFFFFFFFFFFFFFFF

			result = NewUint(x)
			result.Mul(result, NewUint(y))
			result.Mod(result, bigQ)

			require.Equalf(t, MRed(x, MForm(y, q, bredParams), q, mredparams), result.Uint64(), "x = %v, y=%v", x, y)
		}
	})
}

func testMForm(tc *testParams, t *testing.T) {

	t.Run(testString("MForm", tc.ringQ), func(t *testing.T) {

		polWant := tc.uniformSamplerQ.ReadNew()
		polTest := tc.ringQ.NewPoly()

		tc.ringQ.MForm(polWant, polTest)
		tc.ringQ.InvMForm(polTest, polTest)

		require.True(t, tc.ringQ.Equal(polWant, polTest))
	})
}

func testMulScalarBigint(tc *testParams, t *testing.T) {

	t.Run(testString("MulScalarBigint", tc.ringQ), func(t *testing.T) {

		polWant := tc.uniformSamplerQ.ReadNew()
		polTest := polWant.CopyNew()

		rand1 := RandUniform(tc.prng, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)
		rand2 := RandUniform(tc.prng, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)

		scalarBigint := NewUint(rand1)
		scalarBigint.Mul(scalarBigint, NewUint(rand2))

		tc.ringQ.MulScalar(polWant, rand1, polWant)
		tc.ringQ.MulScalar(polWant, rand2, polWant)
		tc.ringQ.MulScalarBigint(polTest, scalarBigint, polTest)

		require.True(t, tc.ringQ.Equal(polWant, polTest))
	})
}

func testExtendBasis(tc *testParams, t *testing.T) {

	N := tc.ringQ.N()

	t.Run(testString("ModUp/QToP", tc.ringQ), func(t *testing.T) {

		basisextender := NewBasisExtender(tc.ringQ, tc.ringP)

		levelQ := tc.ringQ.MaxLevel() - 1
		levelP := tc.ringP.MaxLevel() - 1

		Q := tc.ringQ.ModulusAtLevel[levelQ]

		coeffs := make([]*big.Int, N)
		for i := 0; i < N; i++ {
			coeffs[i] = RandInt(Q)
		}

		PolQHave := tc.ringQ.NewPolyLvl(levelQ)
		PolPTest := tc.ringP.NewPolyLvl(levelP)
		PolPWant := tc.ringP.NewPolyLvl(levelP)

		tc.ringQ.SetCoefficientsBigintLvl(levelQ, coeffs, PolQHave)
		tc.ringP.SetCoefficientsBigintLvl(levelP, coeffs, PolPWant)

		basisextender.ModUpQtoP(levelQ, levelP, PolQHave, PolPTest)
		tc.ringP.Reduce(PolPTest, PolPTest)

		for i := 0; i < levelQ+1; i++ {
			for j := 0; j < N; j++ {
				require.Equal(t, PolPWant.Coeffs[i][j], PolPTest.Coeffs[i][j])
			}
		}
	})

	t.Run(testString("ModUp/PToQ", tc.ringQ), func(t *testing.T) {

		basisextender := NewBasisExtender(tc.ringQ, tc.ringP)

		levelQ := tc.ringQ.MaxLevel() - 1
		levelP := tc.ringP.MaxLevel() - 1

		P := tc.ringP.ModulusAtLevel[levelP]

		coeffs := make([]*big.Int, N)
		for i := 0; i < N; i++ {
			coeffs[i] = RandInt(P)
		}

		PolQWant := tc.ringQ.NewPolyLvl(levelQ)
		PolQTest := tc.ringP.NewPolyLvl(levelQ)
		PolPHave := tc.ringP.NewPolyLvl(levelP)

		tc.ringQ.SetCoefficientsBigintLvl(levelQ, coeffs, PolQWant)
		tc.ringP.SetCoefficientsBigintLvl(levelP, coeffs, PolPHave)

		basisextender.ModUpPtoQ(levelQ, levelP, PolPHave, PolQTest)
		tc.ringQ.Reduce(PolQTest, PolQTest)

		for i := 0; i < levelQ+1; i++ {
			for j := 0; j < N; j++ {
				require.Equal(t, PolQWant.Coeffs[i][j], PolQTest.Coeffs[i][j])
			}
		}
	})

	t.Run(testString("ModDown/QPToQ", tc.ringQ), func(t *testing.T) {

		basisextender := NewBasisExtender(tc.ringQ, tc.ringP)

		levelQ := tc.ringQ.MaxLevel() - 1
		levelP := tc.ringP.MaxLevel() - 1

		Q := tc.ringQ.ModulusAtLevel[levelQ]
		P := tc.ringP.ModulusAtLevel[levelP]

		QP := new(big.Int).Mul(Q, P)

		coeffs := make([]*big.Int, N)
		for i := 0; i < N; i++ {
			coeffs[i] = RandInt(QP)
			coeffs[i].Quo(coeffs[i], NewUint(10))
		}

		coeffsWant := make([]*big.Int, N)
		for i := range coeffs {
			coeffsWant[i] = new(big.Int).Set(coeffs[i])
			coeffsWant[i].Quo(coeffsWant[i], P)
		}

		PolQHave := tc.ringQ.NewPolyLvl(levelQ)
		PolPHave := tc.ringP.NewPolyLvl(levelP)
		PolQWant := tc.ringP.NewPolyLvl(levelQ)

		tc.ringQ.SetCoefficientsBigintLvl(levelQ, coeffs, PolQHave)
		tc.ringP.SetCoefficientsBigintLvl(levelP, coeffs, PolPHave)
		tc.ringQ.SetCoefficientsBigintLvl(levelQ, coeffsWant, PolQWant)

		basisextender.ModDownQPtoQ(levelQ, levelP, PolQHave, PolPHave, PolQHave)
		tc.ringQ.Reduce(PolQHave, PolQHave)

		for i := 0; i < levelQ+1; i++ {
			for j := 0; j < N; j++ {
				require.Equal(t, PolQHave.Coeffs[i][j], PolQWant.Coeffs[i][j])
			}
		}
	})

	t.Run(testString("ModDown/QPToP", tc.ringQ), func(t *testing.T) {

		basisextender := NewBasisExtender(tc.ringQ, tc.ringP)

		levelQ := tc.ringQ.MaxLevel() - 1
		levelP := tc.ringP.MaxLevel() - 1

		Q := tc.ringQ.ModulusAtLevel[levelQ]
		P := tc.ringP.ModulusAtLevel[levelP]

		QP := new(big.Int).Mul(Q, P)

		coeffs := make([]*big.Int, N)
		for i := 0; i < N; i++ {
			coeffs[i] = RandInt(QP)
			coeffs[i].Quo(coeffs[i], NewUint(10))
		}

		coeffsWant := make([]*big.Int, N)
		for i := range coeffs {
			coeffsWant[i] = new(big.Int).Set(coeffs[i])
			coeffsWant[i].Quo(coeffsWant[i], P)
		}

		PolQHave := tc.ringQ.NewPolyLvl(levelQ)
		PolPHave := tc.ringP.NewPolyLvl(levelP)
		PolQWant := tc.ringP.NewPolyLvl(levelQ)

		tc.ringQ.SetCoefficientsBigintLvl(levelQ, coeffs, PolQHave)
		tc.ringP.SetCoefficientsBigintLvl(levelP, coeffs, PolPHave)
		tc.ringQ.SetCoefficientsBigintLvl(levelQ, coeffsWant, PolQWant)

		basisextender.ModDownQPtoQ(levelQ, levelP, PolQHave, PolPHave, PolQHave)
		tc.ringQ.Reduce(PolQHave, PolQHave)

		for i := 0; i < levelQ+1; i++ {
			for j := 0; j < N; j++ {
				require.Equal(t, PolQHave.Coeffs[i][j], PolQWant.Coeffs[i][j])
			}
		}
	})
}

func testMultByMonomial(tc *testParams, t *testing.T) {

	t.Run(testString("MultByMonomial", tc.ringQ), func(t *testing.T) {

		p1 := tc.uniformSamplerQ.ReadNew()

		p3Test := tc.ringQ.NewPoly()
		p3Want := tc.ringQ.NewPoly()

		tc.ringQ.MultByMonomial(p1, 1, p3Test)
		tc.ringQ.MultByMonomial(p3Test, 8, p3Test)

		tc.ringQ.MultByMonomial(p1, 9, p3Want)

		require.Equal(t, p3Want.Coeffs[0][:tc.ringQ.N()], p3Test.Coeffs[0][:tc.ringQ.N()])
	})
}

func testShift(t *testing.T) {

	r, _ := NewRing(16, []uint64{97})
	p1, p2 := r.NewPoly(), r.NewPoly()

	for i := range p1.Coeffs[0] {
		p1.Coeffs[0][i] = uint64(i)
	}

	r.Shift(p1, 3, p2)
	require.Equal(t, p2.Coeffs[0], []uint64{3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2})

}
