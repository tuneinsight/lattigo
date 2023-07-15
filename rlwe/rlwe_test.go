package rlwe

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"runtime"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
	"github.com/tuneinsight/lattigo/v4/utils/buffer"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
	"github.com/tuneinsight/lattigo/v4/utils/structs"
)

var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")

func testString(params Parameters, levelQ, levelP, bpw2 int, opname string) string {
	return fmt.Sprintf("%s/logN=%d/Qi=%d/Pi=%d/Pw2=%d/NTT=%t/RingType=%s",
		opname,
		params.LogN(),
		levelQ+1,
		levelP+1,
		bpw2,
		params.NTTFlag(),
		params.RingType())
}

type DumDum struct {
	A int
	B float64
}

func TestRLWE(t *testing.T) {

	var err error

	defaultParamsLiteral := testParamsLiteral

	if *flagParamString != "" {
		var jsonParams TestParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			t.Fatal(err)
		}
		defaultParamsLiteral = []TestParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, paramsLit := range defaultParamsLiteral[:] {

		for _, NTTFlag := range []bool{true, false}[:] {

			for _, RingType := range []ring.Type{ring.Standard, ring.ConjugateInvariant}[:] {

				paramsLit.NTTFlag = NTTFlag
				paramsLit.RingType = RingType

				var params Parameters
				if params, err = NewParametersFromLiteral(paramsLit.ParametersLiteral); err != nil {
					t.Fatal(err)
				}

				tc, err := NewTestContext(params)
				require.NoError(t, err)

				testParameters(tc, t)
				testKeyGenerator(tc, paramsLit.BaseTwoDecomposition, t)
				testMarshaller(tc, t)
				testWriteAndRead(tc, paramsLit.BaseTwoDecomposition, t)

				for _, level := range []int{0, params.MaxLevel()}[:] {

					for _, testSet := range []func(tc *TestContext, level, bpw2 int, t *testing.T){
						testEncryptor,
						testGadgetProduct,
						testApplyEvaluationKey,
						testAutomorphism,
						testLinearTransform,
					} {
						testSet(tc, level, paramsLit.BaseTwoDecomposition, t)
						runtime.GC()
					}
				}
			}
		}
	}

	testUserDefinedParameters(t)
}

type TestContext struct {
	params Parameters
	kgen   *KeyGenerator
	enc    EncryptorInterface
	dec    *Decryptor
	sk     *SecretKey
	pk     *PublicKey
	eval   *Evaluator
}

func testUserDefinedParameters(t *testing.T) {

	t.Run("Parameters/UnmarshalJSON", func(t *testing.T) {

		var err error
		// checks that ckks.Parameters can be unmarshalled with log-moduli definition without error
		dataWithLogModuli := []byte(`{"LogN":13,"LogQ":[50,50],"LogP":[60]}`)
		var paramsWithLogModuli Parameters
		err = json.Unmarshal(dataWithLogModuli, &paramsWithLogModuli)
		require.Nil(t, err)
		require.Equal(t, 2, paramsWithLogModuli.QCount())
		require.Equal(t, 1, paramsWithLogModuli.PCount())
		require.Equal(t, ring.Standard, paramsWithLogModuli.RingType()) // Omitting the RingType field should result in a standard instance
		require.True(t, paramsWithLogModuli.Xe() == DefaultXe)          // Omitting Xe should result in Default being used
		require.True(t, paramsWithLogModuli.Xs() == DefaultXs)          // Omitting Xs should result in Default being used

		// checks that ckks.Parameters can be unmarshalled with log-moduli definition with empty or omitted P without error
		for _, dataWithLogModuliNoP := range [][]byte{
			[]byte(`{"LogN":13,"LogQ":[50,50],"LogP":[],"RingType": "ConjugateInvariant"}`),
			[]byte(`{"LogN":13,"LogQ":[50,50],"RingType": "ConjugateInvariant"}`),
		} {
			var paramsWithLogModuliNoP Parameters
			err = json.Unmarshal(dataWithLogModuliNoP, &paramsWithLogModuliNoP)
			require.Nil(t, err)
			require.Equal(t, 2, paramsWithLogModuliNoP.QCount())
			require.Equal(t, 0, paramsWithLogModuliNoP.PCount())
			require.Equal(t, ring.ConjugateInvariant, paramsWithLogModuliNoP.RingType())
		}

		// checks that one can provide custom parameters for the secret-key and error distributions
		dataWithCustomSecrets := []byte(`{"LogN":13,"LogQ":[50,50],"LogP":[60],"Xs":{"Type":"Ternary", "H":5462},"Xe":{"Type":"DiscreteGaussian","Sigma":6.4,"Bound":38}}`)
		var paramsWithCustomSecrets Parameters
		err = json.Unmarshal(dataWithCustomSecrets, &paramsWithCustomSecrets)
		require.Nil(t, err)
		require.True(t, paramsWithCustomSecrets.Xe() == ring.DiscreteGaussian{Sigma: 6.4, Bound: 38})
		require.True(t, paramsWithCustomSecrets.Xs() == ring.Ternary{H: 5462})

		var paramsWithBadDist Parameters
		// checks that providing an ambiguous gaussian distribution yields an error
		dataWithBadDist := []byte(`{"LogN":13,"LogQ":[50,50],"LogP":[60],"Xs":{"Type":"DiscreteGaussian", "Sigma":3.2}}`)
		err = json.Unmarshal(dataWithBadDist, &paramsWithBadDist)
		require.NotNil(t, err)
		require.Equal(t, paramsWithBadDist, Parameters{})

		// checks that providing an ambiguous ternary distribution yields an error
		dataWithBadDist = []byte(`{"LogN":13,"LogQ":[50,50],"LogP":[60],"Xs":{"Type":"Ternary", "H":5462,"P":0.3}}`)

		err = json.Unmarshal(dataWithBadDist, &paramsWithBadDist)
		require.NotNil(t, err)
		require.Equal(t, paramsWithBadDist, Parameters{})
	})

}

func NewTestContext(params Parameters) (tc *TestContext, err error) {
	kgen := NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()

	pk, err := kgen.GenPublicKeyNew(sk)
	if err != nil {
		return nil, err
	}

	eval := NewEvaluator(params, nil)

	enc, err := NewEncryptor(params, sk)
	if err != nil {
		return nil, err
	}

	dec, err := NewDecryptor(params, sk)
	if err != nil {
		return nil, err
	}

	return &TestContext{
		params: params,
		kgen:   kgen,
		sk:     sk,
		pk:     pk,
		enc:    enc,
		dec:    dec,
		eval:   eval,
	}, nil
}

func testParameters(tc *TestContext, t *testing.T) {

	params := tc.params

	t.Run(testString(params, params.MaxLevelQ(), params.MaxLevelP(), 0, "ModInvGaloisElement"), func(t *testing.T) {

		N := params.N()
		mask := params.RingQ().NthRoot() - 1

		for i := 1; i < N>>1; i++ {
			galEl := params.GaloisElement(i)
			inv := params.ModInvGaloisElement(galEl)
			res := (inv * galEl) & mask
			require.Equal(t, uint64(1), res)
		}
	})
}

func testKeyGenerator(tc *TestContext, bpw2 int, t *testing.T) {

	params := tc.params
	kgen := tc.kgen
	sk := tc.sk
	pk := tc.pk

	// Checks that the secret-key has exactly params.h non-zero coefficients
	t.Run(testString(params, params.MaxLevelQ(), params.MaxLevelP(), bpw2, "KeyGenerator/GenSecretKey"), func(t *testing.T) {

		switch xs := params.Xs().(type) {
		case ring.Ternary:
			if xs.P != 0 {
				t.Skip("cannot run test for probabilistic ternary distribution")
			}
		default:
			t.Skip("cannot run test for non ternary distribution")
		}

		skINTT := NewSecretKey(params)

		if params.PCount() > 0 {
			params.RingP().AtLevel(sk.LevelP()).INTT(sk.Value.P, skINTT.Value.P)
			for i := range skINTT.Value.P.Coeffs {
				var zeros int
				for j := range skINTT.Value.P.Coeffs[i] {
					if skINTT.Value.P.Coeffs[i][j] == 0 {
						zeros++
					}
				}
				require.Equal(t, params.ringP.N(), zeros+params.XsHammingWeight())
			}
		}

		params.RingQ().AtLevel(sk.LevelQ()).INTT(sk.Value.Q, skINTT.Value.Q)
		for i := range skINTT.Value.Q.Coeffs {
			var zeros int
			for j := range skINTT.Value.Q.Coeffs[i] {
				if skINTT.Value.Q.Coeffs[i][j] == 0 {
					zeros++
				}
			}
			require.Equal(t, params.ringQ.N(), zeros+params.XsHammingWeight())
		}

	})

	// Checks that sum([-as + e, a] + [as])) <= N * 6 * sigma
	t.Run(testString(params, params.MaxLevelQ(), params.MaxLevelP(), bpw2, "KeyGenerator/GenPublicKey"), func(t *testing.T) {

		if params.PCount() > 0 {

			ringQP := params.RingQP()

			zero := ringQP.NewPoly()

			ringQP.MulCoeffsMontgomery(sk.Value, pk.Value[1], zero)
			ringQP.Add(zero, pk.Value[0], zero)
			ringQP.INTT(zero, zero)
			ringQP.IMForm(zero, zero)

			require.GreaterOrEqual(t, math.Log2(params.NoiseFreshSK())+1, params.RingQ().Log2OfStandardDeviation(zero.Q))
			require.GreaterOrEqual(t, math.Log2(params.NoiseFreshSK())+1, params.RingP().Log2OfStandardDeviation(zero.P))
		} else {

			ringQ := params.RingQ()

			zero := ringQ.NewPoly()

			ringQ.MulCoeffsMontgomeryThenAdd(sk.Value.Q, pk.Value[1].Q, zero)
			ringQ.Add(zero, pk.Value[0].Q, zero)
			ringQ.INTT(zero, zero)
			ringQ.IMForm(zero, zero)

			require.GreaterOrEqual(t, math.Log2(params.NoiseFreshSK())+1, params.RingQ().Log2OfStandardDeviation(zero))
		}
	})

	var levelsQ = []int{0}
	var levelsP = []int{0}

	if params.MaxLevelQ() > 0 {
		levelsQ = append(levelsQ, params.MaxLevelQ())
	}

	if params.MaxLevelP() > 0 {
		levelsP = append(levelsP, params.MaxLevelP())
	}

	for _, levelQ := range levelsQ {

		for _, levelP := range levelsP {

			evkParams := EvaluationKeyParameters{LevelQ: levelQ, LevelP: levelP, BaseTwoDecomposition: bpw2}

			// Checks that EvaluationKeys are en encryption under the output key
			// of the RNS decomposition of the input key by
			// 1) Decrypting the RNS decomposed input key
			// 2) Reconstructing the key
			// 3) Checking that the difference with the input key has a small norm
			t.Run(testString(params, levelQ, levelP, bpw2, "KeyGenerator/GenEvaluationKey"), func(t *testing.T) {

				skOut := kgen.GenSecretKeyNew()

				decompRNS := params.DecompRNS(levelQ, levelP)
				decompPW2 := params.DecompPw2(levelQ, levelP, bpw2)

				evk := NewEvaluationKey(params, evkParams)

				// Generates Decomp([-asIn + w*P*sOut + e, a])
				kgen.GenEvaluationKey(sk, skOut, evk)

				require.Equal(t, decompRNS*decompPW2, len(evk.Value)*len(evk.Value[0])) // checks that decomposition size is correct

				require.GreaterOrEqual(t, math.Log2(math.Sqrt(float64(decompRNS))*params.NoiseFreshSK())+1, NoiseEvaluationKey(evk, sk, skOut, params))
			})

			t.Run(testString(params, levelQ, levelP, bpw2, "KeyGenerator/GenRelinearizationKey"), func(t *testing.T) {

				decompRNS := params.DecompRNS(levelQ, levelP)
				decompPW2 := params.DecompPw2(levelQ, levelP, bpw2)

				rlk := NewRelinearizationKey(params, evkParams)

				// Generates Decomp([-asIn + w*P*sOut + e, a])
				kgen.GenRelinearizationKey(sk, rlk)

				require.Equal(t, decompRNS*decompPW2, len(rlk.Value)*len(rlk.Value[0])) // checks that decomposition size is correct

				require.GreaterOrEqual(t, math.Log2(math.Sqrt(float64(decompRNS))*params.NoiseFreshSK())+1, NoiseRelinearizationKey(rlk, sk, params))
			})

			t.Run(testString(params, levelQ, levelP, bpw2, "KeyGenerator/GenGaloisKey"), func(t *testing.T) {

				decompRNS := params.DecompRNS(levelQ, levelP)
				decompPW2 := params.DecompPw2(levelQ, levelP, bpw2)

				gk := NewGaloisKey(params, evkParams)

				// Generates Decomp([-asIn + w*P*sOut + e, a])
				kgen.GenGaloisKey(ring.GaloisGen, sk, gk)

				require.Equal(t, decompRNS*decompPW2, len(gk.Value)*len(gk.Value[0])) // checks that decomposition size is correct

				require.GreaterOrEqual(t, math.Log2(math.Sqrt(float64(decompRNS))*params.NoiseFreshSK())+1, NoiseGaloisKey(gk, sk, params))
			})
		}
	}
}

func testEncryptor(tc *TestContext, level, bpw2 int, t *testing.T) {

	params := tc.params
	kgen := tc.kgen
	sk, pk := tc.sk, tc.pk
	enc := tc.enc
	dec := tc.dec

	t.Run(testString(params, level, params.MaxLevelP(), bpw2, "Encryptor/Encrypt/Pk"), func(t *testing.T) {
		ringQ := params.RingQ().AtLevel(level)

		pt := NewPlaintext(params, level)
		ct := NewCiphertext(params, 1, level)

		encPk, err := enc.WithKey(pk)
		require.NoError(t, err)

		require.NoError(t, encPk.Encrypt(pt, ct))

		dec.Decrypt(ct, pt)

		if pt.IsNTT {
			ringQ.INTT(pt.Value, pt.Value)
		}

		require.GreaterOrEqual(t, math.Log2(params.NoiseFreshPK())+1, ringQ.Log2OfStandardDeviation(pt.Value))
	})

	t.Run(testString(params, level, params.MaxLevelP(), bpw2, "Encryptor/Encrypt/Pk/ShallowCopy"), func(t *testing.T) {
		enc1, err := enc.WithKey(pk)
		require.NoError(t, err)

		enc2 := enc1.ShallowCopy()
		pkEnc1, pkEnc2 := enc1.(*EncryptorPublicKey), enc2.(*EncryptorPublicKey)
		require.True(t, pkEnc1.params.Equal(pkEnc2.params))
		require.True(t, pkEnc1.pk == pkEnc2.pk)
		require.False(t, (pkEnc1.basisextender == pkEnc2.basisextender) && (pkEnc1.basisextender != nil) && (pkEnc2.basisextender != nil))
		require.False(t, pkEnc1.encryptorBuffers == pkEnc2.encryptorBuffers)
		require.False(t, pkEnc1.xsSampler == pkEnc2.xsSampler)
		require.False(t, pkEnc1.xeSampler == pkEnc2.xeSampler)
	})

	t.Run(testString(params, level, params.MaxLevelP(), bpw2, "Encryptor/Encrypt/Sk"), func(t *testing.T) {
		ringQ := params.RingQ().AtLevel(level)

		pt := NewPlaintext(params, level)
		ct := NewCiphertext(params, 1, level)

		enc.Encrypt(pt, ct)
		dec.Decrypt(ct, pt)

		if pt.IsNTT {
			ringQ.INTT(pt.Value, pt.Value)
		}
		require.GreaterOrEqual(t, math.Log2(params.NoiseFreshSK())+1, ringQ.Log2OfStandardDeviation(pt.Value))
	})

	t.Run(testString(params, level, params.MaxLevelP(), bpw2, "Encryptor/Encrypt/Sk/PRNG"), func(t *testing.T) {
		ringQ := params.RingQ().AtLevel(level)

		pt := NewPlaintext(params, level)

		enc, err := NewPRNGEncryptor(params, sk)
		require.NoError(t, err)

		ct := NewCiphertext(params, 1, level)

		prng1, _ := sampling.NewKeyedPRNG([]byte{'a', 'b', 'c'})
		prng2, _ := sampling.NewKeyedPRNG([]byte{'a', 'b', 'c'})

		enc.WithPRNG(prng1).Encrypt(pt, ct)

		samplerQ := ring.NewUniformSampler(prng2, ringQ)

		require.True(t, ringQ.Equal(ct.Value[1], samplerQ.ReadNew()))

		dec.Decrypt(ct, pt)

		if pt.IsNTT {
			ringQ.INTT(pt.Value, pt.Value)
		}

		require.GreaterOrEqual(t, math.Log2(params.NoiseFreshSK())+1, ringQ.Log2OfStandardDeviation(pt.Value))
	})

	t.Run(testString(params, level, params.MaxLevelP(), bpw2, "Encrypt/Sk/ShallowCopy"), func(t *testing.T) {
		enc1, err := NewEncryptor(params, sk)
		require.NoError(t, err)

		enc2 := enc1.ShallowCopy()
		skEnc1, skEnc2 := enc1.(*EncryptorSecretKey), enc2.(*EncryptorSecretKey)
		require.True(t, skEnc1.params.Equal(skEnc2.params))
		require.True(t, skEnc1.sk == skEnc2.sk)
		require.False(t, (skEnc1.basisextender == skEnc2.basisextender) && (skEnc1.basisextender != nil) && (skEnc2.basisextender != nil))
		require.False(t, skEnc1.encryptorBuffers == skEnc2.encryptorBuffers)
		require.False(t, skEnc1.xsSampler == skEnc2.xsSampler)
		require.False(t, skEnc1.xeSampler == skEnc2.xeSampler)
	})

	t.Run(testString(params, level, params.MaxLevelP(), bpw2, "Encrypt/WithKey/Sk->Sk"), func(t *testing.T) {
		sk2 := kgen.GenSecretKeyNew()
		enc1, err := NewEncryptor(params, sk)
		require.NoError(t, err)

		enc2, err := enc1.WithKey(sk2)
		require.NoError(t, err)

		skEnc1, skEnc2 := enc1.(*EncryptorSecretKey), enc2.(*EncryptorSecretKey)
		require.True(t, skEnc1.params.Equal(skEnc2.params))
		require.True(t, skEnc1.sk.Equal(sk))
		require.True(t, skEnc2.sk.Equal(sk2))
		require.True(t, skEnc1.basisextender == skEnc2.basisextender)
		require.True(t, skEnc1.encryptorBuffers == skEnc2.encryptorBuffers)
		require.True(t, skEnc1.xsSampler == skEnc2.xsSampler)
		require.True(t, skEnc1.xeSampler == skEnc2.xeSampler)
	})
}

func testApplyEvaluationKey(tc *TestContext, level, bpw2 int, t *testing.T) {

	params := tc.params
	sk := tc.sk
	kgen := tc.kgen
	eval := tc.eval
	enc := tc.enc
	dec := tc.dec

	var NoiseBound = float64(params.LogN())

	t.Run(testString(params, level, params.MaxLevelP(), bpw2, "Evaluator/ApplyEvaluationKey/SameDegree"), func(t *testing.T) {

		skOut := kgen.GenSecretKeyNew()

		pt := NewPlaintext(params, level)

		ct := NewCiphertext(params, 1, level)

		enc.Encrypt(pt, ct)

		// Test that Dec(KS(Enc(ct, sk), skOut), skOut) has a small norm
		evk, err := kgen.GenEvaluationKeyNew(sk, skOut)
		require.NoError(t, err)

		eval.ApplyEvaluationKey(ct, evk, ct)

		dec, err := NewDecryptor(params, skOut)
		require.NoError(t, err)

		dec.Decrypt(ct, pt)

		ringQ := params.RingQ().AtLevel(level)

		if pt.IsNTT {
			ringQ.INTT(pt.Value, pt.Value)
		}

		require.GreaterOrEqual(t, NoiseBound, ringQ.Log2OfStandardDeviation(pt.Value))
	})

	t.Run(testString(params, level, params.MaxLevelP(), bpw2, "Evaluator/ApplyEvaluationKey/LargeToSmall"), func(t *testing.T) {

		paramsLargeDim := params

		paramsSmallDim, err := NewParametersFromLiteral(ParametersLiteral{
			LogN:     paramsLargeDim.LogN() - 1,
			Q:        paramsLargeDim.Q(),
			P:        []uint64{0x1ffffffff6c80001, 0x1ffffffff6140001}[:paramsLargeDim.PCount()], // some other P to test that the modulus is correctly extended in the keygen
			RingType: paramsLargeDim.RingType(),
		})

		require.Nil(t, err)

		kgenLargeDim := kgen
		skLargeDim := sk
		kgenSmallDim := NewKeyGenerator(paramsSmallDim)
		skSmallDim := kgenSmallDim.GenSecretKeyNew()

		evk, err := kgenLargeDim.GenEvaluationKeyNew(skLargeDim, skSmallDim)
		require.NoError(t, err)

		enc, err := NewEncryptor(paramsLargeDim, skLargeDim)
		require.NoError(t, err)

		ctLargeDim := enc.EncryptZeroNew(level)

		ctSmallDim := NewCiphertext(paramsSmallDim, 1, level)

		// skLarge -> skSmall embeded in N
		eval.ApplyEvaluationKey(ctLargeDim, evk, ctSmallDim)

		// Decrypts with smaller dimension key
		dec, err := NewDecryptor(paramsSmallDim, skSmallDim)
		require.NoError(t, err)

		ptSmallDim := dec.DecryptNew(ctSmallDim)

		ringQSmallDim := paramsSmallDim.RingQ().AtLevel(level)
		if ptSmallDim.IsNTT {
			ringQSmallDim.INTT(ptSmallDim.Value, ptSmallDim.Value)
		}

		require.GreaterOrEqual(t, NoiseBound, ringQSmallDim.Log2OfStandardDeviation(ptSmallDim.Value))
	})

	t.Run(testString(params, level, params.MaxLevelP(), bpw2, "Evaluator/ApplyEvaluationKey/SmallToLarge"), func(t *testing.T) {

		paramsLargeDim := params

		paramsSmallDim, err := NewParametersFromLiteral(ParametersLiteral{
			LogN:     paramsLargeDim.LogN() - 1,
			Q:        paramsLargeDim.Q(),
			P:        []uint64{0x1ffffffff6c80001, 0x1ffffffff6140001}[:paramsLargeDim.PCount()], // some other P to test that the modulus is correctly extended in the keygen
			RingType: paramsLargeDim.RingType(),
		})

		require.Nil(t, err)

		kgenLargeDim := kgen
		skLargeDim := sk
		kgenSmallDim := NewKeyGenerator(paramsSmallDim)
		skSmallDim := kgenSmallDim.GenSecretKeyNew()

		evk, err := kgenLargeDim.GenEvaluationKeyNew(skSmallDim, skLargeDim)
		require.NoError(t, err)

		enc, err := NewEncryptor(paramsSmallDim, skSmallDim)
		require.NoError(t, err)

		ctSmallDim := enc.EncryptZeroNew(level)

		ctLargeDim := NewCiphertext(paramsLargeDim, 1, level)

		eval.ApplyEvaluationKey(ctSmallDim, evk, ctLargeDim)

		ptLargeDim := dec.DecryptNew(ctLargeDim)

		ringQLargeDim := paramsLargeDim.RingQ().AtLevel(level)
		if ptLargeDim.IsNTT {
			ringQLargeDim.INTT(ptLargeDim.Value, ptLargeDim.Value)
		}

		require.GreaterOrEqual(t, NoiseBound, ringQLargeDim.Log2OfStandardDeviation(ptLargeDim.Value))
	})
}

func testGadgetProduct(tc *TestContext, levelQ, bpw2 int, t *testing.T) {

	params := tc.params
	sk := tc.sk
	kgen := tc.kgen
	eval := tc.eval

	ringQ := params.RingQ().AtLevel(levelQ)

	prng, _ := sampling.NewKeyedPRNG([]byte{'a', 'b', 'c'})

	sampler := ring.NewUniformSampler(prng, ringQ)

	var NoiseBound = float64(params.LogN())

	levelsP := []int{0}

	if params.MaxLevelP() > 0 {
		levelsP = append(levelsP, params.MaxLevelP())
	}

	for _, levelP := range levelsP {

		evkParams := EvaluationKeyParameters{LevelQ: levelQ, LevelP: levelP, BaseTwoDecomposition: bpw2}

		t.Run(testString(params, levelQ, levelP, bpw2, "Evaluator/GadgetProduct"), func(t *testing.T) {

			skOut := kgen.GenSecretKeyNew()

			// Generates a random polynomial
			a := sampler.ReadNew()

			// Generate the receiver
			ct := NewCiphertext(params, 1, levelQ)

			evk := NewEvaluationKey(params, evkParams)

			// Generate the evaluationkey [-bs1 + s1, b]
			require.NoError(t, kgen.GenEvaluationKey(sk, skOut, evk))

			// Gadget product: ct = [-cs1 + as0 , c]
			eval.GadgetProduct(levelQ, a, &evk.GadgetCiphertext, ct)

			// pt = as0
			dec, err := NewDecryptor(params, skOut)
			require.NoError(t, err)

			pt := dec.DecryptNew(ct)

			ringQ := params.RingQ().AtLevel(levelQ)

			// pt = as1 - as1 = 0 (+ some noise)
			if !pt.IsNTT {
				ringQ.NTT(pt.Value, pt.Value)
				ringQ.NTT(a, a)
			}

			ringQ.MulCoeffsMontgomeryThenSub(a, sk.Value.Q, pt.Value)
			ringQ.INTT(pt.Value, pt.Value)

			require.GreaterOrEqual(t, NoiseBound, ringQ.Log2OfStandardDeviation(pt.Value))
		})

		t.Run(testString(params, levelQ, levelP, bpw2, "Evaluator/GadgetProductHoisted"), func(t *testing.T) {

			if bpw2 != 0 {
				t.Skip("method is unsupported for BaseTwoDecomposition != 0")
			}

			skOut := kgen.GenSecretKeyNew()

			// Generates a random polynomial
			a := sampler.ReadNew()

			// Generate the receiver
			ct := NewCiphertext(params, 1, levelQ)

			evk := NewEvaluationKey(params, evkParams)

			// Generate the evaluationkey [-bs1 + s1, b]
			kgen.GenEvaluationKey(sk, skOut, evk)

			//Decompose the ciphertext
			eval.DecomposeNTT(levelQ, levelP, levelP+1, a, ct.IsNTT, eval.BuffDecompQP)

			// Gadget product: ct = [-cs1 + as0 , c]
			eval.GadgetProductHoisted(levelQ, eval.BuffDecompQP, &evk.GadgetCiphertext, ct)

			// pt = as0
			dec, err := NewDecryptor(params, skOut)
			require.NoError(t, err)

			pt := dec.DecryptNew(ct)

			ringQ := params.RingQ().AtLevel(levelQ)

			// pt = as1 - as1 = 0 (+ some noise)
			if !pt.IsNTT {
				ringQ.NTT(pt.Value, pt.Value)
				ringQ.NTT(a, a)
			}

			ringQ.MulCoeffsMontgomeryThenSub(a, sk.Value.Q, pt.Value)
			ringQ.INTT(pt.Value, pt.Value)

			require.GreaterOrEqual(t, NoiseBound, ringQ.Log2OfStandardDeviation(pt.Value))
		})
	}
}

func testAutomorphism(tc *TestContext, level, bpw2 int, t *testing.T) {

	params := tc.params
	sk := tc.sk
	kgen := tc.kgen
	eval := tc.eval
	enc := tc.enc
	dec := tc.dec

	var NoiseBound = float64(params.LogN())

	t.Run(testString(params, level, params.MaxLevelP(), bpw2, "Evaluator/Automorphism"), func(t *testing.T) {

		// Generate a plaintext with values up to 2^30
		pt := genPlaintext(params, level, 1<<30)

		// Encrypt
		ct, err := enc.EncryptNew(pt)
		require.NoError(t, err)

		// Chooses a Galois Element (must be coprime with 2N)
		galEl := params.GaloisElement(-1)

		// Generate the GaloisKey
		gk, err := kgen.GenGaloisKeyNew(galEl, sk)
		require.NoError(t, err)

		// Allocate a new EvaluationKeySet and adds the GaloisKey
		evk := NewMemEvaluationKeySet(nil, gk)

		// Evaluate the automorphism
		eval.WithKey(evk).Automorphism(ct, galEl, ct)

		// Apply the same automorphism on the plaintext
		ringQ := params.RingQ().AtLevel(level)

		tmp := ringQ.NewPoly()
		if pt.IsNTT {
			ringQ.AutomorphismNTT(pt.Value, galEl, tmp)
		} else {
			ringQ.Automorphism(pt.Value, galEl, tmp)
		}

		// Decrypt
		dec.Decrypt(ct, pt)

		// Subract the permuted plaintext to the decrypted plaintext
		ringQ.Sub(pt.Value, tmp, pt.Value)

		// Switch out of NTT if required
		if pt.IsNTT {
			ringQ.INTT(pt.Value, pt.Value)
		}

		// Logs the noise
		require.GreaterOrEqual(t, NoiseBound, ringQ.Log2OfStandardDeviation(pt.Value))
	})

	t.Run(testString(params, level, params.MaxLevelP(), bpw2, "Evaluator/AutomorphismHoisted"), func(t *testing.T) {
		// Generate a plaintext with values up to 2^30
		pt := genPlaintext(params, level, 1<<30)

		// Encrypt
		ct, err := enc.EncryptNew(pt)
		require.NoError(t, err)

		// Chooses a Galois Element (must be coprime with 2N)
		galEl := params.GaloisElement(-1)

		// Generate the GaloisKey
		gk, err := kgen.GenGaloisKeyNew(galEl, sk)
		require.NoError(t, err)

		// Allocate a new EvaluationKeySet and adds the GaloisKey
		evk := NewMemEvaluationKeySet(nil, gk)

		//Decompose the ciphertext
		eval.DecomposeNTT(level, params.MaxLevelP(), params.MaxLevelP()+1, ct.Value[1], ct.IsNTT, eval.BuffDecompQP)

		// Evaluate the automorphism
		eval.WithKey(evk).AutomorphismHoisted(level, ct, eval.BuffDecompQP, galEl, ct)

		// Apply the same automorphism on the plaintext
		ringQ := params.RingQ().AtLevel(level)

		tmp := ringQ.NewPoly()
		if pt.IsNTT {
			ringQ.AutomorphismNTT(pt.Value, galEl, tmp)
		} else {
			ringQ.Automorphism(pt.Value, galEl, tmp)
		}

		// Decrypt
		dec.Decrypt(ct, pt)

		// Subract the permuted plaintext to the decrypted plaintext
		ringQ.Sub(pt.Value, tmp, pt.Value)

		// Switch out of NTT if required
		if pt.IsNTT {
			ringQ.INTT(pt.Value, pt.Value)
		}

		// Logs the noise
		require.GreaterOrEqual(t, NoiseBound, ringQ.Log2OfStandardDeviation(pt.Value))
	})

	t.Run(testString(params, level, params.MaxLevelP(), bpw2, "Evaluator/AutomorphismHoistedLazy"), func(t *testing.T) {
		// Generate a plaintext with values up to 2^30
		pt := genPlaintext(params, level, 1<<30)

		// Encrypt
		ct, err := enc.EncryptNew(pt)
		require.NoError(t, err)

		// Chooses a Galois Element (must be coprime with 2N)
		galEl := params.GaloisElement(-1)

		// Generate the GaloisKey
		gk, err := kgen.GenGaloisKeyNew(galEl, sk)
		require.NoError(t, err)

		// Allocate a new EvaluationKeySet and adds the GaloisKey
		evk := NewMemEvaluationKeySet(nil, gk)

		//Decompose the ciphertext
		eval.DecomposeNTT(level, params.MaxLevelP(), params.MaxLevelP()+1, ct.Value[1], ct.IsNTT, eval.BuffDecompQP)

		ctQP := NewOperandQP(params, 1, level, params.MaxLevelP())

		// Evaluate the automorphism
		eval.WithKey(evk).AutomorphismHoistedLazy(level, ct, eval.BuffDecompQP, galEl, ctQP)

		eval.ModDown(level, params.MaxLevelP(), ctQP, ct)

		// Apply the same automorphism on the plaintext
		ringQ := params.RingQ().AtLevel(level)

		tmp := ringQ.NewPoly()
		if pt.IsNTT {
			ringQ.AutomorphismNTT(pt.Value, galEl, tmp)
		} else {
			ringQ.Automorphism(pt.Value, galEl, tmp)
		}

		// Decrypt
		dec.Decrypt(ct, pt)

		// Subract the permuted plaintext to the decrypted plaintext
		ringQ.Sub(pt.Value, tmp, pt.Value)

		// Switch out of NTT if required
		if pt.IsNTT {
			ringQ.INTT(pt.Value, pt.Value)
		}

		// Logs the noise
		require.GreaterOrEqual(t, NoiseBound, ringQ.Log2OfStandardDeviation(pt.Value))
	})
}

func testLinearTransform(tc *TestContext, level, bpw2 int, t *testing.T) {

	params := tc.params
	sk := tc.sk
	kgen := tc.kgen
	eval := tc.eval
	enc := tc.enc
	dec := tc.dec

	t.Run(testString(params, level, params.MaxLevelP(), bpw2, "Evaluator/Expand"), func(t *testing.T) {

		if params.RingType() != ring.Standard {
			t.Skip("Expand not supported for ring.Type = ring.ConjugateInvariant")
		}

		pt := NewPlaintext(params, level)
		ringQ := params.RingQ().AtLevel(level)

		logN := 4
		logGap := 0
		gap := 1 << logGap

		values := make([]uint64, params.N())

		scale := 1 << 22

		for i := 0; i < 1<<logN; i++ { // embeds even coefficients only
			values[i] = uint64(i * scale)
		}

		for i := 0; i < pt.Level()+1; i++ {
			copy(pt.Value.Coeffs[i], values)
		}

		if pt.IsNTT {
			ringQ.NTT(pt.Value, pt.Value)
		}

		ctIn := NewCiphertext(params, 1, level)
		enc.Encrypt(pt, ctIn)

		// GaloisKeys
		var gks, err = kgen.GenGaloisKeysNew(params.GaloisElementsForExpand(logN), sk)
		require.NoError(t, err)

		evk := NewMemEvaluationKeySet(nil, gks...)

		eval := NewEvaluator(params, evk)

		ciphertexts, err := eval.WithKey(evk).Expand(ctIn, logN, logGap)
		require.NoError(t, err)

		Q := ringQ.ModuliChain()

		NoiseBound := float64(params.LogN() - logN)

		for i := range ciphertexts {

			dec.Decrypt(ciphertexts[i], pt)

			if pt.IsNTT {
				ringQ.INTT(pt.Value, pt.Value)
			}

			for j := 0; j < level+1; j++ {
				pt.Value.Coeffs[j][0] = ring.CRed(pt.Value.Coeffs[j][0]+Q[j]-values[i*gap], Q[j])
			}

			// Logs the noise
			require.GreaterOrEqual(t, NoiseBound, ringQ.Log2OfStandardDeviation(pt.Value))
		}
	})

	t.Run(testString(params, level, params.MaxLevelP(), bpw2, "Evaluator/Pack/LogGap=LogN"), func(t *testing.T) {

		if params.RingType() != ring.Standard {
			t.Skip("Pack not supported for ring.Type = ring.ConjugateInvariant")
		}

		pt := NewPlaintext(params, level)
		N := params.N()
		ringQ := tc.params.RingQ().AtLevel(level)
		gap := params.N() / 16

		ptPacked := NewPlaintext(params, level)
		ciphertexts := make(map[int]*Ciphertext)
		slotIndex := make(map[int]bool)
		for i := 0; i < N; i += gap {

			ciphertexts[i] = enc.EncryptZeroNew(level)

			scalar := (1 << 30) + uint64(i)*(1<<20)

			if ciphertexts[i].IsNTT {
				ringQ.AddScalar(ciphertexts[i].Value[0], scalar, ciphertexts[i].Value[0])
			} else {
				for j := 0; j < level+1; j++ {
					ciphertexts[i].Value[0].Coeffs[j][0] = ring.CRed(ciphertexts[i].Value[0].Coeffs[j][0]+scalar, ringQ.SubRings[j].Modulus)
				}
			}

			slotIndex[i] = true

			for j := 0; j < level+1; j++ {
				ptPacked.Value.Coeffs[j][i] = scalar
			}
		}

		// Galois Keys
		galEls, err := params.GaloisElementsForPack(params.LogN())
		require.NoError(t, err)

		gks, err := kgen.GenGaloisKeysNew(galEls, sk)
		require.NoError(t, err)

		evk := NewMemEvaluationKeySet(nil, gks...)

		ct, err := eval.WithKey(evk).Pack(ciphertexts, params.LogN(), false)
		require.NoError(t, err)

		dec.Decrypt(ct, pt)

		if pt.IsNTT {
			ringQ.INTT(pt.Value, pt.Value)
		}

		ringQ.Sub(pt.Value, ptPacked.Value, pt.Value)

		for i := 0; i < N; i++ {
			if i%gap != 0 {
				for j := 0; j < level+1; j++ {
					pt.Value.Coeffs[j][i] = 0
				}
			}
		}

		NoiseBound := 15.0

		// Logs the noise
		require.GreaterOrEqual(t, NoiseBound, ringQ.Log2OfStandardDeviation(pt.Value))
	})

	t.Run(testString(params, level, params.MaxLevelP(), bpw2, "Evaluator/Pack/LogGap=LogN-1"), func(t *testing.T) {

		if params.RingType() != ring.Standard {
			t.Skip("Pack not supported for ring.Type = ring.ConjugateInvariant")
		}

		pt := NewPlaintext(params, level)
		N := params.N()
		ringQ := tc.params.RingQ().AtLevel(level)

		ptPacked := NewPlaintext(params, level)
		ciphertexts := make(map[int]*Ciphertext)
		slotIndex := make(map[int]bool)
		for i := 0; i < N/2; i += params.N() / 16 {

			ciphertexts[i] = enc.EncryptZeroNew(level)

			scalar := (1 << 30) + uint64(i)*(1<<20)

			if ciphertexts[i].IsNTT {
				ringQ.INTT(ciphertexts[i].Value[0], ciphertexts[i].Value[0])
			}

			for j := 0; j < level+1; j++ {
				ciphertexts[i].Value[0].Coeffs[j][0] = ring.CRed(ciphertexts[i].Value[0].Coeffs[j][0]+scalar, ringQ.SubRings[j].Modulus)
				ciphertexts[i].Value[0].Coeffs[j][N/2] = ring.CRed(ciphertexts[i].Value[0].Coeffs[j][N/2]+scalar, ringQ.SubRings[j].Modulus)
			}

			if ciphertexts[i].IsNTT {
				ringQ.NTT(ciphertexts[i].Value[0], ciphertexts[i].Value[0])
			}

			slotIndex[i] = true

			for j := 0; j < level+1; j++ {
				ptPacked.Value.Coeffs[j][i] = scalar
				ptPacked.Value.Coeffs[j][i+N/2] = scalar
			}
		}

		// Galois Keys
		galEls, err := params.GaloisElementsForPack(params.LogN() - 1)
		require.NoError(t, err)

		gks, err := kgen.GenGaloisKeysNew(galEls, sk)
		require.NoError(t, err)

		evk := NewMemEvaluationKeySet(nil, gks...)

		ct, err := eval.WithKey(evk).Pack(ciphertexts, params.LogN()-1, true)
		require.NoError(t, err)

		dec.Decrypt(ct, pt)

		if pt.IsNTT {
			ringQ.INTT(pt.Value, pt.Value)
		}

		ringQ.Sub(pt.Value, ptPacked.Value, pt.Value)

		NoiseBound := 15.0

		// Logs the noise
		require.GreaterOrEqual(t, NoiseBound, ringQ.Log2OfStandardDeviation(pt.Value))
	})

	t.Run(testString(params, level, params.MaxLevelP(), bpw2, "Evaluator/InnerSum"), func(t *testing.T) {

		batch := 5
		n := 7

		ringQ := tc.params.RingQ().AtLevel(level)

		pt := genPlaintext(params, level, 1<<30)
		ptInnerSum := *pt.Value.CopyNew()
		ct, err := enc.EncryptNew(pt)
		require.NoError(t, err)

		// Galois Keys
		gks, err := kgen.GenGaloisKeysNew(params.GaloisElementsForInnerSum(batch, n), sk)
		require.NoError(t, err)

		evk := NewMemEvaluationKeySet(nil, gks...)

		eval.WithKey(evk).InnerSum(ct, batch, n, ct)

		dec.Decrypt(ct, pt)

		if pt.IsNTT {
			ringQ.INTT(pt.Value, pt.Value)
			ringQ.INTT(ptInnerSum, ptInnerSum)
		}

		polyTmp := ringQ.NewPoly()

		// Applies the same circuit (naively) on the plaintext
		polyInnerSum := *ptInnerSum.CopyNew()
		for i := 1; i < n; i++ {
			galEl := params.GaloisElement(i * batch)
			ringQ.Automorphism(ptInnerSum, galEl, polyTmp)
			ringQ.Add(polyInnerSum, polyTmp, polyInnerSum)
		}

		ringQ.Sub(pt.Value, polyInnerSum, pt.Value)

		NoiseBound := float64(params.LogN())

		// Logs the noise
		require.GreaterOrEqual(t, NoiseBound, ringQ.Log2OfStandardDeviation(pt.Value))

	})
}

func genPlaintext(params Parameters, level, max int) (pt *Plaintext) {

	N := params.N()

	step := float64(max) / float64(N)

	pt = NewPlaintext(params, level)

	for i := 0; i < level+1; i++ {
		c := pt.Value.Coeffs[i]
		for j := 0; j < N; j++ {
			c[j] = uint64(float64(j) * step)
		}
	}

	if pt.IsNTT {
		params.RingQ().AtLevel(level).NTT(pt.Value, pt.Value)
	}

	return
}

func testWriteAndRead(tc *TestContext, bpw2 int, t *testing.T) {

	params := tc.params

	sk, pk := tc.sk, tc.pk

	levelQ := params.MaxLevelQ()
	levelP := params.MaxLevelP()

	t.Run(testString(params, levelQ, levelP, bpw2, "WriteAndRead/Operand[ring.Poly]"), func(t *testing.T) {

		prng, _ := sampling.NewPRNG()
		sampler := ring.NewUniformSampler(prng, params.RingQ())

		op := Operand[ring.Poly]{
			Value: structs.Vector[ring.Poly]{
				sampler.ReadNew(),
				sampler.ReadNew(),
			},
			MetaData: &MetaData{
				IsNTT: params.NTTFlag(),
			},
		}

		buffer.RequireSerializerCorrect(t, &op)
	})

	t.Run(testString(params, levelQ, levelP, bpw2, "WriteAndRead/Operand[ringqp.Poly]"), func(t *testing.T) {

		prng, _ := sampling.NewPRNG()
		sampler := ringqp.NewUniformSampler(prng, *params.RingQP())

		op := Operand[ringqp.Poly]{
			Value: structs.Vector[ringqp.Poly]{
				sampler.ReadNew(),
				sampler.ReadNew(),
			},
			MetaData: &MetaData{
				IsNTT: params.NTTFlag(),
			},
		}

		buffer.RequireSerializerCorrect(t, &op)
	})

	t.Run(testString(params, levelQ, levelP, bpw2, "WriteAndRead/Plaintext"), func(t *testing.T) {
		prng, _ := sampling.NewPRNG()
		plaintextWant := NewPlaintext(params, levelQ)
		ring.NewUniformSampler(prng, params.RingQ()).Read(plaintextWant.Value)
		buffer.RequireSerializerCorrect(t, plaintextWant)
	})

	t.Run(testString(params, levelQ, levelP, bpw2, "WriteAndRead/Ciphertext"), func(t *testing.T) {

		prng, _ := sampling.NewPRNG()

		for degree := 0; degree < 4; degree++ {
			t.Run(fmt.Sprintf("degree=%d", degree), func(t *testing.T) {
				buffer.RequireSerializerCorrect(t, NewCiphertextRandom(prng, params, degree, levelQ))
			})
		}
	})

	t.Run(testString(params, levelQ, levelP, bpw2, "WriteAndRead/GadgetCiphertext"), func(t *testing.T) {

		rlk := NewRelinearizationKey(params, EvaluationKeyParameters{LevelQ: levelQ, LevelP: levelP, BaseTwoDecomposition: bpw2})

		tc.kgen.GenRelinearizationKey(tc.sk, rlk)

		buffer.RequireSerializerCorrect(t, &rlk.GadgetCiphertext)
	})

	t.Run(testString(params, levelQ, levelP, bpw2, "WriteAndRead/Sk"), func(t *testing.T) {
		buffer.RequireSerializerCorrect(t, sk)
	})

	t.Run(testString(params, levelQ, levelP, bpw2, "WriteAndRead/Pk"), func(t *testing.T) {
		buffer.RequireSerializerCorrect(t, pk)
	})

	t.Run(testString(params, levelQ, levelP, bpw2, "WriteAndRead/EvaluationKey"), func(t *testing.T) {
		evk, err := tc.kgen.GenEvaluationKeyNew(sk, sk)
		require.NoError(t, err)
		buffer.RequireSerializerCorrect(t, evk)
	})

	t.Run(testString(params, levelQ, levelP, bpw2, "WriteAndRead/RelinearizationKey"), func(t *testing.T) {
		rlk, err := tc.kgen.GenRelinearizationKeyNew(tc.sk)
		require.NoError(t, err)
		buffer.RequireSerializerCorrect(t, rlk)
	})

	t.Run(testString(params, levelQ, levelP, bpw2, "WriteAndRead/GaloisKey"), func(t *testing.T) {
		gk, err := tc.kgen.GenGaloisKeyNew(5, tc.sk)
		require.NoError(t, err)
		buffer.RequireSerializerCorrect(t, gk)
	})

	t.Run(testString(params, levelQ, levelP, bpw2, "WriteAndRead/EvaluationKeySet"), func(t *testing.T) {

		rlk, err := tc.kgen.GenRelinearizationKeyNew(tc.sk)
		require.NoError(t, err)
		galEl := uint64(5)
		gk, err := tc.kgen.GenGaloisKeyNew(galEl, tc.sk)
		require.NoError(t, err)

		buffer.RequireSerializerCorrect(t, &MemEvaluationKeySet{
			Rlk: rlk,
			Gks: map[uint64]*GaloisKey{galEl: gk},
		})
	})

	t.Run(testString(params, levelQ, levelP, bpw2, "WriteAndRead/PowerBasis"), func(t *testing.T) {

		prng, _ := sampling.NewPRNG()

		ct := NewCiphertextRandom(prng, params, 1, levelQ)

		basis := NewPowerBasis(ct, bignum.Chebyshev)

		basis.Value[2] = NewCiphertextRandom(prng, params, 1, levelQ)
		basis.Value[3] = NewCiphertextRandom(prng, params, 2, levelQ)
		basis.Value[4] = NewCiphertextRandom(prng, params, 1, levelQ)
		basis.Value[8] = NewCiphertextRandom(prng, params, 1, levelQ)

		buffer.RequireSerializerCorrect(t, &basis)
	})
}

func testMarshaller(tc *TestContext, t *testing.T) {

	//params := tc.params

	//sk, pk := tc.sk, tc.pk

	/*
		t.Run(testString(params, params.MaxLevel(), "Marshaller/Parameters/Binary"), func(t *testing.T) {
			bytes, err := params.MarshalBinary()

			require.Nil(t, err)
			var p Parameters
			require.Nil(t, p.UnmarshalBinary(bytes))
			require.Equal(t, params, p)
			require.Equal(t, params.RingQ(), p.RingQ())
		})

		t.Run(testString(params, params.MaxLevel(), "Marshaller/Parameters/JSON"), func(t *testing.T) {

			paramsLit := params.ParametersLiteral()

			paramsLit.DefaultScale = NewScale(1 << 45)

			var err error
			params, err = NewParametersFromLiteral(paramsLit)

			require.Nil(t, err)

			data, err := params.MarshalJSON()
			require.Nil(t, err)
			require.NotNil(t, data)

			var p Parameters
			require.Nil(t, p.UnmarshalJSON(data))

			require.Equal(t, params, p)
		})
	*/

	t.Run("Marshaller/MetaData", func(t *testing.T) {
		m := MetaData{PlaintextScale: NewScaleModT(1, 65537), IsNTT: true, IsMontgomery: true}

		data, err := m.MarshalBinary()
		require.Nil(t, err)
		require.NotNil(t, data)

		mHave := &MetaData{}

		require.Nil(t, mHave.UnmarshalBinary(data))

		require.True(t, m.Equal(mHave))
	})
}
