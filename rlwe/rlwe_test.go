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
	"github.com/tuneinsight/lattigo/v4/utils/bignum/polynomial"
	"github.com/tuneinsight/lattigo/v4/utils/buffer"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")

func testString(params Parameters, level int, opname string) string {
	return fmt.Sprintf("%s/logN=%d/Qi=%d/Pi=%d/Bit=%d/NTT=%t/Level=%d/RingType=%s",
		opname,
		params.LogN(),
		params.QCount(),
		params.PCount(),
		params.Pow2Base(),
		params.NTTFlag(),
		level,
		params.RingType())
}

type DumDum struct {
	A int
	B float64
}

func TestRLWE(t *testing.T) {

	var err error

	defaultParamsLiteral := TestParamsLiteral[:]

	if *flagParamString != "" {
		var jsonParams ParametersLiteral
		if err = json.Unmarshal([]byte(*flagParamString), &jsonParams); err != nil {
			t.Fatal(err)
		}
		defaultParamsLiteral = []ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, paramsLit := range defaultParamsLiteral[:] {

		for _, NTTFlag := range []bool{true, false}[:] {

			for _, RingType := range []ring.Type{ring.Standard, ring.ConjugateInvariant}[:] {

				paramsLit.NTTFlag = NTTFlag
				paramsLit.RingType = RingType

				var params Parameters
				if params, err = NewParametersFromLiteral(paramsLit); err != nil {
					t.Fatal(err)
				}

				tc := NewTestContext(params)

				testParameters(tc, t)
				testKeyGenerator(tc, t)
				testMarshaller(tc, t)
				testWriteAndRead(tc, t)

				for _, level := range []int{0, params.MaxLevel()}[:] {

					for _, testSet := range []func(tc *TestContext, level int, t *testing.T){
						testEncryptor,
						testGadgetProduct,
						testApplyEvaluationKey,
						testAutomorphism,
						testLinearTransform,
					} {
						testSet(tc, level, t)
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

func NewTestContext(params Parameters) (tc *TestContext) {
	kgen := NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	pk := kgen.GenPublicKeyNew(sk)
	eval := NewEvaluator(params, nil)

	return &TestContext{
		params: params,
		kgen:   kgen,
		sk:     sk,
		pk:     pk,
		enc:    NewEncryptor(params, sk),
		dec:    NewDecryptor(params, sk),
		eval:   eval,
	}
}

func testParameters(tc *TestContext, t *testing.T) {

	params := tc.params

	t.Run(testString(params, params.MaxLevel(), "ModInvGaloisElement"), func(t *testing.T) {

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

func testKeyGenerator(tc *TestContext, t *testing.T) {

	params := tc.params
	kgen := tc.kgen
	sk := tc.sk
	pk := tc.pk

	// Checks that the secret-key has exactly params.h non-zero coefficients
	t.Run(testString(params, params.MaxLevel(), "KeyGenerator/GenSecretKey"), func(t *testing.T) {

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
	t.Run(testString(params, params.MaxLevel(), "KeyGenerator/GenPublicKey"), func(t *testing.T) {

		if params.PCount() > 0 {

			ringQP := params.RingQP()

			zero := ringQP.NewPoly()

			ringQP.MulCoeffsMontgomery(&sk.Value, &pk.Value[1], zero)
			ringQP.Add(zero, &pk.Value[0], zero)
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

	// Checks that EvaluationKeys are en encryption under the output key
	// of the RNS decomposition of the input key by
	// 1) Decrypting the RNS decomposed input key
	// 2) Reconstructing the key
	// 3) Checking that the difference with the input key has a small norm
	t.Run(testString(params, params.MaxLevel(), "KeyGenerator/GenEvaluationKey"), func(t *testing.T) {

		skOut := kgen.GenSecretKeyNew()
		levelQ, levelP := params.MaxLevelQ(), params.MaxLevelP()
		decompPW2 := params.DecompPw2(levelQ, levelP)
		decompRNS := params.DecompRNS(levelQ, levelP)

		// Generates Decomp([-asIn + w*P*sOut + e, a])
		evk := kgen.GenEvaluationKeyNew(sk, skOut)

		require.Equal(t, decompRNS*decompPW2, len(evk.Value)*len(evk.Value[0])) // checks that decomposition size is correct

		require.True(t, EvaluationKeyIsCorrect(evk, sk, skOut, params, math.Log2(math.Sqrt(float64(decompRNS))*params.NoiseFreshSK())+1))
	})

	t.Run(testString(params, params.MaxLevel(), "KeyGenerator/GenRelinearizationKey"), func(t *testing.T) {

		levelQ, levelP := params.MaxLevelQ(), params.MaxLevelP()
		decompPW2 := params.DecompPw2(levelQ, levelP)
		decompRNS := params.DecompRNS(levelQ, levelP)

		// Generates Decomp([-asIn + w*P*sOut + e, a])
		rlk := kgen.GenRelinearizationKeyNew(sk)

		require.Equal(t, decompRNS*decompPW2, len(rlk.Value)*len(rlk.Value[0])) // checks that decomposition size is correct

		require.True(t, RelinearizationKeyIsCorrect(rlk, sk, params, math.Log2(math.Sqrt(float64(decompRNS))*params.NoiseFreshSK())+1))
	})

	t.Run(testString(params, params.MaxLevel(), "KeyGenerator/GenGaloisKey"), func(t *testing.T) {

		levelQ, levelP := params.MaxLevelQ(), params.MaxLevelP()
		decompPW2 := params.DecompPw2(levelQ, levelP)
		decompRNS := params.DecompRNS(levelQ, levelP)

		// Generates Decomp([-asIn + w*P*sOut + e, a])
		gk := kgen.GenGaloisKeyNew(ring.GaloisGen, sk)

		require.Equal(t, decompRNS*decompPW2, len(gk.Value)*len(gk.Value[0])) // checks that decomposition size is correct

		require.True(t, GaloisKeyIsCorrect(gk, sk, params, math.Log2(math.Sqrt(float64(decompRNS))*params.NoiseFreshSK())+1))
	})
}

func testEncryptor(tc *TestContext, level int, t *testing.T) {

	params := tc.params
	kgen := tc.kgen
	sk, pk := tc.sk, tc.pk
	enc := tc.enc
	dec := tc.dec

	t.Run(testString(params, level, "Encryptor/Encrypt/Pk"), func(t *testing.T) {
		ringQ := params.RingQ().AtLevel(level)

		pt := NewPlaintext(params, level)
		ct := NewCiphertext(params, 1, level)

		enc.WithKey(pk).Encrypt(pt, ct)
		dec.Decrypt(ct, pt)

		if pt.IsNTT {
			ringQ.INTT(pt.Value, pt.Value)
		}

		require.GreaterOrEqual(t, math.Log2(params.NoiseFreshPK())+1, ringQ.Log2OfStandardDeviation(pt.Value))
	})

	t.Run(testString(params, level, "Encryptor/Encrypt/Pk/ShallowCopy"), func(t *testing.T) {
		enc1 := enc.WithKey(pk)
		enc2 := enc1.ShallowCopy()
		pkEnc1, pkEnc2 := enc1.(*EncryptorPublicKey), enc2.(*EncryptorPublicKey)
		require.True(t, pkEnc1.params.Equal(pkEnc2.params))
		require.True(t, pkEnc1.pk == pkEnc2.pk)
		require.False(t, (pkEnc1.basisextender == pkEnc2.basisextender) && (pkEnc1.basisextender != nil) && (pkEnc2.basisextender != nil))
		require.False(t, pkEnc1.encryptorBuffers == pkEnc2.encryptorBuffers)
		require.False(t, pkEnc1.xsSampler == pkEnc2.xsSampler)
		require.False(t, pkEnc1.xeSampler == pkEnc2.xeSampler)
	})

	t.Run(testString(params, level, "Encryptor/Encrypt/Sk"), func(t *testing.T) {
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

	t.Run(testString(params, level, "Encryptor/Encrypt/Sk/PRNG"), func(t *testing.T) {
		ringQ := params.RingQ().AtLevel(level)

		pt := NewPlaintext(params, level)

		enc := NewPRNGEncryptor(params, sk)
		ct := NewCiphertext(params, 1, level)

		prng1, _ := sampling.NewKeyedPRNG([]byte{'a', 'b', 'c'})
		prng2, _ := sampling.NewKeyedPRNG([]byte{'a', 'b', 'c'})

		enc.WithPRNG(prng1).Encrypt(pt, ct)

		samplerQ := ring.NewUniformSampler(prng2, ringQ)

		require.True(t, ringQ.Equal(&ct.Value[1], samplerQ.ReadNew()))

		dec.Decrypt(ct, pt)

		if pt.IsNTT {
			ringQ.INTT(pt.Value, pt.Value)
		}

		require.GreaterOrEqual(t, math.Log2(params.NoiseFreshSK())+1, ringQ.Log2OfStandardDeviation(pt.Value))
	})

	t.Run(testString(params, level, "Encrypt/Sk/ShallowCopy"), func(t *testing.T) {
		enc1 := NewEncryptor(params, sk)
		enc2 := enc1.ShallowCopy()
		skEnc1, skEnc2 := enc1.(*EncryptorSecretKey), enc2.(*EncryptorSecretKey)
		require.True(t, skEnc1.params.Equal(skEnc2.params))
		require.True(t, skEnc1.sk == skEnc2.sk)
		require.False(t, (skEnc1.basisextender == skEnc2.basisextender) && (skEnc1.basisextender != nil) && (skEnc2.basisextender != nil))
		require.False(t, skEnc1.encryptorBuffers == skEnc2.encryptorBuffers)
		require.False(t, skEnc1.xsSampler == skEnc2.xsSampler)
		require.False(t, skEnc1.xeSampler == skEnc2.xeSampler)
	})

	t.Run(testString(params, level, "Encrypt/WithKey/Sk->Sk"), func(t *testing.T) {
		sk2 := kgen.GenSecretKeyNew()
		enc1 := NewEncryptor(params, sk)
		enc2 := enc1.WithKey(sk2)
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

func testApplyEvaluationKey(tc *TestContext, level int, t *testing.T) {

	params := tc.params
	sk := tc.sk
	kgen := tc.kgen
	eval := tc.eval
	enc := tc.enc
	dec := tc.dec

	var NoiseBound = float64(params.LogN())

	t.Run(testString(params, level, "Evaluator/ApplyEvaluationKey/SameDegree"), func(t *testing.T) {

		skOut := kgen.GenSecretKeyNew()

		pt := NewPlaintext(params, level)

		ct := NewCiphertext(params, 1, level)

		enc.Encrypt(pt, ct)

		// Test that Dec(KS(Enc(ct, sk), skOut), skOut) has a small norm
		evk := kgen.GenEvaluationKeyNew(sk, skOut)

		eval.ApplyEvaluationKey(ct, evk, ct)

		NewDecryptor(params, skOut).Decrypt(ct, pt)

		ringQ := params.RingQ().AtLevel(level)

		if pt.IsNTT {
			ringQ.INTT(pt.Value, pt.Value)
		}

		require.GreaterOrEqual(t, NoiseBound, ringQ.Log2OfStandardDeviation(pt.Value))
	})

	t.Run(testString(params, level, "Evaluator/ApplyEvaluationKey/LargeToSmall"), func(t *testing.T) {

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

		evk := kgenLargeDim.GenEvaluationKeyNew(skLargeDim, skSmallDim)

		ctLargeDim := NewEncryptor(paramsLargeDim, skLargeDim).EncryptZeroNew(level)
		ctSmallDim := NewCiphertext(paramsSmallDim, 1, level)

		// skLarge -> skSmall embeded in N
		eval.ApplyEvaluationKey(ctLargeDim, evk, ctSmallDim)

		// Decrypts with smaller dimension key
		ptSmallDim := NewDecryptor(paramsSmallDim, skSmallDim).DecryptNew(ctSmallDim)

		ringQSmallDim := paramsSmallDim.RingQ().AtLevel(level)
		if ptSmallDim.IsNTT {
			ringQSmallDim.INTT(ptSmallDim.Value, ptSmallDim.Value)
		}

		require.GreaterOrEqual(t, NoiseBound, ringQSmallDim.Log2OfStandardDeviation(ptSmallDim.Value))
	})

	t.Run(testString(params, level, "Evaluator/ApplyEvaluationKey/SmallToLarge"), func(t *testing.T) {

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

		evk := kgenLargeDim.GenEvaluationKeyNew(skSmallDim, skLargeDim)

		ctSmallDim := NewEncryptor(paramsSmallDim, skSmallDim).EncryptZeroNew(level)
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

func testGadgetProduct(tc *TestContext, level int, t *testing.T) {

	params := tc.params
	sk := tc.sk
	kgen := tc.kgen
	eval := tc.eval

	ringQ := params.RingQ().AtLevel(level)

	prng, _ := sampling.NewKeyedPRNG([]byte{'a', 'b', 'c'})

	sampler := ring.NewUniformSampler(prng, ringQ)

	var NoiseBound = float64(params.LogN())

	t.Run(testString(params, level, "Evaluator/GadgetProduct"), func(t *testing.T) {

		skOut := kgen.GenSecretKeyNew()

		// Generates a random polynomial
		a := sampler.ReadNew()

		// Generate the receiver
		ct := NewCiphertext(params, 1, level)

		// Generate the evaluationkey [-bs1 + s1, b]
		evk := kgen.GenEvaluationKeyNew(sk, skOut)

		// Gadget product: ct = [-cs1 + as0 , c]
		eval.GadgetProduct(level, a, &evk.GadgetCiphertext, ct)

		// pt = as0
		pt := NewDecryptor(params, skOut).DecryptNew(ct)

		ringQ := params.RingQ().AtLevel(level)

		// pt = as1 - as1 = 0 (+ some noise)
		if !pt.IsNTT {
			ringQ.NTT(pt.Value, pt.Value)
			ringQ.NTT(a, a)
		}

		ringQ.MulCoeffsMontgomeryThenSub(a, sk.Value.Q, pt.Value)
		ringQ.INTT(pt.Value, pt.Value)

		require.GreaterOrEqual(t, NoiseBound, ringQ.Log2OfStandardDeviation(pt.Value))
	})

	t.Run(testString(params, level, "Evaluator/GadgetProductHoisted"), func(t *testing.T) {

		skOut := kgen.GenSecretKeyNew()

		// Generates a random polynomial
		a := sampler.ReadNew()

		// Generate the receiver
		ct := NewCiphertext(params, 1, level)

		// Generate the evaluationkey [-bs1 + s1, b]
		evk := kgen.GenEvaluationKeyNew(sk, skOut)

		//Decompose the ciphertext
		eval.DecomposeNTT(level, params.MaxLevelP(), params.MaxLevelP()+1, a, ct.IsNTT, eval.BuffDecompQP)

		// Gadget product: ct = [-cs1 + as0 , c]
		eval.GadgetProductHoisted(level, eval.BuffDecompQP, &evk.GadgetCiphertext, ct)

		// pt = as0
		pt := NewDecryptor(params, skOut).DecryptNew(ct)

		ringQ := params.RingQ().AtLevel(level)

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

func testAutomorphism(tc *TestContext, level int, t *testing.T) {

	params := tc.params
	sk := tc.sk
	kgen := tc.kgen
	eval := tc.eval
	enc := tc.enc
	dec := tc.dec

	var NoiseBound = float64(params.LogN())

	t.Run(testString(params, level, "Evaluator/Automorphism"), func(t *testing.T) {

		// Generate a plaintext with values up to 2^30
		pt := genPlaintext(params, level, 1<<30)

		// Encrypt
		ct := enc.EncryptNew(pt)

		// Chooses a Galois Element (must be coprime with 2N)
		galEl := params.GaloisElement(-1)

		// Generate the GaloisKey
		gk := kgen.GenGaloisKeyNew(galEl, sk)

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

	t.Run(testString(params, level, "Evaluator/AutomorphismHoisted"), func(t *testing.T) {
		// Generate a plaintext with values up to 2^30
		pt := genPlaintext(params, level, 1<<30)

		// Encrypt
		ct := enc.EncryptNew(pt)

		// Chooses a Galois Element (must be coprime with 2N)
		galEl := params.GaloisElement(-1)

		// Generate the GaloisKey
		gk := kgen.GenGaloisKeyNew(galEl, sk)

		// Allocate a new EvaluationKeySet and adds the GaloisKey
		evk := NewMemEvaluationKeySet(nil, gk)

		//Decompose the ciphertext
		eval.DecomposeNTT(level, params.MaxLevelP(), params.MaxLevelP()+1, &ct.Value[1], ct.IsNTT, eval.BuffDecompQP)

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

	t.Run(testString(params, level, "Evaluator/AutomorphismHoistedLazy"), func(t *testing.T) {
		// Generate a plaintext with values up to 2^30
		pt := genPlaintext(params, level, 1<<30)

		// Encrypt
		ct := enc.EncryptNew(pt)

		// Chooses a Galois Element (must be coprime with 2N)
		galEl := params.GaloisElement(-1)

		// Generate the GaloisKey
		gk := kgen.GenGaloisKeyNew(galEl, sk)

		// Allocate a new EvaluationKeySet and adds the GaloisKey
		evk := NewMemEvaluationKeySet(nil, gk)

		//Decompose the ciphertext
		eval.DecomposeNTT(level, params.MaxLevelP(), params.MaxLevelP()+1, &ct.Value[1], ct.IsNTT, eval.BuffDecompQP)

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

func testLinearTransform(tc *TestContext, level int, t *testing.T) {

	params := tc.params
	sk := tc.sk
	kgen := tc.kgen
	eval := tc.eval
	enc := tc.enc
	dec := tc.dec

	t.Run(testString(params, level, "Evaluator/Expand"), func(t *testing.T) {

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
		var gks = kgen.GenGaloisKeysNew(params.GaloisElementsForExpand(logN), sk)
		evk := NewMemEvaluationKeySet(nil, gks...)

		eval := NewEvaluator(params, evk)

		ciphertexts := eval.WithKey(evk).Expand(ctIn, logN, logGap)

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

	t.Run(testString(params, level, "Evaluator/Pack/LogGap=LogN"), func(t *testing.T) {

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
				ringQ.AddScalar(&ciphertexts[i].Value[0], scalar, &ciphertexts[i].Value[0])
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
		gks := kgen.GenGaloisKeysNew(params.GaloisElementsForPack(params.LogN()), sk)
		evk := NewMemEvaluationKeySet(nil, gks...)

		ct := eval.WithKey(evk).Pack(ciphertexts, params.LogN(), false)

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

	t.Run(testString(params, level, "Evaluator/Pack/LogGap=LogN-1"), func(t *testing.T) {

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
				ringQ.INTT(&ciphertexts[i].Value[0], &ciphertexts[i].Value[0])
			}

			for j := 0; j < level+1; j++ {
				ciphertexts[i].Value[0].Coeffs[j][0] = ring.CRed(ciphertexts[i].Value[0].Coeffs[j][0]+scalar, ringQ.SubRings[j].Modulus)
				ciphertexts[i].Value[0].Coeffs[j][N/2] = ring.CRed(ciphertexts[i].Value[0].Coeffs[j][N/2]+scalar, ringQ.SubRings[j].Modulus)
			}

			if ciphertexts[i].IsNTT {
				ringQ.NTT(&ciphertexts[i].Value[0], &ciphertexts[i].Value[0])
			}

			slotIndex[i] = true

			for j := 0; j < level+1; j++ {
				ptPacked.Value.Coeffs[j][i] = scalar
				ptPacked.Value.Coeffs[j][i+N/2] = scalar
			}
		}

		// Galois Keys
		gks := kgen.GenGaloisKeysNew(params.GaloisElementsForPack(params.LogN()-1), sk)
		evk := NewMemEvaluationKeySet(nil, gks...)

		ct := eval.WithKey(evk).Pack(ciphertexts, params.LogN()-1, true)

		dec.Decrypt(ct, pt)

		if pt.IsNTT {
			ringQ.INTT(pt.Value, pt.Value)
		}

		ringQ.Sub(pt.Value, ptPacked.Value, pt.Value)

		NoiseBound := 15.0

		// Logs the noise
		require.GreaterOrEqual(t, NoiseBound, ringQ.Log2OfStandardDeviation(pt.Value))
	})

	t.Run(testString(params, level, "Evaluator/InnerSum"), func(t *testing.T) {

		batch := 5
		n := 7

		ringQ := tc.params.RingQ().AtLevel(level)

		pt := genPlaintext(params, level, 1<<30)
		ptInnerSum := pt.Value.CopyNew()
		ct := enc.EncryptNew(pt)

		// Galois Keys
		gks := kgen.GenGaloisKeysNew(params.GaloisElementsForInnerSum(batch, n), sk)
		evk := NewMemEvaluationKeySet(nil, gks...)

		eval.WithKey(evk).InnerSum(ct, batch, n, ct)

		dec.Decrypt(ct, pt)

		if pt.IsNTT {
			ringQ.INTT(pt.Value, pt.Value)
			ringQ.INTT(ptInnerSum, ptInnerSum)
		}

		polyTmp := ringQ.NewPoly()

		// Applies the same circuit (naively) on the plaintext
		polyInnerSum := ptInnerSum.CopyNew()
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

func testWriteAndRead(tc *TestContext, t *testing.T) {

	params := tc.params

	sk, pk := tc.sk, tc.pk

	t.Run(testString(params, params.MaxLevel(), "WriteAndRead/OperandQ"), func(t *testing.T) {
		prng, _ := sampling.NewPRNG()
		plaintextWant := NewPlaintext(params, params.MaxLevel())
		ring.NewUniformSampler(prng, params.RingQ()).Read(plaintextWant.Value)
		buffer.RequireSerializerCorrect(t, &plaintextWant.OperandQ)
	})

	t.Run(testString(params, params.MaxLevel(), "WriteAndRead/Plaintext"), func(t *testing.T) {
		prng, _ := sampling.NewPRNG()
		plaintextWant := NewPlaintext(params, params.MaxLevel())
		ring.NewUniformSampler(prng, params.RingQ()).Read(plaintextWant.Value)
		buffer.RequireSerializerCorrect(t, plaintextWant)
	})

	t.Run(testString(params, params.MaxLevel(), "WriteAndRead/Ciphertext"), func(t *testing.T) {

		prng, _ := sampling.NewPRNG()

		for degree := 0; degree < 4; degree++ {
			t.Run(fmt.Sprintf("degree=%d", degree), func(t *testing.T) {
				buffer.RequireSerializerCorrect(t, NewCiphertextRandom(prng, params, degree, params.MaxLevel()))
			})
		}
	})

	t.Run(testString(params, params.MaxLevel(), "WriteAndRead/CiphertextQP"), func(t *testing.T) {
		buffer.RequireSerializerCorrect(t, &OperandQP{Value: tc.pk.Value[:]})
	})

	t.Run(testString(params, params.MaxLevel(), "WriteAndRead/GadgetCiphertext"), func(t *testing.T) {
		buffer.RequireSerializerCorrect(t, &tc.kgen.GenRelinearizationKeyNew(tc.sk).GadgetCiphertext)
	})

	t.Run(testString(params, params.MaxLevel(), "WriteAndRead/Sk"), func(t *testing.T) {
		buffer.RequireSerializerCorrect(t, sk)
	})

	t.Run(testString(params, params.MaxLevel(), "WriteAndRead/Pk"), func(t *testing.T) {
		buffer.RequireSerializerCorrect(t, pk)
	})

	t.Run(testString(params, params.MaxLevel(), "WriteAndRead/EvaluationKey"), func(t *testing.T) {
		buffer.RequireSerializerCorrect(t, tc.kgen.GenEvaluationKeyNew(sk, sk))
	})

	t.Run(testString(params, params.MaxLevel(), "WriteAndRead/RelinearizationKey"), func(t *testing.T) {
		buffer.RequireSerializerCorrect(t, tc.kgen.GenRelinearizationKeyNew(tc.sk))
	})

	t.Run(testString(params, params.MaxLevel(), "WriteAndRead/GaloisKey"), func(t *testing.T) {
		buffer.RequireSerializerCorrect(t, tc.kgen.GenGaloisKeyNew(5, tc.sk))
	})

	t.Run(testString(params, params.MaxLevel(), "WriteAndRead/EvaluationKeySet"), func(t *testing.T) {
		buffer.RequireSerializerCorrect(t, &MemEvaluationKeySet{
			Rlk: tc.kgen.GenRelinearizationKeyNew(tc.sk),
			Gks: map[uint64]*GaloisKey{5: tc.kgen.GenGaloisKeyNew(5, tc.sk)},
		})
	})

	t.Run(testString(params, params.MaxLevel(), "WriteAndRead/PowerBasis"), func(t *testing.T) {

		prng, _ := sampling.NewPRNG()

		ct := NewCiphertextRandom(prng, params, 1, params.MaxLevel())

		basis := NewPowerBasis(ct, polynomial.Chebyshev, nil)

		basis.Value[2] = NewCiphertextRandom(prng, params, 1, params.MaxLevel())
		basis.Value[3] = NewCiphertextRandom(prng, params, 2, params.MaxLevel())
		basis.Value[4] = NewCiphertextRandom(prng, params, 1, params.MaxLevel())
		basis.Value[8] = NewCiphertextRandom(prng, params, 1, params.MaxLevel())

		buffer.RequireSerializerCorrect(t, basis)
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
