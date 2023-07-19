package he

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"runtime"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")

func testString(params rlwe.Parameters, levelQ, levelP, bpw2 int, opname string) string {
	return fmt.Sprintf("%s/logN=%d/Qi=%d/Pi=%d/Pw2=%d/NTT=%t/RingType=%s",
		opname,
		params.LogN(),
		levelQ+1,
		levelP+1,
		bpw2,
		params.NTTFlag(),
		params.RingType())
}

func TestHE(t *testing.T) {
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

				var params rlwe.Parameters
				if params, err = rlwe.NewParametersFromLiteral(paramsLit.ParametersLiteral); err != nil {
					t.Fatal(err)
				}

				tc, err := NewTestContext(params)
				require.NoError(t, err)

				for _, level := range []int{0, params.MaxLevel()}[:] {

					for _, testSet := range []func(tc *TestContext, level, bpw2 int, t *testing.T){
						testLinearTransformation,
					} {
						testSet(tc, level, paramsLit.BaseTwoDecomposition, t)
						runtime.GC()
					}
				}
			}
		}
	}
}

type TestContext struct {
	params rlwe.Parameters
	kgen   *rlwe.KeyGenerator
	enc    *rlwe.Encryptor
	dec    *rlwe.Decryptor
	sk     *rlwe.SecretKey
	pk     *rlwe.PublicKey
	eval   *Evaluator
}

func NewTestContext(params rlwe.Parameters) (tc *TestContext, err error) {
	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()

	pk, err := kgen.GenPublicKeyNew(sk)
	if err != nil {
		return nil, err
	}

	eval := NewEvaluator(params, nil)

	enc, err := rlwe.NewEncryptor(params, sk)
	if err != nil {
		return nil, err
	}

	dec, err := rlwe.NewDecryptor(params, sk)
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

func testLinearTransformation(tc *TestContext, level, bpw2 int, t *testing.T) {

	params := tc.params
	sk := tc.sk
	kgen := tc.kgen
	eval := tc.eval
	enc := tc.enc
	dec := tc.dec

	evkParams := rlwe.EvaluationKeyParameters{LevelQ: level, LevelP: params.MaxLevelP(), BaseTwoDecomposition: bpw2}

	t.Run(testString(params, level, params.MaxLevelP(), bpw2, "Evaluator/Expand"), func(t *testing.T) {

		if params.RingType() != ring.Standard {
			t.Skip("Expand not supported for ring.Type = ring.ConjugateInvariant")
		}

		pt := rlwe.NewPlaintext(params, level)
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

		ctIn := rlwe.NewCiphertext(params, 1, level)
		enc.Encrypt(pt, ctIn)

		// GaloisKeys
		var gks, err = kgen.GenGaloisKeysNew(params.GaloisElementsForExpand(logN), sk, evkParams)
		require.NoError(t, err)

		evk := rlwe.NewMemEvaluationKeySet(nil, gks...)

		eval := NewEvaluator(params, evk)

		ciphertexts, err := eval.WithKey(evk).Expand(ctIn, logN, logGap)
		require.NoError(t, err)

		Q := ringQ.ModuliChain()

		NoiseBound := float64(params.LogN() - logN + bpw2)

		if bpw2 != 0 {
			NoiseBound += float64(level + 5)
		}

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

		pt := rlwe.NewPlaintext(params, level)
		N := params.N()
		ringQ := tc.params.RingQ().AtLevel(level)
		gap := params.N() / 16

		ptPacked := rlwe.NewPlaintext(params, level)
		ciphertexts := make(map[int]*rlwe.Ciphertext)
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

		gks, err := kgen.GenGaloisKeysNew(galEls, sk, evkParams)
		require.NoError(t, err)

		evk := rlwe.NewMemEvaluationKeySet(nil, gks...)

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

		NoiseBound := 15.0 + float64(bpw2)

		if bpw2 != 0 {
			NoiseBound += math.Log2(float64(level)+1.0) + 1.0
		}

		// Logs the noise
		require.GreaterOrEqual(t, NoiseBound, ringQ.Log2OfStandardDeviation(pt.Value))
	})

	t.Run(testString(params, level, params.MaxLevelP(), bpw2, "Evaluator/Pack/LogGap=LogN-1"), func(t *testing.T) {

		if params.RingType() != ring.Standard {
			t.Skip("Pack not supported for ring.Type = ring.ConjugateInvariant")
		}

		pt := rlwe.NewPlaintext(params, level)
		N := params.N()
		ringQ := tc.params.RingQ().AtLevel(level)

		ptPacked := rlwe.NewPlaintext(params, level)
		ciphertexts := make(map[int]*rlwe.Ciphertext)
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

		gks, err := kgen.GenGaloisKeysNew(galEls, sk, evkParams)
		require.NoError(t, err)

		evk := rlwe.NewMemEvaluationKeySet(nil, gks...)

		ct, err := eval.WithKey(evk).Pack(ciphertexts, params.LogN()-1, true)
		require.NoError(t, err)

		dec.Decrypt(ct, pt)

		if pt.IsNTT {
			ringQ.INTT(pt.Value, pt.Value)
		}

		ringQ.Sub(pt.Value, ptPacked.Value, pt.Value)

		NoiseBound := 15.0 + float64(bpw2)

		if bpw2 != 0 {
			NoiseBound += math.Log2(float64(level)+1.0) + 1.0
		}

		// Logs the noise
		require.GreaterOrEqual(t, NoiseBound, ringQ.Log2OfStandardDeviation(pt.Value))
	})

	t.Run(testString(params, level, params.MaxLevelP(), bpw2, "Evaluator/InnerSum"), func(t *testing.T) {

		if params.MaxLevelP() == -1 {
			t.Skip("test requires #P > 0")
		}

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

		evk := rlwe.NewMemEvaluationKeySet(nil, gks...)

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

func genPlaintext(params rlwe.Parameters, level, max int) (pt *rlwe.Plaintext) {

	N := params.N()

	step := float64(max) / float64(N)

	pt = rlwe.NewPlaintext(params, level)

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
