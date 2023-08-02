package circuits

import (
	"encoding/json"
	"fmt"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
	"github.com/tuneinsight/lattigo/v4/utils/buffer"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

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

	defaultParamsLiteral := circuitsTestParamsLiteral

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

				testSerialization(tc, tc.params.MaxLevel(), paramsLit.BaseTwoDecomposition, t)
			}
		}
	}
}

func testSerialization(tc *TestContext, level, bpw2 int, t *testing.T) {

	params := tc.params

	levelQ := level
	levelP := params.MaxLevelP()

	t.Run(testString(params, levelQ, levelP, bpw2, "WriteAndRead/PowerBasis"), func(t *testing.T) {

		prng, _ := sampling.NewPRNG()

		ct := rlwe.NewCiphertextRandom(prng, params, 1, levelQ)

		basis := NewPowerBasis(ct, bignum.Chebyshev)

		basis.Value[2] = rlwe.NewCiphertextRandom(prng, params, 1, levelQ)
		basis.Value[3] = rlwe.NewCiphertextRandom(prng, params, 2, levelQ)
		basis.Value[4] = rlwe.NewCiphertextRandom(prng, params, 1, levelQ)
		basis.Value[8] = rlwe.NewCiphertextRandom(prng, params, 1, levelQ)

		buffer.RequireSerializerCorrect(t, &basis)
	})
}

type TestContext struct {
	params rlwe.Parameters
	kgen   *rlwe.KeyGenerator
	enc    *rlwe.Encryptor
	dec    *rlwe.Decryptor
	sk     *rlwe.SecretKey
	pk     *rlwe.PublicKey
}

func NewTestContext(params rlwe.Parameters) (tc *TestContext, err error) {
	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()

	pk, err := kgen.GenPublicKeyNew(sk)
	if err != nil {
		return nil, err
	}

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
	}, nil
}

type TestParametersLiteral struct {
	BaseTwoDecomposition int
	rlwe.ParametersLiteral
}

var (
	logN = 10
	qi   = []uint64{0x200000440001, 0x7fff80001, 0x800280001, 0x7ffd80001, 0x7ffc80001}
	pj   = []uint64{0x3ffffffb80001, 0x4000000800001}

	circuitsTestParamsLiteral = []TestParametersLiteral{
		// RNS decomposition, no Pw2 decomposition
		{
			BaseTwoDecomposition: 0,

			ParametersLiteral: rlwe.ParametersLiteral{
				LogN:    logN,
				Q:       qi,
				P:       pj,
				NTTFlag: true,
			},
		},
		// RNS decomposition, Pw2 decomposition
		{
			BaseTwoDecomposition: 16,

			ParametersLiteral: rlwe.ParametersLiteral{
				LogN:    logN,
				Q:       qi,
				P:       pj[:1],
				NTTFlag: true,
			},
		},
		// No RNS decomposition, Pw2 decomposition
		{
			BaseTwoDecomposition: 1,

			ParametersLiteral: rlwe.ParametersLiteral{
				LogN:    logN,
				Q:       qi,
				P:       nil,
				NTTFlag: true,
			},
		},
	}
)
