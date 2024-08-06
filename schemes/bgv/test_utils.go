package bgv

import (
	"fmt"
	"math"
	"slices"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

type TestContext struct {
	Params Parameters
	Ecd    *Encoder

	Prng    sampling.PRNG
	Sampler *ring.UniformSampler

	Kgen *rlwe.KeyGenerator
	Sk   *rlwe.SecretKey
	Pk   *rlwe.PublicKey

	Enc *rlwe.Encryptor
	Dec *rlwe.Decryptor

	Evl *Evaluator
}

func NewTestContext(params ParametersLiteral, scaleInvariant bool) *TestContext {
	tc := new(TestContext)

	var err error

	tc.Params, err = NewParametersFromLiteral(params)
	if err != nil {
		panic(err)
	}
	tc.Ecd = NewEncoder(tc.Params)

	tc.Prng, err = sampling.NewPRNG()
	if err != nil {
		panic(err)
	}
	tc.Sampler = ring.NewUniformSampler(tc.Prng, tc.Params.RingT())

	tc.Kgen = rlwe.NewKeyGenerator(tc.Params)
	tc.Sk, tc.Pk = tc.Kgen.GenKeyPairNew()

	tc.Enc = rlwe.NewEncryptor(tc.Params, tc.Pk)
	tc.Dec = rlwe.NewDecryptor(tc.Params, tc.Sk)

	tc.Evl = NewEvaluator(tc.Params, rlwe.NewMemEvaluationKeySet(tc.Kgen.GenRelinearizationKeyNew(tc.Sk)), scaleInvariant)

	return tc
}

func (tc TestContext) String() string {
	return fmt.Sprintf("LogN=%d/logQ=%d/logP=%d/LogSlots=%dx%d/logT=%d/Qi=%d/Pi=%d",
		tc.Params.LogN(),
		int(math.Round(tc.Params.LogQ())),
		int(math.Round(tc.Params.LogP())),
		tc.Params.LogMaxDimensions().Rows,
		tc.Params.LogMaxDimensions().Cols,
		int(math.Round(tc.Params.LogT())),
		tc.Params.QCount(),
		tc.Params.PCount())
}

func VerifyTestVectors(params Parameters, encoder *Encoder, decryptor *rlwe.Decryptor, have interface{}, want []uint64, t *testing.T) {
	values := make([]uint64, params.MaxSlots())

	switch have := have.(type) {
	case *rlwe.Plaintext:
		require.NoError(t, encoder.Decode(have, values))
	case *rlwe.Ciphertext:
		require.NoError(t, encoder.Decode(decryptor.DecryptNew(have), values))
	default:
		t.Error("invalid unsupported test object type")
	}

	require.True(t, slices.Equal(values, want))
}

func NewTestVector(params Parameters, encoder *Encoder, encryptor *rlwe.Encryptor, level int, scale rlwe.Scale) (values []uint64, pt *rlwe.Plaintext, ct *rlwe.Ciphertext) {
	values = make([]uint64, params.MaxSlots())
	for i := range values {
		values[i] = sampling.RandUint64() % params.PlaintextModulus()
	}

	pt = NewPlaintext(params, level)
	pt.Scale = scale
	if err := encoder.Encode(values, pt); err != nil {
		panic(err)
	}
	if encryptor != nil {
		var err error
		ct, err = encryptor.EncryptNew(pt)
		if err != nil {
			panic(err)
		}
	}
	return
}
