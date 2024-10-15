package mod1

import (
	"math"
	"runtime"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

func TestMod1(t *testing.T) {
	var err error

	if runtime.GOARCH == "wasm" {
		t.Skip("skipping homomorphic mod tests for GOARCH=wasm")
	}

	ParametersLiteral := ckks.ParametersLiteral{
		LogN:            10,
		LogQ:            []int{55, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 53},
		LogP:            []int{61, 61, 61, 61, 61},
		Xs:              ring.Ternary{H: 192},
		LogDefaultScale: 45,
	}

	testMod1Marshalling(t)

	var params ckks.Parameters
	if params, err = ckks.NewParametersFromLiteral(ParametersLiteral); err != nil {
		t.Fatal(err)
	}

	for _, testSet := range []func(params ckks.Parameters, t *testing.T){
		testMod1,
	} {
		testSet(params, t)
		runtime.GC()
	}
}

func testMod1Marshalling(t *testing.T) {
	t.Run("Marshalling", func(t *testing.T) {

		evm := ParametersLiteral{
			LevelQ:          12,
			Mod1Type:        SinContinuous,
			LogMessageRatio: 8,
			K:               14,
			Mod1Degree:      127,
			Mod1InvDegree:   7,
			LogScale:        60,
		}

		data, err := evm.MarshalBinary()
		assert.Nil(t, err)

		evmNew := new(ParametersLiteral)
		if err := evmNew.UnmarshalBinary(data); err != nil {
			assert.Nil(t, err)
		}
		assert.Equal(t, evm, *evmNew)
	})
}

func testMod1(params ckks.Parameters, t *testing.T) {

	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	ecd := ckks.NewEncoder(params)
	enc := rlwe.NewEncryptor(params, sk)
	dec := rlwe.NewDecryptor(params, sk)
	eval := ckks.NewEvaluator(params, rlwe.NewMemEvaluationKeySet(kgen.GenRelinearizationKeyNew(sk)))

	t.Run("SineContinuousWithArcSine", func(t *testing.T) {

		evm := ParametersLiteral{
			LevelQ:          12,
			Mod1Type:        SinContinuous,
			LogMessageRatio: 8,
			K:               14,
			Mod1Degree:      127,
			Mod1InvDegree:   7,
			LogScale:        60,
		}

		values, ciphertext := evaluateMod1(evm, params, ecd, enc, dec, eval, t)

		ckks.VerifyTestVectors(params, ecd, dec, values, ciphertext, params.LogDefaultScale(), 0, true, t)
	})

	t.Run("CosDiscrete", func(t *testing.T) {

		evm := ParametersLiteral{
			LevelQ:          12,
			Mod1Type:        CosDiscrete,
			LogMessageRatio: 8,
			K:               12,
			Mod1Degree:      30,
			DoubleAngle:     3,
			LogScale:        60,
		}

		values, ciphertext := evaluateMod1(evm, params, ecd, enc, dec, eval, t)

		ckks.VerifyTestVectors(params, ecd, dec, values, ciphertext, params.LogDefaultScale(), 0, true, t)
	})

	t.Run("CosContinuous", func(t *testing.T) {

		evm := ParametersLiteral{
			LevelQ:          12,
			Mod1Type:        CosContinuous,
			LogMessageRatio: 4,
			K:               325,
			Mod1Degree:      177,
			DoubleAngle:     4,
			LogScale:        60,
		}

		values, ciphertext := evaluateMod1(evm, params, ecd, enc, dec, eval, t)

		ckks.VerifyTestVectors(params, ecd, dec, values, ciphertext, params.LogDefaultScale(), 0, true, t)
	})
}

func evaluateMod1(evm ParametersLiteral, params ckks.Parameters, ecd *ckks.Encoder, enc *rlwe.Encryptor, dec *rlwe.Decryptor, eval *ckks.Evaluator, t *testing.T) ([]float64, *rlwe.Ciphertext) {

	mod1Parameters, err := NewParametersFromLiteral(params, evm)
	require.NoError(t, err)

	values, _, ciphertext := newTestVectorsMod1(params, enc, ecd, mod1Parameters, t)

	// Scale the message to Delta = Q/MessageRatio
	scale := rlwe.NewScale(math.Exp2(math.Round(math.Log2(float64(params.Q()[0]) / mod1Parameters.MessageRatio()))))
	scale = scale.Div(ciphertext.Scale)
	eval.ScaleUp(ciphertext, rlwe.NewScale(math.Round(scale.Float64())), ciphertext)

	// Scale the message up to Sine/MessageRatio
	scale = mod1Parameters.ScalingFactor().Div(ciphertext.Scale)
	scale = scale.Div(rlwe.NewScale(mod1Parameters.MessageRatio()))
	eval.ScaleUp(ciphertext, rlwe.NewScale(math.Round(scale.Float64())), ciphertext)

	// Normalization
	require.NoError(t, eval.Mul(ciphertext, 1/(mod1Parameters.K*mod1Parameters.QDiff), ciphertext))
	require.NoError(t, eval.Rescale(ciphertext, ciphertext))

	// EvalMod
	ciphertext, err = NewEvaluator(eval, polynomial.NewEvaluator(params, eval), mod1Parameters).EvaluateNew(ciphertext)
	require.NoError(t, err)

	// PlaintextCircuit
	for i := range values {
		x := values[i]

		x /= mod1Parameters.MessageRatio()
		x /= mod1Parameters.QDiff
		x = math.Sin(6.28318530717958 * x)

		if evm.Mod1InvDegree > 0 {
			x = math.Asin(x)
		}

		x *= mod1Parameters.MessageRatio()
		x *= mod1Parameters.QDiff
		x /= 6.28318530717958

		values[i] = x
	}

	return values, ciphertext
}

func newTestVectorsMod1(params ckks.Parameters, encryptor *rlwe.Encryptor, encoder *ckks.Encoder, evm Parameters, t *testing.T) (values []float64, plaintext *rlwe.Plaintext, ciphertext *rlwe.Ciphertext) {

	values = make([]float64, params.MaxSlots())

	K := evm.K - 1
	Q := evm.QDiff * evm.MessageRatio()

	for i := range values {
		values[i] = math.Round(sampling.RandFloat64(-K, K))*Q + sampling.RandFloat64(-1, 1)
	}

	values[0] = K*Q + 0.5

	plaintext = ckks.NewPlaintext(params, params.MaxLevel())

	encoder.Encode(values, plaintext)

	if encryptor != nil {
		var err error
		ciphertext, err = encryptor.EncryptNew(plaintext)
		require.NoError(t, err)
	}

	return values, plaintext, ciphertext
}
