package float

import (
	"math"
	"runtime"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

func TestHomomorphicMod(t *testing.T) {
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

	testEvalModMarshalling(t)

	var params ckks.Parameters
	if params, err = ckks.NewParametersFromLiteral(ParametersLiteral); err != nil {
		t.Fatal(err)
	}

	for _, testSet := range []func(params ckks.Parameters, t *testing.T){
		testEvalMod,
	} {
		testSet(params, t)
		runtime.GC()
	}
}

func testEvalModMarshalling(t *testing.T) {
	t.Run("Marshalling", func(t *testing.T) {

		evm := EvalModLiteral{
			LevelStart:      12,
			SineType:        SinContinuous,
			LogMessageRatio: 8,
			K:               14,
			SineDegree:      127,
			ArcSineDegree:   7,
			LogScale:        60,
		}

		data, err := evm.MarshalBinary()
		assert.Nil(t, err)

		evmNew := new(EvalModLiteral)
		if err := evmNew.UnmarshalBinary(data); err != nil {
			assert.Nil(t, err)
		}
		assert.Equal(t, evm, *evmNew)
	})
}

func testEvalMod(params ckks.Parameters, t *testing.T) {

	kgen := ckks.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	ecd := ckks.NewEncoder(params)
	enc := ckks.NewEncryptor(params, sk)
	dec := ckks.NewDecryptor(params, sk)
	eval := ckks.NewEvaluator(params, rlwe.NewMemEvaluationKeySet(kgen.GenRelinearizationKeyNew(sk)))

	modEval := NewHModEvaluator(eval)

	t.Run("SineContinuousWithArcSine", func(t *testing.T) {

		evm := EvalModLiteral{
			LevelStart:      12,
			SineType:        SinContinuous,
			LogMessageRatio: 8,
			K:               14,
			SineDegree:      127,
			ArcSineDegree:   7,
			LogScale:        60,
		}

		values, ciphertext := evaluatexmod1(evm, params, ecd, enc, modEval, t)

		ckks.VerifyTestVectors(params, ecd, dec, values, ciphertext, params.LogDefaultScale(), nil, *printPrecisionStats, t)
	})

	t.Run("CosDiscrete", func(t *testing.T) {

		evm := EvalModLiteral{
			LevelStart:      12,
			SineType:        CosDiscrete,
			LogMessageRatio: 8,
			K:               12,
			SineDegree:      30,
			DoubleAngle:     3,
			LogScale:        60,
		}

		values, ciphertext := evaluatexmod1(evm, params, ecd, enc, modEval, t)

		ckks.VerifyTestVectors(params, ecd, dec, values, ciphertext, params.LogDefaultScale(), nil, *printPrecisionStats, t)
	})

	t.Run("CosContinuous", func(t *testing.T) {

		evm := EvalModLiteral{
			LevelStart:      12,
			SineType:        CosContinuous,
			LogMessageRatio: 4,
			K:               325,
			SineDegree:      177,
			DoubleAngle:     4,
			LogScale:        60,
		}

		values, ciphertext := evaluatexmod1(evm, params, ecd, enc, modEval, t)

		ckks.VerifyTestVectors(params, ecd, dec, values, ciphertext, params.LogDefaultScale(), nil, *printPrecisionStats, t)
	})
}

func evaluatexmod1(evm EvalModLiteral, params ckks.Parameters, ecd *ckks.Encoder, enc *rlwe.Encryptor, eval *HModEvaluator, t *testing.T) ([]float64, *rlwe.Ciphertext) {

	EvalModPoly, err := NewEvalModPolyFromLiteral(params, evm)
	require.NoError(t, err)

	values, _, ciphertext := newTestVectorsEvalMod(params, enc, ecd, EvalModPoly, t)

	// Scale the message to Delta = Q/MessageRatio
	scale := rlwe.NewScale(math.Exp2(math.Round(math.Log2(float64(params.Q()[0]) / EvalModPoly.MessageRatio()))))
	scale = scale.Div(ciphertext.Scale)
	eval.ScaleUp(ciphertext, rlwe.NewScale(math.Round(scale.Float64())), ciphertext)

	// Scale the message up to Sine/MessageRatio
	scale = EvalModPoly.ScalingFactor().Div(ciphertext.Scale)
	scale = scale.Div(rlwe.NewScale(EvalModPoly.MessageRatio()))
	eval.ScaleUp(ciphertext, rlwe.NewScale(math.Round(scale.Float64())), ciphertext)

	// Normalization
	require.NoError(t, eval.Mul(ciphertext, 1/(float64(EvalModPoly.K())*EvalModPoly.QDiff()), ciphertext))
	require.NoError(t, eval.Rescale(ciphertext, ciphertext))

	// EvalMod
	ciphertext, err = eval.EvalModNew(ciphertext, EvalModPoly)
	require.NoError(t, err)

	// PlaintextCircuit
	for i := range values {
		x := values[i]

		x /= EvalModPoly.MessageRatio()
		x /= EvalModPoly.QDiff()
		x = math.Sin(6.28318530717958 * x)

		if evm.ArcSineDegree > 0 {
			x = math.Asin(x)
		}

		x *= EvalModPoly.MessageRatio()
		x *= EvalModPoly.QDiff()
		x /= 6.28318530717958

		values[i] = x
	}

	return values, ciphertext
}

func newTestVectorsEvalMod(params ckks.Parameters, encryptor *rlwe.Encryptor, encoder *ckks.Encoder, evm EvalModPoly, t *testing.T) (values []float64, plaintext *rlwe.Plaintext, ciphertext *rlwe.Ciphertext) {

	logSlots := params.LogMaxDimensions().Cols

	values = make([]float64, 1<<logSlots)

	K := float64(evm.K() - 1)
	Q := float64(params.Q()[0]) / math.Exp2(math.Round(math.Log2(float64(params.Q()[0])))) * evm.MessageRatio()

	for i := uint64(0); i < 1<<logSlots; i++ {
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
