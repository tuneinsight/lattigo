package advanced

import (
	"math"
	"runtime"
	"testing"

	"github.com/stretchr/testify/assert"
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
		LogN: 14,
		Q: []uint64{
			0x80000000080001,   // 55 Q0
			0xffffffffffc0001,  // 60
			0x10000000006e0001, // 60
			0xfffffffff840001,  // 60
			0x1000000000860001, // 60
			0xfffffffff6a0001,  // 60
			0x1000000000980001, // 60
			0xfffffffff5a0001,  // 60
			0x1000000000b00001, // 60
			0x1000000000ce0001, // 60
			0xfffffffff2a0001,  // 60
			0xfffffffff240001,  // 60
			0x1000000000f00001, // 60
			0x200000000e0001,   // 53
		},
		P: []uint64{
			0x1fffffffffe00001, // Pi 61
			0x1fffffffffc80001, // Pi 61
			0x1fffffffffb40001, // Pi 61
			0x1fffffffff500001, // Pi 61
			0x1fffffffff420001, // Pi 61
		},
		H:        192,
		Sigma:    rlwe.DefaultSigma,
		LogSlots: 13,
		LogScale: 45,
	}

	testEvalModMarshalling(t)

	var params ckks.Parameters
	if params, err = ckks.NewParametersFromLiteral(ParametersLiteral); err != nil {
		panic(err)
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
	encoder := ckks.NewEncoder(params)
	encryptor := ckks.NewEncryptor(params, sk)
	decryptor := ckks.NewDecryptor(params, sk)

	evk := rlwe.NewEvaluationKeySet()
	evk.RelinearizationKey = kgen.GenRelinearizationKeyNew(sk)

	eval := NewEvaluator(params, evk)

	t.Run("SineChebyshevWithArcSine", func(t *testing.T) {

		evm := EvalModLiteral{
			LevelStart:      12,
			SineType:        SinContinuous,
			LogMessageRatio: 8,
			K:               14,
			SineDegree:      127,
			ArcSineDegree:   7,
			LogScale:        60,
		}

		EvalModPoly := NewEvalModPolyFromLiteral(params, evm)

		values, _, ciphertext := newTestVectorsEvalMod(params, encryptor, encoder, EvalModPoly, t)

		// Scale the message to Delta = Q/MessageRatio
		scale := rlwe.NewScale(math.Exp2(math.Round(math.Log2(float64(params.Q()[0]) / EvalModPoly.MessageRatio()))))
		scale = scale.Div(ciphertext.Scale)
		eval.ScaleUp(ciphertext, rlwe.NewScale(math.Round(scale.Float64())), ciphertext)

		// Scale the message up to Sine/MessageRatio
		scale = EvalModPoly.ScalingFactor().Div(ciphertext.Scale)
		scale = scale.Div(rlwe.NewScale(EvalModPoly.MessageRatio()))
		eval.ScaleUp(ciphertext, rlwe.NewScale(math.Round(scale.Float64())), ciphertext)

		// Normalization
		eval.MultByConst(ciphertext, 1/(float64(EvalModPoly.K())*EvalModPoly.QDiff()), ciphertext)
		if err := eval.Rescale(ciphertext, params.DefaultScale(), ciphertext); err != nil {
			t.Error(err)
		}

		// EvalMod
		ciphertext = eval.EvalModNew(ciphertext, EvalModPoly)

		// PlaintextCircuit
		//pi2r := 6.283185307179586/complex(math.Exp2(float64(evm.DoubleAngle)), 0)
		for i := range values {
			values[i] -= complex(EvalModPoly.MessageRatio()*EvalModPoly.QDiff()*math.Round(real(values[i])/(EvalModPoly.MessageRatio()/EvalModPoly.QDiff())), 0)
			//values[i] = sin2pi2pi(values[i] / complex(evm.MessageRatio*evm.QDiff(), 0)) * complex(evm.MessageRatio*evm.QDiff(), 0) / 6.283185307179586
		}

		verifyTestVectors(params, encoder, decryptor, values, ciphertext, params.LogSlots(), t)
	})

	t.Run("CosOptimizedChebyshevWithArcSine", func(t *testing.T) {

		evm := EvalModLiteral{
			LevelStart:      12,
			SineType:        CosContinuous,
			LogMessageRatio: 8,
			K:               325,
			SineDegree:      255,
			DoubleAngle:     4,
			LogScale:        60,
		}

		EvalModPoly := NewEvalModPolyFromLiteral(params, evm)

		values, _, ciphertext := newTestVectorsEvalMod(params, encryptor, encoder, EvalModPoly, t)

		// Scale the message to Delta = Q/MessageRatio
		scale := rlwe.NewScale(math.Exp2(math.Round(math.Log2(float64(params.Q()[0]) / EvalModPoly.MessageRatio()))))
		scale = scale.Div(ciphertext.Scale)
		eval.ScaleUp(ciphertext, rlwe.NewScale(math.Round(scale.Float64())), ciphertext)

		// Scale the message up to Sine/MessageRatio
		scale = EvalModPoly.ScalingFactor().Div(ciphertext.Scale)
		scale = scale.Div(rlwe.NewScale(EvalModPoly.MessageRatio()))
		eval.ScaleUp(ciphertext, rlwe.NewScale(math.Round(scale.Float64())), ciphertext)

		// Normalization
		eval.MultByConst(ciphertext, 1/(float64(EvalModPoly.K())*EvalModPoly.QDiff()), ciphertext)
		if err := eval.Rescale(ciphertext, params.DefaultScale(), ciphertext); err != nil {
			t.Error(err)
		}

		// EvalMod
		ciphertext = eval.EvalModNew(ciphertext, EvalModPoly)

		// PlaintextCircuit
		//pi2r := 6.283185307179586/complex(math.Exp2(float64(evm.DoubleAngle)), 0)
		for i := range values {
			values[i] -= complex(EvalModPoly.MessageRatio()*EvalModPoly.QDiff()*math.Round(real(values[i])/(EvalModPoly.MessageRatio()/EvalModPoly.QDiff())), 0)
			//values[i] = sin2pi2pi(values[i] / complex(evm.MessageRatio*evm.QDiff(), 0)) * complex(evm.MessageRatio*evm.QDiff(), 0) / 6.283185307179586
		}

		verifyTestVectors(params, encoder, decryptor, values, ciphertext, params.LogSlots(), t)
	})
}

func newTestVectorsEvalMod(params ckks.Parameters, encryptor rlwe.Encryptor, encoder ckks.Encoder, evm EvalModPoly, t *testing.T) (values []complex128, plaintext *rlwe.Plaintext, ciphertext *rlwe.Ciphertext) {

	logSlots := params.LogSlots()

	values = make([]complex128, 1<<logSlots)

	K := float64(evm.K() - 1)
	Q := float64(params.Q()[0]) / math.Exp2(math.Round(math.Log2(float64(params.Q()[0])))) * evm.MessageRatio()

	for i := uint64(0); i < 1<<logSlots; i++ {
		values[i] = complex(math.Round(sampling.RandFloat64(-K, K))*Q+sampling.RandFloat64(-1, 1), 0)
	}

	values[0] = complex(K*Q+0.5, 0)

	plaintext = ckks.NewPlaintext(params, params.MaxLevel())

	encoder.Encode(values, plaintext, logSlots)

	if encryptor != nil {
		ciphertext = encryptor.EncryptNew(plaintext)
	}

	return values, plaintext, ciphertext
}
