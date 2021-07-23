package ckks

import (
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
	"github.com/stretchr/testify/assert"
	"math"
	"runtime"
	"testing"
)

func TestCKKSAdvancedHomomorphicMod(t *testing.T) {
	var err error

	if runtime.GOARCH == "wasm" {
		t.Skip("skipping homomorphic mod tests for GOARCH=wasm")
	}

	ParametersLiteral := ParametersLiteral{
		LogN:     14,
		LogSlots: 13,
		Scale:    1 << 45,
		Sigma:    rlwe.DefaultSigma,
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
	}

	testEvalModMarshalling(t)

	var params Parameters
	if params, err = NewParametersFromLiteral(ParametersLiteral); err != nil {
		panic(err)
	}

	for _, testSet := range []func(params Parameters, t *testing.T){
		testEvalMod,
	} {
		testSet(params, t)
		runtime.GC()
	}

}

func testEvalModMarshalling(t *testing.T) {
	t.Run("Marshalling", func(t *testing.T) {

		evm := EvalModParameters{
			Q:             0x80000000080001,
			LevelStart:    12,
			SineType:      Sin,
			MessageRatio:  256.0,
			K:             14,
			SineDeg:       127,
			DoubleAngle:   0,
			ArcSineDeg:    7,
			ScalingFactor: 1 << 60,
		}

		data, err := evm.MarshalBinary()
		assert.Nil(t, err)

		evmNew := new(EvalModParameters)
		if err := evmNew.UnmarshalBinary(data); err != nil {
			assert.Nil(t, err)
		}
		assert.Equal(t, evm, *evmNew)

	})

}

func testEvalMod(params Parameters, t *testing.T) {

	kgen := NewKeyGenerator(params)
	sk := kgen.GenSecretKey()
	rlk := kgen.GenRelinearizationKey(sk, 2)
	encoder := NewEncoder(params)
	encryptor := NewEncryptor(params, sk)
	decryptor := NewDecryptor(params, sk)
	eval := NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: nil})

	t.Run("SineChebyshevWithArcSine", func(t *testing.T) {

		evm := EvalModParameters{
			Q:             0x80000000080001,
			LevelStart:    12,
			SineType:      Sin,
			MessageRatio:  256.0,
			K:             14,
			SineDeg:       127,
			DoubleAngle:   0,
			ArcSineDeg:    7,
			ScalingFactor: 1 << 60,
		}

		EvalModPoly := evm.GenPoly()

		values, _, ciphertext := newTestVectorsEvalMod(params, encryptor, encoder, evm, t)

		scale := math.Exp2(math.Round(math.Log2(float64(evm.Q) / evm.MessageRatio)))

		// Scale the message to Delta = Q/MessageRatio
		eval.ScaleUp(ciphertext, math.Round(scale/ciphertext.Scale), ciphertext)

		// Scale the message up to Sine/MessageRatio
		eval.ScaleUp(ciphertext, math.Round((EvalModPoly.ScalingFactor/evm.MessageRatio)/ciphertext.Scale), ciphertext)

		// Normalization
		eval.MultByConst(ciphertext, 1/(float64(evm.K)*evm.QDiff()), ciphertext)
		eval.Rescale(ciphertext, params.Scale(), ciphertext)

		// EvalMod
		ciphertext = eval.EvalMod(ciphertext, EvalModPoly)

		// PlaintextCircuit
		//pi2r := 6.283185307179586/complex(math.Exp2(float64(evm.DoubleAngle)), 0)
		for i := range values {
			values[i] -= complex(evm.MessageRatio*evm.QDiff()*math.Round(real(values[i])/(evm.MessageRatio/evm.QDiff())), 0)
		}

		verifyTestVectors(params, encoder, decryptor, values, ciphertext, params.LogSlots(), 0, t)
	})

	t.Run("CosOptimizedChebyshevWithArcSine", func(t *testing.T) {

		evm := EvalModParameters{
			Q:             0x80000000080001,
			LevelStart:    12,
			SineType:      Cos1,
			MessageRatio:  256.0,
			K:             25,
			SineDeg:       63,
			DoubleAngle:   2,
			ArcSineDeg:    7,
			ScalingFactor: 1 << 60,
		}

		EvalModPoly := evm.GenPoly()

		values, _, ciphertext := newTestVectorsEvalMod(params, encryptor, encoder, evm, t)

		scale := math.Exp2(math.Round(math.Log2(float64(evm.Q) / evm.MessageRatio)))

		// Scale the message to Delta = Q/MessageRatio
		eval.ScaleUp(ciphertext, math.Round(scale/ciphertext.Scale), ciphertext)

		// Scale the message up to Sine/MessageRatio
		eval.ScaleUp(ciphertext, math.Round((EvalModPoly.ScalingFactor/evm.MessageRatio)/ciphertext.Scale), ciphertext)

		// Normalization
		eval.MultByConst(ciphertext, 1/(float64(evm.K)*evm.QDiff()), ciphertext)
		eval.Rescale(ciphertext, params.Scale(), ciphertext)

		// EvalMod
		ciphertext = eval.EvalMod(ciphertext, EvalModPoly)

		// PlaintextCircuit
		//pi2r := 6.283185307179586/complex(math.Exp2(float64(evm.DoubleAngle)), 0)
		for i := range values {
			values[i] -= complex(evm.MessageRatio*evm.QDiff()*math.Round(real(values[i])/(evm.MessageRatio/evm.QDiff())), 0)
		}

		verifyTestVectors(params, encoder, decryptor, values, ciphertext, params.LogSlots(), 0, t)
	})
}

func newTestVectorsEvalMod(params Parameters, encryptor Encryptor, encoder Encoder, evm EvalModParameters, t *testing.T) (values []complex128, plaintext *Plaintext, ciphertext *Ciphertext) {

	logSlots := params.LogSlots()

	values = make([]complex128, 1<<logSlots)

	K := float64(evm.K - 1)
	Q := float64(evm.Q) / math.Exp2(math.Round(math.Log2(float64(evm.Q)))) * evm.MessageRatio

	for i := uint64(0); i < 1<<logSlots; i++ {
		values[i] = complex(math.Round(utils.RandFloat64(-K, K))*Q+utils.RandFloat64(-1, 1), 0)
	}

	values[0] = complex(K*Q+0.5, 0)

	plaintext = NewPlaintext(params, params.MaxLevel(), params.Scale())

	encoder.EncodeNTT(plaintext, values, logSlots)

	if encryptor != nil {
		ciphertext = encryptor.EncryptNew(plaintext)
	}

	return values, plaintext, ciphertext
}
