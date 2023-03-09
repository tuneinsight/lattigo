package bootstrapping

import (
	"flag"
	"fmt"
	"runtime"
	"sync"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

var minPrec float64 = 12.0

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters + secure bootstrapping). Overrides -short and requires -timeout=0.")
var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")

func ParamsToString(params ckks.Parameters, opname string) string {
	return fmt.Sprintf("%slogN=%d/LogSlots=%d/logQP=%f/levels=%d/a=%d/b=%d",
		opname,
		params.LogN(),
		params.LogSlots(),
		params.LogQP(),
		params.MaxLevel()+1,
		params.PCount(),
		params.DecompRNS(params.MaxLevelQ(), params.MaxLevelP()))
}

func TestBootstrapParametersMarshalling(t *testing.T) {

	t.Run("ParametersLiteral", func(t *testing.T) {

		paramsLit := ParametersLiteral{
			CoeffsToSlotsFactorizationDepthAndLogScales: [][]int{{53}, {53}, {53}, {53}},
			SlotsToCoeffsFactorizationDepthAndLogScales: [][]int{{30}, {30, 30}},
			EvalModLogScale:       utils.PointyInt(59),
			EphemeralSecretWeight: utils.PointyInt(1),
			Iterations:            utils.PointyInt(2),
			SineDegree:            utils.PointyInt(32),
			ArcSineDegree:         utils.PointyInt(7),
		}

		data, err := paramsLit.MarshalBinary()
		require.Nil(t, err)

		paramsLitNew := new(ParametersLiteral)
		if err := paramsLitNew.UnmarshalBinary(data); err != nil {
			require.Nil(t, err)
		}
		require.Equal(t, paramsLit, *paramsLitNew)
	})

	t.Run("Parameters", func(t *testing.T) {
		paramSet := DefaultParametersSparse[0]

		_, btpParams, err := NewParametersFromLiteral(paramSet.SchemeParams, paramSet.BootstrappingParams)
		require.Nil(t, err)

		data, err := btpParams.MarshalBinary()
		require.Nil(t, err)
		btpParamsNew := new(Parameters)
		if err := btpParamsNew.UnmarshalBinary(data); err != nil {
			require.Nil(t, err)
		}

		require.Equal(t, btpParams, *btpParamsNew)
	})
}

func TestBootstrap(t *testing.T) {

	if runtime.GOARCH == "wasm" {
		t.Skip("skipping bootstrapping tests for GOARCH=wasm")
	}

	paramSet := DefaultParametersSparse[0]

	ckksParamsLit, btpParams, err := NewParametersFromLiteral(paramSet.SchemeParams, paramSet.BootstrappingParams)
	require.Nil(t, err)

	// Insecure params for fast testing only
	if !*flagLongTest {
		ckksParamsLit.LogN = 13

		// Corrects the message ratio to take into account the smaller number of slots and keep the same precision
		btpParams.EvalModParameters.LogMessageRatio += paramSet.SchemeParams.LogN - 1 - ckksParamsLit.LogN - 1

		ckksParamsLit.LogSlots = ckksParamsLit.LogN - 1
	}

	Xs := ckksParamsLit.Xs

	EphemeralSecretWeight := btpParams.EphemeralSecretWeight

	for _, testSet := range [][]bool{{false, false}, {true, false}, {false, true}, {true, true}} {

		if testSet[0] {
			ckksParamsLit.Xs = &distribution.Ternary{H: EphemeralSecretWeight}
			btpParams.EphemeralSecretWeight = 0
		} else {
			ckksParamsLit.Xs = Xs
			btpParams.EphemeralSecretWeight = EphemeralSecretWeight
		}

		if testSet[1] {
			ckksParamsLit.LogSlots = ckksParamsLit.LogN - 2
		} else {
			ckksParamsLit.LogSlots = ckksParamsLit.LogN - 1
		}

		params, err := ckks.NewParametersFromLiteral(ckksParamsLit)
		if err != nil {
			panic(err)
		}

		testbootstrap(params, testSet[0], btpParams, t)
		runtime.GC()
	}
}

func testbootstrap(params ckks.Parameters, original bool, btpParams Parameters, t *testing.T) {

	btpType := "Encapsulation/"

	if original {
		btpType = "Original/"
	}

	t.Run(ParamsToString(params, "Bootstrapping/FullCircuit/"+btpType), func(t *testing.T) {

		kgen := ckks.NewKeyGenerator(params)
		sk := kgen.GenSecretKeyNew()
		encoder := ckks.NewEncoder(params)
		encryptor := ckks.NewEncryptor(params, sk)
		decryptor := ckks.NewDecryptor(params, sk)

		evk := GenEvaluationKeySetNew(btpParams, params, sk)

		btp, err := NewBootstrapper(params, btpParams, evk)
		if err != nil {
			panic(err)
		}

		values := make([]complex128, 1<<params.LogSlots())
		for i := range values {
			values[i] = sampling.RandComplex128(-1, 1)
		}

		values[0] = complex(0.9238795325112867, 0.3826834323650898)
		values[1] = complex(0.9238795325112867, 0.3826834323650898)
		if 1<<params.LogSlots() > 2 {
			values[2] = complex(0.9238795325112867, 0.3826834323650898)
			values[3] = complex(0.9238795325112867, 0.3826834323650898)
		}

		plaintext := ckks.NewPlaintext(params, 0)
		encoder.Encode(values, plaintext, params.LogSlots())

		n := 2

		ciphertexts := make([]*rlwe.Ciphertext, n)
		bootstrappers := make([]*Bootstrapper, n)
		bootstrappers[0] = btp
		ciphertexts[0] = encryptor.EncryptNew(plaintext)
		for i := 1; i < len(ciphertexts); i++ {
			ciphertexts[i] = encryptor.EncryptNew(plaintext)
			bootstrappers[i] = bootstrappers[0].ShallowCopy()
		}

		var wg sync.WaitGroup
		wg.Add(n)
		for i := range ciphertexts {
			go func(index int) {
				ciphertexts[index] = bootstrappers[index].Bootstrap(ciphertexts[index])
				wg.Done()
			}(i)
		}
		wg.Wait()

		for i := range ciphertexts {
			verifyTestVectors(params, encoder, decryptor, values, ciphertexts[i], params.LogSlots(), t)
		}
	})
}

func verifyTestVectors(params ckks.Parameters, encoder ckks.Encoder, decryptor rlwe.Decryptor, valuesWant []complex128, element interface{}, logSlots int, t *testing.T) {
	precStats := ckks.GetPrecisionStats(params, encoder, decryptor, valuesWant, element, logSlots, nil)
	if *printPrecisionStats {
		t.Log(precStats.String())
	}

	require.GreaterOrEqual(t, precStats.MeanPrecision.Real, minPrec)
	require.GreaterOrEqual(t, precStats.MeanPrecision.Imag, minPrec)
}
