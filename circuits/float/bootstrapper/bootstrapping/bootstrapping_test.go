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

func ParamsToString(params ckks.Parameters, LogSlots int, opname string) string {
	return fmt.Sprintf("%slogN=%d/LogSlots=%d/logQP=%f/levels=%d/a=%d/b=%d",
		opname,
		params.LogN(),
		LogSlots,
		params.LogQP(),
		params.MaxLevel()+1,
		params.PCount(),
		params.BaseRNSDecompositionVectorSize(params.MaxLevelQ(), params.MaxLevelP()))
}

func TestBootstrapParametersMarshalling(t *testing.T) {

	t.Run("ParametersLiteral", func(t *testing.T) {

		paramsLit := ParametersLiteral{
			CoeffsToSlotsFactorizationDepthAndLogScales: [][]int{{53}, {53}, {53}, {53}},
			SlotsToCoeffsFactorizationDepthAndLogScales: [][]int{{30}, {30, 30}},
			EvalModLogScale:       utils.Pointy(59),
			EphemeralSecretWeight: utils.Pointy(1),
			IterationsParameters:  &IterationsParameters{BootstrappingPrecision: []float64{20, 20}, ReservedPrimeBitSize: 20},
			SineDegree:            utils.Pointy(32),
			ArcSineDegree:         utils.Pointy(7),
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

		params, err := ckks.NewParametersFromLiteral(paramSet.SchemeParams)
		require.Nil(t, err)

		btpParams, err := NewParametersFromLiteral(params, paramSet.BootstrappingParams)
		require.Nil(t, err)

		data, err := btpParams.MarshalBinary()
		require.Nil(t, err)
		btpParamsNew := new(Parameters)
		if err := btpParamsNew.UnmarshalBinary(data); err != nil {
			require.Nil(t, err)
		}

		require.True(t, btpParams.Equal(btpParamsNew))
	})
}

func TestBootstrappingWithEncapsulation(t *testing.T) {

	if runtime.GOARCH == "wasm" {
		t.Skip("skipping bootstrapping tests for GOARCH=wasm")
	}

	paramSet := DefaultParametersSparse[0]

	if !*flagLongTest {
		paramSet.SchemeParams.LogN -= 3
	}

	paramSet.BootstrappingParams.LogN = utils.Pointy(paramSet.SchemeParams.LogN)

	for _, LogSlots := range []int{1, paramSet.SchemeParams.LogN - 2, paramSet.SchemeParams.LogN - 1} {
		paramsSetCpy := paramSet

		level := utils.Min(1, len(paramSet.SchemeParams.LogQ))

		paramsSetCpy.SchemeParams.LogQ = paramSet.SchemeParams.LogQ[:level+1]

		paramsSetCpy.BootstrappingParams.LogSlots = &LogSlots

		params, err := ckks.NewParametersFromLiteral(paramsSetCpy.SchemeParams)
		require.NoError(t, err)

		btpParams, err := NewParametersFromLiteral(params, paramsSetCpy.BootstrappingParams)
		require.NoError(t, err)

		// Insecure params for fast testing only
		if !*flagLongTest {
			// Corrects the message ratio to take into account the smaller number of slots and keep the same precision
			btpParams.Mod1ParametersLiteral.LogMessageRatio += utils.Min(utils.Max(15-LogSlots, 0), 8)
		}

		testbootstrap(params, btpParams, level, t)
		runtime.GC()
	}

	testBootstrapHighPrecision(paramSet, t)
}

func TestBootstrappingOriginal(t *testing.T) {

	if runtime.GOARCH == "wasm" {
		t.Skip("skipping bootstrapping tests for GOARCH=wasm")
	}

	paramSet := DefaultParametersDense[0]

	if !*flagLongTest {
		paramSet.SchemeParams.LogN -= 3
	}

	paramSet.BootstrappingParams.LogN = utils.Pointy(paramSet.SchemeParams.LogN)

	for _, LogSlots := range []int{1, paramSet.SchemeParams.LogN - 2, paramSet.SchemeParams.LogN - 1} {

		paramsSetCpy := paramSet

		level := utils.Min(1, len(paramSet.SchemeParams.LogQ))

		paramsSetCpy.SchemeParams.LogQ = paramSet.SchemeParams.LogQ[:level+1]

		paramsSetCpy.BootstrappingParams.LogSlots = &LogSlots

		params, err := ckks.NewParametersFromLiteral(paramsSetCpy.SchemeParams)
		require.NoError(t, err)

		btpParams, err := NewParametersFromLiteral(params, paramsSetCpy.BootstrappingParams)
		require.NoError(t, err)

		// Insecure params for fast testing only
		if !*flagLongTest {
			// Corrects the message ratio to take into account the smaller number of slots and keep the same precision
			btpParams.Mod1ParametersLiteral.LogMessageRatio += utils.Min(utils.Max(15-LogSlots, 0), 8)
		}

		testbootstrap(params, btpParams, level, t)
		runtime.GC()
	}

	testBootstrapHighPrecision(paramSet, t)
}

func testbootstrap(params ckks.Parameters, btpParams Parameters, level int, t *testing.T) {

	t.Run(ParamsToString(params, btpParams.LogMaxSlots(), ""), func(t *testing.T) {

		kgen := ckks.NewKeyGenerator(btpParams.Parameters)
		sk := kgen.GenSecretKeyNew()
		encoder := ckks.NewEncoder(params)

		encryptor := ckks.NewEncryptor(params, sk)
		decryptor := ckks.NewDecryptor(params, sk)

		evk := btpParams.GenEvaluationKeySetNew(sk)

		btp, err := NewBootstrapper(btpParams, evk)
		require.NoError(t, err)

		values := make([]complex128, 1<<btpParams.LogMaxSlots())
		for i := range values {
			values[i] = sampling.RandComplex128(-1, 1)
		}

		values[0] = complex(0.9238795325112867, 0.3826834323650898)
		values[1] = complex(0.9238795325112867, 0.3826834323650898)

		if btpParams.LogMaxSlots() > 1 {
			values[2] = complex(0.9238795325112867, 0.3826834323650898)
			values[3] = complex(0.9238795325112867, 0.3826834323650898)
		}

		plaintext := ckks.NewPlaintext(params, 0)
		plaintext.Scale = params.DefaultScale()
		plaintext.LogDimensions = btpParams.LogMaxDimensions()
		encoder.Encode(values, plaintext)

		n := 1

		ciphertexts := make([]*rlwe.Ciphertext, n)
		bootstrappers := make([]*Bootstrapper, n)
		bootstrappers[0] = btp
		ciphertexts[0], err = encryptor.EncryptNew(plaintext)
		require.NoError(t, err)
		for i := 1; i < len(ciphertexts); i++ {
			ciphertexts[i], err = encryptor.EncryptNew(plaintext)
			require.NoError(t, err)
			bootstrappers[i] = bootstrappers[0].ShallowCopy()
		}

		var wg sync.WaitGroup
		wg.Add(n)
		for i := range ciphertexts {
			go func(index int) {
				var err error
				ciphertexts[index], err = bootstrappers[index].Bootstrap(ciphertexts[index])
				require.NoError(t, err)
				wg.Done()
			}(i)
		}
		wg.Wait()

		for i := range ciphertexts {
			require.True(t, ciphertexts[i].Level() == level)
		}

		for i := range ciphertexts {
			verifyTestVectors(params, encoder, decryptor, values, ciphertexts[i], t)
		}
	})
}

func testBootstrapHighPrecision(paramSet defaultParametersLiteral, t *testing.T) {

	t.Run("HighPrecision", func(t *testing.T) {

		level := utils.Min(4, len(paramSet.SchemeParams.LogQ))

		paramSet.SchemeParams.LogQ = paramSet.SchemeParams.LogQ[:level+1]

		paramSet.BootstrappingParams.IterationsParameters = &IterationsParameters{
			BootstrappingPrecision: []float64{24.5, 24.5, 24.5, 24.5, 24.5},
			ReservedPrimeBitSize:   28,
		}

		params, err := ckks.NewParametersFromLiteral(paramSet.SchemeParams)
		if err != nil {
			panic(err)
		}

		btpParams, err := NewParametersFromLiteral(params, paramSet.BootstrappingParams)

		if err != nil {
			t.Fatal(err)
		}

		// Insecure params for fast testing only
		if !*flagLongTest {
			// Corrects the message ratio to take into account the smaller number of slots and keep the same precision
			btpParams.Mod1ParametersLiteral.LogMessageRatio += utils.Min(utils.Max(16-params.LogN(), 0), 8)
		}

		t.Run(ParamsToString(params, btpParams.LogMaxSlots(), ""), func(t *testing.T) {

			kgen := ckks.NewKeyGenerator(btpParams.Parameters)
			sk := kgen.GenSecretKeyNew()
			encoder := ckks.NewEncoder(params, 164)
			encryptor := ckks.NewEncryptor(params, sk)
			decryptor := ckks.NewDecryptor(params, sk)

			evk := btpParams.GenEvaluationKeySetNew(sk)

			bootstrapper, err := NewBootstrapper(btpParams, evk)
			require.NoError(t, err)

			values := make([]complex128, 1<<btpParams.LogMaxSlots())
			for i := range values {
				values[i] = sampling.RandComplex128(-1, 1)
			}

			values[0] = complex(0.9238795325112867, 0.3826834323650898)
			values[1] = complex(0.9238795325112867, 0.3826834323650898)

			if btpParams.LogMaxSlots() > 1 {
				values[2] = complex(0.9238795325112867, 0.3826834323650898)
				values[3] = complex(0.9238795325112867, 0.3826834323650898)
			}

			plaintext := ckks.NewPlaintext(params, level-1)
			plaintext.Scale = params.DefaultScale()
			for i := 0; i < plaintext.Level(); i++ {
				plaintext.Scale = plaintext.Scale.Mul(rlwe.NewScale(1 << 40))
			}

			plaintext.LogDimensions = btpParams.LogMaxDimensions()
			encoder.Encode(values, plaintext)

			ciphertext, err := encryptor.EncryptNew(plaintext)
			require.NoError(t, err)

			ciphertext, err = bootstrapper.Bootstrap(ciphertext)
			require.NoError(t, err)

			require.True(t, ciphertext.Level() == level)

			verifyTestVectors(params, encoder, decryptor, values, ciphertext, t)
		})

		runtime.GC()
	})
}

func verifyTestVectors(params ckks.Parameters, encoder *ckks.Encoder, decryptor *rlwe.Decryptor, valuesWant, valuesHave interface{}, t *testing.T) {
	precStats := ckks.GetPrecisionStats(params, encoder, decryptor, valuesWant, valuesHave, nil, false)
	if *printPrecisionStats {
		t.Log(precStats.String())
	}

	rf64, _ := precStats.MeanPrecision.Real.Float64()
	if64, _ := precStats.MeanPrecision.Imag.Float64()

	require.GreaterOrEqual(t, rf64, minPrec)
	require.GreaterOrEqual(t, if64, minPrec)
}
