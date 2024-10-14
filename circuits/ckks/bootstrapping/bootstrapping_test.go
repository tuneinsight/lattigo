package bootstrapping

import (
	"flag"
	"math"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters + secure bootstrapping). Overrides -short and requires -timeout=0.")
var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")

var testPrec45 = ckks.ParametersLiteral{
	LogN:            10,
	LogQ:            []int{60, 40},
	LogP:            []int{61},
	LogDefaultScale: 40,
}

func TestBootstrapping(t *testing.T) {

	t.Run("BootstrappingWithoutRingDegreeSwitch", func(t *testing.T) {

		schemeParamsLit := testPrec45
		btpParamsLit := ParametersLiteral{}

		if *flagLongTest {
			schemeParamsLit.LogN = 16
		}

		params, err := ckks.NewParametersFromLiteral(schemeParamsLit)
		require.Nil(t, err)

		btpParamsLit.LogN = utils.Pointy(params.LogN())

		btpParams, err := NewParametersFromLiteral(params, btpParamsLit)
		require.Nil(t, err)

		// Insecure params for fast testing only
		if !*flagLongTest {
			btpParams.SlotsToCoeffsParameters.LogSlots = btpParams.BootstrappingParameters.LogN() - 1
			btpParams.CoeffsToSlotsParameters.LogSlots = btpParams.BootstrappingParameters.LogN() - 1

			// Corrects the message ratio to take into account the smaller number of slots and keep the same precision
			btpParams.Mod1ParametersLiteral.LogMessageRatio += 16 - params.LogN()
		}

		t.Logf("ParamsN2: LogN=%d/LogSlots=%d/LogQP=%f", params.LogN(), params.LogMaxSlots(), params.LogQP())

		sk := rlwe.NewKeyGenerator(btpParams.BootstrappingParameters).GenSecretKeyNew()

		t.Log("Generating Bootstrapping Keys")
		btpKeys, _, err := btpParams.GenEvaluationKeys(sk)
		require.NoError(t, err)

		evaluator, err := NewEvaluator(btpParams, btpKeys)
		require.NoError(t, err)

		ecd := ckks.NewEncoder(params)
		enc := rlwe.NewEncryptor(params, sk)
		dec := rlwe.NewDecryptor(params, sk)

		values := make([]complex128, params.MaxSlots())
		for i := range values {
			values[i] = sampling.RandComplex128(-1, 1)
		}

		values[0] = complex(0.9238795325112867, 0.3826834323650898)
		values[1] = complex(0.9238795325112867, 0.3826834323650898)
		if len(values) > 2 {
			values[2] = complex(0.9238795325112867, 0.3826834323650898)
			values[3] = complex(0.9238795325112867, 0.3826834323650898)
		}

		t.Run("Bootstrapping", func(t *testing.T) {

			plaintext := ckks.NewPlaintext(params, 0)
			ecd.Encode(values, plaintext)

			ctQ0, err := enc.EncryptNew(plaintext)
			require.NoError(t, err)

			// Checks that the input ciphertext is at the level 0
			require.True(t, ctQ0.Level() == 0)

			// Bootstrapps the ciphertext
			ctQL, err := evaluator.Bootstrap(ctQ0)
			require.NoError(t, err)

			// Checks that the output ciphertext is at the max level of paramsN1
			require.True(t, ctQL.Level() == params.MaxLevel())
			require.True(t, ctQL.Scale.Equal(params.DefaultScale()))

			verifyTestVectorsBootstrapping(params, ecd, dec, values, ctQL, t)
		})
	})

	t.Run("BootstrappingWithRingDegreeSwitch", func(t *testing.T) {

		schemeParamsLit := testPrec45
		btpParamsLit := ParametersLiteral{}

		if *flagLongTest {
			schemeParamsLit.LogN = 16
		}

		schemeParamsLit.LogNthRoot = schemeParamsLit.LogN + 1
		schemeParamsLit.LogN--

		params, err := ckks.NewParametersFromLiteral(schemeParamsLit)
		require.Nil(t, err)

		btpParamsLit.LogN = utils.Pointy(params.LogN() + 1)

		btpParams, err := NewParametersFromLiteral(params, btpParamsLit)
		require.Nil(t, err)

		// Insecure params for fast testing only
		if !*flagLongTest {
			btpParams.SlotsToCoeffsParameters.LogSlots = btpParams.BootstrappingParameters.LogN() - 1
			btpParams.CoeffsToSlotsParameters.LogSlots = btpParams.BootstrappingParameters.LogN() - 1

			// Corrects the message ratio to take into account the smaller number of slots and keep the same precision
			btpParams.Mod1ParametersLiteral.LogMessageRatio += 16 - params.LogN()
		}

		t.Logf("Params: LogN=%d/LogSlots=%d/LogQP=%f", btpParams.ResidualParameters.LogN(), btpParams.ResidualParameters.LogMaxSlots(), btpParams.ResidualParameters.LogQP())
		t.Logf("BTPParams: LogN=%d/LogSlots=%d/LogQP=%f", btpParams.BootstrappingParameters.LogN(), btpParams.BootstrappingParameters.LogMaxSlots(), btpParams.BootstrappingParameters.LogQP())

		sk := rlwe.NewKeyGenerator(params).GenSecretKeyNew()

		t.Log("Generating Bootstrapping Keys")
		btpKeys, _, err := btpParams.GenEvaluationKeys(sk)
		require.Nil(t, err)

		evaluator, err := NewEvaluator(btpParams, btpKeys)
		require.Nil(t, err)

		ecd := ckks.NewEncoder(params)
		enc := rlwe.NewEncryptor(params, sk)
		dec := rlwe.NewDecryptor(params, sk)

		values := make([]complex128, params.MaxSlots())
		for i := range values {
			values[i] = sampling.RandComplex128(-1, 1)
		}

		values[0] = complex(0.9238795325112867, 0.3826834323650898)
		values[1] = complex(0.9238795325112867, 0.3826834323650898)
		if len(values) > 2 {
			values[2] = complex(0.9238795325112867, 0.3826834323650898)
			values[3] = complex(0.9238795325112867, 0.3826834323650898)
		}

		t.Run("N1ToN2->Bootstrapping->N2ToN1", func(t *testing.T) {

			plaintext := ckks.NewPlaintext(params, 0)
			ecd.Encode(values, plaintext)

			ctQ0, err := enc.EncryptNew(plaintext)
			require.NoError(t, err)

			// Checks that the input ciphertext is at the level 0
			require.True(t, ctQ0.Level() == 0)

			// Bootstrapps the ciphertext
			ctQL, err := evaluator.Bootstrap(ctQ0)

			if err != nil {
				t.Fatal(err)
			}

			// Checks that the output ciphertext is at the max level of params
			require.True(t, ctQL.Level() == params.MaxLevel())
			require.True(t, ctQL.Scale.Equal(params.DefaultScale()))

			verifyTestVectorsBootstrapping(params, ecd, dec, values, ctQL, t)

		})
	})

	t.Run("BootstrappingPackedWithRingDegreeSwitch", func(t *testing.T) {

		schemeParamsLit := testPrec45
		btpParamsLit := ParametersLiteral{}

		if *flagLongTest {
			schemeParamsLit.LogN = 16
		}

		btpParamsLit.LogN = utils.Pointy(schemeParamsLit.LogN)
		schemeParamsLit.LogNthRoot = schemeParamsLit.LogN + 1
		schemeParamsLit.LogN -= 3

		params, err := ckks.NewParametersFromLiteral(schemeParamsLit)
		require.Nil(t, err)

		btpParams, err := NewParametersFromLiteral(params, btpParamsLit)
		require.Nil(t, err)

		// Insecure params for fast testing only
		if !*flagLongTest {
			btpParams.SlotsToCoeffsParameters.LogSlots = btpParams.BootstrappingParameters.LogN() - 1
			btpParams.CoeffsToSlotsParameters.LogSlots = btpParams.BootstrappingParameters.LogN() - 1

			// Corrects the message ratio to take into account the smaller number of slots and keep the same precision
			btpParams.Mod1ParametersLiteral.LogMessageRatio += 16 - params.LogN()
		}

		t.Logf("Params: LogN=%d/LogSlots=%d/LogQP=%f", btpParams.ResidualParameters.LogN(), btpParams.ResidualParameters.LogMaxSlots(), btpParams.ResidualParameters.LogQP())
		t.Logf("BTPParams: LogN=%d/LogSlots=%d/LogQP=%f", btpParams.BootstrappingParameters.LogN(), btpParams.BootstrappingParameters.LogMaxSlots(), btpParams.BootstrappingParameters.LogQP())

		sk := rlwe.NewKeyGenerator(params).GenSecretKeyNew()

		t.Log("Generating Bootstrapping Keys")
		btpKeys, _, err := btpParams.GenEvaluationKeys(sk)
		require.Nil(t, err)

		evaluator, err := NewEvaluator(btpParams, btpKeys)
		require.Nil(t, err)

		ecd := ckks.NewEncoder(params)
		enc := rlwe.NewEncryptor(params, sk)
		dec := rlwe.NewDecryptor(params, sk)

		values := make([]complex128, params.MaxSlots())
		for i := range values {
			values[i] = sampling.RandComplex128(-1, 1)
		}

		values[0] = complex(0.9238795325112867, 0.3826834323650898)
		values[1] = complex(0.9238795325112867, 0.3826834323650898)
		if len(values) > 2 {
			values[2] = complex(0.9238795325112867, 0.3826834323650898)
			values[3] = complex(0.9238795325112867, 0.3826834323650898)
		}

		pt := ckks.NewPlaintext(params, 0)

		cts := make([]rlwe.Ciphertext, 7)
		for i := range cts {

			require.NoError(t, ecd.Encode(utils.RotateSlice(values, i), pt))

			ct, err := enc.EncryptNew(pt)
			require.NoError(t, err)

			cts[i] = *ct
		}

		if cts, err = evaluator.BootstrapMany(cts); err != nil {
			t.Fatal(err)
		}

		for i := range cts {
			// Checks that the output ciphertext is at the max level of paramsN1
			require.True(t, cts[i].Level() == params.MaxLevel())
			require.True(t, cts[i].Scale.Equal(params.DefaultScale()))

			verifyTestVectorsBootstrapping(params, ecd, dec, utils.RotateSlice(values, i), &cts[i], t)
		}
	})

	t.Run("BootstrappingWithRingTypeSwitch", func(t *testing.T) {

		schemeParamsLit := testPrec45
		schemeParamsLit.RingType = ring.ConjugateInvariant
		btpParamsLit := ParametersLiteral{}

		if *flagLongTest {
			schemeParamsLit.LogN = 16
		}

		btpParamsLit.LogN = utils.Pointy(schemeParamsLit.LogN)
		schemeParamsLit.LogNthRoot = schemeParamsLit.LogN + 1
		schemeParamsLit.LogN--

		params, err := ckks.NewParametersFromLiteral(schemeParamsLit)
		require.Nil(t, err)

		btpParams, err := NewParametersFromLiteral(params, btpParamsLit)
		require.Nil(t, err)

		// Insecure params for fast testing only
		if !*flagLongTest {
			// Corrects the message ratio to take into account the smaller number of slots and keep the same precision
			btpParams.Mod1ParametersLiteral.LogMessageRatio += 16 - params.LogN()
		}

		t.Logf("Params: LogN=%d/LogSlots=%d/LogQP=%f", btpParams.ResidualParameters.LogN(), btpParams.ResidualParameters.LogMaxSlots(), btpParams.ResidualParameters.LogQP())
		t.Logf("BTPParams: LogN=%d/LogSlots=%d/LogQP=%f", btpParams.BootstrappingParameters.LogN(), btpParams.BootstrappingParameters.LogMaxSlots(), btpParams.BootstrappingParameters.LogQP())

		sk := rlwe.NewKeyGenerator(params).GenSecretKeyNew()

		t.Log("Generating Bootstrapping Keys")
		btpKeys, _, err := btpParams.GenEvaluationKeys(sk)
		require.Nil(t, err)

		evaluator, err := NewEvaluator(btpParams, btpKeys)
		require.Nil(t, err)

		ecd := ckks.NewEncoder(params)
		enc := rlwe.NewEncryptor(params, sk)
		dec := rlwe.NewDecryptor(params, sk)

		values := make([]float64, params.MaxSlots())
		for i := range values {
			values[i] = sampling.RandFloat64(-1, 1)
		}

		values[0] = 0.9238795325112867
		values[1] = 0.9238795325112867
		if len(values) > 2 {
			values[2] = 0.9238795325112867
			values[3] = 0.9238795325112867
		}

		t.Run("ConjugateInvariant->Standard->Bootstrapping->Standard->ConjugateInvariant", func(t *testing.T) {

			plaintext := ckks.NewPlaintext(params, 0)
			require.NoError(t, ecd.Encode(values, plaintext))

			ctLeftQ0, err := enc.EncryptNew(plaintext)
			require.NoError(t, err)
			ctRightQ0, err := enc.EncryptNew(plaintext)
			require.NoError(t, err)

			// Checks that the input ciphertext is at the level 0
			require.True(t, ctLeftQ0.Level() == 0)
			require.True(t, ctRightQ0.Level() == 0)

			// Bootstraps the ciphertext
			ctLeftQL, ctRightQL, err := evaluator.EvaluateConjugateInvariant(ctLeftQ0, ctRightQ0)

			require.NoError(t, err)

			// Checks that the output ciphertext is at the max level of paramsN1
			require.True(t, ctLeftQL.Level() == params.MaxLevel())
			require.True(t, ctLeftQL.Scale.Equal(params.DefaultScale()))

			verifyTestVectorsBootstrapping(params, ecd, dec, values, ctLeftQL, t)

			require.True(t, ctRightQL.Level() == params.MaxLevel())
			require.True(t, ctRightQL.Scale.Equal(params.DefaultScale()))
			verifyTestVectorsBootstrapping(params, ecd, dec, values, ctRightQL, t)
		})
	})
}

func verifyTestVectorsBootstrapping(params ckks.Parameters, encoder *ckks.Encoder, decryptor *rlwe.Decryptor, valuesWant, element interface{}, t *testing.T) {
	precStats := ckks.GetPrecisionStats(params, encoder, decryptor, valuesWant, element, 0, false)
	if *printPrecisionStats {
		t.Log(precStats.String())
	}

	rf64 := precStats.AVGLog2Prec.Real
	if64 := precStats.AVGLog2Prec.Imag

	minPrec := math.Log2(params.DefaultScale().Float64()) - float64(params.LogN()+2)
	if minPrec < 0 {
		minPrec = 0
	}

	minPrec -= 10

	require.GreaterOrEqual(t, rf64, minPrec)
	require.GreaterOrEqual(t, if64, minPrec)
}
