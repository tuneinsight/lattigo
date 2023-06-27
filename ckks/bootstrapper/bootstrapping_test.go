package bootstrapper

import (
	"flag"
	"math"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ckks/bootstrapper/bootstrapping"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters + secure bootstrapping). Overrides -short and requires -timeout=0.")
var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")

func TestBootstrapping(t *testing.T) {

	paramSet := bootstrapping.DefaultParametersSparse[0]

	paramsN2Lit, btpParamsN2, err := bootstrapping.NewParametersFromLiteral(paramSet.SchemeParams, paramSet.BootstrappingParams)
	require.Nil(t, err)

	// Insecure params for fast testing only
	if !*flagLongTest {
		paramsN2Lit.LogN = 13
		btpParamsN2.SlotsToCoeffsParameters.LogSlots = paramsN2Lit.LogN - 1
		btpParamsN2.CoeffsToSlotsParameters.LogSlots = paramsN2Lit.LogN - 1

		// Corrects the message ratio to take into account the smaller number of slots and keep the same precision
		btpParamsN2.EvalModParameters.LogMessageRatio += paramSet.SchemeParams.LogN - paramsN2Lit.LogN
	}

	endLevel := len(paramSet.SchemeParams.LogQ) - 1

	require.True(t, endLevel == len(paramsN2Lit.Q)-1-btpParamsN2.Depth()) // Checks the depth of the bootstrapping

	t.Run("BootstrapingWithoutRingDegreeSwitch", func(t *testing.T) {

		paramsN2, err := ckks.NewParametersFromLiteral(paramsN2Lit)
		require.Nil(t, err)

		t.Logf("ParamsN2: LogN=%d/LogSlots=%d/LogQP=%f", paramsN2.LogN(), paramsN2.PlaintextSlots(), paramsN2.LogQP())

		skN2 := ckks.NewKeyGenerator(paramsN2).GenSecretKeyNew()

		t.Log("Generating Bootstrapping Keys")
		btpKeys, err := GenBootstrappingKeys(nil, paramsN2, btpParamsN2, nil, *skN2)
		require.Nil(t, err)

		bootstrapperInterface, err := NewBootstrapper(nil, paramsN2, btpParamsN2, btpKeys)
		require.Nil(t, err)

		bootstrapper := bootstrapperInterface.(*Bootstrapper)

		ecdN2 := ckks.NewEncoder(paramsN2)
		encN2 := ckks.NewEncryptor(paramsN2, skN2)
		decN2 := ckks.NewDecryptor(paramsN2, skN2)

		values := make([]complex128, paramsN2.PlaintextSlots())
		for i := range values {
			values[i] = sampling.RandComplex128(-1, 1)
		}

		values[0] = complex(0.9238795325112867, 0.3826834323650898)
		values[1] = complex(0.9238795325112867, 0.3826834323650898)
		if paramsN2.PlaintextSlots() > 2 {
			values[2] = complex(0.9238795325112867, 0.3826834323650898)
			values[3] = complex(0.9238795325112867, 0.3826834323650898)
		}

		t.Run("Bootstrapping", func(t *testing.T) {

			plaintext := ckks.NewPlaintext(paramsN2, 0)
			ecdN2.Encode(values, plaintext)

			ctN2Q0 := encN2.EncryptNew(plaintext)

			// Checks that the input ciphertext is at the level 0
			require.True(t, ctN2Q0.Level() == 0)

			// Bootstrapps the ciphertext
			ctN2QL := bootstrapper.refreshStandard(ctN2Q0)

			// Checks that the output ciphertext is at the max level of paramsN1
			require.True(t, ctN2QL.Level() == endLevel)
			require.True(t, ctN2QL.PlaintextScale.Equal(paramsN2.PlaintextScale()))

			verifyTestVectorsBootstrapping(paramsN2, ecdN2, decN2, values, ctN2QL, t)
		})
	})

	t.Run("BootstrappingWithRingDegreeSwitch", func(t *testing.T) {

		paramsN2, err := ckks.NewParametersFromLiteral(paramsN2Lit)
		require.Nil(t, err)

		parmasN1Lit := ckks.ParametersLiteral{
			LogN:              paramsN2Lit.LogN - 1,
			Q:                 paramsN2Lit.Q[:endLevel+1],
			P:                 []uint64{0x80000000440001, 0x80000000500001},
			LogPlaintextScale: paramsN2Lit.LogPlaintextScale,
		}

		paramsN1, err := ckks.NewParametersFromLiteral(parmasN1Lit)
		require.Nil(t, err)

		t.Logf("ParamsN2: LogN=%d/LogSlots=%d/LogQP=%f", paramsN2.LogN(), paramsN2.PlaintextLogSlots(), paramsN2.LogQP())
		t.Logf("ParamsN1: LogN=%d/LogSlots=%d/LogQP=%f", paramsN1.LogN(), paramsN1.PlaintextLogSlots(), paramsN1.LogQP())

		skN2 := ckks.NewKeyGenerator(paramsN2).GenSecretKeyNew()
		skN1 := ckks.NewKeyGenerator(paramsN1).GenSecretKeyNew()

		t.Log("Generating Bootstrapping Keys")
		btpKeys, err := GenBootstrappingKeys(&paramsN1, paramsN2, btpParamsN2, skN1, *skN2)
		require.Nil(t, err)

		bootstrapperInterface, err := NewBootstrapper(&paramsN1, paramsN2, btpParamsN2, btpKeys)
		require.Nil(t, err)

		bootstrapper := bootstrapperInterface.(*Bootstrapper)

		ecdN1 := ckks.NewEncoder(paramsN1)
		encN1 := ckks.NewEncryptor(paramsN1, skN1)
		decN1 := ckks.NewDecryptor(paramsN1, skN1)

		values := make([]complex128, paramsN1.PlaintextSlots())
		for i := range values {
			values[i] = sampling.RandComplex128(-1, 1)
		}

		values[0] = complex(0.9238795325112867, 0.3826834323650898)
		values[1] = complex(0.9238795325112867, 0.3826834323650898)
		if paramsN1.PlaintextSlots() > 2 {
			values[2] = complex(0.9238795325112867, 0.3826834323650898)
			values[3] = complex(0.9238795325112867, 0.3826834323650898)
		}

		t.Run("SwitchDimensionN1toN2", func(t *testing.T) {
			ptN1 := ckks.NewPlaintext(paramsN1, 0)
			ecdN1.Encode(values, ptN1)
			ctN1Q0 := encN1.EncryptNew(ptN1)
			ctN2Q0 := ckks.NewCiphertext(paramsN2, ctN1Q0.Degree(), ctN1Q0.Level())
			ctN2Q0.PlaintextScale = ctN1Q0.PlaintextScale
			bootstrapper.ringN1toRingN2(bootstrapper.bootstrapper.Evaluator.Evaluator, ctN1Q0, ctN2Q0)
			verifyTestVectorsBootstrapping(paramsN1, ckks.NewEncoder(paramsN2), ckks.NewDecryptor(paramsN2, skN2), values, ctN2Q0, t)
		})

		t.Run("SwitchDimensionN2toN1", func(t *testing.T) {
			ptN2 := ckks.NewPlaintext(paramsN2, paramsN1.MaxLevel())
			ptN2.PlaintextLogDimensions = paramsN1.PlaintextLogDimensions()

			ecdN2 := ckks.NewEncoder(paramsN2)
			encN2 := ckks.NewEncryptor(paramsN2, skN2)
			ecdN2.Encode(values, ptN2)
			ctN2QL := encN2.EncryptNew(ptN2)

			ctN1QL := ckks.NewCiphertext(paramsN1, ctN2QL.Degree(), ctN2QL.Level())
			ctN1QL.PlaintextScale = ctN2QL.PlaintextScale

			bootstrapper.ringN2toRingN1(bootstrapper.bootstrapper.Evaluator.Evaluator, ctN2QL, ctN1QL)

			verifyTestVectorsBootstrapping(paramsN1, ecdN1, decN1, values, ctN1QL, t)
		})

		t.Run("N1ToN2->Bootstrapping->N2ToN1", func(t *testing.T) {

			plaintext := ckks.NewPlaintext(paramsN1, 0)
			ecdN1.Encode(values, plaintext)

			ctN1Q0 := encN1.EncryptNew(plaintext)

			// Checks that the input ciphertext is at the level 0
			require.True(t, ctN1Q0.Level() == 0)

			// Bootstrapps the ciphertext
			ctN1QL := bootstrapper.refreshRingDegreeSwitch(ctN1Q0)

			// Checks that the output ciphertext is at the max level of paramsN1
			require.True(t, ctN1QL.Level() == paramsN1.MaxLevel())
			require.True(t, ctN1QL.PlaintextScale.Equal(paramsN1.PlaintextScale()))

			verifyTestVectorsBootstrapping(paramsN1, ecdN1, decN1, values, ctN1QL, t)

		})
	})

	t.Run("BootstrappingWithRingTypeSwitch", func(t *testing.T) {

		paramsN2, err := ckks.NewParametersFromLiteral(paramsN2Lit)
		require.Nil(t, err)

		parmasN1Lit := ckks.ParametersLiteral{
			LogN:              paramsN2Lit.LogN - 1,
			Q:                 paramsN2Lit.Q[:endLevel+1],
			P:                 paramsN2Lit.P,
			LogPlaintextScale: paramsN2Lit.LogPlaintextScale,
			RingType:          ring.ConjugateInvariant,
		}

		paramsN1, err := ckks.NewParametersFromLiteral(parmasN1Lit)
		require.Nil(t, err)

		t.Logf("ParamsN2: LogN=%d/LogSlots=%d/LogQP=%f", paramsN2.LogN(), paramsN2.PlaintextLogSlots(), paramsN2.LogQP())
		t.Logf("ParamsN1: LogN=%d/LogSlots=%d/LogQP=%f", paramsN1.LogN(), paramsN1.PlaintextLogSlots(), paramsN1.LogQP())

		skN2 := ckks.NewKeyGenerator(paramsN2).GenSecretKeyNew()
		skN1 := ckks.NewKeyGenerator(paramsN1).GenSecretKeyNew()

		t.Log("Generating Bootstrapping Keys")
		btpKeys, err := GenBootstrappingKeys(&paramsN1, paramsN2, btpParamsN2, skN1, *skN2)
		require.Nil(t, err)

		bootstrapperInterface, err := NewBootstrapper(&paramsN1, paramsN2, btpParamsN2, btpKeys)
		require.Nil(t, err)

		bootstrapper := bootstrapperInterface.(*Bootstrapper)

		ecdN1 := ckks.NewEncoder(paramsN1)
		encN1 := ckks.NewEncryptor(paramsN1, skN1)
		decN1 := ckks.NewDecryptor(paramsN1, skN1)

		values := make([]float64, paramsN1.PlaintextSlots())
		for i := range values {
			values[i] = sampling.RandFloat64(-1, 1)
		}

		values[0] = 0.9238795325112867
		values[1] = 0.9238795325112867
		if paramsN1.PlaintextSlots() > 2 {
			values[2] = 0.9238795325112867
			values[3] = 0.9238795325112867
		}

		t.Run("ConjugateToStandard", func(t *testing.T) {
			ptN1 := ckks.NewPlaintext(paramsN1, 0)
			ecdN1.Encode(values, ptN1)
			ctN1Q0 := encN1.EncryptNew(ptN1)
			ctN2Q0 := ckks.NewCiphertext(paramsN2, ctN1Q0.Degree(), ctN1Q0.Level())
			ctN2Q0.PlaintextScale = ctN1Q0.PlaintextScale
			bootstrapper.ringConjugateToRingStandard(bootstrapper.bootstrapper.Evaluator.Evaluator, ctN1Q0, ctN2Q0)
			verifyTestVectorsBootstrapping(paramsN1, ckks.NewEncoder(paramsN2), ckks.NewDecryptor(paramsN2, skN2), values, ctN2Q0, t)
		})

		t.Run("StandardToConjugate", func(t *testing.T) {
			ptN2 := ckks.NewPlaintext(paramsN2, paramsN1.MaxLevel())
			ecdN2 := ckks.NewEncoder(paramsN2)
			encN2 := ckks.NewEncryptor(paramsN2, skN2)
			ecdN2.Encode(values, ptN2)
			ctN2QL := encN2.EncryptNew(ptN2)
			ctN1QL := ckks.NewCiphertext(paramsN1, ctN2QL.Degree(), ctN2QL.Level())
			ctN1QL.PlaintextScale = ctN2QL.PlaintextScale
			bootstrapper.ringStandardToConjugate(bootstrapper.bootstrapper.Evaluator.Evaluator, ctN2QL, ctN1QL)
			verifyTestVectorsBootstrapping(paramsN1, ecdN1, decN1, values, ctN1QL, t)
		})
		t.Run("ConjugateInvariant->Standard->Bootstrapping->Standard->ConjugateInvariant", func(t *testing.T) {

			plaintext := ckks.NewPlaintext(paramsN1, 0)
			ecdN1.Encode(values, plaintext)

			ctLeftN1Q0 := encN1.EncryptNew(plaintext)
			ctRightN1Q0 := encN1.EncryptNew(plaintext)

			// Checks that the input ciphertext is at the level 0
			require.True(t, ctLeftN1Q0.Level() == 0)
			require.True(t, ctRightN1Q0.Level() == 0)

			// Bootstrapps the ciphertext
			ctLeftN1QL, ctRightN1QL := bootstrapper.refreshConjugateInvariant(ctLeftN1Q0, ctRightN1Q0)

			// Checks that the output ciphertext is at the max level of paramsN1
			require.True(t, ctLeftN1QL.Level() == paramsN1.MaxLevel())
			require.True(t, ctLeftN1QL.PlaintextScale.Equal(paramsN1.PlaintextScale()))

			verifyTestVectorsBootstrapping(paramsN1, ecdN1, decN1, values, ctLeftN1QL, t)

			require.True(t, ctRightN1QL.Level() == paramsN1.MaxLevel())
			require.True(t, ctRightN1QL.PlaintextScale.Equal(paramsN1.PlaintextScale()))
			verifyTestVectorsBootstrapping(paramsN1, ecdN1, decN1, values, ctRightN1QL, t)
		})
	})
}

func verifyTestVectorsBootstrapping(params ckks.Parameters, encoder *ckks.Encoder, decryptor *rlwe.Decryptor, valuesWant, element interface{}, t *testing.T) {
	precStats := ckks.GetPrecisionStats(params, encoder, decryptor, valuesWant, element, nil, false)
	if *printPrecisionStats {
		t.Log(precStats.String())
	}

	rf64, _ := precStats.MeanPrecision.Real.Float64()
	if64, _ := precStats.MeanPrecision.Imag.Float64()

	minPrec := math.Log2(params.PlaintextScale().Float64()) - float64(params.LogN()+2)
	if minPrec < 0 {
		minPrec = 0
	}

	require.GreaterOrEqual(t, rf64, minPrec)
	require.GreaterOrEqual(t, if64, minPrec)
}
