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
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters + secure bootstrapping). Overrides -short and requires -timeout=0.")
var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")

func TestBootstrapping(t *testing.T) {

	paramSet := bootstrapping.DefaultParametersSparse[0]
	paramSet.SchemeParams.LogQ = paramSet.SchemeParams.LogQ[:utils.Min(2, len(paramSet.SchemeParams.LogQ))]

	paramsN2Lit, btpParamsN2, err := bootstrapping.NewParametersFromLiteral(paramSet.SchemeParams, paramSet.BootstrappingParams)
	require.Nil(t, err)

	// Insecure params for fast testing only
	if !*flagLongTest {
		paramsN2Lit.LogN = 13
		btpParamsN2.SlotsToCoeffsParameters.LogSlots = paramsN2Lit.LogN - 1
		btpParamsN2.CoeffsToSlotsParameters.LogSlots = paramsN2Lit.LogN - 1

		// Corrects the message ratio to take into account the smaller number of slots and keep the same precision
		btpParamsN2.EvalModParameters.LogMessageRatio += paramSet.SchemeParams.LogN - paramsN2Lit.LogN + 1

	}

	endLevel := len(paramSet.SchemeParams.LogQ) - 1

	require.True(t, endLevel == len(paramsN2Lit.Q)-1-btpParamsN2.Depth()) // Checks the depth of the bootstrapping

	// Check that the bootstrapper complies to the rlwe.Bootstrapper interface
	var _ rlwe.Bootstrapper = (*Bootstrapper)(nil)

	t.Run("BootstrapingWithoutRingDegreeSwitch", func(t *testing.T) {

		paramsN2, err := ckks.NewParametersFromLiteral(paramsN2Lit)
		require.Nil(t, err)

		t.Logf("ParamsN2: LogN=%d/LogSlots=%d/LogQP=%f", paramsN2.LogN(), paramsN2.PlaintextLogSlots(), paramsN2.LogQP())

		skN2 := ckks.NewKeyGenerator(paramsN2).GenSecretKeyNew()

		t.Log("Generating Bootstrapping Keys")
		btpKeys, err := GenBootstrappingKeys(paramsN2, paramsN2, btpParamsN2, skN2, skN2)
		require.Nil(t, err)

		bootstrapperInterface, err := NewBootstrapper(paramsN2, paramsN2, btpParamsN2, btpKeys)
		require.Nil(t, err)

		bootstrapper := bootstrapperInterface.(*Bootstrapper)

		ecdN2 := ckks.NewEncoder(paramsN2)
		encN2, err := ckks.NewEncryptor(paramsN2, skN2)
		require.NoError(t, err)
		decN2, err := ckks.NewDecryptor(paramsN2, skN2)
		require.NoError(t, err)

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

			ctN2Q0, err := encN2.EncryptNew(plaintext)
			require.NoError(t, err)

			// Checks that the input ciphertext is at the level 0
			require.True(t, ctN2Q0.Level() == 0)

			// Bootstrapps the ciphertext
			ctN2QL, err := bootstrapper.Bootstrap(ctN2Q0)

			if err != nil {
				t.Fatal(err)
			}

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
		btpKeys, err := GenBootstrappingKeys(paramsN1, paramsN2, btpParamsN2, skN1, skN2)
		require.Nil(t, err)

		bootstrapperInterface, err := NewBootstrapper(paramsN1, paramsN2, btpParamsN2, btpKeys)
		require.Nil(t, err)

		bootstrapper := bootstrapperInterface.(*Bootstrapper)

		ecdN1 := ckks.NewEncoder(paramsN1)
		encN1, err := ckks.NewEncryptor(paramsN1, skN1)
		require.NoError(t, err)
		decN1, err := ckks.NewDecryptor(paramsN1, skN1)
		require.NoError(t, err)

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

		t.Run("N1ToN2->Bootstrapping->N2ToN1", func(t *testing.T) {

			plaintext := ckks.NewPlaintext(paramsN1, 0)
			ecdN1.Encode(values, plaintext)

			ctN1Q0, err := encN1.EncryptNew(plaintext)
			require.NoError(t, err)

			// Checks that the input ciphertext is at the level 0
			require.True(t, ctN1Q0.Level() == 0)

			// Bootstrapps the ciphertext
			ctN1QL, err := bootstrapper.Bootstrap(ctN1Q0)

			if err != nil {
				t.Fatal(err)
			}

			// Checks that the output ciphertext is at the max level of paramsN1
			require.True(t, ctN1QL.Level() == paramsN1.MaxLevel())
			require.True(t, ctN1QL.PlaintextScale.Equal(paramsN1.PlaintextScale()))

			verifyTestVectorsBootstrapping(paramsN1, ecdN1, decN1, values, ctN1QL, t)

		})
	})

	t.Run("BootstrappingPackedWithRingDegreeSwitch", func(t *testing.T) {
		paramsN2, err := ckks.NewParametersFromLiteral(paramsN2Lit)
		require.Nil(t, err)

		parmasN1Lit := ckks.ParametersLiteral{
			LogN:              paramsN2Lit.LogN - 5,
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
		btpKeys, err := GenBootstrappingKeys(paramsN1, paramsN2, btpParamsN2, skN1, skN2)
		require.Nil(t, err)

		bootstrapperInterface, err := NewBootstrapper(paramsN1, paramsN2, btpParamsN2, btpKeys)
		require.Nil(t, err)

		bootstrapper := bootstrapperInterface.(*Bootstrapper)

		bootstrapper.skN1 = skN2
		bootstrapper.skN2 = skN2

		ecdN1 := ckks.NewEncoder(paramsN1)
		encN1, err := ckks.NewEncryptor(paramsN1, skN1)
		require.NoError(t, err)
		decN1, err := ckks.NewDecryptor(paramsN1, skN1)
		require.NoError(t, err)

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

		ptN1 := ckks.NewPlaintext(paramsN1, 0)

		cts := make([]rlwe.Ciphertext, 17)
		for i := range cts {

			require.NoError(t, ecdN1.Encode(utils.RotateSlice(values, i), ptN1))

			ct, err := encN1.EncryptNew(ptN1)
			require.NoError(t, err)

			cts[i] = *ct
		}

		if cts, err = bootstrapper.BootstrapMany(cts); err != nil {
			t.Fatal(err)
		}

		for i, ct := range cts {
			// Checks that the output ciphertext is at the max level of paramsN1
			require.True(t, ct.Level() == paramsN1.MaxLevel())
			require.True(t, ct.PlaintextScale.Equal(paramsN1.PlaintextScale()))

			want := utils.RotateSlice(values, i)

			verifyTestVectorsBootstrapping(paramsN1, ecdN1, decN1, want, &ct, t)
		}
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
		btpKeys, err := GenBootstrappingKeys(paramsN1, paramsN2, btpParamsN2, skN1, skN2)
		require.Nil(t, err)

		bootstrapperInterface, err := NewBootstrapper(paramsN1, paramsN2, btpParamsN2, btpKeys)
		require.Nil(t, err)

		bootstrapper := bootstrapperInterface.(*Bootstrapper)

		ecdN1 := ckks.NewEncoder(paramsN1)
		encN1, err := ckks.NewEncryptor(paramsN1, skN1)
		require.NoError(t, err)
		decN1, err := ckks.NewDecryptor(paramsN1, skN1)
		require.NoError(t, err)

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

		t.Run("ConjugateInvariant->Standard->Bootstrapping->Standard->ConjugateInvariant", func(t *testing.T) {

			plaintext := ckks.NewPlaintext(paramsN1, 0)
			require.NoError(t, ecdN1.Encode(values, plaintext))

			ctLeftN1Q0, err := encN1.EncryptNew(plaintext)
			require.NoError(t, err)
			ctRightN1Q0, err := encN1.EncryptNew(plaintext)
			require.NoError(t, err)

			// Checks that the input ciphertext is at the level 0
			require.True(t, ctLeftN1Q0.Level() == 0)
			require.True(t, ctRightN1Q0.Level() == 0)

			// Bootstrapps the ciphertext
			ctLeftN1QL, ctRightN1QL, err := bootstrapper.refreshConjugateInvariant(ctLeftN1Q0, ctRightN1Q0)

			require.NoError(t, err)

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

	minPrec -= 10

	require.GreaterOrEqual(t, rf64, minPrec)
	require.GreaterOrEqual(t, if64, minPrec)
}
