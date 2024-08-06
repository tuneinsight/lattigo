package blindrot

import (
	"fmt"
	"math"
	"runtime"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils"
)

func testString(params rlwe.Parameters, opname string) string {
	return fmt.Sprintf("%slogN=%d/logQ=%f/logP=%f/#Qi=%d/#Pi=%d",
		opname,
		params.LogN(),
		params.LogQ(),
		params.LogP(),
		params.QCount(),
		params.PCount())
}

// TestBlindRotation tests the BlindRotation evaluation.
func TestBlindRotation(t *testing.T) {
	for _, testSet := range []func(t *testing.T){
		testBlindRotation,
	} {
		testSet(t)
		runtime.GC()
	}
}

// Function to evaluate
func sign(x float64) float64 {
	if x > 0 {
		return 1
	} else if x == 0 {
		return 0
	}

	return -1
}

var NTTFlag = true

func testBlindRotation(t *testing.T) {
	var err error

	// RLWE parameters of the BlindRotation
	// N=1024, Q=0x7fff801 -> 131 bit secure
	paramsBR, err := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		LogN:    10,
		Q:       []uint64{0x7fff801},
		NTTFlag: NTTFlag,
	})

	require.NoError(t, err)

	// RLWE parameters of the samples
	// N=512, Q=0x3001 -> 135 bit secure
	paramsLWE, err := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		LogN:    9,
		Q:       []uint64{0x3001},
		NTTFlag: NTTFlag,
	})

	require.NoError(t, err)

	evkParams := rlwe.EvaluationKeyParameters{BaseTwoDecomposition: utils.Pointy(7)}

	require.NoError(t, err)

	t.Run(testString(paramsBR, "BlindRotation/"), func(t *testing.T) {

		// Scale of the RLWE samples
		scaleLWE := float64(paramsLWE.Q()[0]) / 4.0

		// Scale of the test poly
		scaleBR := float64(paramsBR.Q()[0]) / 4.0

		// Number of values samples stored in the RLWE sample
		slots := 16

		// Test poly
		testPoly := InitTestPolynomial(sign, rlwe.NewScale(scaleBR), paramsBR.RingQ(), -1, 1)

		// Index map of which test poly to evaluate on which slot
		testPolyMap := make(map[int]*ring.Poly)
		for i := 0; i < slots; i++ {
			testPolyMap[i] = &testPoly
		}

		// RLWE secret for the samples
		skLWE := rlwe.NewKeyGenerator(paramsLWE).GenSecretKeyNew()

		// RLWE encryptor for the samples
		encryptorLWE := rlwe.NewEncryptor(paramsLWE, skLWE)

		// Values to encrypt in the RLWE sample
		values := make([]float64, slots)
		for i := 0; i < slots; i++ {
			values[i] = -1 + float64(2*i)/float64(slots)
		}

		// Encode multiples values in a single RLWE
		ptLWE := rlwe.NewPlaintext(paramsLWE, paramsLWE.MaxLevel())

		for i := range values {
			if values[i] < 0 {
				ptLWE.Value.Coeffs[0][i] = paramsLWE.Q()[0] - uint64(-values[i]*scaleLWE)
			} else {
				ptLWE.Value.Coeffs[0][i] = uint64(values[i] * scaleLWE)
			}
		}

		if ptLWE.IsNTT {
			paramsLWE.RingQ().NTT(ptLWE.Value, ptLWE.Value)
		}

		// Encrypt the multiples values in a single RLWE
		ctLWE := rlwe.NewCiphertext(paramsLWE, 1, paramsLWE.MaxLevel())
		encryptorLWE.Encrypt(ptLWE, ctLWE)

		// Evaluator for the Blind Rotation evaluation
		eval := NewEvaluator(paramsBR, paramsLWE)

		// Secret of the RGSW ciphertexts encrypting the bits of skLWE
		skBR := rlwe.NewKeyGenerator(paramsBR).GenSecretKeyNew()

		// Collection of RGSW ciphertexts encrypting the bits of skLWE under skBR
		BRK := GenEvaluationKeyNew(paramsBR, skBR, paramsLWE, skLWE, evkParams)

		// Evaluation of BlindRotation(ctLWE)
		// Returns one RLWE sample per slot in ctLWE
		ctsBR, err := eval.Evaluate(ctLWE, testPolyMap, BRK)
		require.NoError(t, err)

		// Decrypts, decodes and compares
		q := paramsBR.Q()[0]
		qHalf := q >> 1
		decryptorBR := rlwe.NewDecryptor(paramsBR, skBR)
		ptBR := rlwe.NewPlaintext(paramsBR, paramsBR.MaxLevel())
		for i := 0; i < slots; i++ {

			decryptorBR.Decrypt(ctsBR[i], ptBR)

			if ptBR.IsNTT {
				paramsBR.RingQ().INTT(ptBR.Value, ptBR.Value)
			}

			c := ptBR.Value.Coeffs[0][0]

			var a float64
			if c >= qHalf {
				a = -float64(q-c) / scaleBR
			} else {
				a = float64(c) / scaleBR
			}

			if values[i] != 0 {
				fmt.Printf("%7.4f - %7.4f - %7.4f\n", math.Round(a*32)/32, math.Round(a*8)/8, values[i])
				//require.Equal(t, sign(values[i]), math.Round(a*8)/8)
			}
		}
	})
}
