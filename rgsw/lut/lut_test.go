package lut

import (
	"fmt"
	"math"
	"runtime"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

func testString(params rlwe.Parameters, opname string) string {
	return fmt.Sprintf("%slogN=%d/logQ=%d/logP=%d/#Qi=%d/#Pi=%d",
		opname,
		params.LogN(),
		params.LogQ(),
		params.LogP(),
		params.QCount(),
		params.PCount())
}

// TestLUT tests the LUT evaluation.
func TestLUT(t *testing.T) {
	for _, testSet := range []func(t *testing.T){
		testLUT,
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

var DefaultNTTFlag = true

func testLUT(t *testing.T) {
	var err error

	// RLWE parameters of the LUT
	// N=1024, Q=0x7fff801 -> 2^131
	paramsLUT, err := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		LogN:           10,
		Q:              []uint64{0x7fff801},
		Pow2Base:       6,
		DefaultNTTFlag: DefaultNTTFlag,
	})

	assert.Nil(t, err)

	// RLWE parameters of the samples
	// N=512, Q=0x3001 -> 2^135
	paramsLWE, err := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		LogN:           9,
		Q:              []uint64{0x3001},
		DefaultNTTFlag: DefaultNTTFlag,
	})

	assert.Nil(t, err)

	t.Run(testString(paramsLUT, "LUT/"), func(t *testing.T) {

		// Scale of the RLWE samples
		scaleLWE := float64(paramsLWE.Q()[0]) / 4.0

		// Scale of the test poly
		scaleLUT := float64(paramsLUT.Q()[0]) / 4.0

		// Number of values samples stored in the RLWE sample
		slots := 16

		// Test poly
		LUTPoly := InitLUT(sign, rlwe.NewScale(scaleLUT), paramsLUT.RingQ(), -1, 1)

		// Index map of which test poly to evaluate on which slot
		lutPolyMap := make(map[int]*ring.Poly)
		for i := 0; i < slots; i++ {
			lutPolyMap[i] = LUTPoly
		}

		// RLWE secret for the samples
		skLWE := rlwe.NewKeyGenerator(paramsLWE).GenSecretKey()

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

		// Evaluator for the LUT evaluation
		eval := NewEvaluator(paramsLUT, paramsLWE, nil)

		// Secret of the RGSW ciphertexts encrypting the bits of skLWE
		skLUT := rlwe.NewKeyGenerator(paramsLUT).GenSecretKey()

		// Collection of RGSW ciphertexts encrypting the bits of skLWE under skLUT
		LUTKEY := GenEvaluationKey(paramsLUT, skLUT, paramsLWE, skLWE)

		// Evaluation of LUT(ctLWE)
		// Returns one RLWE sample per slot in ctLWE
		ctsLUT := eval.Evaluate(ctLWE, lutPolyMap, LUTKEY)

		// Decrypts, decodes and compares
		q := paramsLUT.Q()[0]
		qHalf := q >> 1
		decryptorLUT := rlwe.NewDecryptor(paramsLUT, skLUT)
		ptLUT := rlwe.NewPlaintext(paramsLUT, paramsLUT.MaxLevel())
		for i := 0; i < slots; i++ {

			decryptorLUT.Decrypt(ctsLUT[i], ptLUT)

			if ptLUT.IsNTT {
				paramsLUT.RingQ().InvNTT(ptLUT.Value, ptLUT.Value)
			}

			c := ptLUT.Value.Coeffs[0][0]

			var a float64
			if c >= qHalf {
				a = -float64(q-c) / scaleLUT
			} else {
				a = float64(c) / scaleLUT
			}

			if values[i] != 0 {
				//fmt.Printf("%7.4f - %7.4f - %7.4f\n", math.Round(a*32)/32, math.Round(a*8)/8, values[i])
				assert.Equal(t, sign(values[i]), math.Round(a*8)/8)
			}
		}
	})
}
