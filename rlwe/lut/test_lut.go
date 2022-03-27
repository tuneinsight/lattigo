package lut

import (
	"fmt"
	"github.com/stretchr/testify/assert"
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"math"
	"runtime"
	"testing"
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

// m0 : [0, 1/4]
// m1 : [0, 1/4]
// | 0 0 -> 1 (0/8) -> 2/8
// | 0 1 -> 1 (2/8) -> 2/8
// | 1 0 -> 1 (2/8) -> 2/8
// | 1 1 -> 0 (4/8) -> 0/8
func nandGate(x float64) float64 {
	if x > -1/8.0 && x < 3/8.0 {
		return 2 / 8.0
	}

	return 0
}

func testLUT(t *testing.T) {
	var err error

	// N=1024, Q=0x7fff801 -> 2^131
	paramsLUT, err := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		LogN:     8,
		Q:        []uint64{0x7fff801},
		P:        []uint64{},
		Sigma:    rlwe.DefaultSigma,
		LogBase2: 7,
	})

	assert.Nil(t, err)

	// N=512, Q=0x3001 -> 2^135
	paramsLWE, err := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		LogN:  7,
		Q:     []uint64{0x3001},
		P:     []uint64{},
		Sigma: rlwe.DefaultSigma,
	})

	assert.Nil(t, err)

	t.Run(testString(paramsLUT, "LUT/"), func(t *testing.T) {

		scaleLWE := float64(paramsLWE.Q()[0]) / 4.0
		scaleLUT := float64(paramsLUT.Q()[0]) / 4.0

		slots := 16

		LUTPoly := InitLUT(nandGate, scaleLUT, paramsLUT.RingQ(), -1, 1)

		lutPolyMap := make(map[int]*ring.Poly)
		for i := 0; i < slots; i++ {
			lutPolyMap[i] = LUTPoly
		}

		skLWE := rlwe.NewKeyGenerator(paramsLWE).GenSecretKey()
		encryptorLWE := rlwe.NewEncryptor(paramsLWE, skLWE)

		values := make([]float64, slots)
		for i := 0; i < slots; i++ {
			values[i] = -1 + float64(2*i)/float64(slots)
		}

		ptLWE := rlwe.NewPlaintext(paramsLWE, paramsLWE.MaxLevel())
		for i := range values {
			if values[i] < 0 {
				ptLWE.Value.Coeffs[0][i] = paramsLWE.Q()[0] - uint64(-values[i]*scaleLWE)
			} else {
				ptLWE.Value.Coeffs[0][i] = uint64(values[i] * scaleLWE)
			}
		}
		ctLWE := rlwe.NewCiphertextNTT(paramsLWE, 1, paramsLWE.MaxLevel())
		encryptorLWE.Encrypt(ptLWE, ctLWE)

		eval := NewEvaluator(paramsLUT, paramsLWE, nil)

		skLUT := rlwe.NewKeyGenerator(paramsLUT).GenSecretKey()
		LUTKEY := eval.GenLUTKey(skLUT, skLWE)

		ctsLUT := eval.Evaluate(ctLWE, lutPolyMap, LUTKEY)

		q := paramsLUT.Q()[0]
		qHalf := q >> 1
		decryptorLUT := rlwe.NewDecryptor(paramsLUT, skLUT)
		ptLUT := rlwe.NewPlaintext(paramsLUT, paramsLUT.MaxLevel())
		for i := 0; i < slots; i++ {

			decryptorLUT.Decrypt(ctsLUT[i], ptLUT)

			c := ptLUT.Value.Coeffs[0][i]

			var a float64
			if c >= qHalf {
				a = -float64(q-c) / scaleLUT
			} else {
				a = float64(c) / scaleLUT
			}

			//fmt.Printf("%7.4f - %7.4f - %7.4f\n", math.Round(a*32)/32, math.Round(a*8)/8, values[i])
			assert.Equal(t, nandGate(values[i]), math.Round(a*8)/8)
		}
	})
}
