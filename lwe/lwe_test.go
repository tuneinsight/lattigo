package lwe

import (
	"encoding/json"
	"flag"
	"fmt"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"math"
	"runtime"
	"testing"
)

var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")

var TestParams = []rlwe.ParametersLiteral{rlwe.TestPN12QP109, rlwe.TestPN13QP218}

func testString(params rlwe.Parameters, opname string) string {
	return fmt.Sprintf("%slogN=%d/logQ=%d/logP=%d/#Qi=%d/#Pi=%d",
		opname,
		params.LogN(),
		params.LogQ(),
		params.LogP(),
		params.QCount(),
		params.PCount())
}

func TestLWE(t *testing.T) {
	defaultParams := TestParams // the default test runs for ring degree N=2^12, 2^13, 2^14, 2^15
	if testing.Short() {
		defaultParams = TestParams[:1] // the short test suite runs for ring degree N=2^12, 2^13
	}

	if *flagParamString != "" {
		var jsonParams rlwe.ParametersLiteral
		json.Unmarshal([]byte(*flagParamString), &jsonParams)
		defaultParams = []rlwe.ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, defaultParam := range defaultParams[1:] {

		params, err := rlwe.NewParametersFromLiteral(defaultParam)
		if err != nil {
			panic(err)
		}

		for _, testSet := range []func(params rlwe.Parameters, t *testing.T){
			testRLWEToLWE,
			testLWEToRLWE,
			testManyRLWEToSingleRLWE,
		} {
			testSet(params, t)
			runtime.GC()
		}
	}
}

func testRLWEToLWE(params rlwe.Parameters, t *testing.T) {
	t.Run(testString(params, "RLWEToLWE/"), func(t *testing.T) {
		kgen := rlwe.NewKeyGenerator(params)
		sk := kgen.GenSecretKey()
		encryptor := rlwe.NewEncryptor(params, sk)
		pt := rlwe.NewPlaintext(params, params.MaxLevel())
		ct := rlwe.NewCiphertextNTT(params, 1, params.MaxLevel())
		encryptor.Encrypt(pt, ct)

		skInvNTT := params.RingQ().NewPoly()
		params.RingQ().InvNTT(sk.Value.Q, skInvNTT)

		slotIndex := make(map[int]bool)
		for i := 0; i < params.N(); i++ {
			slotIndex[i] = true
		}

		LWE := RLWEToLWE(ct, params.RingQ(), slotIndex)

		for i := 0; i < params.RingQ().N; i++ {
			if math.Abs(DecryptLWE(LWE[i], params.RingQ(), skInvNTT)) > 19 {
				t.Error()
			}
		}
	})
}

func testLWEToRLWE(params rlwe.Parameters, t *testing.T) {
	t.Run(testString(params, "LWEToRLWE/"), func(t *testing.T) {
		kgen := rlwe.NewKeyGenerator(params)
		sk := kgen.GenSecretKey()
		encryptor := rlwe.NewEncryptor(params, sk)
		decryptor := rlwe.NewDecryptor(params, sk)
		pt := rlwe.NewPlaintext(params, params.MaxLevel())
		ct := rlwe.NewCiphertextNTT(params, 1, params.MaxLevel())
		encryptor.Encrypt(pt, ct)

		skInvNTT := params.RingQ().NewPoly()
		params.RingQ().InvNTT(sk.Value.Q, skInvNTT)

		slotIndex := make(map[int]bool)
		for i := 0; i < params.N(); i++ {
			slotIndex[i] = true
		}

		ctLWE := RLWEToLWE(ct, params.RingQ(), slotIndex)

		DecryptLWE(ctLWE[0], params.RingQ(), skInvNTT)

		handler := NewHandler(params, params, nil)

		ctRLWE := handler.LWEToRLWE(ctLWE)

		for i := 0; i < len(ctRLWE); i++ {
			decryptor.Decrypt(ctRLWE[i], pt)

			for j := 0; j < pt.Level()+1; j++ {

				c := pt.Value.Coeffs[j][0]

				if c >= params.RingQ().Modulus[j]>>1 {
					c = params.RingQ().Modulus[j] - c
				}

				if c > 19 {
					t.Fatal(i, j, c)
				}
			}
		}
	})
}

func testManyRLWEToSingleRLWE(params rlwe.Parameters, t *testing.T) {
	t.Run(testString(params, "ManyToSingleRLWE/"), func(t *testing.T) {
		kgen := rlwe.NewKeyGenerator(params)
		sk := kgen.GenSecretKey()
		encryptor := rlwe.NewEncryptor(params, sk)
		decryptor := rlwe.NewDecryptor(params, sk)
		pt := rlwe.NewPlaintext(params, params.MaxLevel())

		for i := 0; i < pt.Level()+1; i++ {
			for j := 0; j < params.N(); j++ {
				pt.Value.Coeffs[i][j] = (1 << 30) + uint64(j)*(1<<20)
			}
		}

		// Rotation Keys
		rotations := []int{}
		for i := 1; i < params.N(); i <<= 1 {
			rotations = append(rotations, i)
		}

		rtks := kgen.GenRotationKeysForRotations(rotations, true, sk)

		handler := NewHandler(params, params, rtks)

		ct := rlwe.NewCiphertextNTT(params, 1, params.MaxLevel())
		encryptor.Encrypt(pt, ct)

		slotIndex := make(map[int]bool)
		for i := 2; i < 2+params.N()/2; i += 3 {
			slotIndex[i] = true
		}

		ctLWE := RLWEToLWE(ct, params.RingQ(), slotIndex)

		ctRLWE := handler.LWEToRLWE(ctLWE)

		handler.Sk = sk

		ct = handler.MergeRLWE(ctRLWE)

		decryptor.Decrypt(ct, pt)

		bound := uint64(params.N() * params.N())

		for i := 0; i < pt.Level()+1; i++ {

			Q := params.RingQ().Modulus[i]
			QHalf := Q >> 1

			for i, c := range pt.Value.Coeffs[i] {

				if c >= QHalf {
					c = Q - c
				}

				if _, ok := slotIndex[i]; !ok {
					if c > bound {
						t.Fatal(i, c)
					}
				}
			}
		}
	})
}
