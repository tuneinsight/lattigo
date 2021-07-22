package bootstrapping

import (
	"flag"
	"fmt"
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/utils"
	"runtime"
	"testing"
)

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters + secure bootstrapping). Overrides -short and requires -timeout=0.")

func ParamsToString(params ckks.Parameters, opname string) string {
	return fmt.Sprintf("%slogN=%d/LogSlots=%d/logQP=%d/levels=%d/a=%d/b=%d",
		opname,
		params.LogN(),
		params.LogSlots(),
		params.LogQP(),
		params.MaxLevel()+1,
		params.PCount(),
		params.Beta())
}

func TestBootstrap(t *testing.T) {

	if runtime.GOARCH == "wasm" {
		t.Skip("skipping bootstrapping tests for GOARCH=wasm")
	}

	paramSet := 4

	bootstrapParams := DefaultBootstrapParams[paramSet : paramSet+1]

	for paramSet := range bootstrapParams {

		btpParams := bootstrapParams[paramSet]

		// Insecure params for fast testing only
		if !*flagLongTest {
			btpParams.LogN = 14
			btpParams.LogSlots = 13
		}

		// Tests homomorphic modular reduction encoding and bootstrapping on sparse slots
		params, err := btpParams.Params()
		if err != nil {
			panic(err)
		}

		for _, testSet := range []func(params ckks.Parameters, btpParams *BootstrappingParameters, t *testing.T){
			testbootstrap,
		} {
			testSet(params, btpParams, t)
			runtime.GC()
		}

		if !*flagLongTest {
			btpParams.LogSlots = 12
		}

		// Tests homomorphic encoding and bootstrapping on full slots
		params, err = btpParams.Params()
		if err != nil {
			panic(err)
		}

		for _, testSet := range []func(params ckks.Parameters, btpParams *BootstrappingParameters, t *testing.T){
			testbootstrap,
		} {
			testSet(params, btpParams, t)
			runtime.GC()
		}
	}
}

func testbootstrap(params ckks.Parameters, btpParams *BootstrappingParameters, t *testing.T) {

	t.Run(ParamsToString(params, "Bootstrapping/FullCircuit/"), func(t *testing.T) {

		kgen := ckks.NewKeyGenerator(params)
		sk := kgen.GenSecretKeySparse(btpParams.H)
		rlk := kgen.GenRelinearizationKey(sk, 2)
		encoder := ckks.NewEncoder(params)
		encryptor := ckks.NewEncryptor(params, sk)
		decryptor := ckks.NewDecryptor(params, sk)

		rotations := btpParams.RotationsForBootstrapping(params.LogSlots())
		rotkeys := kgen.GenRotationKeysForRotations(rotations, true, sk)
		btpKey := BootstrappingKey{rlk, rotkeys}

		btp, err := NewBootstrapper(params, *btpParams, btpKey)
		if err != nil {
			panic(err)
		}

		values := make([]complex128, 1<<params.LogSlots())
		for i := range values {
			values[i] = utils.RandComplex128(-1, 1)
		}

		values[0] = complex(0.9238795325112867, 0.3826834323650898)
		values[1] = complex(0.9238795325112867, 0.3826834323650898)
		if 1<<params.LogSlots() > 2 {
			values[2] = complex(0.9238795325112867, 0.3826834323650898)
			values[3] = complex(0.9238795325112867, 0.3826834323650898)
		}

		plaintext := ckks.NewPlaintext(params, params.MaxLevel(), params.Scale())
		encoder.Encode(plaintext, values, params.LogSlots())

		ciphertext := encryptor.EncryptNew(plaintext)

		for ciphertext.Level() != 0 {
			btp.DropLevel(ciphertext, 1)
		}

		for i := 0; i < 1; i++ {
			ciphertext = btp.Bootstrapp(ciphertext)
			//btp.SetScale(ciphertext, params.Scale())
			verifyTestVectors(params, encoder, decryptor, values, ciphertext, params.LogSlots(), 0, t)
		}
	})
}

func verifyTestVectors(params ckks.Parameters, encoder ckks.Encoder, decryptor ckks.Decryptor, valuesWant []complex128, element interface{}, logSlots int, bound float64, t *testing.T) {
	precStats := ckks.GetPrecisionStats(params, encoder, decryptor, valuesWant, element, logSlots, bound)
	t.Log(precStats.String())
}
