package bootstrapping

import (
	"flag"
	"fmt"
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
	"github.com/stretchr/testify/assert"
	"runtime"
	"sync"
	"testing"
)

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters + secure bootstrapping). Overrides -short and requires -timeout=0.")
var testBootstrapping = flag.Bool("test-bootstrapping", false, "run the bootstrapping tests (memory intensive)")

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

func TestBootstrapParametersMarshalling(t *testing.T) {
	bootstrapParams := DefaultParameters[0]
	data, err := bootstrapParams.MarshalBinary()
	assert.Nil(t, err)

	bootstrapParamsNew := new(Parameters)
	if err := bootstrapParamsNew.UnmarshalBinary(data); err != nil {
		assert.Nil(t, err)
	}
	assert.Equal(t, bootstrapParams, *bootstrapParamsNew)
}

func TestBootstrap(t *testing.T) {

	if runtime.GOARCH == "wasm" {
		t.Skip("skipping bootstrapping tests for GOARCH=wasm")
	}

	if !*testBootstrapping {
		t.Skip("skipping bootstrapping tests (add -test-bootstrapping to run the bootstrapping tests)")
	}

	paramSet := 0

	ckksParams := DefaultCKKSParameters[paramSet]
	bootstrapParams := DefaultParameters[paramSet]

	// Insecure params for fast testing only
	if !*flagLongTest {
		ckksParams.LogN = 13
		ckksParams.LogSlots = 12
	}

	params, err := ckks.NewParametersFromLiteral(ckksParams)
	if err != nil {
		panic(err)
	}

	for _, testSet := range []func(params ckks.Parameters, btpParams Parameters, t *testing.T){
		testbootstrap,
	} {
		testSet(params, bootstrapParams, t)
		runtime.GC()
	}

	ckksParams.LogSlots--

	params, err = ckks.NewParametersFromLiteral(ckksParams)
	if err != nil {
		panic(err)
	}

	for _, testSet := range []func(params ckks.Parameters, btpParams Parameters, t *testing.T){
		testbootstrap,
	} {
		testSet(params, bootstrapParams, t)
		runtime.GC()
	}
}

func testbootstrap(params ckks.Parameters, btpParams Parameters, t *testing.T) {

	t.Run(ParamsToString(params, "Bootstrapping/FullCircuit/"), func(t *testing.T) {

		kgen := ckks.NewKeyGenerator(params)
		sk := kgen.GenSecretKeySparse(btpParams.H)
		rlk := kgen.GenRelinearizationKey(sk, 2)
		encoder := ckks.NewEncoder(params)
		encryptor := ckks.NewEncryptor(params, sk)
		decryptor := ckks.NewDecryptor(params, sk)

		rotations := btpParams.RotationsForBootstrapping(params.LogN(), params.LogSlots())
		rotkeys := kgen.GenRotationKeysForRotations(rotations, true, sk)

		btp, err := NewBootstrapper(params, btpParams, rlwe.EvaluationKey{Rlk: rlk, Rtks: rotkeys})
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

		plaintext := ckks.NewPlaintext(params, 0, params.Scale())
		encoder.Encode(plaintext, values, params.LogSlots())

		ciphertexts := make([]*ckks.Ciphertext, 2)
		bootstrappers := make([]*Bootstrapper, 2)
		for i := range ciphertexts {
			ciphertexts[i] = encryptor.EncryptNew(plaintext)
			if i == 0 {
				bootstrappers[i] = btp
			} else {
				bootstrappers[i] = bootstrappers[0].ShallowCopy()
			}
		}

		var wg sync.WaitGroup
		wg.Add(2)
		for i := range ciphertexts {
			go func(index int) {
				ciphertexts[index] = bootstrappers[index].Bootstrapp(ciphertexts[index])
				//btp.SetScale(ciphertexts[index], params.Scale())
				wg.Done()
			}(i)
		}
		wg.Wait()

		for i := range ciphertexts {
			verifyTestVectors(params, encoder, decryptor, values, ciphertexts[i], params.LogSlots(), 0, t)
		}
	})
}

func verifyTestVectors(params ckks.Parameters, encoder ckks.Encoder, decryptor ckks.Decryptor, valuesWant []complex128, element interface{}, logSlots int, bound float64, t *testing.T) {
	precStats := ckks.GetPrecisionStats(params, encoder, decryptor, valuesWant, element, logSlots, bound)
	t.Log(precStats.String())
}
