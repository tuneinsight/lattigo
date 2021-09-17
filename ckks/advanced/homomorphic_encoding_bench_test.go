package advanced

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
	"math"
	"runtime"
	"testing"
)

func BenchmarkHomomorphicEncoding(b *testing.B) {
	var err error

	ParametersLiteral := ckks.ParametersLiteral{
		LogN:     14,
		LogSlots: 13,
		Scale:    1 << 45,
		Sigma:    rlwe.DefaultSigma,
		Q: []uint64{
			0x10000000006e0001, // 60 Q0
			0x10000140001,      // 40
			0xffffe80001,       // 40
			0xffffc40001,       // 40
			0x100003e0001,      // 40
			0xffffb20001,       // 40
			0x10000500001,      // 40
			0xffff940001,       // 40
			0xffff8a0001,       // 40
			0xffff820001,       // 40
			0x7fffe60001,       // 39 StC
			0x7fffe40001,       // 39 StC
			0x7fffe00001,       // 39 StC
			0xfffffffff840001,  // 60 Sine (double angle)
			0x1000000000860001, // 60 Sine (double angle)
			0xfffffffff6a0001,  // 60 Sine
			0x1000000000980001, // 60 Sine
			0xfffffffff5a0001,  // 60 Sine
			0x1000000000b00001, // 60 Sine
			0x1000000000ce0001, // 60 Sine
			0xfffffffff2a0001,  // 60 Sine
			0x100000000060001,  // 58 CtS
			0xfffffffff00001,   // 58 CtS
			0xffffffffd80001,   // 58 CtS
			0x1000000002a0001,  // 58 CtS
		},
		P: []uint64{
			0x1fffffffffe00001, // Pi 61
			0x1fffffffffc80001, // Pi 61
			0x1fffffffffb40001, // Pi 61
			0x1fffffffff500001, // Pi 61
			0x1fffffffff420001, // Pi 61
		},
	}

	var params ckks.Parameters
	if params, err = ckks.NewParametersFromLiteral(ParametersLiteral); err != nil {
		panic(err)
	}

	for _, testSet := range []func(params ckks.Parameters, b *testing.B){
		benchCoeffsToSlots,
		benchSlotsToCoeffs,
	} {
		testSet(params, b)
		runtime.GC()
	}

	ParametersLiteral.LogSlots--
	if params, err = ckks.NewParametersFromLiteral(ParametersLiteral); err != nil {
		panic(err)
	}

	for _, testSet := range []func(params ckks.Parameters, b *testing.B){
		benchCoeffsToSlots,
		benchSlotsToCoeffs,
	} {
		testSet(params, b)
		runtime.GC()
	}
}

func benchCoeffsToSlots(params ckks.Parameters, b *testing.B) {

	packing := "FullPacking"
	if params.LogSlots() < params.LogN()-1 {
		packing = "SparsePacking"
	}

	b.Run("CoeffsToSlots/"+packing, func(b *testing.B) {

		CoeffsToSlotsParametersLiteral := EncodingMatrixLiteral{
			LinearTransformType: CoeffsToSlots,
			LevelStart:          params.MaxLevel(),
			BSGSRatio:           16.0,
			BitReversed:         false,
			ScalingFactor: [][]float64{
				{params.QiFloat64(params.MaxLevel() - 3)},
				{params.QiFloat64(params.MaxLevel() - 2)},
				{params.QiFloat64(params.MaxLevel() - 1)},
				{params.QiFloat64(params.MaxLevel() - 0)},
			},
		}

		kgen := ckks.NewKeyGenerator(params)
		sk := kgen.GenSecretKey()
		encoder := ckks.NewEncoder(params)

		n := math.Pow(1.0/float64(2*params.Slots()), 1.0/float64(CoeffsToSlotsParametersLiteral.Depth(true)))

		// Generates the encoding matrices
		CoeffsToSlotMatrices := NewHomomorphicEncodingMatrixFromLiteral(CoeffsToSlotsParametersLiteral, encoder, params.LogN(), params.LogSlots(), complex(n, 0))

		// Gets the rotations indexes for CoeffsToSlots
		rotations := CoeffsToSlotsParametersLiteral.Rotations(params.LogN(), params.LogSlots())

		// Generates the rotation keys
		rotKey := kgen.GenRotationKeysForRotations(rotations, true, sk)

		// Creates an evaluator with the rotation keys
		eval := NewEvaluator(params, rlwe.EvaluationKey{Rlk: nil, Rtks: rotKey})

		prng, _ := utils.NewPRNG()
		ciphertext := ckks.NewCiphertextRandom(prng, params, 1, params.MaxLevel(), params.Scale())

		b.ResetTimer()

		for i := 0; i < b.N; i++ {
			eval.CoeffsToSlotsNew(ciphertext, CoeffsToSlotMatrices)
		}
	})
}

func benchSlotsToCoeffs(params ckks.Parameters, b *testing.B) {

	packing := "FullPacking"
	if params.LogSlots() < params.LogN()-1 {
		packing = "SparsePacking"
	}

	b.Run("SlotsToCoeffs/"+packing, func(b *testing.B) {

		SlotsToCoeffsParametersLiteral := EncodingMatrixLiteral{
			LinearTransformType: SlotsToCoeffs,
			LevelStart:          params.MaxLevel(),
			BSGSRatio:           16.0,
			BitReversed:         false,
			ScalingFactor: [][]float64{
				{params.QiFloat64(params.MaxLevel() - 3)},
				{params.QiFloat64(params.MaxLevel() - 2)},
				{params.QiFloat64(params.MaxLevel() - 1)},
				{params.QiFloat64(params.MaxLevel() - 0)},
			},
		}

		kgen := ckks.NewKeyGenerator(params)
		sk := kgen.GenSecretKey()
		encoder := ckks.NewEncoder(params)

		// Generates the encoding matrices
		SlotsToCoeffsMatrix := NewHomomorphicEncodingMatrixFromLiteral(SlotsToCoeffsParametersLiteral, encoder, params.LogN(), params.LogSlots(), 1.0)

		// Gets the rotations indexes for SlotsToCoeffs
		rotations := SlotsToCoeffsParametersLiteral.Rotations(params.LogN(), params.LogSlots())

		// Generates the rotation keys
		rotKey := kgen.GenRotationKeysForRotations(rotations, true, sk)

		// Creates an evaluator with the rotation keys
		eval := NewEvaluator(params, rlwe.EvaluationKey{Rlk: nil, Rtks: rotKey})

		prng, _ := utils.NewPRNG()
		ct0 := ckks.NewCiphertextRandom(prng, params, 1, params.MaxLevel(), params.Scale())

		var ct1 *ckks.Ciphertext
		if params.LogSlots() == params.LogN()-1 {
			ct1 = ckks.NewCiphertextRandom(prng, params, 1, params.MaxLevel(), params.Scale())
		}

		b.ResetTimer()

		for i := 0; i < b.N; i++ {
			eval.SlotsToCoeffsNew(ct0, ct1, SlotsToCoeffsMatrix)
		}
	})
}
