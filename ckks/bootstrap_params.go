package ckks

import (
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// BootstrappingParameters is a struct for the default bootstrapping parameters
type BootstrappingParameters struct {
	ParametersLiteral
	SlotsToCoeffsParameters EncodingMatricesParameters
	EvalModParameters
	CoeffsToSlotsParameters EncodingMatricesParameters
	H                       int // Hamming weight of the secret key
}

// Params generates a new set of Parameters from the BootstrappingParameters
func (b *BootstrappingParameters) Params() (p Parameters, err error) {
	if p, err = NewParametersFromLiteral(b.ParametersLiteral); err != nil {
		return Parameters{}, err
	}
	return
}

// RotationsForBootstrapping returns the list of rotations performed during the Bootstrapping operation.
func (b *BootstrappingParameters) RotationsForBootstrapping(logSlots int) (rotations []int) {

	// List of the rotation key values to needed for the bootstrapp
	rotations = []int{}

	slots := 1 << logSlots
	dslots := slots
	if logSlots < b.LogN-1 {
		dslots <<= 1
	}

	//SubSum rotation needed X -> Y^slots rotations
	for i := logSlots; i < b.LogN-1; i++ {
		if !utils.IsInSliceInt(1<<i, rotations) {
			rotations = append(rotations, 1<<i)
		}
	}

	indexCtS := computeBootstrappingDFTIndexMap(b.LogN, logSlots, b.CoeffsToSlotsParameters.Depth(false), CoeffsToSlots, b.CoeffsToSlotsParameters.BitReversed)

	// Coeffs to Slots rotations
	for _, pVec := range indexCtS {
		N1 := findbestbabygiantstepsplit(pVec, dslots, b.CoeffsToSlotsParameters.BSGSRatio)
		rotations = addMatrixRotToList(pVec, rotations, N1, slots, false)
	}

	indexStC := computeBootstrappingDFTIndexMap(b.LogN, logSlots, b.SlotsToCoeffsParameters.Depth(false), SlotsToCoeffs, b.SlotsToCoeffsParameters.BitReversed)

	// Slots to Coeffs rotations
	for i, pVec := range indexStC {
		N1 := findbestbabygiantstepsplit(pVec, dslots, b.SlotsToCoeffsParameters.BSGSRatio)
		rotations = addMatrixRotToList(pVec, rotations, N1, slots, logSlots < b.LogN-1 && i == 0)
	}

	return
}

// DefaultBootstrapParams are default bootstrapping params for the bootstrapping.
var DefaultBootstrapParams = []*BootstrappingParameters{

	// SET I
	// 1546
	{
		H: 192,
		ParametersLiteral: ParametersLiteral{
			LogN:     16,
			LogSlots: 15,
			Scale:    1 << 40,
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
		},
		SlotsToCoeffsParameters: EncodingMatricesParameters{
			LinearTransformType: SlotsToCoeffs,
			LevelStart:          12,
			BSGSRatio:           16.0,
			BitReversed:         false,
			ScalingFactor: [][]float64{
				{0x7fffe60001},
				{0x7fffe40001},
				{0x7fffe00001},
			},
		},
		EvalModParameters: EvalModParameters{
			LevelStart:    20,
			SineType:      Cos1,
			MessageRatio:  256.0,
			K:             25,
			SineDeg:       63,
			DoubleAngle:   2,
			ArcSineDeg:    0,
			ScalingFactor: 1 << 60,
		},
		CoeffsToSlotsParameters: EncodingMatricesParameters{
			LinearTransformType: CoeffsToSlots,
			LevelStart:          24,
			BSGSRatio:           16.0,
			BitReversed:         false,
			ScalingFactor: [][]float64{
				{0x100000000060001},
				{0xfffffffff00001},
				{0xffffffffd80001},
				{0x1000000002a0001},
			},
		},
	},

	// SET II
	// 1547
	{
		H: 192,
		ParametersLiteral: ParametersLiteral{
			LogN:     16,
			LogSlots: 15,
			Scale:    1 << 45,
			Sigma:    rlwe.DefaultSigma,
			Q: []uint64{
				0x10000000006e0001, // 60 Q0
				0x2000000a0001,     // 45
				0x2000000e0001,     // 45
				0x1fffffc20001,     // 45
				0x200000440001,     // 45
				0x200000500001,     // 45
				0x3ffffe80001,      //42 StC
				0x3ffffd20001,      //42 StC
				0x3ffffca0001,      //42 StC
				0xffffffffffc0001,  // ArcSine
				0xfffffffff240001,  // ArcSine
				0x1000000000f00001, // ArcSine
				0xfffffffff840001,  // Double angle
				0x1000000000860001, // Double angle
				0xfffffffff6a0001,  // Sine
				0x1000000000980001, // Sine
				0xfffffffff5a0001,  // Sine
				0x1000000000b00001, // Sine
				0x1000000000ce0001, // Sine
				0xfffffffff2a0001,  // Sine
				0x400000000360001,  // 58 CtS
				0x3ffffffffbe0001,  // 58 CtS
				0x400000000660001,  // 58 CtS
				0x4000000008a0001,  // 58 CtS
			},
			P: []uint64{
				0x1fffffffffe00001, // Pi 61
				0x1fffffffffc80001, // Pi 61
				0x1fffffffffb40001, // Pi 61
				0x1fffffffff500001, // Pi 61
			},
		},
		SlotsToCoeffsParameters: EncodingMatricesParameters{
			LinearTransformType: SlotsToCoeffs,
			LevelStart:          8,
			BSGSRatio:           16.0,
			BitReversed:         false,
			ScalingFactor: [][]float64{
				{0x3ffffe80001},
				{0x3ffffd20001},
				{0x3ffffca0001},
			},
		},
		EvalModParameters: EvalModParameters{
			LevelStart:    19,
			SineType:      Cos1,
			MessageRatio:  4.0,
			K:             25,
			SineDeg:       63,
			DoubleAngle:   2,
			ArcSineDeg:    7,
			ScalingFactor: 1 << 60,
		},
		CoeffsToSlotsParameters: EncodingMatricesParameters{
			LinearTransformType: CoeffsToSlots,
			LevelStart:          23,
			BSGSRatio:           16.0,
			BitReversed:         false,
			ScalingFactor: [][]float64{
				{0x400000000360001},
				{0x3ffffffffbe0001},
				{0x400000000660001},
				{0x4000000008a0001},
			},
		},
	},

	// SET III
	// 1553
	{
		H: 192,
		ParametersLiteral: ParametersLiteral{
			LogN:     16,
			LogSlots: 15,
			Scale:    1 << 30,
			Sigma:    rlwe.DefaultSigma,
			Q: []uint64{
				0x80000000080001,   // 55 Q0
				0xffffffffffc0001,  // 60
				0x10000000006e0001, // 60
				0xfffffffff840001,  // 60
				0x1000000000860001, // 60
				0xfffffffff6a0001,  // 60
				0x1000000000980001, // 60
				0xfffffffff5a0001,  // 60
				0x1000000000b00001, // 60 StC  (30)
				0x1000000000ce0001, // 60 StC  (30+30)
				0x80000000440001,   // 55 Sine (double angle)
				0x7fffffffba0001,   // 55 Sine (double angle)
				0x80000000500001,   // 55 Sine
				0x7fffffffaa0001,   // 55 Sine
				0x800000005e0001,   // 55 Sine
				0x7fffffff7e0001,   // 55 Sine
				0x7fffffff380001,   // 55 Sine
				0x80000000ca0001,   // 55 Sine
				0x200000000e0001,   // 53 CtS
				0x20000000140001,   // 53 CtS
				0x20000000280001,   // 53 CtS
				0x1fffffffd80001,   // 53 CtS
			},
			P: []uint64{
				0x1fffffffffe00001, // Pi 61
				0x1fffffffffc80001, // Pi 61
				0x1fffffffffb40001, // Pi 61
				0x1fffffffff500001, // Pi 61
				0x1fffffffff420001, // Pi 61
			},
		},
		SlotsToCoeffsParameters: EncodingMatricesParameters{
			LinearTransformType: SlotsToCoeffs,
			LevelStart:          9,
			BSGSRatio:           16.0,
			BitReversed:         false,
			ScalingFactor: [][]float64{
				{1073741824.0},
				{1073741824.0062866, 1073741824.0062866},
			},
		},
		EvalModParameters: EvalModParameters{
			LevelStart:    17,
			SineType:      Cos1,
			MessageRatio:  256.0,
			K:             25,
			SineDeg:       63,
			DoubleAngle:   2,
			ArcSineDeg:    0,
			ScalingFactor: 1 << 55,
		},
		CoeffsToSlotsParameters: EncodingMatricesParameters{
			LinearTransformType: CoeffsToSlots,
			LevelStart:          21,
			BSGSRatio:           16.0,
			BitReversed:         false,
			ScalingFactor: [][]float64{
				{0x200000000e0001},
				{0x20000000140001},
				{0x20000000280001},
				{0x1fffffffd80001},
			},
		},
	},

	// Set IV
	// 1792
	{
		H: 32768,
		ParametersLiteral: ParametersLiteral{
			LogN:     16,
			LogSlots: 15,
			Scale:    1 << 40,
			Sigma:    rlwe.DefaultSigma,
			Q: []uint64{
				0x4000000120001, // 60 Q0
				0x10000140001,
				0xffffe80001,
				0xffffc40001,
				0x100003e0001,
				0xffffb20001,
				0x10000500001,
				0xffff940001,
				0xffff8a0001,
				0xffff820001,
				0x100000000060001,  // 56 StC (28 + 28)
				0xffa0001,          // 28 StC
				0xffffffffffc0001,  // 60 Sine (double angle)
				0x10000000006e0001, // 60 Sine (double angle)
				0xfffffffff840001,  // 60 Sine (double angle)
				0x1000000000860001, // 60 Sine (double angle)
				0xfffffffff6a0001,  // 60 Sine
				0x1000000000980001, // 60 Sine
				0xfffffffff5a0001,  // 60 Sine
				0x1000000000b00001, // 60 Sine
				0x1000000000ce0001, // 60 Sine
				0xfffffffff2a0001,  // 60 Sine
				0xfffffffff240001,  // 60 Sine
				0x1000000000f00001, // 60 Sine
				0x200000000e0001,   // 53 CtS
				0x20000000140001,   // 53 CtS
				0x20000000280001,   // 53 CtS
				0x1fffffffd80001,   // 53 CtS
			},
			P: []uint64{
				0x1fffffffffe00001, // Pi 61
				0x1fffffffffc80001, // Pi 61
				0x1fffffffffb40001, // Pi 61
				0x1fffffffff500001, // Pi 61
				0x1fffffffff420001, // Pi 61
				0x1fffffffff380001, // Pi 61
			},
		},
		SlotsToCoeffsParameters: EncodingMatricesParameters{
			LinearTransformType: SlotsToCoeffs,
			LevelStart:          11,
			BSGSRatio:           16.0,
			BitReversed:         false,
			ScalingFactor: [][]float64{
				{268435456.0007324, 268435456.0007324},
				{0xffa0001},
			},
		},
		EvalModParameters: EvalModParameters{
			LevelStart:    23,
			SineType:      Cos2,
			MessageRatio:  256.0,
			K:             325,
			SineDeg:       255,
			DoubleAngle:   4,
			ArcSineDeg:    0,
			ScalingFactor: 1 << 60,
		},
		CoeffsToSlotsParameters: EncodingMatricesParameters{
			LinearTransformType: CoeffsToSlots,
			LevelStart:          27,
			BSGSRatio:           16.0,
			BitReversed:         false,
			ScalingFactor: [][]float64{
				{0x200000000e0001},
				{0x20000000140001},
				{0x20000000280001},
				{0x1fffffffd80001},
			},
		},
	},

	// Set V
	// 768
	{
		H: 192,
		ParametersLiteral: ParametersLiteral{
			LogN:     15,
			LogSlots: 14,
			Scale:    1 << 25,
			Sigma:    rlwe.DefaultSigma,
			Q: []uint64{
				0x1fff90001,       // 32 Q0
				0x4000000420001,   // 50
				0x1fc0001,         // 25
				0xffffffffffc0001, // 60 StC (30+30)
				0x4000000120001,   // 50 Sine
				0x40000001b0001,   // 50 Sine
				0x3ffffffdf0001,   // 50 Sine
				0x4000000270001,   // 50 Sine
				0x3ffffffd20001,   // 50 Sine
				0x3ffffffcd0001,   // 50 Sine
				0x4000000350001,   // 50 Sine
				0x3ffffffc70001,   // 50 Sine
				0x1fffffff50001,   // 49 CtS
				0x1ffffffea0001,   // 49 CtS
			},
			P: []uint64{
				0x7fffffffe0001, // 51
				0x8000000110001, // 51
			},
		},
		SlotsToCoeffsParameters: EncodingMatricesParameters{
			LinearTransformType: SlotsToCoeffs,
			LevelStart:          3,
			BSGSRatio:           16.0,
			BitReversed:         false,
			ScalingFactor: [][]float64{
				{1073741823.9998779, 1073741823.9998779},
			},
		},
		EvalModParameters: EvalModParameters{
			LevelStart:    11,
			SineType:      Cos1,
			MessageRatio:  256.0,
			K:             25,
			SineDeg:       63,
			DoubleAngle:   2,
			ArcSineDeg:    0,
			ScalingFactor: 1 << 50,
		},
		CoeffsToSlotsParameters: EncodingMatricesParameters{
			LinearTransformType: CoeffsToSlots,
			LevelStart:          13,
			BSGSRatio:           16.0,
			BitReversed:         false,
			ScalingFactor: [][]float64{
				{0x1fffffff50001},
				{0x1ffffffea0001},
			},
		},
	},
}
