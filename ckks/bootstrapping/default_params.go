package bootstrapping

import (
	"github.com/tuneinsight/lattigo/v3/ckks"
	"github.com/tuneinsight/lattigo/v3/ckks/advanced"
)

type defaultParametersLiteral struct {
	SchemeParams        ckks.ParametersLiteral
	BootstrappingParams Parameters
}

// The parameters provided hereunder are the parameters used in the paper
// Bootstrapping for Approximate Homomorphic Encryption with Negligible
// Failure-Probability by Using Sparse-Secret Encapsulation,
// https://eprint.iacr.org/2022/024

// DefaultParametersSparse is a set of default bootstrapping parameters with H=192 as main secret and H=32 as ephemeral secret.
var DefaultParametersSparse = []defaultParametersLiteral{N16QP1546H192H32, N16QP1547H192H32, N16QP1553H192H32, N15QP768H192H32}

// DefaultParametersDense is a set of default bootstrapping parameters with H=N/2 as main secret and H=32 as ephemeral secret.
var DefaultParametersDense = []defaultParametersLiteral{N16QP1767H32768H32, N16QP1788H32768H32, N16QP1793H32768H32, N15QP880H16384H32}

var (
	// N16QP1546H192H32 is a default bootstrapping parameters for a main secret with H=192 and an ephemeral secret with H=32.
	// Residual Q : []int{60, 40, 40, 40, 40, 40, 40, 40, 40, 40} (420 bits).
	// SlotsToCoeffs Q: []int{39, 39, 39}.
	// EvalMod Q: []int{60, 60, 60, 60, 60, 60, 60, 60}.
	// CoeffsToSlots Q: []int{56, 56, 56, 56}.
	// Precision : 26.6 bits for 2^{15} slots.
	// Failure : 2^{-138.7} for 2^{15} slots.
	N16QP1546H192H32 = defaultParametersLiteral{
		ckks.ParametersLiteral{
			LogN:         16,
			LogQ:         []int{60, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 39, 39, 60, 60, 60, 60, 60, 60, 60, 60, 56, 56, 56, 56},
			LogP:         []int{61, 61, 61, 61, 61},
			H:            192,
			LogSlots:     15,
			DefaultScale: 1 << 40,
		},
		Parameters{
			EphemeralSecretWeight: 32,
			SlotsToCoeffsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.SlotsToCoeffs,
				RepackImag2Real:     true,
				LevelStart:          12,
				BSGSRatio:           2.0,
				ScalingFactor:       [][]float64{{0}, {0}, {0}},
			},
			EvalModParameters: advanced.EvalModLiteral{
				LevelStart:    20,
				SineType:      advanced.Cos1,
				MessageRatio:  256.0,
				K:             16,
				SineDeg:       30,
				DoubleAngle:   3,
				ScalingFactor: 1 << 60,
			},
			CoeffsToSlotsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.CoeffsToSlots,
				RepackImag2Real:     true,
				LevelStart:          24,
				BSGSRatio:           2.0,
				ScalingFactor:       [][]float64{{0}, {0}, {0}, {0}},
			},
		},
	}

	// N16QP1547H192H32 is a default bootstrapping parameters for a main secret with H=192 and an ephemeral secret with H=32.
	// Residual Q : []int{60, 45, 45, 45, 45} (285 bits).
	// SlotsToCoeffs Q: []int{42, 42, 42}.
	// EvalMod Q: []int{60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60}.
	// CoeffsToSlots Q: []int{58, 58, 58, 58}.
	// Precision : 32.1 bits for 2^{15} slots.
	// Failure : 2^{-138.7} for 2^{15} slots.
	N16QP1547H192H32 = defaultParametersLiteral{
		ckks.ParametersLiteral{
			LogN:         16,
			LogQ:         []int{60, 45, 45, 45, 45, 45, 42, 42, 42, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 58, 58, 58, 58},
			LogP:         []int{61, 61, 61, 61},
			H:            192,
			LogSlots:     15,
			DefaultScale: 1 << 45,
		},
		Parameters{
			EphemeralSecretWeight: 32,
			SlotsToCoeffsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.SlotsToCoeffs,
				RepackImag2Real:     true,
				LevelStart:          8,
				BSGSRatio:           2.0,
				ScalingFactor:       [][]float64{{0}, {0}, {0}},
			},
			EvalModParameters: advanced.EvalModLiteral{
				LevelStart:    19,
				SineType:      advanced.Cos1,
				MessageRatio:  4.0,
				K:             16,
				SineDeg:       30,
				DoubleAngle:   3,
				ArcSineDeg:    7,
				ScalingFactor: 1 << 60,
			},
			CoeffsToSlotsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.CoeffsToSlots,
				RepackImag2Real:     true,
				LevelStart:          23,
				BSGSRatio:           2.0,
				ScalingFactor:       [][]float64{{0}, {0}, {0}, {0}},
			},
		},
	}

	// N16QP1553H192H32 is a default bootstrapping parameters for a main secret with H=192 and an ephemeral secret with H=32.
	// Residual Q : []int{55, 60, 60, 60, 60, 60, 60, 60, 60} (505 bits).
	// SlotsToCoeffs Q: []int{60, 60}.
	// EvalMod Q: []int{55, 55, 55, 55, 55, 55, 55, 55}.
	// CoeffsToSlots Q: []int{53, 53, 53, 53}.
	// Precision : 19.1 bits for 2^{15} slots.
	// Failure : 2^{-138.7} for 2^{15} slots.
	N16QP1553H192H32 = defaultParametersLiteral{
		ckks.ParametersLiteral{
			LogN:         16,
			LogQ:         []int{55, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 55, 55, 55, 55, 55, 55, 55, 53, 53, 53, 53},
			LogP:         []int{61, 61, 61, 61, 61},
			H:            192,
			LogSlots:     15,
			DefaultScale: 1 << 30,
		},
		Parameters{
			EphemeralSecretWeight: 32,
			SlotsToCoeffsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.SlotsToCoeffs,
				RepackImag2Real:     true,
				LevelStart:          9,
				BSGSRatio:           2.0,
				ScalingFactor:       [][]float64{{0}, {0, 0}},
			},
			EvalModParameters: advanced.EvalModLiteral{
				LevelStart:    17,
				SineType:      advanced.Cos1,
				MessageRatio:  256.0,
				K:             16,
				SineDeg:       30,
				DoubleAngle:   3,
				ScalingFactor: 1 << 55,
			},
			CoeffsToSlotsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.CoeffsToSlots,
				RepackImag2Real:     true,
				LevelStart:          21,
				BSGSRatio:           2.0,
				ScalingFactor:       [][]float64{{0}, {0}, {0}, {0}},
			},
		},
	}

	// N15QP768H192H32 is a default bootstrapping parameters for a main secret with H=192 and an ephemeral secret with H=32.
	// Residual Q : []int{32, 50, 25} (110 bits).
	// SlotsToCoeffs Q: []int{60}.
	// EvalMod Q: []int{50, 50, 50, 50, 50, 50, 50, 50}.
	// CoeffsToSlots Q: []int{49, 49}.
	// Precision : 15.4 bits for 2^{14} slots.
	// Failure : 2^{-139.7} for 2^{14} slots.
	N15QP768H192H32 = defaultParametersLiteral{
		ckks.ParametersLiteral{
			LogN:         15,
			LogQ:         []int{33, 50, 25, 60, 50, 50, 50, 50, 50, 50, 50, 50, 49, 49},
			LogP:         []int{51, 51},
			H:            192,
			LogSlots:     14,
			DefaultScale: 1 << 25,
		},
		Parameters{
			EphemeralSecretWeight: 32,
			SlotsToCoeffsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.SlotsToCoeffs,
				RepackImag2Real:     true,
				LevelStart:          3,
				BSGSRatio:           2.0,
				ScalingFactor:       [][]float64{{0, 0}},
			},
			EvalModParameters: advanced.EvalModLiteral{
				LevelStart:    11,
				SineType:      advanced.Cos1,
				MessageRatio:  256.0,
				K:             16,
				SineDeg:       30,
				DoubleAngle:   3,
				ScalingFactor: 1 << 50,
			},
			CoeffsToSlotsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.CoeffsToSlots,
				RepackImag2Real:     true,
				LevelStart:          13,
				BSGSRatio:           2.0,
				ScalingFactor:       [][]float64{{0}, {0}},
			},
		},
	}

	// N16QP1767H32768H32 is a default bootstrapping parameters for a main secret with H=32768 and an ephemeral secret with H=32.
	// Residual Q : []int{60, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40} (580 bits).
	// SlotsToCoeffs Q: []int{39, 39, 39}.
	// EvalMod Q: []int{60, 60, 60, 60, 60, 60, 60, 60, 60}.
	// CoeffsToSlots Q: []int{56, 56, 56, 56}.
	// Precision : 23.0 bits for 2^{15} slots.
	// Failure : 2^{-138.7} for 2^{15} slots.
	N16QP1767H32768H32 = defaultParametersLiteral{
		ckks.ParametersLiteral{
			LogN:         16,
			LogQ:         []int{60, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 39, 39, 60, 60, 60, 60, 60, 60, 60, 60, 60, 56, 56, 56, 56},
			LogP:         []int{61, 61, 61, 61, 61, 61},
			H:            32768,
			LogSlots:     15,
			DefaultScale: 1 << 40,
		},
		Parameters{
			EphemeralSecretWeight: 32,
			SlotsToCoeffsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.SlotsToCoeffs,
				RepackImag2Real:     true,
				LevelStart:          16,
				BSGSRatio:           2.0,
				ScalingFactor:       [][]float64{{0}, {0}, {0}},
			},
			EvalModParameters: advanced.EvalModLiteral{
				LevelStart:    24,
				SineType:      advanced.Cos1,
				MessageRatio:  256.0,
				K:             16,
				SineDeg:       30,
				DoubleAngle:   3,
				ScalingFactor: 1 << 60,
			},
			CoeffsToSlotsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.CoeffsToSlots,
				RepackImag2Real:     true,
				LevelStart:          28,
				BSGSRatio:           2.0,
				ScalingFactor:       [][]float64{{0}, {0}, {0}, {0}},
			},
		},
	}

	// N16QP1788H32768H32 is a default bootstrapping parameters for a main secret with H=32768 and an ephemeral secret with H=32.
	// Residual Q : []int{60, 45, 45, 45, 45, 45, 45, 45, 45, 45} (465 bits).
	// SlotsToCoeffs Q: []int{42, 42, 42}.
	// EvalMod Q: []int{60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60}.
	// CoeffsToSlots Q: []int{58, 58, 58, 58}.
	// Precision : 29.0 bits for 2^{15} slots.
	// Failure : 2^{-138.7} for 2^{15} slots.
	N16QP1788H32768H32 = defaultParametersLiteral{
		ckks.ParametersLiteral{
			LogN:         16,
			LogQ:         []int{60, 45, 45, 45, 45, 45, 45, 45, 45, 45, 42, 42, 42, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 58, 58, 58, 58},
			LogP:         []int{61, 61, 61, 61, 61},
			H:            32768,
			LogSlots:     15,
			DefaultScale: 1 << 45,
		},
		Parameters{
			EphemeralSecretWeight: 32,
			SlotsToCoeffsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.SlotsToCoeffs,
				RepackImag2Real:     true,
				LevelStart:          12,
				BSGSRatio:           2.0,
				ScalingFactor:       [][]float64{{0}, {0}, {0}},
			},
			EvalModParameters: advanced.EvalModLiteral{
				LevelStart:    23,
				SineType:      advanced.Cos1,
				MessageRatio:  4.0,
				K:             16,
				SineDeg:       30,
				DoubleAngle:   3,
				ArcSineDeg:    7,
				ScalingFactor: 1 << 60,
			},
			CoeffsToSlotsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.CoeffsToSlots,
				RepackImag2Real:     true,
				LevelStart:          27,
				BSGSRatio:           2.0,
				ScalingFactor:       [][]float64{{0}, {0}, {0}, {0}},
			},
		},
	}

	// N16QP1793H32768H32 is a default bootstrapping parameters for a main secret with H=32768 and an ephemeral secret with H=32.
	// Residual Q : []int{55, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 30} 745 bits.
	// SlotsToCoeffs Q: []int{30, 60}.
	// EvalMod Q: []int{55, 55, 55, 55, 55, 55, 55, 55}.
	// CoeffsToSlots Q: []int{53, 53, 53, 53}.
	// Precision : 17.8 bits for 2^{15} slots.
	// Failure : 2^{-138.7} for 2^{15} slots.
	N16QP1793H32768H32 = defaultParametersLiteral{
		ckks.ParametersLiteral{
			LogN:         16,
			LogQ:         []int{55, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 55, 55, 55, 55, 55, 55, 55, 55, 53, 53, 53, 53},
			LogP:         []int{61, 61, 61, 61, 61},
			H:            32768,
			LogSlots:     15,
			DefaultScale: 1 << 30,
		},
		Parameters{
			EphemeralSecretWeight: 32,
			SlotsToCoeffsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.SlotsToCoeffs,
				RepackImag2Real:     true,
				LevelStart:          13,
				BSGSRatio:           2.0,
				ScalingFactor:       [][]float64{{0}, {0, 0}},
			},
			EvalModParameters: advanced.EvalModLiteral{
				LevelStart:    21,
				SineType:      advanced.Cos1,
				MessageRatio:  256.0,
				K:             16,
				SineDeg:       30,
				DoubleAngle:   3,
				ScalingFactor: 1 << 55,
			},
			CoeffsToSlotsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.CoeffsToSlots,
				RepackImag2Real:     true,
				LevelStart:          25,
				BSGSRatio:           2.0,
				ScalingFactor:       [][]float64{{0}, {0}, {0}, {0}},
			},
		},
	}

	// N15QP880H16384H32 is a default bootstrapping parameters for a main secret with H=16384 and an ephemeral secret with H=32.
	// Residual Q : []int{40, 31, 31, 31, 31} (166 bits).
	// SlotsToCoeffs Q: []int{60}.
	// EvalMod Q: []int{55, 55, 55, 55, 55, 55, 55, 55}.
	// CoeffsToSlots Q: []int{52, 52}.
	// Precision : 17.3 bits for 2^{14} slots.
	// Failure : 2^{-139.7} for 2^{14} slots.
	N15QP880H16384H32 = defaultParametersLiteral{
		ckks.ParametersLiteral{
			LogN:         15,
			LogQ:         []int{40, 31, 31, 31, 31, 60, 55, 55, 55, 55, 55, 55, 55, 55, 52, 52},
			LogP:         []int{56, 56},
			H:            16384,
			LogSlots:     14,
			DefaultScale: 1 << 31,
		},
		Parameters{
			EphemeralSecretWeight: 32,
			SlotsToCoeffsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.SlotsToCoeffs,
				RepackImag2Real:     true,
				LevelStart:          5,
				BSGSRatio:           2.0,
				ScalingFactor:       [][]float64{{0, 0}},
			},
			EvalModParameters: advanced.EvalModLiteral{
				LevelStart:    13,
				SineType:      advanced.Cos1,
				MessageRatio:  256.0,
				K:             16,
				SineDeg:       30,
				DoubleAngle:   3,
				ScalingFactor: 1 << 55,
			},
			CoeffsToSlotsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.CoeffsToSlots,
				RepackImag2Real:     true,
				LevelStart:          15,
				BSGSRatio:           2.0,
				ScalingFactor:       [][]float64{{0}, {0}},
			},
		},
	}
)
