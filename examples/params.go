package examples

import (
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/he/heint"
	"github.com/tuneinsight/lattigo/v5/ring"
)

var (

	// HEIntParamsN12QP109 is an example parameter set for the `heint` package logN=12 and logQP=109.
	// These parameters expect the user to use the regular tensoring (i.e. Evaluator.Mul) followed
	// by the rescaling (i.e. Evaluator.Rescale).
	HEIntParamsN12QP109 = heint.ParametersLiteral{
		LogN:             12,
		LogQ:             []int{39, 31},
		LogP:             []int{39},
		PlaintextModulus: 0x10001,
	}

	// HEIntParamsN13QP218 is an example parameter set for the `heint` package with logN=13 and logQP=218.
	// These parameters expect the user to use the regular tensoring (i.e. Evaluator.Mul) followed
	// by the rescaling (i.e. Evaluator.Rescale).
	HEIntParamsN13QP218 = heint.ParametersLiteral{
		LogN:             13,
		LogQ:             []int{42, 33, 33, 33, 33},
		LogP:             []int{44},
		PlaintextModulus: 0x10001,
	}

	// HEIntParamsN14QP438 is an example parameter set for the `heint` package with logN=14 and logQP=438.
	// These parameters expect the user to use the regular tensoring (i.e. Evaluator.Mul) followed
	// by the rescaling (i.e. Evaluator.Rescale).
	HEIntParamsN14QP438 = heint.ParametersLiteral{
		LogN:             14,
		LogQ:             []int{44, 34, 34, 34, 34, 34, 34, 34, 34, 34},
		LogP:             []int{44, 44},
		PlaintextModulus: 0x10001,
	}

	// HEIntParamsN15QP880 is an example parameter set for the `heint` package with logN=15 and logQP=881.
	// These parameters expect the user to use the regular tensoring (i.e. Evaluator.Mul) followed
	// by the rescaling (i.e. Evaluator.Rescale).
	HEIntParamsN15QP880 = heint.ParametersLiteral{
		LogN:             15,
		LogQ:             []int{47, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34},
		LogP:             []int{47, 47, 47, 47},
		PlaintextModulus: 0x10001,
	}

	// HEIntScaleInvariantParamsN12QP109 is an example parameter set for the `heint` package logN=12 and logQP=109.
	// These parameters expect the user to use the scale invariant tensoring (i.e. Evaluator.MulScaleInvariant).
	HEIntScaleInvariantParamsN12QP109 = heint.ParametersLiteral{
		LogN:             12,
		LogQ:             []int{39, 39},
		LogP:             []int{31},
		PlaintextModulus: 0x10001,
	}

	// HEIntScaleInvariantParamsN13QP218 is an example parameter set for the `heint` package with logN=13 and logQP=218.
	// These parameters expect the user to use the scale invariant tensoring (i.e. Evaluator.MulScaleInvariant).
	HEIntScaleInvariantParamsN13QP218 = heint.ParametersLiteral{
		LogN:             13,
		LogQ:             []int{55, 54, 54},
		LogP:             []int{55},
		PlaintextModulus: 0x10001,
	}

	// HEIntScaleInvariantParamsN14QP438 is an example parameter set for the `heint` package with logN=14 and logQP=438.
	// These parameters expect the user to use the scale invariant tensoring (i.e. Evaluator.MulScaleInvariant).
	HEIntScaleInvariantParamsN14QP438 = heint.ParametersLiteral{
		LogN:             14,
		LogQ:             []int{55, 55, 55, 54, 54, 54},
		LogP:             []int{56, 55},
		PlaintextModulus: 0x10001,
	}

	// HEIntScaleInvariantParamsN15QP880 is an example parameter set for the `heint` package with logN=15 and logQP=881.
	// These parameters expect the user to use the scale invariant tensoring (i.e. Evaluator.MulScaleInvariant).
	HEIntScaleInvariantParamsN15QP880 = heint.ParametersLiteral{
		LogN:             15,
		LogQ:             []int{60, 60, 59, 58, 58, 58, 58, 58, 58, 58, 58, 58},
		LogP:             []int{60, 60, 60},
		PlaintextModulus: 0x10001,
	}

	// HEFloatComplexParamsN12QP109 is an example parameter set for the `hefloat` package with logN=12 and logQP=109.
	// These parameters instantiate `hefloat` over the complex field with N/2 SIMD slots.
	HEFloatComplexParamsN12QP109 = hefloat.ParametersLiteral{
		LogN:            12,
		LogQ:            []int{38, 32},
		LogP:            []int{39},
		LogDefaultScale: 32,
	}

	// HEFloatComplexParamsN13QP218 is an example parameter set for the `hefloat` package with logN=13 and logQP=218.
	// These parameters instantiate `hefloat` over the complex field with N/2 SIMD slots.
	HEFloatComplexParamsN13QP218 = hefloat.ParametersLiteral{
		LogN:            13,
		LogQ:            []int{33, 30, 30, 30, 30, 30},
		LogP:            []int{35},
		LogDefaultScale: 30,
	}
	// HEFloatComplexParamsN14QP438 is an example parameter set for the `hefloat` package with logN=14 and logQP=438.
	// These parameters instantiate `hefloat` over the complex field with N/2 SIMD slots.
	HEFloatComplexParamsN14QP438 = hefloat.ParametersLiteral{
		LogN:            14,
		LogQ:            []int{45, 34, 34, 34, 34, 34, 34, 34, 34, 34},
		LogP:            []int{44, 43},
		LogDefaultScale: 34,
	}

	// HEFloatComplexParamsN15QP880 is an example parameter set for the `hefloat` package with logN=15 and logQP=881.
	// These parameters instantiate `hefloat` over the complex field with N/2 SIMD slots.
	HEFloatComplexParamsN15QP881 = hefloat.ParametersLiteral{
		LogN:            15,
		LogQ:            []int{51, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40},
		LogP:            []int{50, 50, 50},
		LogDefaultScale: 40,
	}
	// HEFloatComplexParamsPN16QP1761 is an example parameter set for the `hefloat` package with logN=16 and logQP = 1761.
	// These parameters instantiate `hefloat` over the complex field with N/2 SIMD slots.
	HEFloatComplexParamsPN16QP1761 = hefloat.ParametersLiteral{
		LogN:            16,
		LogQ:            []int{56, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
		LogP:            []int{55, 55, 55, 55},
		LogDefaultScale: 45,
	}

	// HEFloatRealParamsN12QP109 is an example parameter set for the `hefloat` package with conjugate-invariant CKKS and logN=12 and logQP=109.
	// These parameters instantiate `hefloat` over the real field with N SIMD slots.
	HEFloatRealParamsN12QP109 = hefloat.ParametersLiteral{
		LogN:            12,
		LogQ:            []int{38, 32},
		LogP:            []int{39},
		LogDefaultScale: 32,
		RingType:        ring.ConjugateInvariant,
	}

	// HEFloatRealParamsN13QP218 is an example parameter set for the `hefloat` package with conjugate-invariant CKKS and logN=13 and logQP=218
	// These parameters instantiate `hefloat` over the real field with N SIMD slots.
	HEFloatRealParamsN13QP218 = hefloat.ParametersLiteral{
		LogN:            13,
		LogQ:            []int{33, 30, 30, 30, 30, 30},
		LogP:            []int{35},
		LogDefaultScale: 30,
		RingType:        ring.ConjugateInvariant,
	}

	// HEFloatRealParamsN14QP438 is an example parameter set for the `hefloat` package with logN=14 and logQP=438.
	// These parameters instantiate `hefloat` over the real field with N SIMD slots.
	HEFloatRealParamsN14QP438 = hefloat.ParametersLiteral{
		LogN:            14,
		LogQ:            []int{46, 34, 34, 34, 34, 34, 34, 34, 34, 34},
		LogP:            []int{43, 43},
		LogDefaultScale: 34,
		RingType:        ring.ConjugateInvariant,
	}

	// HEFloatRealParamsN15QP880 is an example parameter set for the `hefloat` package with logN=15 and logQP=881.
	// These parameters instantiate `hefloat` over the real field with N SIMD slots.
	HEFloatRealParamsN15QP881 = hefloat.ParametersLiteral{
		LogN:            15,
		LogQ:            []int{51, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40},
		LogP:            []int{50, 50, 50},
		LogDefaultScale: 40,
		RingType:        ring.ConjugateInvariant,
	}

	// HEFloatRealParamsPN16QP1761 is an example parameter set for the `hefloat` package with logN=16 and logQP = 1761
	// These parameters instantiate `hefloat` over the real field with N SIMD slots.
	HEFloatRealParamsPN16QP1761 = hefloat.ParametersLiteral{
		LogN:            16,
		LogQ:            []int{56, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
		LogP:            []int{55, 55, 55, 55},
		LogDefaultScale: 45,
		RingType:        ring.ConjugateInvariant,
	}
)

var HEIntParams = []heint.ParametersLiteral{HEIntParamsN12QP109, HEIntParamsN13QP218, HEIntParamsN14QP438, HEIntParamsN15QP880}

var HEIntScaleInvariantParams = []heint.ParametersLiteral{HEIntScaleInvariantParamsN12QP109, HEIntScaleInvariantParamsN13QP218, HEIntScaleInvariantParamsN14QP438, HEIntScaleInvariantParamsN15QP880}

var HEFloatComplexParams = []hefloat.ParametersLiteral{HEFloatComplexParamsN12QP109, HEFloatComplexParamsN13QP218, HEFloatComplexParamsN14QP438, HEFloatComplexParamsN15QP881, HEFloatComplexParamsPN16QP1761}

var HEFloatRealParams = []hefloat.ParametersLiteral{HEFloatRealParamsN12QP109, HEFloatRealParamsN13QP218, HEFloatRealParamsN14QP438, HEFloatRealParamsN15QP881, HEFloatRealParamsPN16QP1761}
