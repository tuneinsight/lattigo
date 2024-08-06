package examples

import (
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/bgv"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

var (

	// BGVParamsN12QP109 is an example parameter set for the `bgv` package logN=12 and logQP=109.
	// These parameters expect the user to use the regular tensoring (i.e. Evaluator.Mul) followed
	// by the rescaling (i.e. Evaluator.Rescale).
	BGVParamsN12QP109 = bgv.ParametersLiteral{
		LogN:             12,
		LogQ:             []int{39, 31},
		LogP:             []int{39},
		PlaintextModulus: 0x10001,
	}

	// BGVParamsN13QP218 is an example parameter set for the `bgv` package with logN=13 and logQP=218.
	// These parameters expect the user to use the regular tensoring (i.e. Evaluator.Mul) followed
	// by the rescaling (i.e. Evaluator.Rescale).
	BGVParamsN13QP218 = bgv.ParametersLiteral{
		LogN:             13,
		LogQ:             []int{42, 33, 33, 33, 33},
		LogP:             []int{44},
		PlaintextModulus: 0x10001,
	}

	// BGVParamsN14QP438 is an example parameter set for the `bgv` package with logN=14 and logQP=438.
	// These parameters expect the user to use the regular tensoring (i.e. Evaluator.Mul) followed
	// by the rescaling (i.e. Evaluator.Rescale).
	BGVParamsN14QP438 = bgv.ParametersLiteral{
		LogN:             14,
		LogQ:             []int{44, 34, 34, 34, 34, 34, 34, 34, 34, 34},
		LogP:             []int{44, 44},
		PlaintextModulus: 0x10001,
	}

	// BGVParamsN15QP880 is an example parameter set for the `bgv` package with logN=15 and logQP=881.
	// These parameters expect the user to use the regular tensoring (i.e. Evaluator.Mul) followed
	// by the rescaling (i.e. Evaluator.Rescale).
	BGVParamsN15QP880 = bgv.ParametersLiteral{
		LogN:             15,
		LogQ:             []int{47, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34},
		LogP:             []int{47, 47, 47, 47},
		PlaintextModulus: 0x10001,
	}

	// BGVScaleInvariantParamsN12QP109 is an example parameter set for the `bgv` package logN=12 and logQP=109.
	// These parameters expect the user to use the scale invariant tensoring (i.e. Evaluator.MulScaleInvariant).
	BGVScaleInvariantParamsN12QP109 = bgv.ParametersLiteral{
		LogN:             12,
		LogQ:             []int{39, 39},
		LogP:             []int{31},
		PlaintextModulus: 0x10001,
	}

	// BGVScaleInvariantParamsN13QP218 is an example parameter set for the `bgv` package with logN=13 and logQP=218.
	// These parameters expect the user to use the scale invariant tensoring (i.e. Evaluator.MulScaleInvariant).
	BGVScaleInvariantParamsN13QP218 = bgv.ParametersLiteral{
		LogN:             13,
		LogQ:             []int{55, 54, 54},
		LogP:             []int{55},
		PlaintextModulus: 0x10001,
	}

	// BGVScaleInvariantParamsN14QP438 is an example parameter set for the `bgv` package with logN=14 and logQP=438.
	// These parameters expect the user to use the scale invariant tensoring (i.e. Evaluator.MulScaleInvariant).
	BGVScaleInvariantParamsN14QP438 = bgv.ParametersLiteral{
		LogN:             14,
		LogQ:             []int{55, 55, 55, 54, 54, 54},
		LogP:             []int{56, 55},
		PlaintextModulus: 0x10001,
	}

	// BGVScaleInvariantParamsN15QP880 is an example parameter set for the `bgv` package with logN=15 and logQP=881.
	// These parameters expect the user to use the scale invariant tensoring (i.e. Evaluator.MulScaleInvariant).
	BGVScaleInvariantParamsN15QP880 = bgv.ParametersLiteral{
		LogN:             15,
		LogQ:             []int{60, 60, 59, 58, 58, 58, 58, 58, 58, 58, 58, 58},
		LogP:             []int{60, 60, 60},
		PlaintextModulus: 0x10001,
	}

	// CKKSComplexParamsN12QP109 is an example parameter set for the `ckks` package with logN=12 and logQP=109.
	// These parameters instantiate `ckks` over the complex field with N/2 SIMD slots.
	CKKSComplexParamsN12QP109 = ckks.ParametersLiteral{
		LogN:            12,
		LogQ:            []int{38, 32},
		LogP:            []int{39},
		LogDefaultScale: 32,
	}

	// CKKSComplexParamsN13QP218 is an example parameter set for the `ckks` package with logN=13 and logQP=218.
	// These parameters instantiate `ckks` over the complex field with N/2 SIMD slots.
	CKKSComplexParamsN13QP218 = ckks.ParametersLiteral{
		LogN:            13,
		LogQ:            []int{33, 30, 30, 30, 30, 30},
		LogP:            []int{35},
		LogDefaultScale: 30,
	}
	// CKKSComplexParamsN14QP438 is an example parameter set for the `ckks` package with logN=14 and logQP=438.
	// These parameters instantiate `ckks` over the complex field with N/2 SIMD slots.
	CKKSComplexParamsN14QP438 = ckks.ParametersLiteral{
		LogN:            14,
		LogQ:            []int{45, 34, 34, 34, 34, 34, 34, 34, 34, 34},
		LogP:            []int{44, 43},
		LogDefaultScale: 34,
	}

	// CKKSComplexParamsN15QP880 is an example parameter set for the `ckks` package with logN=15 and logQP=881.
	// These parameters instantiate `ckks` over the complex field with N/2 SIMD slots.
	CKKSComplexParamsN15QP881 = ckks.ParametersLiteral{
		LogN:            15,
		LogQ:            []int{51, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40},
		LogP:            []int{50, 50, 50},
		LogDefaultScale: 40,
	}
	// CKKSComplexParamsPN16QP1761 is an example parameter set for the `ckks` package with logN=16 and logQP = 1761.
	// These parameters instantiate `ckks` over the complex field with N/2 SIMD slots.
	CKKSComplexParamsPN16QP1761 = ckks.ParametersLiteral{
		LogN:            16,
		LogQ:            []int{56, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
		LogP:            []int{55, 55, 55, 55},
		LogDefaultScale: 45,
	}

	// CKKSRealParamsN12QP109 is an example parameter set for the `ckks` package with conjugate-invariant CKKS and logN=12 and logQP=109.
	// These parameters instantiate `ckks` over the real field with N SIMD slots.
	CKKSRealParamsN12QP109 = ckks.ParametersLiteral{
		LogN:            12,
		LogQ:            []int{38, 32},
		LogP:            []int{39},
		LogDefaultScale: 32,
		RingType:        ring.ConjugateInvariant,
	}

	// CKKSRealParamsN13QP218 is an example parameter set for the `ckks` package with conjugate-invariant CKKS and logN=13 and logQP=218
	// These parameters instantiate `ckks` over the real field with N SIMD slots.
	CKKSRealParamsN13QP218 = ckks.ParametersLiteral{
		LogN:            13,
		LogQ:            []int{33, 30, 30, 30, 30, 30},
		LogP:            []int{35},
		LogDefaultScale: 30,
		RingType:        ring.ConjugateInvariant,
	}

	// CKKSRealParamsN14QP438 is an example parameter set for the `ckks` package with logN=14 and logQP=438.
	// These parameters instantiate `ckks` over the real field with N SIMD slots.
	CKKSRealParamsN14QP438 = ckks.ParametersLiteral{
		LogN:            14,
		LogQ:            []int{46, 34, 34, 34, 34, 34, 34, 34, 34, 34},
		LogP:            []int{43, 43},
		LogDefaultScale: 34,
		RingType:        ring.ConjugateInvariant,
	}

	// CKKSRealParamsN15QP880 is an example parameter set for the `ckks` package with logN=15 and logQP=881.
	// These parameters instantiate `ckks` over the real field with N SIMD slots.
	CKKSRealParamsN15QP881 = ckks.ParametersLiteral{
		LogN:            15,
		LogQ:            []int{51, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40},
		LogP:            []int{50, 50, 50},
		LogDefaultScale: 40,
		RingType:        ring.ConjugateInvariant,
	}

	// CKKSRealParamsPN16QP1761 is an example parameter set for the `ckks` package with logN=16 and logQP = 1761
	// These parameters instantiate `ckks` over the real field with N SIMD slots.
	CKKSRealParamsPN16QP1761 = ckks.ParametersLiteral{
		LogN:            16,
		LogQ:            []int{56, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
		LogP:            []int{55, 55, 55, 55},
		LogDefaultScale: 45,
		RingType:        ring.ConjugateInvariant,
	}
)

var BGVParams = []bgv.ParametersLiteral{BGVParamsN12QP109, BGVParamsN13QP218, BGVParamsN14QP438, BGVParamsN15QP880}

var BGVScaleInvariantParams = []bgv.ParametersLiteral{BGVScaleInvariantParamsN12QP109, BGVScaleInvariantParamsN13QP218, BGVScaleInvariantParamsN14QP438, BGVScaleInvariantParamsN15QP880}

var CKKSComplexParams = []ckks.ParametersLiteral{CKKSComplexParamsN12QP109, CKKSComplexParamsN13QP218, CKKSComplexParamsN14QP438, CKKSComplexParamsN15QP881, CKKSComplexParamsPN16QP1761}

var CKKSRealParams = []ckks.ParametersLiteral{CKKSRealParamsN12QP109, CKKSRealParamsN13QP218, CKKSRealParamsN14QP438, CKKSRealParamsN15QP881, CKKSRealParamsPN16QP1761}
