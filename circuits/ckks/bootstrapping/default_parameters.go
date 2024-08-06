package bootstrapping

import (
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils"
)

type defaultParametersLiteral struct {
	SchemeParams        ckks.ParametersLiteral
	BootstrappingParams ParametersLiteral
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
			LogN:            16,
			LogQ:            []int{60, 40, 40, 40, 40, 40, 40, 40, 40, 40},
			LogP:            []int{61, 61, 61, 61, 61},
			Xs:              ring.Ternary{H: 192},
			LogDefaultScale: 40,
		},
		ParametersLiteral{},
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
			LogN:            16,
			LogQ:            []int{60, 45, 45, 45, 45, 45},
			LogP:            []int{61, 61, 61, 61},
			Xs:              ring.Ternary{H: 192},
			LogDefaultScale: 45,
		},
		ParametersLiteral{
			SlotsToCoeffsFactorizationDepthAndLogScales: [][]int{{42}, {42}, {42}},
			CoeffsToSlotsFactorizationDepthAndLogScales: [][]int{{58}, {58}, {58}, {58}},
			LogMessageRatio: utils.Pointy(2),
			Mod1InvDegree:   utils.Pointy(7),
		},
	}

	// N16QP1553H192H32 is a default bootstrapping parameters for a main secret with H=192 and an ephemeral secret with H=32.
	// Residual Q : []int{55, 60, 60, 60, 60, 60, 60, 60, 30} (505 bits).
	// SlotsToCoeffs Q: []int{30, {30, 30}}.
	// EvalMod Q: []int{55, 55, 55, 55, 55, 55, 55, 55}.
	// CoeffsToSlots Q: []int{53, 53, 53, 53}.
	// Precision : 19.1 bits for 2^{15} slots.
	// Failure : 2^{-138.7} for 2^{15} slots.
	N16QP1553H192H32 = defaultParametersLiteral{
		ckks.ParametersLiteral{
			LogN:            16,
			LogQ:            []int{55, 60, 60, 60, 60, 60, 60, 60},
			LogP:            []int{61, 61, 61, 61, 61},
			Xs:              ring.Ternary{H: 192},
			LogDefaultScale: 30,
		},
		ParametersLiteral{
			SlotsToCoeffsFactorizationDepthAndLogScales: [][]int{{30}, {30, 30}},
			CoeffsToSlotsFactorizationDepthAndLogScales: [][]int{{53}, {53}, {53}, {53}},
			EvalModLogScale: utils.Pointy(55),
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
			LogN:            15,
			LogQ:            []int{33, 50, 25},
			LogP:            []int{51, 51},
			Xs:              ring.Ternary{H: 192},
			LogDefaultScale: 25,
		},
		ParametersLiteral{
			SlotsToCoeffsFactorizationDepthAndLogScales: [][]int{{30, 30}},
			CoeffsToSlotsFactorizationDepthAndLogScales: [][]int{{49}, {49}},
			EvalModLogScale: utils.Pointy(50),
		},
	}

	// N16QP1767H32768H32 is a default bootstrapping parameters for a main secret with H=32768 and an ephemeral secret with H=32.
	// Residual Q : []int{60, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40} (580 bits).
	// SlotsToCoeffs Q: []int{39, 39, 39}.
	// EvalMod Q: []int{60, 60, 60, 60, 60, 60, 60, 60, 60}.
	// CoeffsToSlots Q: []int{56, 56, 56, 56}.
	// Precision : 23.8 bits for 2^{15} slots.
	// Failure : 2^{-138.7} for 2^{15} slots.
	N16QP1767H32768H32 = defaultParametersLiteral{
		ckks.ParametersLiteral{
			LogN:            16,
			LogQ:            []int{60, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40},
			LogP:            []int{61, 61, 61, 61, 61, 61},
			Xs:              ring.Ternary{H: 32768},
			LogDefaultScale: 40,
		},
		ParametersLiteral{},
	}

	// N16QP1788H32768H32 is a default bootstrapping parameters for a main secret with H=32768 and an ephemeral secret with H=32.
	// Residual Q : []int{60, 45, 45, 45, 45, 45, 45, 45, 45, 45} (465 bits).
	// SlotsToCoeffs Q: []int{42, 42, 42}.
	// EvalMod Q: []int{60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60}.
	// CoeffsToSlots Q: []int{58, 58, 58, 58}.
	// Precision : 29.8 bits for 2^{15} slots.
	// Failure : 2^{-138.7} for 2^{15} slots.
	N16QP1788H32768H32 = defaultParametersLiteral{
		ckks.ParametersLiteral{
			LogN:            16,
			LogQ:            []int{60, 45, 45, 45, 45, 45, 45, 45, 45, 45},
			LogP:            []int{61, 61, 61, 61, 61},
			Xs:              ring.Ternary{H: 32768},
			LogDefaultScale: 45,
		},
		ParametersLiteral{
			SlotsToCoeffsFactorizationDepthAndLogScales: [][]int{{42}, {42}, {42}},
			CoeffsToSlotsFactorizationDepthAndLogScales: [][]int{{58}, {58}, {58}, {58}},
			LogMessageRatio: utils.Pointy(2),
			Mod1InvDegree:   utils.Pointy(7),
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
			LogN:            16,
			LogQ:            []int{55, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 30},
			LogP:            []int{61, 61, 61, 61, 61},
			Xs:              ring.Ternary{H: 32768},
			LogDefaultScale: 30,
		},
		ParametersLiteral{
			SlotsToCoeffsFactorizationDepthAndLogScales: [][]int{{30}, {30, 30}},
			CoeffsToSlotsFactorizationDepthAndLogScales: [][]int{{53}, {53}, {53}, {53}},
			EvalModLogScale: utils.Pointy(55),
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
			LogN:            15,
			LogQ:            []int{40, 31, 31, 31, 31},
			LogP:            []int{56, 56},
			Xs:              ring.Ternary{H: 16384},
			LogDefaultScale: 31,
		},
		ParametersLiteral{
			SlotsToCoeffsFactorizationDepthAndLogScales: [][]int{{30, 30}},
			CoeffsToSlotsFactorizationDepthAndLogScales: [][]int{{52}, {52}},
			EvalModLogScale: utils.Pointy(55),
		},
	}
)
