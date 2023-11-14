package examples

import (
	"github.com/tuneinsight/lattigo/v4/he/hefloat"
	"github.com/tuneinsight/lattigo/v4/he/heint"
	"github.com/tuneinsight/lattigo/v4/ring"
)

var (

	// HEIntParamsN12QP109 is an example parameter set for the `heint` package logN=12 and logQP=109
	HEIntParamsN12QP109 = heint.ParametersLiteral{
		LogN:             12,
		Q:                []uint64{0x7ffffec001, 0x8000016001}, // 39 + 39 bits
		P:                []uint64{0x40002001},                 // 30 bits
		PlaintextModulus: 65537,
	}

	// HEIntParamsN13QP218 is an example parameter set for the `heint` package with logN=13 and logQP=218
	HEIntParamsN13QP218 = heint.ParametersLiteral{
		LogN:             13,
		Q:                []uint64{0x3fffffffef8001, 0x4000000011c001, 0x40000000120001}, // 54 + 54 + 54 bits
		P:                []uint64{0x7ffffffffb4001},                                     // 55 bits
		PlaintextModulus: 65537,
	}

	// HEIntParamsN14QP438 is an example parameter set for the `heint` package with logN=14 and logQP=438
	HEIntParamsN14QP438 = heint.ParametersLiteral{
		LogN: 14,
		Q: []uint64{0x100000000060001, 0x80000000068001, 0x80000000080001,
			0x3fffffffef8001, 0x40000000120001, 0x3fffffffeb8001}, // 56 + 55 + 55 + 54 + 54 + 54 bits
		P:                []uint64{0x80000000130001, 0x7fffffffe90001}, // 55 + 55 bits
		PlaintextModulus: 65537,
	}

	// HEIntParamsN15QP880 is an example parameter set for the `heint` package with logN=15 and logQP=880
	HEIntParamsN15QP880 = heint.ParametersLiteral{
		LogN: 15,
		Q: []uint64{0x7ffffffffe70001, 0x7ffffffffe10001, 0x7ffffffffcc0001, // 59 + 59 + 59 bits
			0x400000000270001, 0x400000000350001, 0x400000000360001, // 58 + 58 + 58 bits
			0x3ffffffffc10001, 0x3ffffffffbe0001, 0x3ffffffffbd0001, // 58 + 58 + 58 bits
			0x4000000004d0001, 0x400000000570001, 0x400000000660001}, // 58 + 58 + 58 bits
		P:                []uint64{0xffffffffffc0001, 0x10000000001d0001, 0x10000000006e0001}, // 60 + 60 + 60 bits
		PlaintextModulus: 65537,
	}

	// HEFloatParamsN12QP109 is an example parameter set for the `hefloat` package with logN=12 and logQP=109
	HEFloatParamsN12QP109 = hefloat.ParametersLiteral{
		LogN:            12,
		Q:               []uint64{0x200000e001, 0x100006001}, // 37 + 32},
		P:               []uint64{0x3ffffea001},              // 38
		LogDefaultScale: 32,
	}

	// HEFloatParamsN13QP218 is an example parameter set for the `hefloat` package with logN=13 and logQP=218
	HEFloatParamsN13QP218 = hefloat.ParametersLiteral{
		LogN: 13,
		Q: []uint64{0x1fffec001, // 33 + 5 x 30
			0x3fff4001,
			0x3ffe8001,
			0x40020001,
			0x40038001,
			0x3ffc0001},
		P:               []uint64{0x800004001}, // 35
		LogDefaultScale: 30,
	}
	// HEFloatParamsN14QP438 is an example parameter set for the `hefloat` package with logN=14 and logQP=438
	HEFloatParamsN14QP438 = hefloat.ParametersLiteral{
		LogN: 14,
		Q: []uint64{0x200000008001, 0x400018001, // 45 + 9 x 34
			0x3fffd0001, 0x400060001,
			0x400068001, 0x3fff90001,
			0x400080001, 0x4000a8001,
			0x400108001, 0x3ffeb8001},
		P:               []uint64{0x7fffffd8001, 0x7fffffc8001}, // 43, 43
		LogDefaultScale: 34,
	}

	// HEFloatParamsN15QP880 is an example parameter set for the `hefloat` package with logN=15 and logQP=880
	HEFloatParamsN15QP880 = hefloat.ParametersLiteral{
		LogN: 15,
		Q: []uint64{0x4000000120001, 0x10000140001, 0xffffe80001, // 50 + 17 x 40
			0x10000290001, 0xffffc40001, 0x100003e0001,
			0x10000470001, 0x100004b0001, 0xffffb20001,
			0x10000500001, 0x10000650001, 0xffff940001,
			0xffff8a0001, 0xffff820001, 0xffff780001,
			0x10000890001, 0xffff750001, 0x10000960001},
		P:               []uint64{0x40000001b0001, 0x3ffffffdf0001, 0x4000000270001}, // 50, 50, 50
		LogDefaultScale: 40,
	}
	// HEFloatParamsPN16QP1761 is an example parameter set for the `hefloat` package with logN=16 and logQP = 1761
	HEFloatParamsPN16QP1761 = hefloat.ParametersLiteral{
		LogN: 16,
		Q: []uint64{0x80000000080001, 0x2000000a0001, 0x2000000e0001, 0x1fffffc20001, // 55 + 33 x 45
			0x200000440001, 0x200000500001, 0x200000620001, 0x1fffff980001,
			0x2000006a0001, 0x1fffff7e0001, 0x200000860001, 0x200000a60001,
			0x200000aa0001, 0x200000b20001, 0x200000c80001, 0x1fffff360001,
			0x200000e20001, 0x1fffff060001, 0x200000fe0001, 0x1ffffede0001,
			0x1ffffeca0001, 0x1ffffeb40001, 0x200001520001, 0x1ffffe760001,
			0x2000019a0001, 0x1ffffe640001, 0x200001a00001, 0x1ffffe520001,
			0x200001e80001, 0x1ffffe0c0001, 0x1ffffdee0001, 0x200002480001,
			0x1ffffdb60001, 0x200002560001},
		P:               []uint64{0x80000000440001, 0x7fffffffba0001, 0x80000000500001, 0x7fffffffaa0001}, // 4 x 55
		LogDefaultScale: 45,
	}

	// HEFloatCIParamsN12QP109 is an example parameter set for the `hefloat` package with conjugate-invariant CKKS and logN=12 and logQP=109
	HEFloatCIParamsN12QP109 = hefloat.ParametersLiteral{
		LogN:            12,
		Q:               []uint64{0x1ffffe0001, 0x100014001}, // 37 + 32
		P:               []uint64{0x4000038001},              // 38
		RingType:        ring.ConjugateInvariant,
		LogDefaultScale: 32,
	}

	// HEFloatCIParamsN13QP218 is an example parameter set for the `hefloat` package with conjugate-invariant CKKS and logN=13 and logQP=218
	HEFloatCIParamsN13QP218 = hefloat.ParametersLiteral{
		LogN: 13,
		Q: []uint64{0x200038001, // 33 + 5 x 30
			0x3ffe8001,
			0x40020001,
			0x40038001,
			0x3ffc0001,
			0x40080001},
		P:               []uint64{0x800008001}, // 35
		RingType:        ring.ConjugateInvariant,
		LogDefaultScale: 30,
	}
)

var HEIntParams = []heint.ParametersLiteral{HEIntParamsN12QP109, HEIntParamsN13QP218, HEIntParamsN14QP438, HEIntParamsN15QP880}

var HEFloatParams = []hefloat.ParametersLiteral{HEFloatParamsN12QP109, HEFloatParamsN13QP218, HEFloatParamsN14QP438, HEFloatParamsN15QP880, HEFloatParamsPN16QP1761, HEFloatCIParamsN12QP109, HEFloatCIParamsN13QP218}
