package drlwe

import (
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

var (
	logN = 10
	qi   = []uint64{0x200000440001, 0x7fff80001, 0x800280001, 0x7ffd80001, 0x7ffc80001}
	pj   = []uint64{0x3ffffffb80001, 0x4000000800001}

	testBitDecomp16P1 = rlwe.ParametersLiteral{
		LogN:     logN,
		Q:        qi,
		Pow2Base: 16,
		P:        pj[:1],
		NTTFlag:  true,
	}

	testBitDecomp0P2 = rlwe.ParametersLiteral{
		LogN:    logN,
		Q:       qi,
		P:       pj,
		NTTFlag: true,
	}

	testParamsLiteral = []rlwe.ParametersLiteral{testBitDecomp16P1, testBitDecomp0P2}
)
