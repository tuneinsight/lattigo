package rlwe

type TestParametersLiteral struct {
	BaseTwoDecomposition int
	ParametersLiteral
}

var (
	logN = 10
	qi   = []uint64{0x200000440001, 0x7fff80001, 0x800280001, 0x7ffd80001, 0x7ffc80001}
	pj   = []uint64{0x3ffffffb80001, 0x4000000800001}

	testParamsLiteral = []TestParametersLiteral{
		{
			BaseTwoDecomposition: 16,

			ParametersLiteral: ParametersLiteral{
				LogN:    logN,
				Q:       qi,
				P:       pj[:1],
				NTTFlag: true,
			},
		},

		{
			BaseTwoDecomposition: 0,

			ParametersLiteral: ParametersLiteral{
				LogN:    logN,
				Q:       qi,
				P:       pj,
				NTTFlag: true,
			},
		},
	}
)
