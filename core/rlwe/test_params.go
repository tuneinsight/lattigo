package rlwe

type TestParametersLiteral struct {
	BaseTwoDecomposition int
	Compressed           bool
	ParametersLiteral
}

var (
	logN = 10
	qi   = []uint64{0x200000440001, 0x7fff80001, 0x800280001, 0x7ffd80001, 0x7ffc80001}
	pj   = []uint64{0x3ffffffb80001, 0x4000000800001}

	// testInsecure are insecure parameters used for the sole purpose of fast testing.
	testInsecure = []TestParametersLiteral{
		// RNS decomposition, no Pw2 decomposition
		{
			BaseTwoDecomposition: 0,

			ParametersLiteral: ParametersLiteral{
				LogN:    logN,
				Q:       qi,
				P:       pj,
				NTTFlag: true,
			},
		},
		// RNS decomposition, Pw2 decomposition
		{
			BaseTwoDecomposition: 16,

			ParametersLiteral: ParametersLiteral{
				LogN:    logN,
				Q:       qi,
				P:       pj[:1],
				NTTFlag: true,
			},
		},
		// No RNS decomposition, Pw2 decomposition
		{
			BaseTwoDecomposition: 2,

			ParametersLiteral: ParametersLiteral{
				LogN:    logN,
				Q:       qi,
				P:       nil,
				NTTFlag: true,
			},
		},
	}
)
