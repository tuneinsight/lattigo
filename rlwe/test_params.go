package rlwe

var (
	LogN = 13
	Q    = []uint64{0x200000440001, 0x7fff80001, 0x800280001, 0x7ffd80001, 0x7ffc80001}
	P    = []uint64{0x3ffffffb80001, 0x4000000800001}

	TESTBITDECOMP16P1 = ParametersLiteral{
		LogN:     LogN,
		Q:        Q,
		Pow2Base: 16,
		P:        P[:1],
		NTTFlag:  true,
	}

	TESTBITDECOMP0P2 = ParametersLiteral{
		LogN:    LogN,
		Q:       Q,
		P:       P,
		NTTFlag: true,
	}

	TestParamsLiteral = []ParametersLiteral{TESTBITDECOMP16P1, TESTBITDECOMP0P2}
)
