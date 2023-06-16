package ckks

var (
	// ExampleParameters128BitLogN14LogQP438 is an example parameters set with logN=14, logQP=435
	// offering 128-bit of security.
	ExampleParameters128BitLogN14LogQP438 = ParametersLiteral{
		LogN: 14,
		Q: []uint64{
			0x80000000080001, // 55
			0x2000000a0001, // 45
			0x2000000e0001, // 45
			0x2000001d0001, // 45
			0x1fffffcf0001, // 45
			0x1fffffc20001, // 45
			0x200000440001, // 45
		},
		P: []uint64{
			0x80000000130001, // 55
			0x7fffffffe90001, // 55
		},
		LogPlaintextScale: 45,
	}

)
