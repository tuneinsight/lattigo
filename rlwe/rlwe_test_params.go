package rlwe

var (
	// TestPN12QP109 is a set of default parameters with logN=12 and logQP=109
	TestPN12QP109 = ParametersLiteral{
		LogN:  12,
		Q:     []uint64{0x7ffffec001, 0x40002001}, // 39 + 39 bits
		P:     []uint64{0x8000016001},             // 30 bits
		Sigma: DefaultSigma,
	}
	// TestPN13QP218 is a set of default parameters with logN=13 and logQP=218
	TestPN13QP218 = ParametersLiteral{
		LogN:  13,
		Q:     []uint64{0x3fffffffef8001, 0x4000000011c001, 0x40000000120001}, // 54 + 54 + 54 bits
		P:     []uint64{0x7ffffffffb4001},                                     // 55 bits
		Sigma: DefaultSigma,
	}

	// TestPN14QP438 is a set of default parameters with logN=14 and logQP=438
	TestPN14QP438 = ParametersLiteral{
		LogN: 14,
		Q: []uint64{0x100000000060001, 0x80000000068001, 0x80000000080001,
			0x3fffffffef8001, 0x40000000120001, 0x3fffffffeb8001}, // 56 + 55 + 55 + 54 + 54 + 54 bits
		P:     []uint64{0x80000000130001, 0x7fffffffe90001}, // 55 + 55 bits
		Sigma: DefaultSigma,
	}

	// TestPN15QP880 is a set of default parameters with logN=15 and logQP=880
	TestPN15QP880 = ParametersLiteral{
		LogN: 15,
		Q: []uint64{0x7ffffffffe70001, 0x7ffffffffe10001, 0x7ffffffffcc0001, // 59 + 59 + 59 bits
			0x400000000270001, 0x400000000350001, 0x400000000360001, // 58 + 58 + 58 bits
			0x3ffffffffc10001, 0x3ffffffffbe0001, 0x3ffffffffbd0001, // 58 + 58 + 58 bits
			0x4000000004d0001, 0x400000000570001, 0x400000000660001}, // 58 + 58 + 58 bits
		P:     []uint64{0xffffffffffc0001, 0x10000000001d0001, 0x10000000006e0001}, // 60 + 60 + 60 bits
		Sigma: DefaultSigma,
	}
	// TestPN16QP240 is a set of default parameters with logN=16 and logQP=240
	TestPN16QP240 = ParametersLiteral{
		LogN:  16,
		LogQ:  []int{60, 60, 60}, // 58 + 58 + 58 bits
		LogP:  []int{60},         // 60 + 60 + 60 bits
		Sigma: DefaultSigma,
	}

	// TestPN17QP360 is a set of default parameters with logN=17 and logQP=360
	TestPN17QP360 = ParametersLiteral{
		LogN:  17,
		LogQ:  []int{60, 60, 60, 60},
		LogP:  []int{60, 60},
		Sigma: DefaultSigma,
	}
)
