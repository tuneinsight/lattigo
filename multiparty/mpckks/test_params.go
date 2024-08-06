package mpckks

import "github.com/tuneinsight/lattigo/v6/schemes/ckks"

var (

	// testInsecurePrec45 are insecure parameters used for the sole purpose of fast testing.
	testInsecurePrec45 = ckks.ParametersLiteral{
		LogN: 10,
		Q: []uint64{
			0x80000000080001,
			0x2000000a0001,
			0x2000000e0001,
			0x2000001d0001,
			0x1fffffcf0001,
			0x1fffffc20001,
			0x200000440001,
		},
		P: []uint64{
			0x80000000130001,
			0x7fffffffe90001,
		},
		LogDefaultScale: 45,
	}

	// testInsecurePrec90 are insecure parameters used for the sole purpose of fast testing.
	testInsecurePrec90 = ckks.ParametersLiteral{
		LogN: 10,
		Q: []uint64{
			0x80000000080001,
			0x80000000440001,
			0x2000000a0001,
			0x2000000e0001,
			0x1fffffc20001,
			0x200000440001,
			0x200000500001,
			0x200000620001,
			0x1fffff980001,
			0x2000006a0001,
			0x1fffff7e0001,
			0x200000860001,
		},
		P: []uint64{
			0xffffffffffc0001,
			0x10000000006e0001,
		},
		LogDefaultScale: 90,
	}

	testParamsLiteral = []ckks.ParametersLiteral{testInsecurePrec45, testInsecurePrec90}
)
