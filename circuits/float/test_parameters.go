package float

import (
	"github.com/tuneinsight/lattigo/v4/ckks"
)

var (
	testPrec45 = ckks.ParametersLiteral{
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

	testPrec90 = ckks.ParametersLiteral{
		LogN: 10,
		Q: []uint64{
			0x80000000080001,
			0x80000000440001,
			0x1fffffff9001, // 44.99999999882438
			0x200000008001, // 45.00000000134366
			0x1ffffffe7001, // 44.99999999580125
			0x20000001c001, // 45.00000000470269
			0x1ffffffe1001, // 44.99999999479353
			0x1ffffffce001, // 44.99999999160245
			0x200000041001, // 45.00000001091691
			0x200000046001, // 45.00000001175667
			0x200000053001, // 45.00000001394004
			0x1ffffffab001, // 44.99999998572414
			0x1ffffffa7001, // 44.99999998505233
			0x1ffffffa2001, // 44.99999998421257
		},
		P: []uint64{
			0xffffffffffc0001,
			0x10000000006e0001,
		},
		LogDefaultScale: 90,
	}

	testParametersLiteralFloat = []ckks.ParametersLiteral{testPrec45, testPrec90}
)
