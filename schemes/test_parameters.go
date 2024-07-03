package schemes

import (
	"github.com/tuneinsight/lattigo/v5/schemes/bgv"
	"github.com/tuneinsight/lattigo/v5/schemes/ckks"
)

var (
	// BgvTestInsecure are insecure parameters used for the sole purpose of fast testing.
	BgvTestInsecure = bgv.ParametersLiteral{
		LogN: 10,
		Q:    []uint64{0x3fffffa8001, 0x1000090001, 0x10000c8001, 0x10000f0001, 0xffff00001},
		P:    []uint64{0x7fffffd8001},
	}

	BgvTestPlaintextModulus = []uint64{0x101, 0xffc001}

	BgvTestParams = []bgv.ParametersLiteral{BgvTestInsecure}
)

var (
	// BgvTestInsecurePrec45 are insecure parameters used for the sole purpose of fast testing.
	CkksTestInsecurePrec45 = ckks.ParametersLiteral{
		LogN:            10,
		LogQ:            []int{55, 45, 45, 45, 45, 45, 45},
		LogP:            []int{60},
		LogDefaultScale: 45,
	}

	// BgvTestInsecurePrec90 are insecure parameters used for the sole purpose of fast testing.
	CkksTestInsecurePrec90 = ckks.ParametersLiteral{
		LogN:            10,
		LogQ:            []int{55, 55, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
		LogP:            []int{60, 60},
		LogDefaultScale: 90,
	}

	CkksTestParametersLiteral = []ckks.ParametersLiteral{CkksTestInsecurePrec45, CkksTestInsecurePrec90}
)
