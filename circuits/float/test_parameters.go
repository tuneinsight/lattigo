package float

import (
	"github.com/tuneinsight/lattigo/v4/ckks"
)

var (
	TestPrec45 = ckks.ParametersLiteral{
		LogN:            10,
		LogQ:            []int{55, 45, 45, 45, 45, 45, 45},
		LogP:            []int{60},
		LogDefaultScale: 45,
	}

	TestPrec90 = ckks.ParametersLiteral{
		LogN:            10,
		LogQ:            []int{55, 55, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
		LogP:            []int{60, 60},
		LogDefaultScale: 90,
	}

	TestParametersLiteral = []ckks.ParametersLiteral{TestPrec45, TestPrec90}
)
