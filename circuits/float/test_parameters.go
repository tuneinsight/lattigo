package float

import (
	"github.com/tuneinsight/lattigo/v4/ckks"
)

var (
	testPrec45 = ckks.ParametersLiteral{
		LogN:            10,
		LogQ:            []int{55, 45, 45, 45, 45, 45, 45},
		LogP:            []int{60},
		LogDefaultScale: 45,
	}

	testPrec90 = ckks.ParametersLiteral{
		LogN:            10,
		LogQ:            []int{55, 55, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
		LogP:            []int{60, 60},
		LogDefaultScale: 90,
	}

	testParametersLiteralFloat = []ckks.ParametersLiteral{testPrec45, testPrec90}
)
