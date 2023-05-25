package polynomial

import (
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

type PolynomialBSGS struct {
	MetaData
	Coeffs [][]*bignum.Complex
}

func OptimalSplit(logDegree int) (logSplit int) {
	logSplit = logDegree >> 1
	a := (1 << logSplit) + (1 << (logDegree - logSplit)) + logDegree - logSplit - 3
	b := (1 << (logSplit + 1)) + (1 << (logDegree - logSplit - 1)) + logDegree - logSplit - 4
	if a > b {
		logSplit++
	}

	return
}
