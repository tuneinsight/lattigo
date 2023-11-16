package bignum

import (
	"math/big"
)

// Interval is a struct storing information about interval
// for a polynomial approximation.
// Nodes: the number of points used for the interpolation.
// [A, B]: the domain of the interpolation.
type Interval struct {
	Nodes int
	A, B  big.Float
}
