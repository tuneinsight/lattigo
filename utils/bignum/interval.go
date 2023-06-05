package bignum

import(
	"math/big"
)

type Interval struct {
	Nodes int
	A, B *big.Float
}