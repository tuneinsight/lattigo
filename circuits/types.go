package circuits

import (
	"github.com/tuneinsight/lattigo/v4/bgv"
	"github.com/tuneinsight/lattigo/v4/ckks"
)

type Numeric interface {
	ckks.Float | bgv.Integer
}

type Float interface {
	ckks.Float
}

type Integer interface {
	bgv.Integer
}
