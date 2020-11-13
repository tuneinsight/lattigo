package ring

import (
	"github.com/ldsec/lattigo/v2/utils"
)

const precision = uint64(56)

type baseSampler struct {
	prng     utils.PRNG
	baseRing *Ring
}
