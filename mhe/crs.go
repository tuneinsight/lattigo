package mhe

import (
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
)

// CRS is an interface for Common Reference Strings.
// CRSs are PRNGs for which the read bits are the same for
// all parties.
type CRS interface {
	sampling.PRNG
}
