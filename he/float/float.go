// Package float implements Homomorphic Encryption for encrypted arithmetic over floating point numbers.
package float

import (
	"github.com/tuneinsight/lattigo/v4/ckks"
)

type Float interface {
	ckks.Float
}
