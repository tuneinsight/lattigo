// Package integer implements Homomorphic Encryption for encrypted modular arithmetic over the integers.
package integer

import (
	"github.com/tuneinsight/lattigo/v4/bgv"
)

type Integer interface {
	bgv.Integer
}
