// Package integer implements advanced homomorphic circuit for encrypted arithmetic modular arithmetic with integers.
package integer

import (
	"github.com/tuneinsight/lattigo/v4/bgv"
)

type Integer interface {
	bgv.Integer
}
