package bgv

import (
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
)

// NewLinearTransformation allocates a new LinearTransformation with zero values according to the parameters specified by the LinearTranfromationParameters.
func NewLinearTransformation[T int64 | uint64](params rlwe.ParametersInterface, lt rlwe.LinearTranfromationParameters[T]) rlwe.LinearTransformation {
	return rlwe.NewLinearTransformation(params, lt)
}

// EncodeLinearTransformation encodes a linear transformation on a pre-allocated linear transformation.
// The method will return an error if the non-zero diagonals between the pre-allocated linear transformation and the parameters of the linear transformation to encode do not match.
func EncodeLinearTransformation[T int64 | uint64](allocated rlwe.LinearTransformation, params rlwe.LinearTranfromationParameters[T], ecd *Encoder) (err error) {
	return rlwe.EncodeLinearTransformation[T](allocated, params, &encoder[T, ringqp.Poly]{ecd})
}
