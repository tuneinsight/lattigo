package circuits

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v4/bgv"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// EncoderInterface defines a set of common and scheme agnostic method provided by an Encoder struct.
type EncoderInterface[T Numeric, U *ring.Poly | ringqp.Poly | *rlwe.Plaintext] interface {
	Encode(values []T, metaData *rlwe.MetaData, output U) (err error)
}

// EncodeIntegerLinearTransformation encodes a linear transformation on a pre-allocated linear transformation.
// The method will return an error if the non-zero diagonals between the pre-allocated linear transformation and the parameters of the linear transformation to encode do not match.
func EncodeIntegerLinearTransformation[T int64 | uint64](params LinearTransformationParameters, ecd *bgv.Encoder, diagonals Diagonals[T], allocated LinearTransformation) (err error) {
	return EncodeLinearTransformation[T](params, &intEncoder[T, ringqp.Poly]{ecd}, diagonals, allocated)
}

type intEncoder[T int64 | uint64, U ring.Poly | ringqp.Poly | *rlwe.Plaintext] struct {
	*bgv.Encoder
}

func (e intEncoder[T, U]) Encode(values []T, metadata *rlwe.MetaData, output U) (err error) {
	return e.Embed(values, false, metadata, output)
}

// EncodeFloatLinearTransformation encodes a linear transformation on a pre-allocated linear transformation.
// The method will return an error if the non-zero diagonals between the pre-allocated linear transformation and the parameters of the linear transformation to encode do not match.
func EncodeFloatLinearTransformation[T Float](params LinearTransformationParameters, ecd *ckks.Encoder, diagonals Diagonals[T], allocated LinearTransformation) (err error) {
	return EncodeLinearTransformation[T](params, &floatEncoder[T, ringqp.Poly]{ecd}, diagonals, allocated)
}

type floatEncoder[T float64 | complex128 | *big.Float | *bignum.Complex, U ring.Poly | ringqp.Poly | *rlwe.Plaintext] struct {
	*ckks.Encoder
}

func (e *floatEncoder[T, U]) Encode(values []T, metadata *rlwe.MetaData, output U) (err error) {
	return e.Encoder.Embed(values, metadata, output)
}
