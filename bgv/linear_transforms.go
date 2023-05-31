package bgv

import (
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
)

// NewLinearTransform allocates a new LinearTransform with zero plaintexts at the specified level.
// If LogBSGSRatio < 0, the LinearTransform is set to not use the BSGS approach.
func NewLinearTransform(params Parameters, nonZeroDiags []int, level int, scale rlwe.Scale, LogBSGSRatio int) rlwe.LinearTransform {
	return rlwe.NewLinearTransform(params, nonZeroDiags, level, scale, params.MaxLogSlots(), LogBSGSRatio)
}

func EncodeLinearTransform[T int64 | uint64](LT rlwe.LinearTransform, diagonals map[int][]T, ecd *Encoder) (err error) {
	return rlwe.EncodeLinearTransform[T](LT, diagonals, &encoder[T, ringqp.Poly]{ecd})
}

func GenLinearTransform[T int64 | uint64](diagonals map[int][]T, ecd *Encoder, level int, scale rlwe.Scale) (LT rlwe.LinearTransform, err error) {
	return rlwe.GenLinearTransform[T](diagonals, &encoder[T, ringqp.Poly]{ecd}, level, scale, ecd.Parameters().MaxLogSlots())
}

func GenLinearTransformBSGS[T int64 | uint64](diagonals map[int][]T, ecd *Encoder, level int, scale rlwe.Scale, LogBSGSRatio int) (LT rlwe.LinearTransform, err error) {
	return rlwe.GenLinearTransformBSGS[T](diagonals, &encoder[T, ringqp.Poly]{ecd}, level, scale, ecd.Parameters().MaxLogSlots(), LogBSGSRatio)
}
