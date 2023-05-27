package bgv

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// LinearTransformEncoder is a struct complying
// to the rlwe.LinearTransformEncoder interface.
type LinearTransformEncoder struct {
	*Encoder
	buf       *ring.Poly
	values    []uint64
	diagonals map[int][]uint64
}

// NewLinearTransformEncoder creates a new LinearTransformEncoder, which implements the rlwe.LinearTransformEncoder interface.
// See rlwe.LinearTransformEncoder for additional informations.
func NewLinearTransformEncoder(ecd *Encoder, diagonals map[int][]uint64) rlwe.LinearTransformEncoder {
	return LinearTransformEncoder{
		Encoder:   ecd,
		buf:       ecd.Parameters().RingT().NewPoly(),
		values:    make([]uint64, ecd.Parameters().N()),
		diagonals: diagonals,
	}
}

// Parameters returns the rlwe.Parametrs of the underlying LinearTransformEncoder.
func (l LinearTransformEncoder) Parameters() rlwe.Parameters {
	return l.Encoder.Parameters().Parameters
}

// NonZeroDiagonals retursn the list of non-zero diagonales of the matrix
// representing the linear transformation.
func (l LinearTransformEncoder) NonZeroDiagonals() []int {
	return utils.GetKeys(l.diagonals)
}

// EncodeLinearTransformDiagonal encodes the i-th non-zero diagonal  of size at most 2^{LogSlots} rotated by `rot` positions
// to the left of the internaly stored matrix at the given Scale on the outut ringqp.Poly.
func (l LinearTransformEncoder) EncodeLinearTransformDiagonal(i, rot int, scale rlwe.Scale, logSlots int, output ringqp.Poly) (err error) {

	ecd := l.Encoder
	buf := l.buf
	slots := 1 << logSlots

	levelQ, levelP := output.LevelQ(), output.LevelP()
	ringQP := ecd.Parameters().RingQP().AtLevel(levelQ, levelP)

	rot &= (slots - 1)

	// manages inputs that have rotation between 0 and slots-1 or between -slots/2 and slots/2-1
	v, ok := l.diagonals[i]
	if !ok {
		if v, ok = l.diagonals[i-slots]; !ok {
			return fmt.Errorf("cannot EncodeLinearTransformDiagonalNaive: diagonal [%d] doesn't exist", i)
		}
	}

	var values []uint64
	if rot != 0 {

		values = l.values

		if len(v) > slots {
			utils.RotateSliceAllocFree(v[slots:], rot, values[slots:])
		}

		utils.RotateSliceAllocFree(v[:slots], rot, values[:slots])

	} else {
		values = v
	}

	l.EncodeRingT(values, scale, buf)

	l.RingT2Q(levelQ, false, buf, output.Q)
	l.RingT2Q(levelP, false, buf, output.P)

	ringQP.NTT(&output, &output)
	ringQP.MForm(&output, &output)

	return
}
