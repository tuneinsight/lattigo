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

func NewLinearTransformEncoder(ecd *Encoder, diagonals map[int][]uint64) rlwe.LinearTransformEncoder {
	return LinearTransformEncoder{
		Encoder:   ecd,
		buf:       ecd.Parameters().RingT().NewPoly(),
		values:    make([]uint64, ecd.Parameters().N()),
		diagonals: diagonals,
	}
}

func (l LinearTransformEncoder) Parameters() rlwe.Parameters {
	return l.Encoder.Parameters().Parameters
}

// NonZeroDiagonals returns the list of non-zero diagonals.
func (l LinearTransformEncoder) NonZeroDiagonals() []int {
	return utils.GetKeys(l.diagonals)
}

// EncodeLinearTransformDiagonalNaive encodes the i-th non-zero diagonal of the internaly stored matrix at the given scale on the outut polynomial.
func (l LinearTransformEncoder) EncodeLinearTransformDiagonalNaive(i int, scale rlwe.Scale, logslots int, output ringqp.Poly) (err error) {

	ecd := l.Encoder
	buf := l.buf
	levelQ, levelP := output.LevelQ(), output.LevelP()
	ringQP := ecd.Parameters().RingQP().AtLevel(levelQ, levelP)

	if diag, ok := l.diagonals[i]; ok {
		l.EncodeRingT(diag, scale, buf)
		l.RingT2Q(levelQ, false, buf, output.Q)
		l.RingT2Q(levelP, false, buf, output.P)
		ringQP.NTT(&output, &output)
		ringQP.MForm(&output, &output)
	} else {
		return fmt.Errorf("cannot EncodeLinearTransformDiagonalNaive: diagonal [%d] doesn't exist", i)
	}

	return
}

// EncodeLinearTransformDiagonalBSGS encodes the i-th non-zero diagonal of the internaly stored matrix at the given scale on the outut polynomial.
func (l LinearTransformEncoder) EncodeLinearTransformDiagonalBSGS(i, rot int, scale rlwe.Scale, logSlots int, output ringqp.Poly) (err error) {

	ecd := l.Encoder
	buf := l.buf
	slots := 1 << logSlots
	values := l.values
	levelQ, levelP := output.LevelQ(), output.LevelP()
	ringQP := ecd.Parameters().RingQP().AtLevel(levelQ, levelP)

	// manages inputs that have rotation between 0 and slots-1 or between -slots/2 and slots/2-1
	v, ok := l.diagonals[i]
	if !ok {
		v = l.diagonals[i-slots]
	}

	if len(v) > slots {
		rotateAndCopyInplace(values[slots:], v[slots:], rot)
	}

	rotateAndCopyInplace(values[:slots], v, rot)

	l.EncodeRingT(values, scale, buf)

	l.RingT2Q(levelQ, false, buf, output.Q)
	l.RingT2Q(levelP, false, buf, output.P)

	ringQP.NTT(&output, &output)
	ringQP.MForm(&output, &output)

	return
}

func rotateAndCopyInplace(values, v []uint64, rot int) {
	n := len(values)
	if len(v) > rot {
		copy(values[:n-rot], v[rot:])
		copy(values[n-rot:], v[:rot])
	} else {
		copy(values[n-rot:], v)
	}
}
