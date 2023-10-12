package float

import (
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/circuits"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

type simEvaluator struct {
	params                     ckks.Parameters
	levelsConsumedPerRescaling int
}

func (d simEvaluator) PolynomialDepth(degree int) int {
	return d.levelsConsumedPerRescaling * (bits.Len64(uint64(degree)) - 1)
}

// Rescale rescales the target circuits.SimOperand n times and returns it.
func (d simEvaluator) Rescale(op0 *circuits.SimOperand) {
	for i := 0; i < d.levelsConsumedPerRescaling; i++ {
		op0.Scale = op0.Scale.Div(rlwe.NewScale(d.params.Q()[op0.Level]))
		op0.Level--
	}
}

// Mul multiplies two circuits.SimOperand, stores the result the target circuits.SimOperand and returns the result.
func (d simEvaluator) MulNew(op0, op1 *circuits.SimOperand) (opOut *circuits.SimOperand) {
	opOut = new(circuits.SimOperand)
	opOut.Level = utils.Min(op0.Level, op1.Level)
	opOut.Scale = op0.Scale.Mul(op1.Scale)
	return
}

func (d simEvaluator) UpdateLevelAndScaleBabyStep(lead bool, tLevelOld int, tScaleOld rlwe.Scale) (tLevelNew int, tScaleNew rlwe.Scale) {

	tLevelNew = tLevelOld
	tScaleNew = tScaleOld

	if lead {
		for i := 0; i < d.levelsConsumedPerRescaling; i++ {
			tScaleNew = tScaleNew.Mul(rlwe.NewScale(d.params.Q()[tLevelNew-i]))
		}
	}

	return
}

func (d simEvaluator) UpdateLevelAndScaleGiantStep(lead bool, tLevelOld int, tScaleOld, xPowScale rlwe.Scale) (tLevelNew int, tScaleNew rlwe.Scale) {

	Q := d.params.Q()

	var qi *big.Int
	if lead {
		qi = bignum.NewInt(Q[tLevelOld])
		for i := 1; i < d.levelsConsumedPerRescaling; i++ {
			qi.Mul(qi, bignum.NewInt(Q[tLevelOld-i]))
		}
	} else {
		qi = bignum.NewInt(Q[tLevelOld+d.levelsConsumedPerRescaling])
		for i := 1; i < d.levelsConsumedPerRescaling; i++ {
			qi.Mul(qi, bignum.NewInt(Q[tLevelOld+d.levelsConsumedPerRescaling-i]))
		}
	}

	tLevelNew = tLevelOld + d.levelsConsumedPerRescaling
	tScaleNew = tScaleOld.Mul(rlwe.NewScale(qi))
	tScaleNew = tScaleNew.Div(xPowScale)

	return
}

func (d simEvaluator) GetPolynmialDepth(degree int) int {
	return d.levelsConsumedPerRescaling * (bits.Len64(uint64(degree)) - 1)
}
