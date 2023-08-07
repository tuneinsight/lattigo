package integer

import (
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/bgv"
	"github.com/tuneinsight/lattigo/v4/circuits"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
)

type simIntegerPolynomialEvaluator struct {
	params             bgv.Parameters
	InvariantTensoring bool
}

func (d simIntegerPolynomialEvaluator) PolynomialDepth(degree int) int {
	if d.InvariantTensoring {
		return 0
	}
	return bits.Len64(uint64(degree)) - 1
}

// Rescale rescales the target circuits.SimOperand n times and returns it.
func (d simIntegerPolynomialEvaluator) Rescale(op0 *circuits.SimOperand) {
	if !d.InvariantTensoring {
		op0.Scale = op0.Scale.Div(rlwe.NewScale(d.params.Q()[op0.Level]))
		op0.Level--
	}
}

// Mul multiplies two circuits.SimOperand, stores the result the taret circuits.SimOperand and returns the result.
func (d simIntegerPolynomialEvaluator) MulNew(op0, op1 *circuits.SimOperand) (opOut *circuits.SimOperand) {
	opOut = new(circuits.SimOperand)
	opOut.Level = utils.Min(op0.Level, op1.Level)

	if d.InvariantTensoring {
		opOut.Scale = bgv.MulScaleInvariant(d.params, op0.Scale, op1.Scale, opOut.Level)
	} else {
		opOut.Scale = op0.Scale.Mul(op1.Scale)
	}

	return
}

func (d simIntegerPolynomialEvaluator) UpdateLevelAndScaleBabyStep(lead bool, tLevelOld int, tScaleOld rlwe.Scale) (tLevelNew int, tScaleNew rlwe.Scale) {
	tLevelNew = tLevelOld
	tScaleNew = tScaleOld
	if !d.InvariantTensoring && lead {
		tScaleNew = tScaleOld.Mul(d.params.NewScale(d.params.Q()[tLevelOld]))
	}
	return
}

func (d simIntegerPolynomialEvaluator) UpdateLevelAndScaleGiantStep(lead bool, tLevelOld int, tScaleOld, xPowScale rlwe.Scale) (tLevelNew int, tScaleNew rlwe.Scale) {

	Q := d.params.Q()

	tLevelNew = tLevelOld
	tScaleNew = tScaleOld.Div(xPowScale)

	// tScaleNew = targetScale*currentQi/XPow.Scale
	if !d.InvariantTensoring {

		var currentQi uint64
		if lead {
			currentQi = Q[tLevelNew]
		} else {
			currentQi = Q[tLevelNew+1]
		}

		tScaleNew = tScaleNew.Mul(d.params.NewScale(currentQi))

	} else {

		T := d.params.PlaintextModulus()

		// -Q mod T
		qModTNeg := new(big.Int).Mod(d.params.RingQ().ModulusAtLevel[tLevelNew], new(big.Int).SetUint64(T)).Uint64()
		qModTNeg = T - qModTNeg
		tScaleNew = tScaleNew.Mul(d.params.NewScale(qModTNeg))
	}

	if !d.InvariantTensoring {
		tLevelNew++
	}

	return
}
