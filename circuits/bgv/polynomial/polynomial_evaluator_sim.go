package polynomial

import (
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v6/circuits/common/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/bgv"
	"github.com/tuneinsight/lattigo/v6/utils"
)

// simEvaluator is a struct used to pre-computed the scaling
// factors of the polynomial coefficients used by the inlined
// polynomial evaluation by running the polynomial evaluation
// with dummy operands.
// This struct implements the interface [polynomial.SimEvaluator].
type simEvaluator struct {
	params             bgv.Parameters
	InvariantTensoring bool
}

// PolynomialDepth returns the depth of the polynomial.
func (d simEvaluator) PolynomialDepth(degree int) int {
	if d.InvariantTensoring {
		return 0
	}
	return bits.Len64(uint64(degree)) - 1
}

// Rescale rescales the target polynomial.SimOperand n times and returns it.
func (d simEvaluator) Rescale(op0 *polynomial.SimOperand) {
	if !d.InvariantTensoring {
		op0.Scale = op0.Scale.Div(rlwe.NewScale(d.params.Q()[op0.Level]))
		op0.Level--
	}
}

// MulNew multiplies two polynomial.SimOperand, stores the result the target [polynomial.SimOperand] and returns the result.
func (d simEvaluator) MulNew(op0, op1 *polynomial.SimOperand) (opOut *polynomial.SimOperand) {
	opOut = new(polynomial.SimOperand)
	opOut.Level = utils.Min(op0.Level, op1.Level)

	if d.InvariantTensoring {
		opOut.Scale = bgv.MulScaleInvariant(d.params, op0.Scale, op1.Scale, opOut.Level)
	} else {
		opOut.Scale = op0.Scale.Mul(op1.Scale)
	}

	return
}

// UpdateLevelAndScaleBabyStep returns the updated level and scale for a baby-step.
func (d simEvaluator) UpdateLevelAndScaleBabyStep(lead bool, tLevelOld int, tScaleOld rlwe.Scale) (tLevelNew int, tScaleNew rlwe.Scale) {
	tLevelNew = tLevelOld
	tScaleNew = tScaleOld
	if !d.InvariantTensoring && lead {
		tScaleNew = tScaleOld.Mul(d.params.NewScale(d.params.Q()[tLevelOld]))
	}
	return
}

// UpdateLevelAndScaleGiantStep returns the updated level and scale for a giant-step.
func (d simEvaluator) UpdateLevelAndScaleGiantStep(lead bool, tLevelOld int, tScaleOld, xPowScale rlwe.Scale) (tLevelNew int, tScaleNew rlwe.Scale) {

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
