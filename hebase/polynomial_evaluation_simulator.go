package hebase

import (
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// DummyOperand is a dummy operand
// that only stores the level and the scale.
type DummyOperand struct {
	Level          int
	PlaintextScale rlwe.Scale
}

type DummyEvaluator interface {
	MulNew(op0, op1 *DummyOperand) *DummyOperand
	Rescale(op0 *DummyOperand)
	PolynomialDepth(degree int) int
	UpdateLevelAndScaleGiantStep(lead bool, tLevelOld int, tScaleOld, xPowScale rlwe.Scale) (tLevelNew int, tScaleNew rlwe.Scale)
	UpdateLevelAndScaleBabyStep(lead bool, tLevelOld int, tScaleOld rlwe.Scale) (tLevelNew int, tScaleNew rlwe.Scale)
}

// DummyPowerBasis is a map storing powers of DummyOperands indexed by their power.
type DummyPowerBasis map[int]*DummyOperand

// GenPower populates the target DummyPowerBasis with the nth power.
func (d DummyPowerBasis) GenPower(params rlwe.GetRLWEParameters, n int, eval DummyEvaluator) {

	if n < 2 {
		return
	}

	a, b := SplitDegree(n)

	d.GenPower(params, a, eval)
	d.GenPower(params, b, eval)

	d[n] = eval.MulNew(d[a], d[b])
	eval.Rescale(d[n])
}
