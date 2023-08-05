package circuits

import (
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/bgv"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

type SimOperand struct {
	Level int
	Scale rlwe.Scale
}

type SimEvaluator interface {
	MulNew(op0, op1 *SimOperand) *SimOperand
	Rescale(op0 *SimOperand)
	PolynomialDepth(degree int) int
	UpdateLevelAndScaleGiantStep(lead bool, tLevelOld int, tScaleOld, xPowScale rlwe.Scale) (tLevelNew int, tScaleNew rlwe.Scale)
	UpdateLevelAndScaleBabyStep(lead bool, tLevelOld int, tScaleOld rlwe.Scale) (tLevelNew int, tScaleNew rlwe.Scale)
}

// SimPowerBasis is a map storing powers of SimOperands indexed by their power.
type SimPowerBasis map[int]*SimOperand

// GenPower populates the target SimPowerBasis with the nth power.
func (d SimPowerBasis) GenPower(params rlwe.ParameterProvider, n int, eval SimEvaluator) {

	if n < 2 {
		return
	}

	a, b := SplitDegree(n)

	d.GenPower(params, a, eval)
	d.GenPower(params, b, eval)

	d[n] = eval.MulNew(d[a], d[b])
	eval.Rescale(d[n])
}

type floatSimEvaluator struct {
	params                      ckks.Parameters
	levelsConsummedPerRescaling int
}

func (d floatSimEvaluator) PolynomialDepth(degree int) int {
	return d.levelsConsummedPerRescaling * (bits.Len64(uint64(degree)) - 1)
}

// Rescale rescales the target SimOperand n times and returns it.
func (d floatSimEvaluator) Rescale(op0 *SimOperand) {
	for i := 0; i < d.levelsConsummedPerRescaling; i++ {
		op0.Scale = op0.Scale.Div(rlwe.NewScale(d.params.Q()[op0.Level]))
		op0.Level--
	}
}

// Mul multiplies two SimOperand, stores the result the taret SimOperand and returns the result.
func (d floatSimEvaluator) MulNew(op0, op1 *SimOperand) (opOut *SimOperand) {
	opOut = new(SimOperand)
	opOut.Level = utils.Min(op0.Level, op1.Level)
	opOut.Scale = op0.Scale.Mul(op1.Scale)
	return
}

func (d floatSimEvaluator) UpdateLevelAndScaleBabyStep(lead bool, tLevelOld int, tScaleOld rlwe.Scale) (tLevelNew int, tScaleNew rlwe.Scale) {

	tLevelNew = tLevelOld
	tScaleNew = tScaleOld

	if lead {
		for i := 0; i < d.levelsConsummedPerRescaling; i++ {
			tScaleNew = tScaleNew.Mul(rlwe.NewScale(d.params.Q()[tLevelNew-i]))
		}
	}

	return
}

func (d floatSimEvaluator) UpdateLevelAndScaleGiantStep(lead bool, tLevelOld int, tScaleOld, xPowScale rlwe.Scale) (tLevelNew int, tScaleNew rlwe.Scale) {

	Q := d.params.Q()

	var qi *big.Int
	if lead {
		qi = bignum.NewInt(Q[tLevelOld])
		for i := 1; i < d.levelsConsummedPerRescaling; i++ {
			qi.Mul(qi, bignum.NewInt(Q[tLevelOld-i]))
		}
	} else {
		qi = bignum.NewInt(Q[tLevelOld+d.levelsConsummedPerRescaling])
		for i := 1; i < d.levelsConsummedPerRescaling; i++ {
			qi.Mul(qi, bignum.NewInt(Q[tLevelOld+d.levelsConsummedPerRescaling-i]))
		}
	}

	tLevelNew = tLevelOld + d.levelsConsummedPerRescaling
	tScaleNew = tScaleOld.Mul(rlwe.NewScale(qi))
	tScaleNew = tScaleNew.Div(xPowScale)

	return
}

func (d floatSimEvaluator) GetPolynmialDepth(degree int) int {
	return d.levelsConsummedPerRescaling * (bits.Len64(uint64(degree)) - 1)
}

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

// Rescale rescales the target SimOperand n times and returns it.
func (d simIntegerPolynomialEvaluator) Rescale(op0 *SimOperand) {
	if !d.InvariantTensoring {
		op0.Scale = op0.Scale.Div(rlwe.NewScale(d.params.Q()[op0.Level]))
		op0.Level--
	}
}

// Mul multiplies two SimOperand, stores the result the taret SimOperand and returns the result.
func (d simIntegerPolynomialEvaluator) MulNew(op0, op1 *SimOperand) (opOut *SimOperand) {
	opOut = new(SimOperand)
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
