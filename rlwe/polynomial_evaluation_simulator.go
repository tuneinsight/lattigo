package rlwe

import (
	"github.com/tuneinsight/lattigo/v4/utils"
)

// DummyOperand is a dummy operand
// that only stores the level and the scale.
type DummyOperand struct {
	Level int
	Scale Scale
}

// Rescale rescales the target DummyOperand n times and returns it.
func (d *DummyOperand) Rescale(params Parameters, n int) *DummyOperand {
	for i := 0; i < n; i++ {
		d.Scale = d.Scale.Div(NewScale(params.Q()[d.Level]))
		d.Level--
	}
	return d
}

// Mul multiplies two DummyOperand, stores the result the taret DummyOperand and returns the result.
func (d *DummyOperand) Mul(a, b *DummyOperand) *DummyOperand {
	d.Level = utils.Min(a.Level, b.Level)
	d.Scale = a.Scale.Mul(b.Scale)
	return d
}

// DummyPowerBasis is a map storing powers of DummyOperands indexed by their power.
type DummyPowerBasis map[int]*DummyOperand

// GenPower populates the target DummyPowerBasis with the nth power.
func (d DummyPowerBasis) GenPower(params Parameters, n, nbModuliPerRescale int) {

	if n < 2 {
		return
	}

	a, b := SplitDegree(n)

	d.GenPower(params, a, nbModuliPerRescale)
	d.GenPower(params, b, nbModuliPerRescale)

	d[n] = new(DummyOperand).Mul(d[a], d[b])
	d[n].Rescale(params, nbModuliPerRescale)
}
