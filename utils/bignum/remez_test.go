package bignum

import (
	"math/big"
	"testing"
)

func TestRemez(t *testing.T) {
	sigmoid := func(x *big.Float) (y *big.Float) {
		z := new(big.Float).Set(x)
		z.Neg(z)
		z = Exp(z)
		z.Add(z, NewFloat(1, x.Prec()))
		y = NewFloat(1, x.Prec())
		y.Quo(y, z)
		return
	}

	prec := uint(96)

	scanStep := NewFloat(2, prec)
	scanStep.Quo(scanStep, NewFloat(1000, prec))

	intervals := []Interval{
		{A: *NewFloat(-6, prec), B: *NewFloat(-5, prec), Nodes: 4},
		{A: *NewFloat(-3, prec), B: *NewFloat(-2, prec), Nodes: 4},
		{A: *NewFloat(-1, prec), B: *NewFloat(1, prec), Nodes: 4},
		{A: *NewFloat(2, prec), B: *NewFloat(3, prec), Nodes: 4},
		{A: *NewFloat(5, prec), B: *NewFloat(6, prec), Nodes: 4},
	}

	params := RemezParameters{
		Function:        sigmoid,
		Basis:           Chebyshev,
		Intervals:       intervals,
		ScanStep:        scanStep,
		Prec:            prec,
		OptimalScanStep: true,
	}

	r := NewRemez(params)
	r.Approximate(200, 1e-15)
	r.ShowCoeffs(50)
	r.ShowError(50)
}
