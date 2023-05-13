package approximation

import (
	"math/big"
	"testing"

	"github.com/tuneinsight/lattigo/v4/utils/bignum"
	"github.com/tuneinsight/lattigo/v4/utils/bignum/polynomial"
)

func TestRemez(t *testing.T) {
	sigmoid := func(x *big.Float) (y *big.Float) {
		z := new(big.Float).Set(x)
		z.Neg(z)
		z = bignum.Exp(z)
		z.Add(z, bignum.NewFloat(1, x.Prec()))
		y = bignum.NewFloat(1, x.Prec())
		y.Quo(y, z)
		return
	}

	prec := uint(96)

	scanStep := bignum.NewFloat(2, prec)
	scanStep.Quo(scanStep, bignum.NewFloat(1000, prec))

	intervals := []bignum.Interval{
		{A: bignum.NewFloat(-6, prec), B: bignum.NewFloat(-5, prec), Nodes: 4},
		{A: bignum.NewFloat(-3, prec), B: bignum.NewFloat(-2, prec), Nodes: 4},
		{A: bignum.NewFloat(-1, prec), B: bignum.NewFloat(1, prec), Nodes: 4},
		{A: bignum.NewFloat(2, prec), B: bignum.NewFloat(3, prec), Nodes: 4},
		{A: bignum.NewFloat(5, prec), B: bignum.NewFloat(6, prec), Nodes: 4},
	}

	params := RemezParameters{
		Function:        sigmoid,
		Basis:           polynomial.Chebyshev,
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
