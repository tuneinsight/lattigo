// Package ringqp is implements a wrapper for both the ringQ and ringP.
package ringqp

import (
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// Ring is a structure that implements the operation in the ring R_QP.
// This type is simply a union type between the two Ring types representing
// R_Q and R_P.
type Ring struct {
	RingQ, RingP *ring.Ring
}

func (r Ring) N() int {
	if r.RingQ != nil {
		return r.RingQ.N()
	}

	if r.RingP != nil {
		return r.RingP.N()
	}

	return 0
}

// AtLevel returns a shallow copy of the target ring configured to
// carry on operations at the specified levels.
func (r Ring) AtLevel(levelQ, levelP int) Ring {

	var ringQ, ringP *ring.Ring

	if levelQ > -1 && r.RingQ != nil {
		ringQ = r.RingQ.AtLevel(levelQ)
	}

	if levelP > -1 && r.RingP != nil {
		ringP = r.RingP.AtLevel(levelP)
	}

	return Ring{
		RingQ: ringQ,
		RingP: ringP,
	}
}

// PolyToBigintCentered reconstructs p1 and returns the result in an array of Int.
// Coefficients are centered around Q/2
// gap defines coefficients X^{i*gap} that will be reconstructed.
// For example, if gap = 1, then all coefficients are reconstructed, while
// if gap = 2 then only coefficients X^{2*i} are reconstructed.
func (r Ring) PolyToBigintCentered(p1 Poly, gap int, coeffsBigint []*big.Int) {

	LevelQ := r.LevelQ()
	LevelP := r.LevelP()

	crtReconstructionQ := make([]*big.Int, LevelQ+1)
	crtReconstructionP := make([]*big.Int, LevelP+1)

	tmp := new(big.Int)

	modulusBigint := new(big.Int).SetUint64(1)

	if LevelQ > -1 {
		modulusBigint.Mul(modulusBigint, r.RingQ.ModulusAtLevel[LevelQ])
	}

	if LevelP > -1 {
		modulusBigint.Mul(modulusBigint, r.RingP.ModulusAtLevel[LevelP])
	}

	// Q
	if LevelQ > -1 {
		var QiB = new(big.Int)
		for i, table := range r.RingQ.SubRings[:LevelQ+1] {
			QiB.SetUint64(table.Modulus)
			crtReconstructionQ[i] = new(big.Int).Quo(modulusBigint, QiB)
			tmp.ModInverse(crtReconstructionQ[i], QiB)
			tmp.Mod(tmp, QiB)
			crtReconstructionQ[i].Mul(crtReconstructionQ[i], tmp)
		}
	}

	// P
	if LevelP > -1 {
		var PiB = new(big.Int)
		for i, table := range r.RingP.SubRings[:LevelP+1] {
			PiB.SetUint64(table.Modulus)
			crtReconstructionP[i] = new(big.Int).Quo(modulusBigint, PiB)
			tmp.ModInverse(crtReconstructionP[i], PiB)
			tmp.Mod(tmp, PiB)
			crtReconstructionP[i].Mul(crtReconstructionP[i], tmp)
		}
	}

	modulusBigintHalf := new(big.Int)
	modulusBigintHalf.Rsh(modulusBigint, 1)

	N := r.N()

	var sign int
	for i, j := 0, 0; j < N; i, j = i+1, j+gap {

		tmp.SetUint64(0)
		coeffsBigint[i].SetUint64(0)

		if LevelQ > -1 {
			for k := 0; k < LevelQ+1; k++ {
				coeffsBigint[i].Add(coeffsBigint[i], tmp.Mul(bignum.NewInt(p1.Q.Coeffs[k][j]), crtReconstructionQ[k]))
			}
		}

		if LevelP > -1 {
			for k := 0; k < LevelP+1; k++ {
				coeffsBigint[i].Add(coeffsBigint[i], tmp.Mul(bignum.NewInt(p1.P.Coeffs[k][j]), crtReconstructionP[k]))
			}
		}

		coeffsBigint[i].Mod(coeffsBigint[i], modulusBigint)

		// Centers the coefficients
		sign = coeffsBigint[i].Cmp(modulusBigintHalf)

		if sign == 1 || sign == 0 {
			coeffsBigint[i].Sub(coeffsBigint[i], modulusBigint)
		}
	}
}

// Log2OfStandardDeviation returns base 2 logarithm of the standard deviation of the coefficients
// of the polynomial.
func (r Ring) Log2OfStandardDeviation(poly Poly) (std float64) {

	N := r.N()

	prec := uint(128)

	coeffs := make([]*big.Int, N)

	for i := 0; i < N; i++ {
		coeffs[i] = new(big.Int)
	}

	r.PolyToBigintCentered(poly, 1, coeffs)

	mean := bignum.NewFloat(0, prec)
	tmp := bignum.NewFloat(0, prec)

	for i := 0; i < N; i++ {
		mean.Add(mean, tmp.SetInt(coeffs[i]))
	}

	mean.Quo(mean, bignum.NewFloat(float64(N), prec))

	stdFloat := bignum.NewFloat(0, prec)

	for i := 0; i < N; i++ {
		tmp.SetInt(coeffs[i])
		tmp.Sub(tmp, mean)
		tmp.Mul(tmp, tmp)
		stdFloat.Add(stdFloat, tmp)
	}

	stdFloat.Quo(stdFloat, bignum.NewFloat(float64(N-1), prec))

	stdFloat.Sqrt(stdFloat)

	stdF64, _ := stdFloat.Float64()

	return math.Log2(stdF64)
}

// LevelQ returns the level at which the target
// ring operates for the modulus Q.
func (r Ring) LevelQ() int {
	if r.RingQ != nil {
		return r.RingQ.Level()
	}

	return -1
}

// LevelP returns the level at which the target
// ring operates for the modulus P.
func (r Ring) LevelP() int {
	if r.RingP != nil {
		return r.RingP.Level()
	}

	return -1
}

func (r Ring) Equal(p1, p2 Poly) (v bool) {
	v = true
	if r.RingQ != nil {
		v = v && r.RingQ.Equal(p1.Q, p2.Q)
	}

	if r.RingP != nil {
		v = v && r.RingP.Equal(p1.P, p2.P)
	}

	return
}

// NewPoly creates a new polynomial with all coefficients set to 0.
func (r Ring) NewPoly() Poly {
	var Q, P ring.Poly
	if r.RingQ != nil {
		Q = r.RingQ.NewPoly()
	}

	if r.RingP != nil {
		P = r.RingP.NewPoly()
	}
	return Poly{Q, P}
}
