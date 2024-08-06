package factorization

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

// Weierstrass is an elliptic curve y^2 = x^3 + ax + b mod N.
type Weierstrass struct {
	A, B, N *big.Int
}

// Point represents an elliptic curve point in standard coordinates.
type Point struct {
	X, Y *big.Int
}

// Add adds two Weierstrass points together with respect
// to the underlying Weierstrass curve.
// This method does not check if the points lie on
// the underlying curve.
func (w *Weierstrass) Add(P, Q Point) Point {

	tmp := new(big.Int)

	xR, yR := new(big.Int), new(big.Int)

	if P.X.Cmp(tmp.SetUint64(0)) == 0 && P.Y.Cmp(tmp.SetUint64(1)) == 0 {

		return Point{xR.Set(Q.X), yR.Set(Q.Y)}
	}

	if Q.X.Cmp(tmp.SetUint64(0)) == 0 && Q.Y.Cmp(tmp.SetUint64(1)) == 0 {
		return Point{xR.Set(P.X), yR.Set(P.Y)}
	}

	xP, yP := P.X, P.Y
	xQ, yQ := Q.X, Q.Y

	N := w.N

	if xP.Cmp(xQ) == 0 && yP.Cmp(new(big.Int).Sub(N, yQ)) == 0 {
		return Point{xR.SetUint64(0), yR.SetUint64(0)}
	}

	S := new(big.Int) // slope

	if xP != xQ {

		// S = (yQ-yP)/(xQ-xP)
		S.Sub(yQ, yP)
		tmp.Sub(xQ, xP)
		tmp.ModInverse(tmp, N)
		S.Mul(S, tmp)
		S.Mod(S, N)

	} else {

		// S = (3*(xP^2) + a)/(2*yP)
		S.Mul(xP, xP)
		S.Mod(S, N)
		S.Mul(S, new(big.Int).SetUint64(3))
		S.Add(S, w.A)
		S.Mod(S, N)
		tmp.Add(yP, yP)
		tmp.ModInverse(tmp, N)
		S.Mul(S, tmp)
		S.Mod(S, N)
	}

	// s^2 - xP - xQ
	xR.Mul(S, S)
	xR.Mod(xR, N)
	xR.Sub(xR, xP)
	xR.Sub(xR, xQ)
	xR.Mod(xR, N)

	// s*(xP-xR)-yP
	yR.Sub(xP, xR)
	yR.Mul(yR, S)
	yR.Mod(yR, N)
	yR.Sub(yR, yP)
	yR.Mod(yR, N)

	return Point{X: xR, Y: yR}
}

// NewRandomWeierstrassCurve generates a new random Weierstrass curve modulo N,
// along with a random point that lies on the curve.
func NewRandomWeierstrassCurve(N *big.Int) (Weierstrass, Point) {

	var A, B, xG, yG *big.Int
	for {

		// Select random values for A, xG and yG
		A = sampling.RandInt(N)
		xG = sampling.RandInt(N)
		yG = sampling.RandInt(N)

		// Deduces B from Y^2 = X^3 + A * X + B evaluated at point (xG, yG)
		yGpow2 := new(big.Int).Mul(yG, yG)
		yGpow2.Mod(yGpow2, N)

		xGpow3 := new(big.Int).Mul(xG, xG)
		xGpow3.Mod(xGpow3, N)
		xGpow3.Sub(xGpow3, A)
		xGpow3.Mul(xGpow3, xG)
		xGpow3.Mod(xGpow3, N)

		B = new(big.Int).Sub(yGpow2, xGpow3) // B = yG^2 - xG*(xG^2 - A)
		B.Mod(B, N)

		// Checks that 4A^3 + 27B^2 != 0
		fourACube := new(big.Int).Add(A, A)
		fourACube.Mul(fourACube, fourACube)
		fourACube.Mod(fourACube, N)
		fourACube.Mul(fourACube, A)

		twentySevenBSquare := new(big.Int).Mul(B, B)
		twentySevenBSquare.Mod(twentySevenBSquare, N)
		twentySevenBSquare.Mul(twentySevenBSquare, new(big.Int).SetUint64(27))
		twentySevenBSquare.Mod(twentySevenBSquare, N)

		jInvariantQuotient := new(big.Int).Add(fourACube, twentySevenBSquare)
		jInvariantQuotient.Mod(jInvariantQuotient, N)

		if jInvariantQuotient.Cmp(new(big.Int).SetUint64(0)) != 0 && new(big.Int).GCD(nil, nil, N, jInvariantQuotient).Cmp(new(big.Int).SetUint64(1)) == 0 {
			return Weierstrass{
				A: A,
				B: B,
				N: N,
			}, Point{X: xG, Y: yG}
		}
	}
}
