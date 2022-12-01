package ring

import (
	"math"
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/utils"
)

// FactorizeECM finds a factor of N
// using ECM factorization.
func FactorizeECM(N uint64) (gcd uint64) {
	ecm := NewECM(N)

	for {

		P := ecm.G

		// !B * P
		for i := 1; i < ecm.B; i++ {
			if P, gcd = ecm.checkThenMul(uint64(i), P); gcd != 1 {
				return
			}
		}

		ecm.NewCurve()
	}
}

// ECM is a struct to store the necessary
// values for ECM factorization
type ECM struct {
	Weierstrass
	N uint64
	B int
	C uint64
	G Point
}

// NewECM instantiates a new ECM factorization for the number N.
func NewECM(N uint64) ECM {

	curve, G := NewRandomWeierstrassCurve(N)

	return ECM{
		N:           N,
		B:           int(math.Exp(math.Sqrt(2*math.Log(float64(N))*math.Log(math.Log(float64(N))))) + 0.5),
		C:           N + 1 + uint64(math.Sqrt(float64(N))+0.5),
		Weierstrass: curve,
		G:           G,
	}
}

func (ecm *ECM) NewCurve() {
	ecm.Weierstrass, ecm.G = NewRandomWeierstrassCurve(ecm.N)
}

func (ecm *ECM) checkThenAdd(P, Q Point) (S Point, gcd uint64) {

	N := ecm.N
	if P.X == Q.X && P.Y == Q.Y {
		if gcd = utils.GCD(P.Y<<1, N); gcd != 1 {
			return
		}
	} else {
		if gcd = utils.GCD(Q.X+N-P.X, N); gcd != 1 {
			return
		}
	}

	return ecm.Weierstrass.Add(P, Q), 1
}

func (ecm *ECM) checkThenMul(k uint64, P Point) (Q Point, gcd uint64) {

	Q = Point{0, 1}

	for k > 0 {
		if k&1 == 1 {
			if Q, gcd = ecm.checkThenAdd(P, Q); gcd != 1 {
				return
			}
		}

		if P, gcd = ecm.checkThenAdd(P, P); gcd != 1 {
			return
		}

		k >>= 1
	}
	return
}

// Weierstrass is an elliptic curve y^2 = x^3 + ax + b mod N.
type Weierstrass struct {
	A, B, N   uint64
	bredParam []uint64
}

// Add adds two Weierstrass points together with respect
// to the underlying Weierstrass curve.
// This method does not check if the point are lying on
// the underlying curve.
func (w *Weierstrass) Add(P, Q Point) Point {
	if P.X == 0 && P.Y == 1 {
		return Q
	}

	if Q.X == 0 && Q.Y == 1 {
		return P
	}

	xP, yP := P.X, P.Y
	xQ, yQ := Q.X, Q.Y

	N := w.N

	if xP == xQ && yP == N-yQ {
		return Point{0, 0}
	}

	bredParam := w.bredParam

	var S uint64
	if xP != xQ {
		S = BRed(CRed(yQ+N-yP, N), ModInv(CRed(xQ+N-xP, N), N), N, bredParam) // (yQ-yP)/(xQ-xP)
	} else {
		S = BRed(3*BRed(xP, xP, N, bredParam)+w.A, ModInv(2*yP, N), N, bredParam) // (3*(xP^2) + a)/(2*yP)
	}

	xR := BRedAdd(BRed(S, S, N, bredParam)+N-xP+N-xQ, N, bredParam) // s^2 - xP - xQ
	yR := CRed(BRed(S, xP+N-xR, N, bredParam)+N-yP, N)              // s*(xP-xR)-yP

	return Point{X: xR, Y: yR}
}

func ModInv(x, p uint64) (xInv uint64) {
	return new(big.Int).ModInverse(new(big.Int).SetUint64(x), new(big.Int).SetUint64(p)).Uint64()
}

// NewRandomWeierstrassCurve generates a new random Weierstrass curve modulo N,
// along with a random point that lies on the curve.
func NewRandomWeierstrassCurve(N uint64) (Weierstrass, Point) {

	bredParam := BRedParams(N)

	mask := uint64(1<<bits.Len64(N)) - 1

	var A, B, xG, yG uint64
	for {

		// Select random values for A, xG and yG
		A = utils.RandUint64() & mask
		for A >= N {
			A = utils.RandUint64() & mask
		}

		xG = utils.RandUint64() & mask
		for xG >= N {
			xG = utils.RandUint64() & mask
		}

		yG = utils.RandUint64() & mask
		for yG >= N {
			yG = utils.RandUint64() & mask
		}

		// Deduces B from Y^2 = X^3 + A * X + B evaluated at point (xG, yG)
		yGpow2 := BRed(yG, yG, N, bredParam)
		xGpow3 := BRed(BRed(xG, xG, N, bredParam), xG, N, bredParam)
		AxG := BRed(A, xG, N, bredParam)
		B = BRedAdd(yGpow2+(N-xGpow3)+(N-AxG), N, bredParam) // B = yG^2 - xG^3 - A * xG

		// Checks that 4A^3 + 27B^2 != 0
		fourACube := (BRed(BRed(A, A, N, bredParam), A, N, bredParam) << 2) % N
		twentySevenBSquare := (27 * BRed(B, B, N, bredParam)) % N
		jInvariantQuotient := (fourACube + twentySevenBSquare) % N
		if jInvariantQuotient != 0 && utils.GCD(N, jInvariantQuotient) == 1 {
			return Weierstrass{
				A:         A,
				B:         B,
				N:         N,
				bredParam: bredParam,
			}, Point{X: xG, Y: yG}
		}
	}
}

// Point represent an elliptic curve point.
type Point struct {
	X, Y uint64
}
