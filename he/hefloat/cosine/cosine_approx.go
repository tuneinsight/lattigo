// Package cosine method is the Go implementation of the polynomial-approximation algorithm by Han and Ki in
//
//	"Better Bootstrapping for Approximate Homomorphic Encryption", <https://epring.iacr.org/2019/688O>.
//
// The algorithm was originally implemented in C++, available at
//
//	https://github.com/DohyeongKi/better-homomorphic-sine-evaluation
package cosine

import (
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v5/utils/bignum"
)

const (
	EncodingPrecision = uint(256)
)

var (
	log2TwoPi = math.Log2(2 * math.Pi)
	aQuarter  = bignum.NewFloat(0.25, EncodingPrecision)
	pi        = bignum.Pi(EncodingPrecision)
)

// ApproximateCos computes a polynomial approximation of degree "degree" in Chebyshev basis of the function
// cos(2*pi*x/2^"scnum") in the range -"K" to "K"
// The nodes of the Chebyshev approximation are are located from -dev to +dev at each integer value between -K and -K
func ApproximateCos(K, degree int, dev float64, scnum int) []*big.Float {

	// Gets the list of degree per interval and the total degree
	deg, totdeg := genDegrees(degree, K, dev)

	// Generates the nodes for each interval, updates the total degree if needed
	nodes, y := genNodes(deg, dev, totdeg, K, scnum)

	// Solves the linear system and returns the coefficients
	return solve(totdeg, K, scnum, nodes, y)[:totdeg]
}

// y = cos(2 * pi * (x - 0.25)/r)
func cos2PiXMinusQuarterOverR(x, r *big.Float) (y *big.Float) {
	//y = 2 * pi
	y = bignum.NewFloat(2.0, EncodingPrecision)
	y.Mul(y, pi)

	// x = (x - 0.25)/r
	x.Sub(x, aQuarter)
	x.Quo(x, r)

	// y = 2 * pi * (x - 0.25)/r
	y.Mul(y, x)

	// y = cos(2 * pi * (x - 0.25)/r)
	return bignum.Cos(y)
}

func log2(x float64) float64 {
	return math.Log2(x)
}

func abs(x float64) float64 {
	return math.Abs(x)
}

func maxIndex(array []float64) (maxind int) {
	max := array[0]
	for i := 1; i < len(array); i++ {
		if array[i] > max {
			maxind = i
			max = array[i]
		}
	}

	return
}

// genDegrees returns the optimal list of nodes for each of the 0 <= i < K intervals [i +/- dev]
// such that the sum of the nodes of all intervals is equal to degree.
func genDegrees(degree, K int, dev float64) ([]int, int) {

	var degbdd = degree + 1

	var totdeg = 2*K - 1

	var err = 1.0 / dev

	var deg = make([]int, K)
	for i := 0; i < K; i++ {
		deg[i] = 1
	}

	var bdd = make([]float64, K)
	var temp = float64(0)
	for i := 1; i <= (2*K - 1); i++ {
		temp -= log2(float64(i))
	}
	temp += (2*float64(K) - 1) * log2TwoPi
	temp += log2(err)

	for i := 0; i < K; i++ {
		bdd[i] = temp
		for j := 1; j <= K-1-i; j++ {
			bdd[i] += log2(float64(j) + err)
		}
		for j := 1; j <= K-1+i; j++ {
			bdd[i] += log2(float64(j) + err)
		}
	}

	var maxiter = 200
	var iter int

	for iter = 0; iter < maxiter; iter++ {
		if totdeg >= degbdd {
			break
		}
		var maxi = maxIndex(bdd)

		if maxi != 0 {
			if totdeg+2 > degbdd {
				break
			}

			for i := 0; i < K; i++ {
				bdd[i] -= log2(float64(totdeg + 1))
				bdd[i] -= log2(float64(totdeg + 2))
				bdd[i] += 2.0 * log2TwoPi

				if i != maxi {
					bdd[i] += log2(abs(float64(i-maxi)) + err)
					bdd[i] += log2(float64(i+maxi) + err)
				} else {
					bdd[i] += log2(err) - 1.0
					bdd[i] += log2(2.0*float64(i) + err)
				}
			}

			totdeg += 2
		} else {
			bdd[0] -= log2(float64(totdeg + 1))
			bdd[0] += log2(err) - 1.0
			bdd[0] += log2TwoPi
			for i := 1; i < K; i++ {
				bdd[i] -= log2(float64(totdeg + 1))
				bdd[i] += log2TwoPi
				bdd[i] += log2(float64(i) + err)
			}

			totdeg++
		}

		deg[maxi]++
	}

	return deg, totdeg
}

func genNodes(deg []int, dev float64, totdeg, K, scnum int) ([]*big.Float, []*big.Float) {

	var scfac = bignum.NewFloat(1<<scnum, EncodingPrecision)

	// Interval [i+e, i-e] with e = 1/dev
	var intersize = bignum.NewFloat(1.0/dev, EncodingPrecision)

	// ===================
	// Allocates the nodes
	// ===================

	// nodes
	var nodes = make([]*big.Float, totdeg)
	for i := range nodes {
		nodes[i] = bignum.NewFloat(0, EncodingPrecision)
	}

	// Ensures deg[0] is even
	var cnt int
	if deg[0]%2 != 0 {
		cnt++
	}

	// Loops over the intervals
	// [-K +/- nodes] U ... U [-1 +/- nodes] U [+/- nodes]
	tmp := new(big.Float)
	for i := K - 1; i > 0; i-- {

		twodegi := bignum.NewFloat(2*deg[i], EncodingPrecision)
		iF := bignum.NewFloat(i, EncodingPrecision)

		// For each node in the interval
		for j := 0; j < deg[i]; j++ {

			tmp.Mul(pi, new(big.Float).SetInt64(int64((2 * j))))
			tmp.Quo(tmp, twodegi)
			tmp.Mul(bignum.Cos(tmp), intersize)

			//   i + cos(pi * (2j-1) / (2*deg[i])) * (1/intersize)
			nodes[cnt].Add(iF, tmp)
			cnt++

			//  -i - cos(pi * (2j-1) / (2*deg[i])) * (1/intersize)
			nodes[cnt].Neg(nodes[cnt-1])
			cnt++
		}
	}

	// Center interval
	// [+/- nodes]
	twodegi := new(big.Float).SetInt64(int64(2 * deg[0]))
	for j := 0; j < deg[0]/2; j++ {

		tmp.Mul(pi, new(big.Float).SetInt64(int64((2 * j))))
		tmp.Quo(tmp, twodegi)
		tmp.Mul(bignum.Cos(tmp), intersize)

		// 0 + cos(pi * (2j-1) / (2*deg[i])) * (1/intersize)
		nodes[cnt].Set(tmp)
		cnt++

		// 0 - cos(pi * (2j-1) / (2*deg[i])) * (1/intersize)
		nodes[cnt].Neg(nodes[cnt-1])
		cnt++
	}

	// Evaluates the nodes y[i] = f(nodes[i])
	var y = make([]*big.Float, totdeg)
	for i := 0; i < totdeg; i++ {
		// y[i] = cos(2*pi*(nodes[i]-0.25)/r)
		y[i] = cos2PiXMinusQuarterOverR(nodes[i], scfac)
	}

	return nodes, y
}

func solve(totdeg, K, scnum int, nodes, y []*big.Float) []*big.Float {

	// 2^r
	scfac := bignum.NewFloat(float64(int(1<<scnum)), EncodingPrecision)

	tmp := new(big.Float)

	//=========================
	// Solves the linear system
	//=========================

	for j := 1; j < totdeg; j++ {
		for i := 0; i < totdeg-j; i++ {

			// y[i] = y[i+1] - y[i]
			y[i].Sub(y[i+1], y[i])

			// y[i] = (y[i+1] - y[i])/(nodes[i+j] - nodes[i])
			tmp.Sub(nodes[i+j], nodes[i])
			y[i].Quo(y[i], tmp)
		}
	}

	totdeg++

	var x = make([]*big.Float, totdeg)
	for i := 0; i < totdeg; i++ {
		// x[i] = K
		x[i] = bignum.NewFloat(float64(K), EncodingPrecision)

		// x[i] = K/r
		x[i].Quo(x[i], scfac)

		// x[i] = (K/r) * cos(PI * i/(totdeg-1))
		tmp.Mul(new(big.Float).SetInt64(int64(i)), pi)
		tmp.Quo(tmp, new(big.Float).SetInt64(int64(totdeg-1)))
		x[i].Mul(x[i], bignum.Cos(tmp))
	}

	var p = make([]*big.Float, totdeg)
	for i := 0; i < totdeg; i++ {
		p[i] = new(big.Float).Set(y[0])
		for j := 1; j < totdeg-1; j++ {
			tmp.Sub(x[i], nodes[j])
			p[i].Mul(p[i], tmp)
			p[i].Add(p[i], y[j])
		}
	}

	var T = make([][]*big.Float, totdeg)
	for i := 0; i < totdeg; i++ {
		T[i] = make([]*big.Float, totdeg)
	}

	KBig := new(big.Float).SetInt64(int64(K))

	// Constructs the totdeg x totdeg linear system using x
	for i := 0; i < totdeg; i++ {

		T[i][0] = bignum.NewFloat(1.0, EncodingPrecision)

		T[i][1] = new(big.Float).Set(x[i])

		tmp.Quo(KBig, scfac)

		T[i][1].Quo(T[i][1], tmp)

		for j := 2; j < totdeg; j++ {

			T[i][j] = bignum.NewFloat(2.0, EncodingPrecision)

			tmp.Quo(KBig, scfac)
			tmp.Quo(x[i], tmp)
			T[i][j].Mul(T[i][j], tmp)
			T[i][j].Mul(T[i][j], T[i][j-1])
			T[i][j].Sub(T[i][j], T[i][j-2])
		}
	}

	var maxabs = new(big.Float)
	var maxindex int
	for i := 0; i < totdeg-1; i++ {

		// Solves by max value of the current i-th column
		// which minimizes numerical errors and avoids
		// division by zero.

		// Finds the max index in the i-th column
		// of the triangular matrix
		maxabs.Abs(T[i][i])
		maxindex = i
		for j := i + 1; j < totdeg; j++ {
			tmp.Abs(T[j][i])
			if tmp.Cmp(maxabs) == 1 {
				maxabs.Abs(T[j][i])
				maxindex = j
			}
		}

		// Swaps the row with the max index with the current row
		if i != maxindex {
			// swaps T[maxindex][j] with T[i][j]
			for j := i; j < totdeg; j++ {
				T[maxindex][j], T[i][j] = T[i][j], T[maxindex][j]
			}

			// Does the same for the target vector
			p[maxindex], p[i] = p[i], p[maxindex]
		}

		// SOLVES THE LINEAR SYSTEM
		for j := i + 1; j < totdeg; j++ {
			T[i][j].Quo(T[i][j], T[i][i])
		}

		p[i].Quo(p[i], T[i][i])
		T[i][i] = bignum.NewFloat(1.0, EncodingPrecision)

		for j := i + 1; j < totdeg; j++ {
			tmp.Mul(T[j][i], p[i])
			p[j].Sub(p[j], tmp)
			for l := i + 1; l < totdeg; l++ {
				tmp.Mul(T[j][i], T[i][l])
				T[j][l].Sub(T[j][l], tmp)
			}
			T[j][i] = bignum.NewFloat(0.0, EncodingPrecision)
		}
	}

	// GETS THE COEFFICIENTS
	var c = make([]*big.Float, totdeg)
	c[totdeg-1] = p[totdeg-1]
	for i := totdeg - 2; i >= 0; i-- {
		c[i] = new(big.Float).Set(p[i])
		for j := i + 1; j < totdeg; j++ {
			tmp.Mul(T[i][j], c[j])
			c[i].Sub(c[i], tmp)
		}
	}

	return c
}
