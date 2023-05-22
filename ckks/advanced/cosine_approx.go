package advanced

// This is the Go implementation of the approximation polynomial algorithm from Han and Ki in
//    "Better Bootstrapping for Approximate Homomorphic Encryption", <https://epring.iacr.org/2019/688O>.
// The algorithm was originally implemented in C++, available at
//    https://github.com/DohyeongKi/better-homomorphic-sine-evaluation

import (
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

const (
	defaultPrecision = uint(512)
)

var (
	log2TwoPi = math.Log2(2 * math.Pi)
	aQuarter  = bignum.NewFloat(0.25, defaultPrecision)
	pi        = bignum.Pi(defaultPrecision)
)

// ApproximateCos computes a polynomial approximation of degree "degree" in Chevyshev basis of the function
// cos(2*pi*x/2^"scnum") in the range -"K" to "K"
// The nodes of the Chevyshev approximation are are located from -dev to +dev at each integer value between -K and -K
func ApproximateCos(K, degree int, dev float64, scnum int) []*big.Float {

	var scfac = bignum.NewFloat(float64(int(1<<scnum)), defaultPrecision)

	deg, totdeg := genDegrees(degree, K, dev)

	x, p, c, totdeg := genNodes(deg, dev, totdeg, K, scnum)

	tmp := new(big.Float)

	var T = make([][]*big.Float, totdeg)
	for i := 0; i < totdeg; i++ {
		T[i] = make([]*big.Float, totdeg)
	}

	KBig := new(big.Float).SetInt64(int64(K))

	for i := 0; i < totdeg; i++ {

		T[i][0] = bignum.NewFloat(1.0, defaultPrecision)

		T[i][1] = new(big.Float).Set(x[i])

		tmp.Quo(KBig, scfac)

		T[i][1].Quo(T[i][1], tmp)

		for j := 2; j < totdeg; j++ {

			T[i][j] = bignum.NewFloat(2.0, defaultPrecision)

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
		maxabs.Abs(T[i][i])
		maxindex = i
		for j := i + 1; j < totdeg; j++ {
			tmp.Abs(T[j][i])
			if tmp.Cmp(maxabs) == 1 {
				maxabs.Abs(T[j][i])
				maxindex = j
			}
		}

		if i != maxindex {
			for j := i; j < totdeg; j++ {
				tmp.Set(T[maxindex][j])
				T[maxindex][j].Set(T[i][j])
				T[i][j].Set(tmp)
			}

			tmp.Set(p[maxindex])
			p[maxindex].Set(p[i])
			p[i].Set(tmp)
		}

		for j := i + 1; j < totdeg; j++ {
			T[i][j].Quo(T[i][j], T[i][i])
		}

		p[i].Quo(p[i], T[i][i])
		T[i][i] = bignum.NewFloat(1.0, defaultPrecision)

		for j := i + 1; j < totdeg; j++ {
			tmp.Mul(T[j][i], p[i])
			p[j].Sub(p[j], tmp)
			for l := i + 1; l < totdeg; l++ {
				tmp.Mul(T[j][i], T[i][l])
				T[j][l].Sub(T[j][l], tmp)
			}
			T[j][i] = bignum.NewFloat(0.0, defaultPrecision)
		}
	}

	c[totdeg-1] = p[totdeg-1]
	for i := totdeg - 2; i >= 0; i-- {
		c[i] = new(big.Float).Set(p[i])
		for j := i + 1; j < totdeg; j++ {
			tmp.Mul(T[i][j], c[j])
			c[i].Sub(c[i], tmp)
		}
	}

	return c[:totdeg-1]
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

	/*
		fmt.Println("==============================================")
		fmt.Println("==Degree Searching Result=====================")
		fmt.Println("==============================================")
		if iter == maxiter{
			fmt.Println("More Iterations Needed")
		}else{
			fmt.Println("Degree of Polynomial :", totdeg-1)
			fmt.Println("Degree :", deg)
		}
		fmt.Println("==============================================")
	*/

	return deg, totdeg
}

func genNodes(deg []int, dev float64, totdeg, K, scnum int) ([]*big.Float, []*big.Float, []*big.Float, int) {

	var scfac = bignum.NewFloat(1<<scnum, defaultPrecision)

	var intersize = bignum.NewFloat(1.0/dev, defaultPrecision)

	var z = make([]*big.Float, totdeg)

	for i := range z {
		z[i] = bignum.NewFloat(0, defaultPrecision)
	}

	var cnt int
	if deg[0]%2 != 0 {
		cnt++
	}

	tmp := new(big.Float)

	for i := K - 1; i > 0; i-- {
		for j := 1; j <= deg[i]; j++ {

			tmp.Mul(pi, new(big.Float).SetInt64(int64((2*j - 1))))
			tmp.Quo(tmp, new(big.Float).SetInt64(int64(2*deg[i])))
			tmp = bignum.Cos(tmp)
			tmp.Mul(tmp, intersize)

			z[cnt].Add(new(big.Float).SetInt64(int64(i)), tmp)
			cnt++

			z[cnt].Sub(new(big.Float).SetInt64(int64(-i)), tmp)
			cnt++
		}
	}

	for j := 1; j <= deg[0]/2; j++ {

		tmp.Mul(pi, new(big.Float).SetInt64(int64((2*j - 1))))
		tmp.Quo(tmp, new(big.Float).SetInt64(int64(2*deg[j])))
		tmp = bignum.Cos(tmp)
		tmp.Mul(tmp, intersize)

		z[cnt].Add(z[cnt], tmp)
		cnt++

		z[cnt].Sub(z[cnt], tmp)
		cnt++
	}

	// cos(2*pi*(x-0.25)/r)
	var d = make([]*big.Float, totdeg)
	for i := 0; i < totdeg; i++ {
		d[i] = cos2PiXMinusQuarterOverR(z[i], scfac)
	}

	for j := 1; j < totdeg; j++ {
		for l := 0; l < totdeg-j; l++ {

			// d[l] = d[l+1] - d[l]
			d[l].Sub(d[l+1], d[l])

			// d[l] = (d[l+1] - d[l])/(z[l+j] - z[l])
			tmp.Sub(z[l+j], z[l])
			d[l].Quo(d[l], tmp)
		}
	}

	totdeg++

	var x = make([]*big.Float, totdeg)
	for i := 0; i < totdeg; i++ {
		// x[i] = K
		x[i] = bignum.NewFloat(float64(K), defaultPrecision)

		// x[i] = K/r
		x[i].Quo(x[i], scfac)

		// x[i] = (K/r) * cos(PI * i/(totdeg-1))
		tmp.Mul(new(big.Float).SetInt64(int64(i)), pi)
		tmp.Quo(tmp, new(big.Float).SetInt64(int64(totdeg-1)))
		x[i].Mul(x[i], bignum.Cos(tmp))
	}

	var c = make([]*big.Float, totdeg)
	var p = make([]*big.Float, totdeg)
	for i := 0; i < totdeg; i++ {
		p[i] = new(big.Float).Set(d[0])
		for j := 1; j < totdeg-1; j++ {
			tmp.Sub(x[i], z[j])
			p[i].Mul(p[i], tmp)
			p[i].Add(p[i], d[j])
		}
	}

	return x, p, c, totdeg
}
