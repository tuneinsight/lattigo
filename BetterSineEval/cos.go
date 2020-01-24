package main

import (
    "fmt"
    "math/big"
    "math"
)

var pi = "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491298336733624406566430860213949463952247371907021798609437027705392171762931767523846748184676694051320005681271452635608277857713427577896091736371787214684409012249534301465495853710507922796892589235420199561121290219608640344181598136297747713099605187072113499999983729780499510597317328160963185950244594553469083026425223082533446850352619311881710100031378387528865875332083814206171776691473035982534904287554687311595628638823537875937519577818577805321712268066130019278766111959092164201989"
var M_PI = 3.141592653589793238462643383279502884

func NewFloat(x float64) (y *big.Float) {
    y = new(big.Float)
    y.SetPrec(1000) // decimal precision
    y.SetFloat64(x)
    return
}

func maxIndex(array []float64) (max_ind int){
	max := array[0]
	for i := 1 ; i < len(array); i++ {
		if array[i] > max {
			max_ind = i
			max = array[i]
		}
	}

	return
}

// Arbitrary precision computation of Cos(x)
// Iterative process with an error of ~10^{−0.60206*k} after k iterations.
// ref : Johansson, B. Tomas, An elementary algorithm to evaluatetrigonometric functions to high precision, 2018
func Cos(x *big.Float)(cosx *big.Float){
    tmp := new(big.Float)
    
    k := 1000 // number of iterations
    t := NewFloat(0.5)
    half := new(big.Float).Copy(t)
    
    for i := 1 ; i < k-1; i++ {
        t.Mul(t, half)
    }
    
    s := new(big.Float).Mul(x, t)
    s.Mul(s,x)
    s.Mul(s,t)
    
    four := NewFloat(4.0)
    
    for i := 1 ; i < k; i++ {
        tmp.Sub(four, s)
        s.Mul(s, tmp)
    }
    
    cosx = new(big.Float).Quo(s, NewFloat(2.0))
    cosx.Sub(NewFloat(1.0), cosx)
    return

}

func log2(x float64) (float64){
	return math.Log2(x)
}

func abs(x float64) (float64){
	return math.Abs(x)
}

func main() {

	var PI = new(big.Float)
	PI.SetPrec(1000)
	PI.SetString(pi)

	var K = int(12)

	var deg_bdd = 60+1

	var deg = make([]int, K)
	for i := 0 ; i < K; i++ {
		deg[i] = 1
	}
	var tot_deg = 2*K-1

	var dev = 10
	var err = 1.0/float64(int(1<<dev))

	var sc_num = 1
	var sc_fac = NewFloat(float64(int(1<<sc_num)))

	var bdd = make([]float64, K)
	var temp = float64(0)
	for i := 1 ; i <= (2*K-1); i++ {
		temp -= log2(float64(i))
	}
	temp += (2*float64(K)-1)*log2(2*M_PI)
	temp += log2(err)

	for i := 0 ; i < K ; i++ {
		bdd[i] = temp
		for j := 1 ; j <= K-1-i; j++ {
			bdd[i] += log2(float64(j) + err)
		}
		for j := 1 ; j <= K-1+i; j++ {
			bdd[i] += log2(float64(j) + err)
		}
	}


	var max_iter = 200
	var iter int

	for iter = 0 ; iter < max_iter; iter++ {
		if tot_deg >= deg_bdd{
			break
		}
		var maxi = maxIndex(bdd)

		if maxi != 0 {
			if tot_deg + 2 > deg_bdd{
				break
			}

			for i := 0 ; i < K; i++ {
				bdd[i] -= log2(float64(tot_deg+1))
				bdd[i] -= log2(float64(tot_deg+2))
				bdd[i] += 2.0*log2(2.0*M_PI)

				if i != maxi {
					bdd[i] += log2(abs(float64(i-maxi)) + err)
					bdd[i] += log2(float64(i+maxi) + err)
				}else{
					bdd[i] += log2(err)-1.0
					bdd[i] += log2(2.0*float64(i) + err)
				}
			}

			tot_deg += 2
		} else {
			bdd[0] -= log2(float64(tot_deg+1))
			bdd[0] += log2(err)-1.0
			bdd[0] += log2(2.0*M_PI)
			for i := 1 ; i < K; i++ {
				bdd[i] -= log2(float64(tot_deg+1))
				bdd[i] += log2(2.0*M_PI)
				bdd[i] += log2(float64(i) + err)
			}

			tot_deg += 1
		}

		deg[maxi] += 1
	}

	fmt.Println("==============================================")
	fmt.Println("==Degree Searching Result=====================")
	fmt.Println("==============================================")
	if iter == max_iter{
		fmt.Println("More Iterations Needed")
	}else{
		fmt.Println("Degree of Polynomial :", tot_deg-1)
		fmt.Println("Degree :", deg)
	}
	fmt.Println("==============================================")


	var inter_size = NewFloat(1.0 / float64(int(1<<dev)))

	_=inter_size

	var z = make([]*big.Float, tot_deg)
	var cnt int
	if deg[0]%2 != 0 {
		z[cnt] = NewFloat(0)
		cnt++
	}

	var tmp *big.Float

	for i := K-1; i > 0 ; i-- {
		for j := 1 ; j <= deg[i]; j++ {

			tmp = NewFloat(float64(2*j-1))
			tmp.Mul(tmp, PI)
			tmp.Quo(tmp, NewFloat(float64(2*deg[i])))
			tmp = Cos(tmp)
			tmp.Mul(tmp, inter_size)

			z[cnt] = NewFloat(float64(i))
			z[cnt].Add(z[cnt], tmp)
			cnt++

			z[cnt] = NewFloat(float64(-i))
			z[cnt].Sub(z[cnt], tmp)
			cnt++

		}
	}

	for j := 1; j <= deg[0]/2; j++ {
		tmp = NewFloat(float64(2*j-1))
		tmp.Mul(tmp, PI)
		tmp.Quo(tmp, NewFloat(float64(2*deg[0])))
		tmp = Cos(tmp)
		tmp.Mul(tmp, inter_size)

		z[cnt] = new(big.Float).Add(z[cnt], tmp)
		cnt++

		z[cnt] = new(big.Float).Sub(z[cnt], tmp)
		cnt++
	}

	

	// cos(2*pi*(x-0.25)/r)
	var d = make([]*big.Float, tot_deg)
	for i := 0; i < tot_deg; i++ {

		d[i] = NewFloat(2.0)
		d[i].Mul(d[i], PI)

		z[i].Sub(z[i], NewFloat(0.25))
		z[i].Quo(z[i], sc_fac)

		d[i].Mul(d[i], z[i])
		d[i] = Cos(d[i])

		
	}

	for j := 1 ; j < tot_deg; j++ {
		for l := 0; l < tot_deg-j; l++{
			d[l].Sub(d[l+1], d[l])
			tmp.Sub(z[l+j], z[l])
			d[l].Quo(d[l], tmp)
		}
	}

	tot_deg += 1

	var x = make([]*big.Float, tot_deg)
	for i := 0; i < tot_deg; i++ {
		x[i] = NewFloat(float64(K))
		x[i].Quo(x[i], sc_fac)
		tmp.Mul(NewFloat(float64(i)), PI)
		tmp.Quo(tmp, NewFloat(float64(tot_deg-1)))
		x[i].Mul(x[i], Cos(tmp))
	}

	var c = make([]*big.Float, tot_deg)
	var p = make([]*big.Float, tot_deg)
	for i := 0 ; i < tot_deg; i++ {
		p[i] = new(big.Float).Copy(d[0])
		for j := 1 ; j < tot_deg-1; j++ {
			tmp.Sub(x[i], z[j])
			p[i].Mul(p[i], tmp)
			p[i].Add(p[i], d[j])
		}
	}

	

	

	var T = make([][]*big.Float, tot_deg)
	for i := 0 ; i < tot_deg; i++ {
		T[i] = make([]*big.Float, tot_deg)
	}

	for i := 0 ; i < tot_deg; i++ {

		T[i][0] = NewFloat(1.0)

		T[i][1] = new(big.Float).Copy(x[i])

		tmp.Quo(NewFloat(float64(K)), sc_fac)

		T[i][1].Quo(T[i][1], tmp)

		for j := 2 ; j < tot_deg; j++ {

			T[i][j] = NewFloat(2.0)

			tmp.Quo(NewFloat(float64(K)), sc_fac)
			tmp.Quo(x[i], tmp)
			T[i][j].Mul(T[i][j], tmp)
			T[i][j].Mul(T[i][j], T[i][j-1])
			T[i][j].Sub(T[i][j], T[i][j-2])
		}

		
	}

	var max_abs = new(big.Float)
	var max_index int
	for i := 0; i < tot_deg-1; i++ {
		max_abs.Abs(T[i][i])
		max_index = i
		for j := i+1; j < tot_deg; j++ {
			tmp.Abs(T[j][i])
			if tmp.Cmp(max_abs) == 1 {
				max_abs.Abs(T[j][i])
				max_index = j
			}
		}

		if i != max_index {
			for j := i ; j < tot_deg; j++ {
				tmp.Copy(T[max_index][j])
				T[max_index][j].Set(T[i][j])
				T[i][j].Set(tmp)
			}

			tmp.Set(p[max_index])
			p[max_index].Set(p[i])
			p[i].Set(tmp)
		}

		for j := i+1; j < tot_deg; j++ {
			T[i][j].Quo(T[i][j], T[i][i])
		}

		p[i].Quo(p[i], T[i][i])
		T[i][i] = NewFloat(1.0)

		for j := i+1; j < tot_deg; j++ {
			tmp.Mul(T[j][i], p[i])
			p[j].Sub(p[j], tmp)
			for l := i+1; l<tot_deg; l++ {
				tmp.Mul(T[j][i], T[i][l])
				T[j][l].Sub(T[j][l], tmp)
			}
			T[j][i] = NewFloat(0.0)
		}
	}

	c[tot_deg-1] = p[tot_deg-1]
	for i := tot_deg-2; i>=0; i--{
		c[i] = new(big.Float)
		c[i].Copy(p[i])
		for j := i+1; j < tot_deg; j++ {
			tmp.Mul(T[i][j], c[j])
			c[i].Sub(c[i], tmp)
		}
	}

	tot_deg -= 1

	res := make([]float64, tot_deg)
	fmt.Printf("[")
	for i := 0 ; i < tot_deg; i++ {
		res[i], _ = c[i].Float64()
		fmt.Printf("%.20f, ", res[i])
	}
	fmt.Printf("]\n")

	
	

}