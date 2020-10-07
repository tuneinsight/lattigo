package ring

import (
	//"fmt"
	"math/big"
	"testing"
	//"time"
)

var pi = "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491298336733624406566430860213949463952247371907021798609437027705392171762931767523846748184676694051320005681271452635608277857713427577896091736371787214684409012249534301465495853710507922796892589235420199561121290219608640344181598136297747713099605187072113499999983729780499510597317328160963185950244594553469083026425223082533446850352619311881710100031378387528865875332083814206171776691473035982534904287554687311595628638823537875937519577818577805321712268066130019278766111959092164201989"

func Test_Complex128_Mul(t *testing.T) {

	prec := uint64(128)

	var PI = new(big.Float)
	PI.SetPrec(uint(prec))
	PI.SetString(pi)

	var PIHalf = new(big.Float)
	PIHalf.SetPrec(uint(prec))
	PIHalf.SetString(pi)
	PIHalf.Quo(PIHalf, NewFloat(2.0, prec))

	cMul := NewComplexMultiplier()

	LogN := 13
	N := uint64(1 << LogN)

	m := uint64(2 << LogN)

	rotGroup := make([]uint64, m>>1)
	fivePows := uint64(1)
	for i := uint64(0); i < m>>2; i++ {
		rotGroup[i] = fivePows
		fivePows *= 5
		fivePows &= (m - 1)
	}

	//time0 := time.Now()
	var angle *big.Float
	roots := make([]*Complex, m+1)
	for i := uint64(0); i < m; i++ {
		angle = NewFloat(2, prec)
		angle.Mul(angle, PI)
		angle.Mul(angle, NewFloat(float64(i), prec))
		angle.Quo(angle, NewFloat(float64(m), prec))

		real := Cos(angle)

		angle.Sub(PIHalf, angle)

		imag := Cos(angle)

		roots[i] = NewComplex(real, imag)
	}
	//time1 := time.Now()

	//fmt.Printf("Computed roots in %s sec \n", time1.Sub(time0))

	//fmt.Println(roots[1][0])

	values := make([]*Complex, N)

	for i := range values {
		values[i] = NewComplex(NewFloat(1.0, prec), NewFloat(0.0, prec))
	}

	//time0 = time.Now()
	var lenh, lenq, gap, idx uint64
	u := NewComplex(nil, nil)
	v := NewComplex(nil, nil)

	for len := N; len >= 1; len >>= 1 {
		for i := uint64(0); i < N; i += len {

			lenh = len >> 1
			lenq = len << 2

			gap = m / lenq

			for j := uint64(0); j < lenh; j++ {

				idx = (lenq - (rotGroup[j] % lenq)) * gap

				u.Add(values[i+j], values[i+j+lenh])
				v.Sub(values[i+j], values[i+j+lenh])

				cMul.Mul(v, roots[idx], v)

				values[i+j].Set(u)
				values[i+j+lenh].Set(v)

				//fmt.Println(values[0][0].Prec())
			}
		}
	}

	NBig := NewFloat(float64(N), prec)
	for i := range values {
		values[i][0].Quo(values[i][0], NBig)
		values[i][1].Quo(values[i][1], NBig)
	}

	//time1 = time.Now()

	//fmt.Printf("Computed FFT for LogSlots = %d in %s sec \n", N, time1.Sub(time0))

	//fmt.Println(values[0][0])
	//fmt.Println(values[0][1])
}
