package lwe

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"math/big"
)

func normalizeInv(x, a, b float64) (y float64) {
	return (x*(b-a) + b + a) / 2.0
}

func scaleUp(value float64, scale float64, Q uint64) (res uint64) {

	var isNegative bool
	var xFlo *big.Float
	var xInt *big.Int

	isNegative = false
	if value < 0 {
		isNegative = true
		xFlo = big.NewFloat(-scale * value)
	} else {
		xFlo = big.NewFloat(scale * value)
	}

	xFlo.Add(xFlo, big.NewFloat(0.5))

	xInt = new(big.Int)
	xFlo.Int(xInt)
	xInt.Mod(xInt, ring.NewUint(Q))

	res = xInt.Uint64()

	if isNegative {
		res = Q - res
	}

	return
}

// DecryptLWE decrypts an LWE sample
func DecryptLWE(ct *Ciphertext, ringQ *ring.Ring, skMont *ring.Poly) float64 {

	level := ct.Level()

	pol := ringQ.NewPolyLvl(ct.Level())
	for i := 0; i < level+1; i++ {
		copy(pol.Coeffs[i], ct.Value[i][1:])
	}

	ringQ.MulCoeffsMontgomeryLvl(level, pol, skMont, pol)

	a := make([]uint64, level+1)

	for i := 0; i < level+1; i++ {
		qi := ringQ.Modulus[i]
		tmp := pol.Coeffs[i]
		a[i] = ct.Value[i][0]
		for j := 0; j < ringQ.N; j++ {
			a[i] = ring.CRed(a[i]+tmp[j], qi)
		}
	}

	crtReconstruction := make([]*big.Int, level+1)

	QiB := new(big.Int)
	tmp := new(big.Int)
	modulusBigint := ring.NewUint(1)

	for i := 0; i < level+1; i++ {

		qi := ringQ.Modulus[i]
		QiB.SetUint64(qi)

		modulusBigint.Mul(modulusBigint, QiB)

		crtReconstruction[i] = new(big.Int)
		crtReconstruction[i].Quo(ringQ.ModulusBigint, QiB)
		tmp.ModInverse(crtReconstruction[i], QiB)
		tmp.Mod(tmp, QiB)
		crtReconstruction[i].Mul(crtReconstruction[i], tmp)
	}

	tmp.SetUint64(0)
	coeffsBigint := ring.NewUint(0)

	modulusBigintHalf := new(big.Int)
	modulusBigintHalf.Rsh(modulusBigint, 1)

	var sign int
	for i := 0; i < level+1; i++ {
		coeffsBigint.Add(coeffsBigint, tmp.Mul(ring.NewUint(a[i]), crtReconstruction[i]))
	}

	coeffsBigint.Mod(coeffsBigint, modulusBigint)

	// Centers the coefficients
	sign = coeffsBigint.Cmp(modulusBigintHalf)

	if sign == 1 || sign == 0 {
		coeffsBigint.Sub(coeffsBigint, modulusBigint)
	}

	flo := new(big.Float)
	flo.SetInt(coeffsBigint)
	flo64, _ := flo.Float64()

	return flo64
}
