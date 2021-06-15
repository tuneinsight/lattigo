package dckks

import (
	"github.com/ldsec/lattigo/v2/ring"
	"math/big"
)

func extendBasisSmallNormAndCenter(ringQ, ringP *ring.Ring, polQ, polP *ring.Poly) {
	var coeff, Q, QHalf, sign uint64
	Q = ringQ.Modulus[0]
	QHalf = Q >> 1

	for j := 0; j < ringQ.N; j++ {

		coeff = polQ.Coeffs[0][j]

		sign = 1
		if coeff > QHalf {
			coeff = Q - coeff
			sign = 0
		}

		for i, pi := range ringP.Modulus {
			polP.Coeffs[i][j] = (coeff * sign) | (pi-coeff)*(sign^1)
		}
	}
}

func centerAndExtendBasisLargeNorm(minLevel, maxLevel int, ringQ *ring.Ring, polIn *ring.Poly, vBigint []*big.Int, polOut *ring.Poly) {

	// Bigint modulus at minimum level
	QminLevel := ring.NewUint(ringQ.Modulus[0])
	for i := 1; i < minLevel+1; i++ {
		QminLevel.Mul(QminLevel, ring.NewUint(ringQ.Modulus[i]))
	}

	QminLevelHalf := new(big.Int).Rsh(QminLevel, 1)

	// Sets the coefficients outside of the NTT domain
	ringQ.InvNTTLvl(minLevel, polIn, polOut)

	// Sets the coefficients outside of the RNS domain
	ringQ.PolyToBigintLvl(minLevel, polOut, vBigint)

	// Centers the coefficients
	var sign int
	for i := 0; i < ringQ.N; i++ {
		sign = vBigint[i].Cmp(QminLevelHalf)
		if sign == 1 || sign == 0 {
			vBigint[i].Sub(vBigint[i], QminLevel)
		}
	}

	// Sets the coefficients back to the RNS domain
	ringQ.SetCoefficientsBigintLvl(maxLevel, vBigint, polOut)

	// Sets the coefficients back to the RNS NTT domain
	ringQ.NTTLvl(maxLevel, polOut, polOut)
}
