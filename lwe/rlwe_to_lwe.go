package lwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// RLWEToLWESingle extract the first coefficient of the input RLWE and returns it
// as a LWE sample.
func RLWEToLWESingle(ct *rlwe.Ciphertext, ringQ *ring.Ring) (lwe *Ciphertext) {

	level := ct.Level()

	c0 := ringQ.NewPolyLvl(level)
	c1 := ringQ.NewPolyLvl(level)
	acc := ringQ.NewPolyLvl(level)

	ringQ.InvNTTLvl(level, ct.Value[0], c0)
	ringQ.InvNTTLvl(level, ct.Value[1], c1)

	// Copy coefficients multiplied by X^{N-1} in reverse order:
	// a_{0} -a_{N-1} -a2_{N-2} ... -a_{1}
	for i, qi := range ringQ.Modulus[:level+1] {
		tmp0 := acc.Coeffs[i]
		tmp1 := c1.Coeffs[i]
		tmp0[0] = tmp1[0]
		for j := 1; j < ringQ.N; j++ {
			tmp0[j] = qi - tmp1[ringQ.N-j]
		}
	}

	N := ringQ.N

	lwe = NewCiphertext(N, level)

	for j := 0; j < level+1; j++ {
		lwe.Value[j][0] = c0.Coeffs[j][0]
		copy(lwe.Value[j][1:], acc.Coeffs[j])
	}

	return
}

// RLWEToLWE extracts all LWE samples from a RLWE ciphertext.
func RLWEToLWE(ct *rlwe.Ciphertext, ringQ *ring.Ring, logSlots int) (LWE []*Ciphertext) {

	n := 1 << logSlots

	LWE = make([]*Ciphertext, n)

	level := ct.Level()

	c0 := ringQ.NewPolyLvl(level)
	c1 := ringQ.NewPolyLvl(level)
	acc := ringQ.NewPolyLvl(level)

	ringQ.InvNTTLvl(level, ct.Value[0], c0)
	ringQ.InvNTTLvl(level, ct.Value[1], c1)

	// Copy coefficients multiplied by X^{N-1} in reverse order:
	// a_{0} -a_{N-1} -a2_{N-2} ... -a_{1}
	for i, qi := range ringQ.Modulus[:level+1] {
		tmp0 := acc.Coeffs[i]
		tmp1 := c1.Coeffs[i]
		tmp0[0] = tmp1[0]
		for j := 1; j < ringQ.N; j++ {
			tmp0[j] = qi - tmp1[ringQ.N-j]
		}
	}

	N := ringQ.N
	gap := N / n

	// Real values
	for i, idx := 0, 0; i < n; i, idx = i+1, idx+gap {

		LWE[i] = NewCiphertext(N, level)

		for j := 0; j < level+1; j++ {
			LWE[i].Value[j][0] = c0.Coeffs[j][idx]
			copy(LWE[i].Value[j][1:], acc.Coeffs[j])
		}

		// Multiplies the accumulator by X^{N/(2*slots)}
		MulBySmallMonomial(ringQ, acc, gap)
	}

	return
}

//MulBySmallMonomial multiplies pol by x^n
func MulBySmallMonomial(ringQ *ring.Ring, pol *ring.Poly, n int) {
	for i, qi := range ringQ.Modulus[:pol.Level()+1] {
		pol.Coeffs[i] = append(pol.Coeffs[i][ringQ.N-n:], pol.Coeffs[i][:ringQ.N-n]...)
		tmp := pol.Coeffs[i]
		for j := 0; j < n; j++ {
			tmp[j] = qi - tmp[j]
		}
	}
}
