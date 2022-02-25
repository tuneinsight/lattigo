package lwe

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
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
func RLWEToLWE(ct *rlwe.Ciphertext, ringQ *ring.Ring, slotIndex map[int]bool) (LWE map[int]*Ciphertext) {

	LWE = make(map[int]*Ciphertext)

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

	var prevIndex int
	for index := 0; index < ringQ.N; index++ {

		if _, ok := slotIndex[index]; ok {

			// Multiplies the accumulator by X^{N/(2*slots)}
			MulBySmallMonomial(ringQ, acc, index-prevIndex)
			prevIndex = index

			LWE[index] = NewCiphertext(ringQ.N, level)

			for j := 0; j < level+1; j++ {
				LWE[index].Value[j][0] = c0.Coeffs[j][index]
				copy(LWE[index].Value[j][1:], acc.Coeffs[j])
			}
		}
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
