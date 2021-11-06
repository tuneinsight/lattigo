package lwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"math/big"
)

// LWEToRLWE transforms a set of LWE samples into their respective RLWE ciphertext such that decrypt(RLWE)[0] = decrypt(LWE)
func LWEToRLWE(ctLWE []*Ciphertext, params rlwe.Parameters) (ctRLWE []*rlwe.Ciphertext) {

	level := ctLWE[0].Level()

	ringQ := params.RingQ()
	acc := ringQ.NewPolyLvl(level)
	ctRLWE = make([]*rlwe.Ciphertext, len(ctLWE))
	for i := 0; i < len(ctLWE); i++ {

		// Alocates ciphertext
		ctRLWE[i] = rlwe.NewCiphertextNTT(params, 1, level)

		for u := 0; u < level+1; u++ {

			ctRLWE[i].Value[0].Coeffs[u][0] = ctLWE[i].Value[u][0]

			// Copy coefficients multiplied by X^{N-1} in reverse order:
			// a_{0} -a_{N-1} -a2_{N-2} ... -a_{1}
			tmp0, tmp1 := acc.Coeffs[u], ctLWE[i].Value[u][1:]
			tmp0[0] = tmp1[0]
			for k := 1; k < ringQ.N; k++ {
				tmp0[k] = ringQ.Modulus[u] - tmp1[ringQ.N-k]
			}

			copy(ctRLWE[i].Value[1].Coeffs[u], acc.Coeffs[u])
		}

		// Switches to NTT domain
		ringQ.NTTLvl(level, ctRLWE[i].Value[0], ctRLWE[i].Value[0])
		ringQ.NTTLvl(level, ctRLWE[i].Value[1], ctRLWE[i].Value[1])
	}

	return
}

func ExtractLWEFromRLWESingle(ct *rlwe.Ciphertext, ringQ *ring.Ring) (lwe *Ciphertext) {

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

// RLWEToLWE extracts LWE samples from a RLWE sample
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
