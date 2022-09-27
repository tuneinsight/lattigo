package bfv

import (
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v3/rlwe"
)

// Noise decrypts a ciphertext and returns the log2
// of the standard deviation, minimum and maximum norm of the noise
// assuming the decryption is correct.
// This function is used for testing/profiling/evaluation purposes
func Noise(cthave *Ciphertext, dec Decryptor) (std, min, max float64) {

	level := cthave.Level()

	params := dec.(*decryptor).params
	ringQ := params.RingQ()

	ecd := NewEncoder(params).(*encoder)

	pt := &Plaintext{Plaintext: &rlwe.Plaintext{Value: ecd.tmpPoly}}

	dec.Decrypt(cthave, pt)

	ecd.ScaleDown(pt, ecd.tmpPtRt)
	ecd.ScaleUp(ecd.tmpPtRt, pt)

	ringQ.SubLvl(level, cthave.Value[0], pt.Value, cthave.Value[0])

	dec.Decrypt(cthave, pt)

	bigintCoeffs := make([]*big.Int, ringQ.N)
	ringQ.PolyToBigintLvl(level, pt.Value, 1, bigintCoeffs)

	Q := new(big.Int).SetUint64(1)
	for i := 0; i < level+1; i++ {
		Q.Mul(Q, new(big.Int).SetUint64(ringQ.Modulus[i]))
	}

	center(bigintCoeffs, Q)
	stdErr, minErr, maxErr := errorStats(bigintCoeffs)
	return math.Log2(stdErr), math.Log2(minErr), math.Log2(maxErr)
}

func errorStats(vec []*big.Int) (float64, float64, float64) {

	vecfloat := make([]*big.Float, len(vec))
	minErr := new(big.Float).SetFloat64(0)
	maxErr := new(big.Float).SetFloat64(0)
	tmp := new(big.Float)
	minErr.SetInt(vec[0])
	minErr.Abs(minErr)
	for i := range vec {
		vecfloat[i] = new(big.Float)
		vecfloat[i].SetInt(vec[i])

		tmp.Abs(vecfloat[i])

		if minErr.Cmp(tmp) == 1 {
			minErr.Set(tmp)
		}

		if maxErr.Cmp(tmp) == -1 {
			maxErr.Set(tmp)
		}
	}

	n := new(big.Float).SetFloat64(float64(len(vec)))

	mean := new(big.Float).SetFloat64(0)

	for _, c := range vecfloat {
		mean.Add(mean, c)
	}

	mean.Quo(mean, n)

	err := new(big.Float).SetFloat64(0)
	for _, c := range vecfloat {
		tmp.Sub(c, mean)
		tmp.Mul(tmp, tmp)
		err.Add(err, tmp)
	}

	err.Quo(err, n)
	err.Sqrt(err)

	x, _ := err.Float64()
	y, _ := minErr.Float64()
	z, _ := maxErr.Float64()

	return x, y, z

}

func center(coeffs []*big.Int, Q *big.Int) {
	qHalf := new(big.Int)
	qHalf.Set(Q)
	qHalf.Rsh(qHalf, 1)
	var sign int
	for i := range coeffs {
		coeffs[i].Mod(coeffs[i], Q)
		sign = coeffs[i].Cmp(qHalf)
		if sign == 1 || sign == 0 {
			coeffs[i].Sub(coeffs[i], Q)
		}
	}
}
