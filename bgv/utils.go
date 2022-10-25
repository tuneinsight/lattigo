package bgv

import (
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// Norm returns the log2 of the standard deviation, minimum and maximum absolute norm of
// the decrypted ciphertext, before the decoding (i.e. including the error).
func Norm(ct *rlwe.Ciphertext, dec Decryptor) (std, min, max float64) {

	params := dec.(*decryptor).params

	coeffsBigint := make([]*big.Int, params.N())
	for i := range coeffsBigint {
		coeffsBigint[i] = new(big.Int)
	}

	buffQ := dec.(*decryptor).buffQ
	pt := rlwe.NewPlaintextNTTAtLevelFromPoly(ct.Level(), buffQ)
	dec.(*decryptor).Decryptor.Decrypt(ct.El(), pt)
	params.RingQ().InvNTTLvl(ct.Level(), buffQ, buffQ)
	params.RingQ().PolyToBigintCenteredLvl(ct.Level(), buffQ, 1, coeffsBigint)

	return errorStats(coeffsBigint)
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

	return math.Log2(x), math.Log2(y), math.Log2(z)
}
