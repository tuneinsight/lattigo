package rlwe

import (
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// PublicKeyIsCorrect returns true if pk is a correct RLWE public-key for secret-key sk and parameters params.
func PublicKeyIsCorrect(pk *PublicKey, sk *SecretKey, params Parameters, log2Bound int) bool {

	pk = pk.CopyNew()

	levelQ, levelP := params.QCount()-1, params.PCount()-1
	ringQP := params.RingQP().AtLevel(levelQ, levelP)

	// [-as + e] + [as]
	ringQP.MulCoeffsMontgomeryAndAdd(sk.Value, pk.Value[1], pk.Value[0])
	ringQP.InvNTT(pk.Value[0], pk.Value[0])
	ringQP.InvMForm(pk.Value[0], pk.Value[0])

	if log2Bound <= ringQP.RingQ.Log2OfInnerSum(pk.Value[0].Q) {
		return false
	}

	if ringQP.RingP != nil && log2Bound <= ringQP.RingP.Log2OfInnerSum(pk.Value[0].P) {
		return false
	}

	return true
}

// RelinearizationKeyIsCorrect returns true if swk is a correct RLWE relinearization-key for secret-key sk and parameters params.
func RelinearizationKeyIsCorrect(rlk *SwitchingKey, skIdeal *SecretKey, params Parameters, log2Bound int) bool {
	levelQ, levelP := params.QCount()-1, params.PCount()-1
	skIn := skIdeal.CopyNew()
	skOut := skIdeal.CopyNew()
	params.RingQP().AtLevel(levelQ, levelP).MulCoeffsMontgomery(skIn.Value, skIn.Value, skIn.Value)
	return SwitchingKeyIsCorrect(rlk, skIn, skOut, params, log2Bound)
}

// RotationKeyIsCorrect returns true if swk is a correct RLWE switching-key for galois element galEl, secret-key sk and parameters params.
func RotationKeyIsCorrect(swk *SwitchingKey, galEl uint64, skIdeal *SecretKey, params Parameters, log2Bound int) bool {
	swk = swk.CopyNew()
	skIn := skIdeal.CopyNew()
	skOut := skIdeal.CopyNew()
	galElInv := ring.ModExp(galEl, uint64(4*params.N()-1), uint64(4*params.N()))
	ringQ, ringP := params.RingQ(), params.RingP()

	ringQ.PermuteNTT(skIdeal.Value.Q, galElInv, skOut.Value.Q)
	if ringP != nil {
		ringP.PermuteNTT(skIdeal.Value.P, galElInv, skOut.Value.P)
	}

	return SwitchingKeyIsCorrect(swk, skIn, skOut, params, log2Bound)
}

// SwitchingKeyIsCorrect returns true if swk is a correct RLWE switching-key for input key skIn, output key skOut and parameters params.
func SwitchingKeyIsCorrect(swk *SwitchingKey, skIn, skOut *SecretKey, params Parameters, log2Bound int) bool {
	swk = swk.CopyNew()
	skIn = skIn.CopyNew()
	skOut = skOut.CopyNew()
	levelQ, levelP := params.QCount()-1, params.PCount()-1
	ringQP := params.RingQP().AtLevel(levelQ, levelP)
	ringQ, ringP := ringQP.RingQ, ringQP.RingP
	decompPw2 := params.DecompPw2(levelQ, levelP)

	// Decrypts
	// [-asIn + w*P*sOut + e, a] + [asIn]
	for i := range swk.Value {
		for j := range swk.Value[i] {
			ringQP.MulCoeffsMontgomeryAndAdd(swk.Value[i][j].Value[1], skOut.Value, swk.Value[i][j].Value[0])
		}
	}

	// Sums all bases together (equivalent to multiplying with CRT decomposition of 1)
	// sum([1]_w * [RNS*PW2*P*sOut + e]) = PWw*P*sOut + sum(e)
	for i := range swk.Value { // RNS decomp
		if i > 0 {
			for j := range swk.Value[i] { // PW2 decomp
				ringQP.Add(swk.Value[0][j].Value[0], swk.Value[i][j].Value[0], swk.Value[0][j].Value[0])
			}
		}
	}

	if levelP != -1 {
		// sOut * P
		ringQ.MulScalarBigint(skIn.Value.Q, ringP.Modulus(), skIn.Value.Q)
	}

	for i := 0; i < decompPw2; i++ {

		// P*s^i + sum(e) - P*s^i = sum(e)
		ringQ.Sub(swk.Value[0][i].Value[0].Q, skIn.Value.Q, swk.Value[0][i].Value[0].Q)

		// Checks that the error is below the bound
		// Worst error bound is N * floor(6*sigma) * #Keys
		ringQP.InvNTT(swk.Value[0][i].Value[0], swk.Value[0][i].Value[0])
		ringQP.InvMForm(swk.Value[0][i].Value[0], swk.Value[0][i].Value[0])

		// Worst bound of inner sum
		// N*#Keys*(N * #Parties * floor(sigma*6) + #Parties * floor(sigma*6) + N * #Parties  +  #Parties * floor(6*sigma))

		if log2Bound < ringQ.Log2OfInnerSum(swk.Value[0][i].Value[0].Q) {
			return false
		}

		if levelP != -1 {
			if log2Bound < ringP.Log2OfInnerSum(swk.Value[0][i].Value[0].P) {
				return false
			}
		}

		// sOut * P * PW2
		ringQ.MulScalar(skIn.Value.Q, 1<<params.Pow2Base(), skIn.Value.Q)
	}

	return true
}

// Norm returns the log2 of the standard deviation, minimum and maximum absolute norm of
// the decrypted Ciphertext, before the decoding (i.e. including the error).
func Norm(ct *Ciphertext, dec Decryptor) (std, min, max float64) {

	params := dec.(*decryptor).params

	coeffsBigint := make([]*big.Int, params.N())
	for i := range coeffsBigint {
		coeffsBigint[i] = new(big.Int)
	}

	pt := NewPlaintext(params, ct.Level())

	dec.Decrypt(ct, pt)

	ringQ := params.RingQ().AtLevel(ct.Level())

	if pt.IsNTT {
		ringQ.InvNTT(pt.Value, pt.Value)
	}

	ringQ.PolyToBigintCentered(pt.Value, 1, coeffsBigint)

	return NormStats(coeffsBigint)
}

func NormStats(vec []*big.Int) (float64, float64, float64) {

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
