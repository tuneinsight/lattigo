package rlwe

import (
	"math"
	"math/big"
	"slices"

	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils"
)

// NoisePublicKey returns the log2 of the standard deviation of the input [PublicKey] with respect to the given [Secret] and parameters.
func NoisePublicKey(pk *PublicKey, sk *SecretKey, params Parameters) float64 {

	pk = pk.CopyNew()

	ringQP := params.RingQP().AtLevel(pk.LevelQ(), pk.LevelP())

	// [-as + e] + [as]
	ringQP.MulCoeffsMontgomeryThenAdd(sk.Value, pk.Value[1], pk.Value[0])
	ringQP.INTT(pk.Value[0], pk.Value[0])
	ringQP.IMForm(pk.Value[0], pk.Value[0])

	return ringQP.Log2OfStandardDeviation(pk.Value[0])
}

// NoiseRelinearizationKey the log2 of the standard deviation of the noise of the input [RelinearizationKey] with respect to the given [Secret] and paramters.
func NoiseRelinearizationKey(rlk *RelinearizationKey, sk *SecretKey, params Parameters) float64 {
	sk2 := sk.CopyNew()
	params.RingQP().AtLevel(rlk.LevelQ(), rlk.LevelP()).MulCoeffsMontgomery(sk2.Value, sk2.Value, sk2.Value)
	return NoiseEvaluationKey(&rlk.EvaluationKey, sk2, sk, params)
}

// NoiseGaloisKey the log2 of the standard deviation of the noise of the input [GaloisKey] key with respect to the given [SecretKey] and paramters.
func NoiseGaloisKey(gk *GaloisKey, sk *SecretKey, params Parameters) float64 {

	skIn := sk.CopyNew()
	skOut := sk.CopyNew()

	nthRoot := params.RingQ().NthRoot()

	galElInv := ring.ModExp(gk.GaloisElement, nthRoot-1, nthRoot)

	params.RingQP().AtLevel(gk.LevelQ(), gk.LevelP()).AutomorphismNTT(sk.Value, galElInv, skOut.Value)

	return NoiseEvaluationKey(&gk.EvaluationKey, skIn, skOut, params)
}

// NoiseGadgetCiphertext returns the log2 of the standard deviation of the noise of the input [GadgetCiphertext] with respect to the given [Plaintext], [SecretKey] and [Parameters].
// The polynomial pt is expected to be in the NTT and Montgomery domain.
func NoiseGadgetCiphertext(gct *GadgetCiphertext, pt ring.Poly, sk *SecretKey, params Parameters) float64 {

	gct = gct.CopyNew()
	pt = *pt.CopyNew()
	levelQ, levelP := gct.LevelQ(), gct.LevelP()
	ringQP := params.RingQP().AtLevel(levelQ, levelP)
	ringQ, ringP := ringQP.RingQ, ringQP.RingP
	BaseTwoDecompositionVectorSize := slices.Min(gct.BaseTwoDecompositionVectorSize()) // required else the check becomes very complicated

	// Decrypts
	// [-asIn + w*P*sOut + e, a] + [asIn]
	for i := range gct.Value {
		for j := range gct.Value[i] {
			ringQP.MulCoeffsMontgomeryThenAdd(gct.Value[i][j][1], sk.Value, gct.Value[i][j][0])
		}
	}

	// Sums all bases together (equivalent to multiplying with CRT decomposition of 1)
	// sum([1]_w * [RNS*PW2*P*sOut + e]) = PWw*P*sOut + sum(e)
	for i := range gct.Value { // RNS decomp
		if i > 0 {
			for j := range gct.Value[i] { // PW2 decomp
				ringQP.Add(gct.Value[0][j][0], gct.Value[i][j][0], gct.Value[0][j][0])
			}
		}
	}

	if levelP != -1 {
		// sOut * P
		ringQ.MulScalarBigint(pt, ringP.ModulusAtLevel[levelP], pt)
	}

	var maxLog2Std float64

	for i := 0; i < BaseTwoDecompositionVectorSize; i++ {

		// P*s^i + sum(e) - P*s^i = sum(e)
		ringQ.Sub(gct.Value[0][i][0].Q, pt, gct.Value[0][i][0].Q)

		// Checks that the error is below the bound
		// Worst error bound is N * floor(6*sigma) * #Keys
		ringQP.INTT(gct.Value[0][i][0], gct.Value[0][i][0])
		ringQP.IMForm(gct.Value[0][i][0], gct.Value[0][i][0])

		maxLog2Std = utils.Max(maxLog2Std, ringQP.Log2OfStandardDeviation(gct.Value[0][i][0]))

		// sOut * P * PW2
		ringQ.MulScalar(pt, 1<<gct.BaseTwoDecomposition, pt)
	}

	return maxLog2Std
}

// NoiseEvaluationKey the log2 of the standard deviation of the noise of the input [GaloisKey] with respect to the given [SecretKey] and [Parameters].
func NoiseEvaluationKey(evk *EvaluationKey, skIn, skOut *SecretKey, params Parameters) float64 {
	return NoiseGadgetCiphertext(&evk.GadgetCiphertext, skIn.Value.Q, skOut, params)
}

// Norm returns the log2 of the standard deviation, minimum and maximum absolute norm of
// the decrypted [Ciphertext], before the decoding (i.e. including the error).
func Norm(ct *Ciphertext, dec *Decryptor) (std, min, max float64) {

	params := dec.params

	coeffsBigint := make([]*big.Int, params.N())
	for i := range coeffsBigint {
		coeffsBigint[i] = new(big.Int)
	}

	pt := NewPlaintext(params, ct.Level())

	dec.Decrypt(ct, pt)

	ringQ := params.RingQ().AtLevel(ct.Level())

	if pt.IsNTT {
		ringQ.INTT(pt.Value, pt.Value)
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

// NTTSparseAndMontgomery takes a polynomial Z[Y] outside of the NTT domain and maps it to a polynomial Z[X] in the NTT domain where Y = X^(gap).
// This method is used to accelerate the NTT of polynomials that encode sparse polynomials.
func NTTSparseAndMontgomery(r *ring.Ring, metadata *MetaData, pol ring.Poly) {

	if 1<<metadata.LogDimensions.Cols == r.NthRoot()>>2 {

		if metadata.IsNTT {
			r.NTT(pol, pol)
		}

		if metadata.IsMontgomery {
			r.MForm(pol, pol)
		}

	} else {

		var n int
		var NTT func(p1, p2 []uint64, N int, Q, QInv uint64, BRedConstant [2]uint64, nttPsi []uint64)
		switch r.Type() {
		case ring.Standard:
			n = 2 << metadata.LogDimensions.Cols
			NTT = ring.NTTStandard
		case ring.ConjugateInvariant:
			n = 1 << metadata.LogDimensions.Cols
			NTT = ring.NTTConjugateInvariant
		}

		N := r.N()
		gap := N / n
		for i, s := range r.SubRings[:r.Level()+1] {

			coeffs := pol.Coeffs[i]

			if metadata.IsMontgomery {
				s.MForm(coeffs[:n], coeffs[:n])
			}

			if metadata.IsNTT {
				// NTT in dimension n but with roots of N
				// This is a small hack to perform at reduced cost an NTT of dimension N on a vector in Y = X^{N/n}, i.e. sparse polynomials.
				NTT(coeffs[:n], coeffs[:n], n, s.Modulus, s.MRedConstant, s.BRedConstant, s.RootsForward)

				// Maps NTT in dimension n to NTT in dimension N
				for j := n - 1; j >= 0; j-- {
					c := coeffs[j]
					for w := 0; w < gap; w++ {
						coeffs[j*gap+w] = c
					}
				}
			} else {
				for j := n - 1; j >= 0; j-- {
					coeffs[j*gap] = coeffs[j]
					for j := 1; j < gap; j++ {
						coeffs[j*gap-j] = 0
					}
				}
			}
		}
	}
}

// ExtendBasisSmallNormAndCenterNTTMontgomery extends a small-norm polynomial polQ in R_Q to a polynomial
// polP in R_P.
// This method can be used to extend from Q0 to QL.
// Input and output are in the NTT and Montgomery domain.
func ExtendBasisSmallNormAndCenterNTTMontgomery(rQ, rP *ring.Ring, polQ, buff, polP ring.Poly) {
	rQ = rQ.AtLevel(0)

	levelP := rP.Level()

	// Switches Q[0] out of the NTT and Montgomery domain.
	rQ.INTT(polQ, buff)
	rQ.IMForm(buff, buff)

	// Reconstruct P from Q
	Q := rQ.SubRings[0].Modulus
	QHalf := Q >> 1

	P := rP.ModuliChain()
	N := rQ.N()

	var sign uint64
	for j := 0; j < N; j++ {

		coeff := buff.Coeffs[0][j]

		sign = 1
		if coeff > QHalf {
			coeff = Q - coeff
			sign = 0
		}

		for i := 0; i < levelP+1; i++ {
			polP.Coeffs[i][j] = (coeff * sign) | (P[i]-coeff)*(sign^1)
		}
	}

	rP.NTT(polP, polP)
	rP.MForm(polP, polP)
}
