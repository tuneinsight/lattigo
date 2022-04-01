package bfv

import (
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// Scaler is an interface that rescales polynomial coefficients by a fraction t/Q.
type Scaler interface {
	// DivByQOverTRoundedLvl returns p1 scaled by a factor t/Q and mod t on the receiver p2.
	DivByQOverTRoundedLvl(level int, p1, p2 *ring.Poly)
}

// RNSScaler implements the Scaler interface by performing a scaling by t/Q in the RNS domain.
// This implementation of the Scaler interface is preferred over the SimpleScaler implementation.
type RNSScaler struct {
	ringQ, ringT *ring.Ring

	polypoolQ *ring.Poly
	polypoolT *ring.Poly

	qHalf     []*big.Int // (q-1)/2
	qHalfModT []uint64   // (q-1)/2 mod t
	qInv      []uint64   //(q mod t)^-1 mod t

	paramsQP []ring.ModupParams

	tDividesQ bool
}

// NewRNSScaler creates a new SimpleScaler from t, the modulus under which the reconstruction is returned, the Ring in which the polynomial to reconstruct is represented.
func NewRNSScaler(ringQ *ring.Ring, T uint64) (rnss *RNSScaler) {

	if utils.IsInSliceUint64(T, ringQ.Modulus) && ringQ.Modulus[0] != T {
		panic("T must be Q[0] if T|Q")
	}

	rnss = new(RNSScaler)

	rnss.ringQ = ringQ

	rnss.polypoolQ = ringQ.NewPoly()

	rnss.ringT = new(ring.Ring)
	rnss.ringT.N = ringQ.N
	rnss.ringT.Modulus = []uint64{T}
	rnss.ringT.MredParams = []uint64{ring.MRedParams(T)}
	rnss.polypoolT = rnss.ringT.NewPoly()

	rnss.tDividesQ = T == ringQ.Modulus[0]

	if !rnss.tDividesQ {

		rnss.qHalf = make([]*big.Int, len(ringQ.Modulus))
		rnss.qInv = make([]uint64, len(ringQ.Modulus))
		rnss.qHalfModT = make([]uint64, len(ringQ.Modulus))
		rnss.paramsQP = make([]ring.ModupParams, len(ringQ.Modulus))

		bigQ := new(big.Int).SetUint64(1)
		tmp := new(big.Int)
		bredParams := ring.BRedParams(T)
		TBig := ring.NewUint(T)
		for i := range ringQ.Modulus {
			rnss.paramsQP[i] = ring.GenModUpParams(ringQ.Modulus[:i+1], rnss.ringT.Modulus)

			bigQ.Mul(bigQ, ring.NewUint(ringQ.Modulus[i]))

			rnss.qInv[i] = tmp.Mod(bigQ, TBig).Uint64()
			rnss.qInv[i] = ring.ModExp(rnss.qInv[i], T-2, T)
			rnss.qInv[i] = ring.MForm(rnss.qInv[i], T, bredParams)

			rnss.qHalf[i] = new(big.Int).Set(bigQ)
			rnss.qHalf[i].Rsh(rnss.qHalf[i], 1)

			rnss.qHalfModT[i] = tmp.Mod(rnss.qHalf[i], TBig).Uint64()
		}
	}

	return
}

// DivByQOverTRoundedLvl returns p1 scaled by a factor t/Q and mod t on the receiver p2.
func (rnss *RNSScaler) DivByQOverTRoundedLvl(level int, p1Q, p2T *ring.Poly) {

	ringQ := rnss.ringQ

	if level > 0 {
		if rnss.tDividesQ {

			ringQ.DivRoundByLastModulusManyLvl(level, level, p1Q, rnss.polypoolQ, p2T)

		} else {

			ringT := rnss.ringT
			T := ringT.Modulus[0]
			tInv := ringT.MredParams[0]
			p2tmp := p2T.Coeffs[0]
			p3tmp := rnss.polypoolT.Coeffs[0]
			qInv := T - rnss.qInv[level]
			qHalfModT := T - rnss.qHalfModT[level]

			// Multiply P_{Q} by t and extend the basis from P_{Q} to t*(P_{Q}||P_{t})
			// Since the coefficients of P_{t} are multiplied by t, they are all zero,
			// hence the basis extension can be omitted
			ringQ.MulScalarLvl(level, p1Q, T, rnss.polypoolQ)

			// Center t*P_{Q} around (Q-1)/2 to round instead of floor during the division
			ringQ.AddScalarBigintLvl(level, rnss.polypoolQ, rnss.qHalf[level], rnss.polypoolQ)

			// Extend the basis of (t*P_{Q} + (Q-1)/2) to (t*P_{t} + (Q-1)/2)
			ring.ModUpExact(rnss.polypoolQ.Coeffs[:level+1], rnss.polypoolT.Coeffs, ringQ, ringT, rnss.paramsQP[level])

			// Compute [Q^{-1} * (t*P_{t} -   (t*P_{Q} - ((Q-1)/2 mod t)))] mod t which returns round(t/Q * P_{Q}) mod t
			ring.AddScalarNoModAndMulScalarMontgomeryVec(p3tmp, p2tmp, qHalfModT, qInv, T, tInv)
		}
	} else {
		qOverT := float64(ringQ.Modulus[0]) / float64(rnss.ringT.Modulus[0])
		tmp0, tmp1 := p2T.Coeffs[0], p1Q.Coeffs[0]
		for i := 0; i < ringQ.N; i++ {
			tmp0[i] = uint64(float64(tmp1[i])/qOverT + 0.5)
		}
	}
}

// ScaleUpTCoprimeWithQVec takes a Poly pIn in ringT, scales its coefficients up by (Q/T) mod Q, and writes the result in a
// Poly pOut in ringQ.
func ScaleUpTCoprimeWithQVec(ringQ, ringT *ring.Ring, rescaleParams, tmp []uint64, pIn, pOut *ring.Poly) {

	qModTmontgomery := ring.MForm(new(big.Int).Mod(ringQ.ModulusBigint, ringT.ModulusBigint).Uint64(), ringT.Modulus[0], ringT.BredParams[0])

	t := ringT.Modulus[0]
	tHalf := t >> 1
	tInv := ringT.MredParams[0]

	// (x * Q + T/2) mod T
	ring.MulScalarMontgomeryAndAddScalarVec(pIn.Coeffs[0], tmp, tHalf, qModTmontgomery, t, tInv)

	// (x * T^-1 - T/2) mod Qi
	for i := 0; i < len(pOut.Coeffs); i++ {
		p0tmp := tmp
		p1tmp := pOut.Coeffs[i]
		qi := ringQ.Modulus[i]
		bredParams := ringQ.BredParams[i]
		mredParams := ringQ.MredParams[i]
		rescaleParams := qi - rescaleParams[i]
		tHalfNegQi := qi - ring.BRedAdd(tHalf, qi, bredParams)

		ring.AddScalarNoModAndMulScalarMontgomeryVec(p0tmp, p1tmp, tHalfNegQi, rescaleParams, qi, mredParams)
	}
}

// ScaleUpTDividesQVec takes a Poly pIn in ringT, scales its coefficients up by (Q/T) mod Q, and writes the result in a
// Poly pOut in ringQ.
// T is in this case assumed to be the first prime in the moduli chain.
func ScaleUpTDividesQVec(ringQ *ring.Ring, pIn, pOut *ring.Poly) {
	// Q/T mod T
	tmp := new(big.Int).Set(ringQ.ModulusBigint)
	tmp.Quo(tmp, new(big.Int).SetUint64(ringQ.Modulus[0]))
	QOverTMont := tmp.Mod(tmp, new(big.Int).SetUint64(ringQ.Modulus[0])).Uint64()

	// pOut = Q/T * pIn
	ring.MForm(QOverTMont, ringQ.Modulus[0], ringQ.BredParams[0])
	ring.MulScalarMontgomeryVec(pIn.Coeffs[0], pOut.Coeffs[0], QOverTMont, ringQ.Modulus[0], ringQ.MredParams[0])
}

// DecryptAndPrintError decrypts a ciphertext and prints the log2 of the error.
func DecryptAndPrintError(ptWant *Plaintext, cthave *Ciphertext, ringQ *ring.Ring, decryptor Decryptor) {

	level := cthave.Level()

	ringQ.SubLvl(level, cthave.Value[0], ptWant.Value, cthave.Value[0])
	plaintext := decryptor.DecryptNew(cthave)

	bigintCoeffs := make([]*big.Int, ringQ.N)
	ringQ.PolyToBigintLvl(level, plaintext.Value, 1, bigintCoeffs)

	Q := new(big.Int).SetUint64(1)
	for i := 0; i < level+1; i++ {
		Q.Mul(Q, new(big.Int).SetUint64(ringQ.Modulus[i]))
	}

	center(bigintCoeffs, Q)
	stdErr, minErr, maxErr := errorStats(bigintCoeffs)
	fmt.Printf("STD : %f - Min : %f - Max : %f\n", math.Log2(stdErr), math.Log2(minErr), math.Log2(maxErr))
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
