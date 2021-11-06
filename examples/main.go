package main

import (
	"fmt"
	"github.com/ldsec/lattigo/v2/lwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
	"time"
	//"math"
	"math/big"
)

// LogN of the ring degree of the used parameters
var LogN = 15

func main() {
	LUT()
}

func ScaleUpExact(value float64, scale float64, Q uint64) (res uint64) {

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

func InitLUT(scale float64, ringQ *ring.Ring) (F *ring.Poly) {
	F = ringQ.NewPoly()
	Q := ringQ.Modulus[0]

	// Zero
	F.Coeffs[0][0] = ScaleUpExact(0, scale, Q)

	// Negatives
	for i := 1; i < (ringQ.N>>1)+1; i++ {
		F.Coeffs[0][i] = ScaleUpExact(-float64(1), scale, Q)
	}

	// Positives
	for i := (ringQ.N >> 1) + 1; i < ringQ.N; i++ {
		F.Coeffs[0][i] = ScaleUpExact(-float64(1), scale, Q)
	}

	ringQ.NTT(F, F)

	return
}

func LUT() {
	var params, paramsLWE rlwe.Parameters
	var err error
	if params, err = rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		LogN:     11,
		Q:        []uint64{0x3fffffffef8001}, // 56
		P:        []uint64{0x80000000068001}, // 60 bits
		Sigma:    rlwe.DefaultSigma,
		RingType: rlwe.RingStandard,
	}); err != nil {
		panic(err)
	}

	if paramsLWE, err = rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		LogN:     10,
		Q:        []uint64{0x3fffffffef8001}, // 56
		P:        []uint64{},
		Sigma:    rlwe.DefaultSigma,
		RingType: rlwe.RingStandard,
	}); err != nil {
		panic(err)
	}

	ringQ := params.RingQ()

	values := make([]float64, paramsLWE.N())
	pt := rlwe.NewPlaintext(paramsLWE, paramsLWE.MaxLevel())
	pt.Value.IsNTT = false
	var scale float64 = 1 << 48
	for i := 0; i < paramsLWE.N(); i++ {
		values[i] = utils.RandFloat64(-16, 16)
		if values[i] < 0 {
			pt.Value.Coeffs[0][i] = paramsLWE.Q()[0] - uint64(values[i]*-scale)
		} else {
			pt.Value.Coeffs[0][i] = uint64(values[i] * scale)
		}
	}

	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKey()

	LUTScale := float64(int(1<<30)) * scale * float64(2*params.N()) / float64(paramsLWE.Q()[0])

	kgenLWE := rlwe.NewKeyGenerator(paramsLWE)
	skLWE := kgenLWE.GenSecretKey()

	encryptorLWE := rlwe.NewEncryptor(paramsLWE, skLWE)
	ciphertextLWE := rlwe.NewCiphertextNTT(paramsLWE, 1, pt.Level())
	encryptorLWE.Encrypt(pt, ciphertextLWE)

	handler := lwe.NewHandler(params, paramsLWE, nil)

	fmt.Printf("Encrypting bits of skLWE in RGSW...\n ")
	now := time.Now()
	LUTKEY := handler.GenLUTKey(sk, skLWE)
	fmt.Printf("Done (%s)\n", time.Since(now))

	fmt.Printf("Generating LUT... ")
	now = time.Now()
	LUTPoly := InitLUT(LUTScale, ringQ)
	fmt.Printf("Done (%s)\n", time.Since(now))

	fmt.Printf("Evaluating LUT... ")
	now = time.Now()
	ciphertexts := handler.ExtractAndEvaluateLUT(ciphertextLWE, LUTPoly, LUTKEY)
	fmt.Printf("Done (%s)\n", time.Since(now))

	for i := range ciphertexts[:] {
		fmt.Printf("%8.4f -> ", values[i])
		DecryptAndCenter(ciphertexts[i].Value[0], ciphertexts[i].Value[1], sk.Value.Q, ringQ, false, LUTScale)
	}
}

func RescaleLWE(lwe *lwe.Ciphertext, twoN, Q uint64) {
	for i := range lwe.Value[0] {
		lwe.Value[0][i] = CenterMulAndDivRound(lwe.Value[0][i], twoN, Q)
	}
}

func CenterMulAndDivRound(a, twoN, Q uint64) uint64 {

	aBig := ring.NewUint(a)
	QBig := ring.NewUint(Q)
	aBig.Mul(aBig, ring.NewUint(twoN))
	ring.DivRound(aBig, QBig, aBig)
	a = aBig.Uint64()
	return a & (twoN - 1)
}

func PrintPoly(pol *ring.Poly, scale float64, Q uint64) {
	fmt.Printf("[")
	for _, c := range pol.Coeffs[0][:1] {
		if c > Q>>1 {
			fmt.Printf("%8.4f, ", float64(int(c)-int(Q))/scale)
		} else {
			fmt.Printf("%8.4f, ", float64(int(c))/scale)
		}
	}
	fmt.Printf("]\n")
}

func DecryptAndCenter(b, a, sk *ring.Poly, ringQ *ring.Ring, mForm bool, scale float64) {
	pt := ringQ.NewPoly()
	ringQ.MulCoeffsMontgomery(a, sk, pt)
	ringQ.Add(pt, b, pt)
	ringQ.InvNTT(pt, pt)
	if mForm {
		ringQ.InvMForm(pt, pt)
	}

	Q := ringQ.Modulus[0]
	fmt.Printf("[")
	for _, c := range pt.Coeffs[0][:1] {
		if c > Q>>1 {
			fmt.Printf("%8.4f, ", (float64(c)-float64(Q))/scale)
		} else {
			fmt.Printf("%8.4f, ", float64(c)/scale)
		}
	}
	fmt.Printf("]\n")
}

func AddRGSW(rgsw *rlwe.RGSWCiphertext, ringQP *rlwe.RingQP, res *rlwe.RGSWCiphertext) {

	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	for i := range rgsw.Value {
		ringQ.Add(res.Value[i][0][0].Q, rgsw.Value[i][0][0].Q, res.Value[i][0][0].Q)
		ringP.Add(res.Value[i][0][0].P, rgsw.Value[i][0][0].P, res.Value[i][0][0].P)

		ringQ.Add(res.Value[i][0][1].Q, rgsw.Value[i][0][1].Q, res.Value[i][0][1].Q)
		ringP.Add(res.Value[i][0][1].P, rgsw.Value[i][0][1].P, res.Value[i][0][1].P)

		ringQ.Add(res.Value[i][1][0].Q, rgsw.Value[i][1][0].Q, res.Value[i][1][0].Q)
		ringP.Add(res.Value[i][1][0].P, rgsw.Value[i][1][0].P, res.Value[i][1][0].P)

		ringQ.Add(res.Value[i][1][1].Q, rgsw.Value[i][1][1].Q, res.Value[i][1][1].Q)
		ringP.Add(res.Value[i][1][1].P, rgsw.Value[i][1][1].P, res.Value[i][1][1].P)
	}
}

func MulRGSWByXPowAlphaMinusOne(rgsw *rlwe.RGSWCiphertext, alpha uint64, powXMinusOne []rlwe.PolyQP, ringQP *rlwe.RingQP, res *rlwe.RGSWCiphertext) {

	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	for i := range rgsw.Value {
		ringQ.MulCoeffsMontgomery(rgsw.Value[i][0][0].Q, powXMinusOne[alpha].Q, res.Value[i][0][0].Q)
		ringP.MulCoeffsMontgomery(rgsw.Value[i][0][0].P, powXMinusOne[alpha].P, res.Value[i][0][0].P)

		ringQ.MulCoeffsMontgomery(rgsw.Value[i][0][1].Q, powXMinusOne[alpha].Q, res.Value[i][0][1].Q)
		ringP.MulCoeffsMontgomery(rgsw.Value[i][0][1].P, powXMinusOne[alpha].P, res.Value[i][0][1].P)

		ringQ.MulCoeffsMontgomery(rgsw.Value[i][1][0].Q, powXMinusOne[alpha].Q, res.Value[i][1][0].Q)
		ringP.MulCoeffsMontgomery(rgsw.Value[i][1][0].P, powXMinusOne[alpha].P, res.Value[i][1][0].P)

		ringQ.MulCoeffsMontgomery(rgsw.Value[i][1][1].Q, powXMinusOne[alpha].Q, res.Value[i][1][1].Q)
		ringP.MulCoeffsMontgomery(rgsw.Value[i][1][1].P, powXMinusOne[alpha].P, res.Value[i][1][1].P)
	}
}

func MulRGSWByXPowAlphaMinusOneAndAdd(rgsw *rlwe.RGSWCiphertext, alpha uint64, powXMinusOne []rlwe.PolyQP, ringQP *rlwe.RingQP, res *rlwe.RGSWCiphertext) {

	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	for i := range rgsw.Value {
		ringQ.MulCoeffsMontgomeryAndAdd(rgsw.Value[i][0][0].Q, powXMinusOne[alpha].Q, res.Value[i][0][0].Q)
		ringP.MulCoeffsMontgomeryAndAdd(rgsw.Value[i][0][0].P, powXMinusOne[alpha].P, res.Value[i][0][0].P)

		ringQ.MulCoeffsMontgomeryAndAdd(rgsw.Value[i][0][1].Q, powXMinusOne[alpha].Q, res.Value[i][0][1].Q)
		ringP.MulCoeffsMontgomeryAndAdd(rgsw.Value[i][0][1].P, powXMinusOne[alpha].P, res.Value[i][0][1].P)

		ringQ.MulCoeffsMontgomeryAndAdd(rgsw.Value[i][1][0].Q, powXMinusOne[alpha].Q, res.Value[i][1][0].Q)
		ringP.MulCoeffsMontgomeryAndAdd(rgsw.Value[i][1][0].P, powXMinusOne[alpha].P, res.Value[i][1][0].P)

		ringQ.MulCoeffsMontgomeryAndAdd(rgsw.Value[i][1][1].Q, powXMinusOne[alpha].Q, res.Value[i][1][1].Q)
		ringP.MulCoeffsMontgomeryAndAdd(rgsw.Value[i][1][1].P, powXMinusOne[alpha].P, res.Value[i][1][1].P)
	}
}
