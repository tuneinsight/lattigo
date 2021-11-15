package main

import (
	"fmt"
	"github.com/ldsec/lattigo/v2/lwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/utils"
	"time"
	"math"
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

// InitLUT takes a function g, and creates an LUT table for the function between the intervals a, b.
// When evaluating the LUT, the input x is assumed to be normalized to (2*x - a - b)/(b-a) (rescaled to be between -1 a 1)
func InitLUT(g func(x float64) (y float64), scale float64, ringQ *ring.Ring, a, b float64) (F *ring.Poly) {
	F = ringQ.NewPoly()
	Q := ringQ.Modulus[0]

	// Zero
	interval := 2.0 / float64(ringQ.N)

	// Negatives
	for i := 0; i < (ringQ.N>>1)+1; i++ {
		F.Coeffs[0][i] = ScaleUpExact(g(normalizeInv(-interval * float64(i), a, b)), scale, Q)
	}

	// Positives
	for i := (ringQ.N>>1)+1; i < ringQ.N; i++ {
		F.Coeffs[0][i] = ScaleUpExact(-g(normalizeInv(interval * float64(ringQ.N - i), a, b)), scale, Q)
	}

	ringQ.NTT(F, F)

	return
}

func normalize(x, a, b float64) (y float64){
	return (2* x - a - b)/(b-a)
}

func normalizeInv(x, a, b float64) (y float64){
	return (x * (b-a) + b + a)/2.0
}

func sign(x float64) (y float64){
	if x > 0{
		return 1
	}else if x < 0{
		return -1
	}else{
		return 0
	}
}

func sigmoid(x float64) (y float64){
	return 1.0/(math.Exp(-x)+1)
}

func identity(x float64) (y float64){
	return  x
}

func relu(x float64) (y float64){
	if x < 0{
		return 0
	}else{
		return x
	}
}

func LUT() {
	var params rlwe.Parameters
	var err error
	if params, err = rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		LogN:     11,
		Q:        []uint64{0x4000000120001}, // 27 bits
		P:        []uint64{0x80000000440001}, // 27 bits
		Sigma:    rlwe.DefaultSigma,
		RingType: rlwe.RingStandard,
	}); err != nil {
		panic(err)
	}

	var paramsLWE ckks.Parameters
	if paramsLWE, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:     10,
		LogSlots: 9,
		Q:        []uint64{0x8007001, 0x8008001}, // 27 bits
		P:        []uint64{},
		Scale:1<<20,
		Sigma:    rlwe.DefaultSigma,
		RingType: rlwe.RingStandard,
	}); err != nil {
		panic(err)
	}

	encoder := ckks.NewEncoder(paramsLWE)

	ringQ := params.RingQ()

	
	a, b := -4.0, 4.0
	
	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKey()

	LUTScale := float64(1<<40)

	fmt.Println(LUTScale)

	kgenLWE := ckks.NewKeyGenerator(paramsLWE)
	skLWE := kgenLWE.GenSecretKey()
	eval := ckks.NewEvaluator(paramsLWE, rlwe.EvaluationKey{nil, nil})

	encryptorLWE := ckks.NewEncryptor(paramsLWE, skLWE)

	handler := lwe.NewHandler(params, paramsLWE.Parameters, nil)

	fmt.Printf("Encrypting bits of skLWE in RGSW...\n ")
	now := time.Now()
	LUTKEY := handler.GenLUTKey(sk, skLWE)
	fmt.Printf("Done (%s)\n", time.Since(now))

	values := make([]float64, paramsLWE.N())
	valuesMinusAMinusB := make([]float64, paramsLWE.N())
	for i := 0; i < paramsLWE.N(); i++ {
		values[i] = utils.RandFloat64(a, b)
		valuesMinusAMinusB[i] = (-a-b)/(b-a)
	}

	pt := ckks.NewPlaintext(paramsLWE, paramsLWE.MaxLevel(), paramsLWE.Scale())
	encoder.EncodeCoeffsNTT(values, pt)

	ptMinusAminusB := ckks.NewPlaintext(paramsLWE, paramsLWE.MaxLevel(), (paramsLWE.QiFloat64(0) / 4.0) * paramsLWE.QiFloat64(1))
	encoder.EncodeCoeffsNTT(valuesMinusAMinusB, ptMinusAminusB)

	ciphertextLWE := encryptorLWE.EncryptNew(pt)

	fmt.Println(ciphertextLWE.Value[0].IsNTT)
	

	fmt.Printf("Generating LUT... ")
	now = time.Now()
	LUTPoly := InitLUT(relu, LUTScale, ringQ, a, b)
	fmt.Printf("Done (%s)\n", time.Since(now))

	fmt.Printf("Evaluating LUT... ")
	now = time.Now()

	fmt.Println()
	DecryptAndCenter(8, ciphertextLWE.Value[0], ciphertextLWE.Value[1], skLWE.Value.Q, paramsLWE.RingQ(), false, ciphertextLWE.Scale)

	// Change of basis and changes the scale from Delta to Q/4
	diffScale := paramsLWE.QiFloat64(0) / (4.0 * ciphertextLWE.Scale)
	eval.MultByConst(ciphertextLWE, (2.0/(b-a)) * diffScale , ciphertextLWE)
	ciphertextLWE.Scale =  (paramsLWE.QiFloat64(0) / 4.0) * paramsLWE.QiFloat64(1)
	eval.Add(ciphertextLWE, ptMinusAminusB, ciphertextLWE)
	eval.Rescale(ciphertextLWE, paramsLWE.Scale(), ciphertextLWE)
	
	DecryptAndCenter(8, ciphertextLWE.Value[0], ciphertextLWE.Value[1], skLWE.Value.Q, paramsLWE.RingQ(), false, ciphertextLWE.Scale)

	ciphertexts := handler.ExtractAndEvaluateLUT(ciphertextLWE.Ciphertext, LUTPoly, LUTKEY)
	fmt.Printf("Done (%s)\n", time.Since(now))

	for i := range ciphertexts[:16] {
		fmt.Printf("%8.4f -> ", values[i])
		DecryptAndCenter(1, ciphertexts[i].Value[0], ciphertexts[i].Value[1], sk.Value.Q, ringQ, false, LUTScale)
	}
}

func AddConst(p *ring.Poly, scalar, scale float64, ringQ *ring.Ring){

	fmt.Println(scale)

	ringQ.InvNTTLvl(p.Level(), p, p)
	fmt.Println(scalar)
	for i := 0; i < p.Level()+1; i++{
		qi := ringQ.Modulus[i]
		ring.AddScalarVec(p.Coeffs[i], p.Coeffs[i], ScaleUpExact(scalar, scale, qi), qi)
	}
	ringQ.NTTLvl(p.Level(), p, p)
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

func DecryptAndCenter(n int, b, a, sk *ring.Poly, ringQ *ring.Ring, mForm bool, scale float64) {

	pt := ringQ.NewPolyLvl(0)
	ringQ.MulCoeffsMontgomeryLvl(0, a, sk, pt)
	ringQ.AddLvl(0, pt, b, pt)
	ringQ.InvNTTLvl(0, pt, pt)
	if mForm {
		ringQ.InvMFormLvl(0, pt, pt)
	}

	Q := ringQ.Modulus[0]
	fmt.Printf("[")
	for _, c := range pt.Coeffs[0][:n] {
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
