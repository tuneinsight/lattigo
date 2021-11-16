package main

import (
	"fmt"
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/lwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
	"math"
	"time"
)

func main() {
	LUT()
}

func sign(x float64) (y float64) {
	if x > 0 {
		return 1
	} else if x < 0 {
		return -1
	} else {
		return 0
	}
}

func sigmoid(x float64) (y float64) {
	return 1.0 / (math.Exp(-x) + 1)
}

func identity(x float64) (y float64) {
	return x
}

func relu(x float64) (y float64) {
	if x < 0 {
		return 0
	} else {
		return x
	}
}

func LUT() {
	var params rlwe.Parameters
	var err error
	if params, err = rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		LogN:     11,
		Q:        []uint64{0x4000000120001},  // 27 bits
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
		Scale:    1 << 20,
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

	LUTScale := float64(1 << 40)

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
		valuesMinusAMinusB[i] = (-a - b) / (b - a)
	}

	pt := ckks.NewPlaintext(paramsLWE, paramsLWE.MaxLevel(), paramsLWE.Scale())
	encoder.EncodeCoeffsNTT(values, pt)

	ptMinusAminusB := ckks.NewPlaintext(paramsLWE, paramsLWE.MaxLevel(), (paramsLWE.QiFloat64(0)/4.0)*paramsLWE.QiFloat64(1))
	encoder.EncodeCoeffsNTT(valuesMinusAMinusB, ptMinusAminusB)

	ciphertextLWE := encryptorLWE.EncryptNew(pt)

	fmt.Println(ciphertextLWE.Value[0].IsNTT)

	fmt.Printf("Generating LUT... ")
	now = time.Now()
	LUTPoly := lwe.InitLUT(relu, LUTScale, ringQ, a, b)
	fmt.Printf("Done (%s)\n", time.Since(now))

	fmt.Printf("Evaluating LUT... ")
	now = time.Now()

	fmt.Println()
	DecryptAndCenter(8, ciphertextLWE.Value[0], ciphertextLWE.Value[1], skLWE.Value.Q, paramsLWE.RingQ(), false, ciphertextLWE.Scale)

	// Change of basis and changes the scale from Delta to Q/4
	diffScale := paramsLWE.QiFloat64(0) / (4.0 * ciphertextLWE.Scale)
	eval.MultByConst(ciphertextLWE, (2.0/(b-a))*diffScale, ciphertextLWE)
	ciphertextLWE.Scale = (paramsLWE.QiFloat64(0) / 4.0) * paramsLWE.QiFloat64(1)
	eval.Add(ciphertextLWE, ptMinusAminusB, ciphertextLWE)
	eval.Rescale(ciphertextLWE, paramsLWE.Scale(), ciphertextLWE)

	DecryptAndCenter(8, ciphertextLWE.Value[0], ciphertextLWE.Value[1], skLWE.Value.Q, paramsLWE.RingQ(), false, ciphertextLWE.Scale)

	lutPolyMap := make(map[int]*ring.Poly)

	for i := 0; i < 32; i++ {
		lutPolyMap[i] = LUTPoly
	}

	lutPolyMap[64] = LUTPoly

	ciphertexts := handler.ExtractAndEvaluateLUT(ciphertextLWE.Ciphertext, lutPolyMap, LUTKEY)
	fmt.Printf("Done (%s)\n", time.Since(now))

	for i := range ciphertexts {
		fmt.Printf("Slot %4d : %8.4f -> ", i, values[i])
		DecryptAndCenter(1, ciphertexts[i].Value[0], ciphertexts[i].Value[1], sk.Value.Q, ringQ, false, LUTScale)
	}
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
