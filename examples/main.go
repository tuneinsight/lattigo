package main

import (
	"fmt"
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
	"runtime"
	"time"
	//"math"
	"math/big"
)

// LogN of the ring degree of the used parameters
var LogN = 15

func main() {
	LUT()
	//SimulatedComplex()
	//Complex()
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
		F.Coeffs[0][i] = ScaleUpExact(-float64(-i), scale, Q)
	}

	// Positives
	for i := (ringQ.N>>1)+1; i < ringQ.N; i++ {
		F.Coeffs[0][i] = ScaleUpExact(-float64(ringQ.N-i), scale, Q)
	}

	ringQ.NTT(F, F)

	return
}

func LUT() {
	var params rlwe.Parameters
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

	ringQ := params.RingQ()

	ks := rlwe.NewKeySwitcher(params)
	_ = ks

	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKey()
	encryptor := rlwe.NewEncryptor(params, sk)
	//decryptor := rlwe.NewDecryptor(params, sk)

	fmt.Printf("Encrypting bits of skLWE in RGSW... ")
	now := time.Now()

	prng, _ := utils.NewPRNG()

	LWEDegree := 1 << 10
	LUTScale := float64(int(1 << 30))
	mask := uint64(2*params.N() - 1)

	skLWE := make([]uint64, LWEDegree)
	for i := range skLWE {

		si := ring.RandInt32(prng, 3)

		for si == 3 {
			si = ring.RandInt32(prng, 3)
		}

		skLWE[i] = (si - 1) & mask
	}

	plaintextRGSWOne := rlwe.NewPlaintext(params, params.MaxLevel())
	plaintextRGSWOne.Value.IsNTT = true
	plaintextRGSWZer := rlwe.NewPlaintext(params, params.MaxLevel())
	plaintextRGSWZer.Value.IsNTT = true

	for i := 0; i < params.N(); i++ {
		plaintextRGSWOne.Value.Coeffs[0][i] = 1
	}

	EncOneRGSW := rlwe.NewCiphertextRGSWNTT(params, params.MaxLevel())

	encryptor.EncryptRGSW(plaintextRGSWOne, EncOneRGSW)

	skRGSWPos := make([]*rlwe.RGSWCiphertext, params.N())
	skRGSWNeg := make([]*rlwe.RGSWCiphertext, params.N())

	for i := 0; i < len(skLWE); i++ {

		skRGSWPos[i] = rlwe.NewCiphertextRGSWNTT(params, params.MaxLevel())
		skRGSWNeg[i] = rlwe.NewCiphertextRGSWNTT(params, params.MaxLevel())

		si := skLWE[i]

		if si == 1 {
			encryptor.EncryptRGSW(plaintextRGSWOne, skRGSWPos[i])
			encryptor.EncryptRGSW(plaintextRGSWZer, skRGSWNeg[i])
		} else if si == mask {
			encryptor.EncryptRGSW(plaintextRGSWZer, skRGSWPos[i])
			encryptor.EncryptRGSW(plaintextRGSWOne, skRGSWNeg[i])
		} else {
			encryptor.EncryptRGSW(plaintextRGSWZer, skRGSWPos[i])
			encryptor.EncryptRGSW(plaintextRGSWZer, skRGSWNeg[i])
		}
	}

	fmt.Printf("Done (%s)\n", time.Since(now))

	fmt.Printf("Generating LUT... ")
	now = time.Now()
	LUTPoly := InitLUT(LUTScale, ringQ)
	fmt.Printf("Done (%s)\n", time.Since(now))

	fmt.Printf("Generating Table for NTT(X^alpha - 1)... ")
	now = time.Now()
	powXMinusOne := make([]rlwe.PolyQP, 2*params.N())
	powX := make([]*ring.Poly, 2*params.N())
	Q := params.Q()[0]
	P := params.P()[0]
	ringQP := params.RingQP()
	ringP := params.RingP()
	// Stores NTT(X^i - 1)
	for i := 0; i < params.N(); i++ {
		powX[i] = ringQ.NewPoly()
		powXMinusOne[i] = ringQP.NewPoly()

		powXMinusOne[i].Q.Coeffs[0][i] = 1
		powXMinusOne[i].P.Coeffs[0][i] = 1

		ringQ.NTT(powXMinusOne[i].Q, powX[i])
		ringQ.MForm(powX[i], powX[i])

		powXMinusOne[i].Q.Coeffs[0][0] += Q - 1
		powXMinusOne[i].P.Coeffs[0][0] += P - 1

		ringQ.NTT(powXMinusOne[i].Q, powXMinusOne[i].Q)
		ringQ.MForm(powXMinusOne[i].Q, powXMinusOne[i].Q)
		ringP.NTT(powXMinusOne[i].P, powXMinusOne[i].P)
		ringP.MForm(powXMinusOne[i].P, powXMinusOne[i].P)
	}

	for i, j := params.N(), 0; i < 2*params.N(); i, j = i+1, j+1 {
		powX[i] = ringQ.NewPoly()
		powXMinusOne[i] = ringQP.NewPoly()

		powXMinusOne[i].Q.Coeffs[0][j] = Q - 1
		powXMinusOne[i].P.Coeffs[0][j] = P - 1

		ringQ.NTT(powXMinusOne[i].Q, powX[i])
		ringQ.MForm(powX[i], powX[i])

		powXMinusOne[i].Q.Coeffs[0][0] += Q - 1
		powXMinusOne[i].P.Coeffs[0][0] += P - 1

		ringQ.NTT(powXMinusOne[i].Q, powXMinusOne[i].Q)
		ringQ.MForm(powXMinusOne[i].Q, powXMinusOne[i].Q)
		ringP.NTT(powXMinusOne[i].P, powXMinusOne[i].P)
		ringP.MForm(powXMinusOne[i].P, powXMinusOne[i].P)
	}
	fmt.Printf("Done (%s)\n", time.Since(now))

	// Generates some LWE

	m := 1024.0

	a := make([]uint64, LWEDegree)
	for i := range a {
		a[i] = ring.RandInt32(prng, mask)
	}

	var b uint64
	for i := 0; i < LWEDegree; i++ {
		b -= a[i] * skLWE[i]
	}

	b += ScaleUpExact(m, 1.0, uint64(2*LWEDegree))
	b &= mask

	fmt.Printf("Evaluating LUT... ")

	acc := rlwe.NewCiphertextNTT(params, 1, params.MaxLevel())
	ringQ.MulCoeffsMontgomery(LUTPoly, powXMinusOne[b].Q, acc.Value[0])
	ringQ.Add(acc.Value[0], LUTPoly, acc.Value[0])

	tmpRGSW := rlwe.NewCiphertextRGSWNTT(params, params.MaxLevel())

	now = time.Now()
	for i := 0; i < LWEDegree; i++ {
		MulRGSWByXPowAlphaMinusOne(skRGSWPos[i], a[i], powXMinusOne, ringQP, tmpRGSW)
		MulRGSWByXPowAlphaMinusOneAndAdd(skRGSWNeg[i], -a[i]&mask, powXMinusOne, ringQP, tmpRGSW)
		AddRGSW(EncOneRGSW, ringQP, tmpRGSW)
		ks.MulRGSW(acc, tmpRGSW, acc)
	}
	fmt.Printf("Done (%s)\n", time.Since(now))

	tmp := ringQ.NewPoly()
	ringQ.MulCoeffsMontgomery(LUTPoly, powX[b], tmp)

	for i := 0; i < LWEDegree; i++ {
		if skLWE[i] == 1 {
			b += a[i] & mask
			ringQ.MulCoeffsMontgomery(tmp, powX[+a[i]&mask], tmp)
		} else if skLWE[i] == mask {
			b += -a[i] & mask
			ringQ.MulCoeffsMontgomery(tmp, powX[-a[i]&mask], tmp)
		}
	}

	ringQ.InvNTT(tmp, tmp)
	fmt.Printf("Want : ")
	PrintPoly(tmp, LUTScale, Q)

	fmt.Printf("Have : ")
	DecryptAndCenter(acc.Value[0], acc.Value[1], sk.Value.Q, ringQ, false, LUTScale)

}

func PrintPoly(pol *ring.Poly, scale float64, Q uint64) {
	fmt.Printf("[")
	for _, c := range pol.Coeffs[0][:4] {
		if c > Q>>1 {
			fmt.Printf("%0.4f, ", float64(int(c)-int(Q))/scale)
		} else {
			fmt.Printf("%0.4f, ", float64(int(c))/scale)
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
	for _, c := range pt.Coeffs[0][:4] {
		if c > Q>>1 {
			fmt.Printf("%0.4f, ", (float64(c)-float64(Q))/scale)
		} else {
			fmt.Printf("%0.4f, ", float64(c)/scale)
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

// SimulatedComplex implements complex arithmetic using RCKKS
func SimulatedComplex() {
	// Schemes parameters are created from scratch
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:     LogN,
		LogQ:     []int{55, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40},
		LogP:     []int{45, 45, 45},
		Sigma:    rlwe.DefaultSigma,
		LogSlots: LogN,
		Scale:    float64(1 << 40),
		RingType: rlwe.RingConjugateInvariant,
	})
	if err != nil {
		panic(err)
	}

	kgen := ckks.NewKeyGenerator(params)
	sk := kgen.GenSecretKey()
	encryptor := ckks.NewEncryptor(params, sk)
	decryptor := ckks.NewDecryptor(params, sk)
	encoder := ckks.NewEncoder(params)

	rlk := kgen.GenRelinearizationKey(sk, 2)

	valuesReal := make([]complex128, params.Slots())
	valuesImag := make([]complex128, params.Slots())
	for i := range valuesReal {
		valuesReal[i] = complex(1, 0)
		valuesImag[i] = complex(1, 0)
	}

	ptReal := ckks.NewPlaintext(params, params.MaxLevel(), params.Scale())
	ptImag := ckks.NewPlaintext(params, params.MaxLevel(), params.Scale())

	encoder.Encode(ptReal, valuesReal, params.LogSlots())
	encoder.Encode(ptImag, valuesImag, params.LogSlots())

	ctA := Ciphertext{encryptor.EncryptNew(ptReal), encryptor.EncryptNew(ptImag)}
	ctB := Ciphertext{encryptor.EncryptNew(ptReal), encryptor.EncryptNew(ptImag)}

	diagMatrixReal := make(map[int][]complex128)
	diagMatrixImag := make(map[int][]complex128)

	//diagMatrixReal[-15] = make([]complex128, params.Slots())
	//diagMatrixReal[-4] = make([]complex128, params.Slots())
	diagMatrixReal[-1] = make([]complex128, params.Slots())
	diagMatrixReal[0] = make([]complex128, params.Slots())
	diagMatrixReal[1] = make([]complex128, params.Slots())
	//diagMatrixReal[2] = make([]complex128, params.Slots())
	//diagMatrixReal[3] = make([]complex128, params.Slots())
	//diagMatrixReal[4] = make([]complex128, params.Slots())
	//diagMatrixReal[15] = make([]complex128, params.Slots())

	//diagMatrixImag[-15] = make([]complex128, params.Slots())
	//diagMatrixImag[-4] = make([]complex128, params.Slots())
	diagMatrixImag[-1] = make([]complex128, params.Slots())
	diagMatrixImag[0] = make([]complex128, params.Slots())
	diagMatrixImag[1] = make([]complex128, params.Slots())
	//diagMatrixImag[2] = make([]complex128, params.Slots())
	//diagMatrixImag[3] = make([]complex128, params.Slots())
	//diagMatrixImag[4] = make([]complex128, params.Slots())
	//diagMatrixImag[15] = make([]complex128, params.Slots())

	for i := 0; i < params.Slots(); i++ {
		//diagMatrixReal[-15][i] = complex(1, 0)
		//diagMatrixReal[-4][i] = complex(1, 0)
		diagMatrixReal[-1][i] = complex(1, 0)
		diagMatrixReal[0][i] = complex(1, 0)
		diagMatrixReal[1][i] = complex(1, 0)
		//diagMatrixReal[2][i] = complex(1, 0)
		//diagMatrixReal[3][i] = complex(1, 0)
		//diagMatrixReal[4][i] = complex(1, 0)
		//diagMatrixReal[15][i] = complex(1, 0)

		//diagMatrixImag[-15][i] = complex(1, 0)
		//diagMatrixImag[-4][i] = complex(1, 0)
		diagMatrixImag[-1][i] = complex(1, 0)
		diagMatrixImag[0][i] = complex(1, 0)
		diagMatrixImag[1][i] = complex(1, 0)
		//diagMatrixImag[2][i] = complex(1, 0)
		//diagMatrixImag[3][i] = complex(1, 0)
		//diagMatrixImag[4][i] = complex(1, 0)
		//diagMatrixImag[15][i] = complex(1, 0)
	}

	ptDiagReal := encoder.EncodeDiagMatrixBSGSAtLvl(params.MaxLevel(), diagMatrixReal, params.QiFloat64(params.MaxLevel()), 4.0, params.LogSlots())
	ptDiagImag := encoder.EncodeDiagMatrixBSGSAtLvl(params.MaxLevel(), diagMatrixImag, params.QiFloat64(params.MaxLevel()), 4.0, params.LogSlots())

	diagMat := PtDiagMatrix{
		N1:       ptDiagReal.N1,
		Level:    ptDiagReal.Level,
		LogSlots: ptDiagReal.LogSlots,
		Scale:    ptDiagReal.Scale,
		VecReal:  ptDiagReal.Vec,
		VecImag:  ptDiagImag.Vec}

	rots := params.RotationsForDiagMatrixMult(ptDiagReal)

	rotKey := kgen.GenRotationKeysForRotations(rots, false, sk)

	eval := NewEvaluator(params, rlk, rotKey)

	var ctC Ciphertext
	var tot time.Duration
	for i := 0; i < 100; i++ {
		now := time.Now()
		ctC = eval.MulRelinNew(ctA, ctB)
		eval.Rescale(ctC.Real, eval.params.Scale(), ctC.Real)
		eval.Rescale(ctC.Imag, eval.params.Scale(), ctC.Imag)
		tot += time.Since(now)
	}
	fmt.Println("RCKKS - LogN :", params.LogN())
	fmt.Printf("Done : %s\n", tot/100.)

	//ctC := ctA

	PoolDecompRealQP := make([]rlwe.PolyQP, params.Beta())
	PoolDecompImagQP := make([]rlwe.PolyQP, params.Beta())

	for i := range PoolDecompImagQP {
		PoolDecompRealQP[i] = params.RingQP().NewPoly()
		PoolDecompImagQP[i] = params.RingQP().NewPoly()
	}

	ks := eval.GetKeySwitcher()

	ks.DecomposeNTT(ctC.Level(), params.PCount()-1, params.PCount(), ctC.Real.Value[1], PoolDecompRealQP)
	ks.DecomposeNTT(ctC.Level(), params.PCount()-1, params.PCount(), ctC.Imag.Value[1], PoolDecompImagQP)

	eval.MultiplyByDiagMatrixBSGS(ctC, diagMat, PoolDecompRealQP, PoolDecompImagQP, ctC)

	fmt.Println(ctC.Real.Scale)

	vReal := encoder.DecodePublic(decryptor.DecryptNew(ctC.Real), params.LogSlots(), 0)
	vImag := encoder.DecodePublic(decryptor.DecryptNew(ctC.Imag), params.LogSlots(), 0)

	for i := range vReal[:4] {
		fmt.Println(i, real(vReal[i]), real(vImag[i]))
	}
}

// Complex implement complex arithmetic using CKKS
func Complex() {
	// Schemes parameters are created from scratch
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:     LogN + 1,
		LogQ:     []int{55, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40},
		LogP:     []int{45, 45, 45},
		Sigma:    rlwe.DefaultSigma,
		LogSlots: 15,
		Scale:    float64(1 << 40),
		RingType: rlwe.RingStandard,
	})
	if err != nil {
		panic(err)
	}

	kgen := ckks.NewKeyGenerator(params)
	sk := kgen.GenSecretKey()
	encryptor := ckks.NewEncryptor(params, sk)
	decryptor := ckks.NewDecryptor(params, sk)
	encoder := ckks.NewEncoder(params)

	rlk := kgen.GenRelinearizationKey(sk, 2)
	eval := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: nil})

	values := make([]complex128, params.Slots())
	for i := range values {
		values[i] = complex(0.9238795325112867, 0.3826834323650898)
	}

	plaintext := ckks.NewPlaintext(params, params.MaxLevel(), params.Scale())
	encoder.Encode(plaintext, values, params.LogSlots())
	ct := encryptor.EncryptNew(plaintext)

	var A *ckks.Ciphertext
	var tot time.Duration
	for i := 0; i < 100; i++ {
		now := time.Now()
		A = eval.MulRelinNew(ct, ct)
		eval.Rescale(A, params.Scale(), A)
		tot += time.Since(now)
	}
	fmt.Println("CKKS - LogN :", params.LogN())
	fmt.Printf("Done : %s\n", tot/100.)

	v := encoder.DecodePublic(decryptor.DecryptNew(A), params.LogSlots(), 0)

	for i := range v[:4] {
		fmt.Println(i, v[i])
	}
}

// Evaluator is a custom evaluator for this example
type Evaluator struct {
	ckks.Evaluator
	params          ckks.Parameters
	rlk             *rlwe.RelinearizationKey
	rtks            *rlwe.RotationKeySet
	permuteNTTIndex map[uint64][]uint64
	ctxpool         Ciphertext
}

func (eval *Evaluator) permuteNTTIndexesForKey(rtks *rlwe.RotationKeySet) *map[uint64][]uint64 {
	if rtks == nil {
		return &map[uint64][]uint64{}
	}
	permuteNTTIndex := make(map[uint64][]uint64, len(rtks.Keys))
	for galEl := range rtks.Keys {
		permuteNTTIndex[galEl] = eval.params.RingQ().PermuteNTTIndex(galEl)
	}
	return &permuteNTTIndex
}

// NewEvaluator creates a new custom Evaluator
func NewEvaluator(params ckks.Parameters, rlk *rlwe.RelinearizationKey, rtks *rlwe.RotationKeySet) (eval Evaluator) {

	eval = Evaluator{}
	eval.Evaluator = ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: rtks})
	eval.params = params
	eval.ctxpool = NewCiphertext(params, 2, params.MaxLevel(), 0)
	eval.rlk = rlk
	eval.rtks = rtks

	if eval.rtks != nil {
		eval.permuteNTTIndex = *eval.permuteNTTIndexesForKey(rtks)
	}
	return
}

// Ciphertext is a custom ciphertext for complex arithmetic with RCKKS.
type Ciphertext struct {
	Real *ckks.Ciphertext
	Imag *ckks.Ciphertext
}

// NewCiphertext creates a new Ciphertext
func NewCiphertext(params ckks.Parameters, degree, level int, scale float64) Ciphertext {
	return Ciphertext{Real: ckks.NewCiphertext(params, degree, level, scale), Imag: ckks.NewCiphertext(params, degree, level, scale)}
}

// Plaintext is a custom plaintext for complex arithmetic with RCKKS
type Plaintext struct {
	Real *ckks.Plaintext
	Imag *ckks.Plaintext
}

// MulRelinNew mul + relinearize
func (eval *Evaluator) MulRelinNew(A, B Ciphertext) (C Ciphertext) {
	level := utils.MinInt(utils.MinInt(A.Real.Level(), A.Imag.Level()), utils.MinInt(B.Real.Level(), B.Imag.Level()))
	C = Ciphertext{ckks.NewCiphertext(eval.params, 1, level, 0), ckks.NewCiphertext(eval.params, 1, level, 0)}
	eval.mulRelin(A, B, C, true)
	return
}

// MulRelin mul + relinearize
func (eval *Evaluator) MulRelin(A, B, C Ciphertext) {
	eval.mulRelin(A, B, C, true)
	return
}

// MulNew mul without relinearize
func (eval *Evaluator) MulNew(A, B Ciphertext) (C Ciphertext) {
	level := utils.MinInt(utils.MinInt(A.Real.Level(), A.Imag.Level()), utils.MinInt(B.Real.Level(), B.Imag.Level()))
	C = Ciphertext{ckks.NewCiphertext(eval.params, 2, level, 0), ckks.NewCiphertext(eval.params, 2, level, 0)}
	eval.mulRelin(A, B, C, false)
	return
}

// Mul mul without relinearize
func (eval *Evaluator) Mul(A, B, C Ciphertext) {
	eval.mulRelin(A, B, C, false)
	return
}

// Operand interface for the Evaluator elements
type Operand interface {
	El() *rlwe.Ciphertext
	Degree() int
	Level() int
}

// MulPlain multiplied ct with pt
func (eval *Evaluator) MulPlain(A Ciphertext, B Plaintext, C Ciphertext) {

}

func (eval *Evaluator) mulRelin(A, B, C Ciphertext, relin bool) {

	if A.Real.Degree() > 1 || A.Imag.Degree() > 1 || B.Real.Degree() > 1 || B.Imag.Degree() > 1 {
		panic("cannot MulRelin: input elements must be degree 0 or 1")
	}

	level := utils.MinInt(utils.MinInt(A.Real.Level(), A.Imag.Level()), utils.MinInt(B.Real.Level(), B.Imag.Level()))

	ringQ := eval.params.RingQ()

	C.Real.Scale = A.Real.Scale * B.Real.Scale
	C.Imag.Scale = A.Imag.Scale * B.Imag.Scale

	var c00, c01 *ring.Poly

	ks := eval.GetKeySwitcher()
	pool0, pool1, pool2, pool3 := ks.Pool[1].Q, ks.Pool[2].Q, ks.Pool[3].Q, ks.Pool[4].Q

	if A.Real.Degree()+A.Imag.Degree() == 2 && B.Real.Degree()+B.Imag.Degree() == 2 {

		var c2Real, c2Imag *ring.Poly
		if relin == false {
			if C.Real.Degree() < 2 {
				C.Real.El().Resize(eval.params.Parameters, 2)
			}
			if C.Real.Degree() < 2 {
				C.Real.El().Resize(eval.params.Parameters, 2)
			}

			c2Real = C.Real.Value[2]
			c2Imag = C.Imag.Value[2]
		} else {
			c2Real = pool0
			c2Imag = pool1
		}

		c00 = eval.PoolQMul()[0]
		c01 = eval.PoolQMul()[1]

		// (a + b) * (c + d)
		if A.Real != B.Real && A.Imag != B.Imag {

			// Mont(a)
			ringQ.MFormConstantLvl(level, A.Real.Value[0], c00)
			ringQ.MFormConstantLvl(level, A.Real.Value[1], c01)

			// ac
			ringQ.MulCoeffsMontgomeryConstantLvl(level, c00, B.Real.Value[0], C.Real.Value[0])    // r0 : 2q
			ringQ.MulCoeffsMontgomeryConstantLvl(level, c00, B.Real.Value[1], C.Real.Value[1])    // r1 : 2q
			ringQ.MulCoeffsMontgomeryAndAddNoModLvl(level, c01, B.Real.Value[0], C.Real.Value[1]) // r1 : 3q
			ringQ.MulCoeffsMontgomeryConstantLvl(level, c01, B.Real.Value[1], c2Real)             // r2 : 2q

			// ad
			ringQ.MulCoeffsMontgomeryConstantLvl(level, c00, B.Imag.Value[0], C.Imag.Value[0])    // i0 : 2q
			ringQ.MulCoeffsMontgomeryConstantLvl(level, c00, B.Imag.Value[1], C.Imag.Value[1])    // i1 : 2q
			ringQ.MulCoeffsMontgomeryAndAddNoModLvl(level, c01, B.Imag.Value[0], C.Imag.Value[1]) // i1 : 3q
			ringQ.MulCoeffsMontgomeryConstantLvl(level, c01, B.Imag.Value[1], c2Imag)             // i2 : 2q

			// Mont(b)
			ringQ.MFormConstantLvl(level, A.Imag.Value[0], c00)
			ringQ.MFormConstantLvl(level, A.Imag.Value[1], c01)

			// ac - bd
			ringQ.MulCoeffsMontgomeryAndSubNoModLvl(level, c00, B.Imag.Value[0], C.Real.Value[0]) // r0 : 3q
			ringQ.MulCoeffsMontgomeryAndSubNoModLvl(level, c00, B.Imag.Value[1], C.Real.Value[1]) // r1 : 4q
			ringQ.MulCoeffsMontgomeryAndSubNoModLvl(level, c01, B.Imag.Value[0], C.Real.Value[1]) // r1 : 5q
			ringQ.MulCoeffsMontgomeryAndSubNoModLvl(level, c01, B.Imag.Value[1], c2Real)          // r2 : 3q

			// ad + bc
			ringQ.MulCoeffsMontgomeryAndAddNoModLvl(level, c00, B.Real.Value[0], C.Imag.Value[0]) // i0 : 3q
			ringQ.MulCoeffsMontgomeryAndAddNoModLvl(level, c00, B.Real.Value[1], C.Imag.Value[1]) // i1 : 4q
			ringQ.MulCoeffsMontgomeryAndAddNoModLvl(level, c01, B.Real.Value[0], C.Imag.Value[1]) // i1 : 5q
			ringQ.MulCoeffsMontgomeryAndAddNoModLvl(level, c01, B.Real.Value[1], c2Imag)          // i2 : 3q
		}

		// (a + b) * (a + c)
		if A.Real == B.Real && A.Imag != B.Imag {
			// Mont(a)
			ringQ.MFormConstantLvl(level, A.Real.Value[0], c00)
			ringQ.MFormConstantLvl(level, A.Real.Value[1], c01)

			// a * a
			ringQ.MulCoeffsMontgomeryConstantLvl(level, c00, B.Real.Value[0], C.Real.Value[0])
			ringQ.MulCoeffsMontgomeryConstantLvl(level, c00, B.Real.Value[1], C.Real.Value[1])
			ringQ.MulCoeffsMontgomeryAndAddNoModLvl(level, c01, B.Real.Value[0], C.Real.Value[1])
			ringQ.MulCoeffsMontgomeryConstantLvl(level, c01, B.Real.Value[1], c2Real)

			// b + c
			ringQ.AddNoModLvl(level, A.Imag.Value[0], B.Imag.Value[0], pool2)
			ringQ.AddNoModLvl(level, A.Imag.Value[1], B.Imag.Value[1], pool3)

			// a * (b + c)
			ringQ.MulCoeffsMontgomeryConstantLvl(level, c00, pool2, C.Imag.Value[0])
			ringQ.MulCoeffsMontgomeryConstantLvl(level, c00, pool3, C.Imag.Value[1])
			ringQ.MulCoeffsMontgomeryAndAddNoModLvl(level, c01, pool2, C.Imag.Value[1])
			ringQ.MulCoeffsMontgomeryConstantLvl(level, c01, pool3, c2Imag)

			// Mont(b)
			ringQ.MFormConstantLvl(level, A.Imag.Value[0], c00)
			ringQ.MFormConstantLvl(level, A.Imag.Value[1], c01)

			// a * a - b * c
			ringQ.MulCoeffsMontgomeryAndSubNoModLvl(level, c00, B.Imag.Value[0], C.Real.Value[0]) // r0 : 3q
			ringQ.MulCoeffsMontgomeryAndSubNoModLvl(level, c00, B.Imag.Value[1], C.Real.Value[1]) // r1 : 4q
			ringQ.MulCoeffsMontgomeryAndSubNoModLvl(level, c01, B.Imag.Value[0], C.Real.Value[1]) // r1 : 5q
			ringQ.MulCoeffsMontgomeryAndSubNoModLvl(level, c01, B.Imag.Value[1], c2Real)          // r2 : 3q
		}

		// (a + b) * (c + b)
		if A.Real != B.Real && A.Imag == B.Imag {

			// Mont(b)
			ringQ.MFormConstantLvl(level, A.Imag.Value[0], c00)
			ringQ.MFormConstantLvl(level, A.Imag.Value[1], c01)

			// b * b
			ringQ.MulCoeffsMontgomeryConstantAndNegLvl(level, c00, B.Imag.Value[0], C.Real.Value[0]) // r0 : 3q
			ringQ.MulCoeffsMontgomeryConstantAndNegLvl(level, c00, B.Imag.Value[1], C.Real.Value[1]) // r1 : 4q
			ringQ.MulCoeffsMontgomeryAndSubNoModLvl(level, c01, B.Imag.Value[0], C.Real.Value[1])    // r1 : 5q
			ringQ.MulCoeffsMontgomeryConstantAndNegLvl(level, c01, B.Imag.Value[1], c2Real)          // r2 : 3q

			// a + c
			ringQ.AddNoModLvl(level, A.Real.Value[0], B.Real.Value[0], pool2)
			ringQ.AddNoModLvl(level, A.Real.Value[1], B.Real.Value[1], pool3)

			// b * (a + c)
			ringQ.MulCoeffsMontgomeryConstantLvl(level, c00, pool2, C.Imag.Value[0])
			ringQ.MulCoeffsMontgomeryConstantLvl(level, c00, pool3, C.Imag.Value[1])
			ringQ.MulCoeffsMontgomeryAndAddNoModLvl(level, c01, pool2, C.Imag.Value[1])
			ringQ.MulCoeffsMontgomeryConstantLvl(level, c01, pool3, c2Imag)

			// Mont(a)
			ringQ.MFormConstantLvl(level, A.Real.Value[0], c00)
			ringQ.MFormConstantLvl(level, A.Real.Value[1], c01)

			// a * c - b * b
			ringQ.MulCoeffsMontgomeryAndAddNoModLvl(level, c00, B.Real.Value[0], C.Real.Value[0]) // r0 : 3q
			ringQ.MulCoeffsMontgomeryAndAddNoModLvl(level, c00, B.Real.Value[1], C.Real.Value[1]) // r1 : 4q
			ringQ.MulCoeffsMontgomeryAndAddNoModLvl(level, c01, B.Real.Value[0], C.Real.Value[1]) // r1 : 5q
			ringQ.MulCoeffsMontgomeryAndAddNoModLvl(level, c01, B.Real.Value[1], c2Real)          // r2 : 3q
		}

		// (a + b) * (a + b)
		if A.Real == B.Real && A.Imag == B.Imag {

			// Mont(a)
			ringQ.MFormConstantLvl(level, A.Real.Value[0], c00)
			ringQ.MFormConstantLvl(level, A.Real.Value[1], c01)

			// aa
			ringQ.MulCoeffsMontgomeryConstantLvl(level, c00, B.Real.Value[0], C.Real.Value[0])    // r0 : 2q
			ringQ.MulCoeffsMontgomeryConstantLvl(level, c00, B.Real.Value[1], C.Real.Value[1])    // r1 : 2q
			ringQ.MulCoeffsMontgomeryAndAddNoModLvl(level, c01, B.Real.Value[0], C.Real.Value[1]) // r1 : 3q
			ringQ.MulCoeffsMontgomeryConstantLvl(level, c01, B.Real.Value[1], c2Real)             // r2 : 2q

			// ab
			ringQ.MulCoeffsMontgomeryConstantLvl(level, c00, B.Imag.Value[0], C.Imag.Value[0])    // i0 : 2q
			ringQ.MulCoeffsMontgomeryConstantLvl(level, c00, B.Imag.Value[1], C.Imag.Value[1])    // i1 : 2q
			ringQ.MulCoeffsMontgomeryAndAddNoModLvl(level, c01, B.Imag.Value[0], C.Imag.Value[1]) // i1 : 3q
			ringQ.MulCoeffsMontgomeryConstantLvl(level, c01, B.Imag.Value[1], c2Imag)             // i2 : 2q

			// 2ab
			ringQ.AddNoModLvl(level, C.Imag.Value[0], C.Imag.Value[0], C.Imag.Value[0])
			ringQ.AddNoModLvl(level, C.Imag.Value[1], C.Imag.Value[1], C.Imag.Value[1])
			ringQ.AddNoModLvl(level, c2Imag, c2Imag, c2Imag)

			// Mont(b)
			ringQ.MFormConstantLvl(level, A.Imag.Value[0], c00)
			ringQ.MFormConstantLvl(level, A.Imag.Value[1], c01)

			// aa - bb
			ringQ.MulCoeffsMontgomeryAndSubNoModLvl(level, c00, B.Imag.Value[0], C.Real.Value[0])
			ringQ.MulCoeffsMontgomeryAndSubNoModLvl(level, c00, B.Imag.Value[1], C.Real.Value[1])
			ringQ.MulCoeffsMontgomeryAndSubNoModLvl(level, c01, B.Imag.Value[0], C.Real.Value[1])
			ringQ.MulCoeffsMontgomeryAndSubNoModLvl(level, c01, B.Imag.Value[1], c2Real)
		}

		ringQ.ReduceConstantLvl(level, C.Real.Value[0], C.Real.Value[0])
		ringQ.ReduceConstantLvl(level, C.Real.Value[1], C.Real.Value[1])
		ringQ.ReduceConstantLvl(level, c2Real, c2Real)
		ringQ.ReduceConstantLvl(level, C.Imag.Value[0], C.Imag.Value[0])
		ringQ.ReduceConstantLvl(level, C.Imag.Value[1], C.Imag.Value[1])
		ringQ.ReduceConstantLvl(level, c2Imag, c2Imag)

		if relin {
			c2Real.IsNTT = true
			ks.SwitchKeysInPlace(level, c2Real, eval.rlk.Keys[0], pool2, pool3)
			ringQ.AddLvl(level, C.Real.Value[0], pool2, C.Real.Value[0])
			ringQ.AddLvl(level, C.Real.Value[1], pool3, C.Real.Value[1])
			c2Imag.IsNTT = true
			ks.SwitchKeysInPlace(level, c2Imag, eval.rlk.Keys[0], pool2, pool3)
			ringQ.AddLvl(level, C.Imag.Value[0], pool2, C.Imag.Value[0])
			ringQ.AddLvl(level, C.Imag.Value[1], pool3, C.Imag.Value[1])
		}
	}
}

// Level returns the level of a custom ciphertext for RCKKS
func (ct *Ciphertext) Level() int {
	return utils.MinInt(ct.Real.Level(), ct.Imag.Level())
}

// PtDiagMatrix custom plaintext matrix for RCKKS with complex arithmetic
type PtDiagMatrix struct {
	N1       int
	Level    int
	LogSlots int
	Scale    float64
	VecReal  map[int]rlwe.PolyQP
	VecImag  map[int]rlwe.PolyQP
}

// MultiplyByDiagMatrixBSGS multiplies the ciphertext "ctIn" by the plaintext matrix "matrix" and returns the result on the ciphertext
// "ctOut". Memory pools for the decomposed ciphertext PoolDecompQ, PoolDecompP must be provided, those are list of poly of ringQ and ringP
// respectively, each of size params.Beta().
// The BSGS approach is used (double hoisting with baby-step giant-step), which is faster than MultiplyByDiagMatrix
// for matrix with more than a few non-zero diagonals and uses much less keys.
func (eval *Evaluator) MultiplyByDiagMatrixBSGS(ctIn Ciphertext, matrix PtDiagMatrix, PoolDecompRealQP, PoolDecompImagQP []rlwe.PolyQP, ctOut Ciphertext) {

	ks := eval.GetKeySwitcher()

	// N1*N2 = N
	N1 := matrix.N1

	ringQ := eval.params.RingQ()
	ringP := eval.params.RingP()
	ringQP := eval.params.RingQP()

	levelQ := utils.MinInt(ctOut.Level(), utils.MinInt(ctIn.Level(), matrix.Level))
	levelP := len(ringP.Modulus) - 1

	QiOverF := eval.params.QiOverflowMargin(levelQ)
	PiOverF := eval.params.PiOverflowMargin(levelP)

	// Computes the rotations indexes of the non-zero rows of the diagonalized DFT matrix for the baby-step giang-step algorithm

	index, rotations := ckks.BsgsIndex(matrix.VecReal, 1<<matrix.LogSlots, matrix.N1)

	var ctInTmpReal0, ctInTmpReal1 *ring.Poly
	var ctInTmpImag0, ctInTmpImag1 *ring.Poly
	if ctIn == ctOut {
		ring.CopyValuesLvl(levelQ, ctIn.Real.Value[0], eval.ctxpool.Real.Value[0])
		ring.CopyValuesLvl(levelQ, ctIn.Real.Value[1], eval.ctxpool.Real.Value[1])
		ring.CopyValuesLvl(levelQ, ctIn.Imag.Value[0], eval.ctxpool.Imag.Value[0])
		ring.CopyValuesLvl(levelQ, ctIn.Imag.Value[1], eval.ctxpool.Imag.Value[1])
		ctInTmpReal0, ctInTmpReal1 = eval.ctxpool.Real.Value[0], eval.ctxpool.Real.Value[1]
		ctInTmpImag0, ctInTmpImag1 = eval.ctxpool.Imag.Value[0], eval.ctxpool.Imag.Value[1]
	} else {
		ctInTmpReal0, ctInTmpReal1 = ctIn.Real.Value[0], ctIn.Real.Value[1]
		ctInTmpImag0, ctInTmpImag1 = ctIn.Imag.Value[0], ctIn.Imag.Value[1]
	}

	// Pre-rotates ciphertext for the baby-step giant-step algorithm, does not divide by P yet
	ctInRotRealQP := eval.RotateHoistedNoModDownNew(levelQ, rotations, ctInTmpReal0, PoolDecompRealQP)
	ctInRotImagQP := eval.RotateHoistedNoModDownNew(levelQ, rotations, ctInTmpImag0, PoolDecompImagQP)

	// Accumulator inner loop
	tmpReal0QP := ringQP.NewPolyLvl(levelQ, levelP)
	tmpReal1QP := ringQP.NewPolyLvl(levelQ, levelP)
	tmpImag0QP := ringQP.NewPolyLvl(levelQ, levelP)
	tmpImag1QP := ringQP.NewPolyLvl(levelQ, levelP)

	// Accumulator outer loop
	c0RealQP := ringQP.NewPolyLvl(levelQ, levelP)
	c1RealQP := ringQP.NewPolyLvl(levelQ, levelP)
	c0ImagQP := ringQP.NewPolyLvl(levelQ, levelP)
	c1ImagQP := ringQP.NewPolyLvl(levelQ, levelP)

	// Result in QP
	c0OutRealQP := rlwe.PolyQP{Q: ctOut.Real.Value[0], P: ringP.NewPolyLvl(levelP)}
	c1OutRealQP := rlwe.PolyQP{Q: ctOut.Real.Value[1], P: ringP.NewPolyLvl(levelP)}
	c0OutImagQP := rlwe.PolyQP{Q: ctOut.Imag.Value[0], P: ringP.NewPolyLvl(levelP)}
	c1OutImagQP := rlwe.PolyQP{Q: ctOut.Imag.Value[1], P: ringP.NewPolyLvl(levelP)}

	ringQ.MulScalarBigintLvl(levelQ, ctInTmpReal0, ringP.ModulusBigint, ctInTmpReal0) // P*c0_real
	ringQ.MulScalarBigintLvl(levelQ, ctInTmpReal1, ringP.ModulusBigint, ctInTmpReal1) // P*c1_real
	ringQ.MulScalarBigintLvl(levelQ, ctInTmpImag0, ringP.ModulusBigint, ctInTmpImag0) // P*c0_imag
	ringQ.MulScalarBigintLvl(levelQ, ctInTmpImag1, ringP.ModulusBigint, ctInTmpImag1) // P*c1_imag

	// OUTER LOOP
	var cnt0 int
	for j := range index {
		// INNER LOOP
		var cnt1 int
		for _, i := range index[j] {
			if i == 0 {
				if cnt1 == 0 {
					ringQ.MulCoeffsMontgomeryLvl(levelQ, matrix.VecReal[N1*j].Q, ctInTmpReal0, tmpReal0QP.Q)
					ringQ.MulCoeffsMontgomeryLvl(levelQ, matrix.VecReal[N1*j].Q, ctInTmpReal1, tmpReal1QP.Q)
					ringQ.MulCoeffsMontgomeryAndSubLvl(levelQ, matrix.VecImag[N1*j].Q, ctInTmpImag0, tmpReal0QP.Q)
					ringQ.MulCoeffsMontgomeryAndSubLvl(levelQ, matrix.VecImag[N1*j].Q, ctInTmpImag1, tmpReal1QP.Q)

					ringQ.MulCoeffsMontgomeryLvl(levelQ, matrix.VecReal[N1*j].Q, ctInTmpImag0, tmpImag0QP.Q)
					ringQ.MulCoeffsMontgomeryLvl(levelQ, matrix.VecReal[N1*j].Q, ctInTmpImag1, tmpImag1QP.Q)
					ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.VecImag[N1*j].Q, ctInTmpReal0, tmpImag0QP.Q)
					ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.VecImag[N1*j].Q, ctInTmpReal1, tmpImag1QP.Q)

					tmpReal0QP.P.Zero()
					tmpReal1QP.P.Zero()
					tmpImag0QP.P.Zero()
					tmpImag1QP.P.Zero()
				} else {
					ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.VecReal[N1*j].Q, ctInTmpReal0, tmpReal0QP.Q)
					ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.VecReal[N1*j].Q, ctInTmpReal1, tmpReal1QP.Q)
					ringQ.MulCoeffsMontgomeryAndSubLvl(levelQ, matrix.VecImag[N1*j].Q, ctInTmpImag0, tmpReal0QP.Q)
					ringQ.MulCoeffsMontgomeryAndSubLvl(levelQ, matrix.VecImag[N1*j].Q, ctInTmpImag1, tmpReal1QP.Q)

					ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.VecReal[N1*j].Q, ctInTmpImag0, tmpImag0QP.Q)
					ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.VecReal[N1*j].Q, ctInTmpImag1, tmpImag1QP.Q)
					ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.VecImag[N1*j].Q, ctInTmpReal0, tmpImag0QP.Q)
					ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.VecImag[N1*j].Q, ctInTmpReal1, tmpImag1QP.Q)
				}
			} else {
				if cnt1 == 0 {
					ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, matrix.VecReal[N1*j+i], ctInRotRealQP[i][0], tmpReal0QP)
					ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, matrix.VecReal[N1*j+i], ctInRotRealQP[i][1], tmpReal1QP)
					ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, matrix.VecImag[N1*j+i], ctInRotImagQP[i][0], tmpReal0QP)
					ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, matrix.VecImag[N1*j+i], ctInRotImagQP[i][1], tmpReal1QP)

					ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, matrix.VecReal[N1*j+i], ctInRotImagQP[i][0], tmpImag0QP)
					ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, matrix.VecReal[N1*j+i], ctInRotImagQP[i][1], tmpImag1QP)
					ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, matrix.VecImag[N1*j+i], ctInRotRealQP[i][0], tmpImag0QP)
					ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, matrix.VecImag[N1*j+i], ctInRotRealQP[i][1], tmpImag1QP)
				} else {
					ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, matrix.VecReal[N1*j+i], ctInRotRealQP[i][0], tmpReal0QP)
					ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, matrix.VecReal[N1*j+i], ctInRotRealQP[i][1], tmpReal1QP)
					ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, matrix.VecImag[N1*j+i], ctInRotImagQP[i][0], tmpReal0QP)
					ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, matrix.VecImag[N1*j+i], ctInRotImagQP[i][1], tmpReal1QP)

					ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, matrix.VecReal[N1*j+i], ctInRotImagQP[i][0], tmpImag0QP)
					ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, matrix.VecReal[N1*j+i], ctInRotImagQP[i][1], tmpImag1QP)
					ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, matrix.VecImag[N1*j+i], ctInRotRealQP[i][0], tmpImag0QP)
					ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, matrix.VecImag[N1*j+i], ctInRotRealQP[i][1], tmpImag1QP)
				}
			}

			if cnt1%(QiOverF>>1) == (QiOverF>>1)-1 {
				ringQ.ReduceLvl(levelQ, tmpReal0QP.Q, tmpReal0QP.Q)
				ringQ.ReduceLvl(levelQ, tmpReal1QP.Q, tmpReal1QP.Q)
				ringQ.ReduceLvl(levelQ, tmpImag0QP.Q, tmpImag0QP.Q)
				ringQ.ReduceLvl(levelQ, tmpImag1QP.Q, tmpImag1QP.Q)
			}

			if cnt1%(PiOverF>>1) == (PiOverF>>1)-1 {
				ringP.ReduceLvl(levelP, tmpReal0QP.P, tmpReal0QP.P)
				ringP.ReduceLvl(levelP, tmpReal1QP.P, tmpReal1QP.P)
				ringP.ReduceLvl(levelP, tmpImag0QP.P, tmpImag0QP.P)
				ringP.ReduceLvl(levelP, tmpImag1QP.P, tmpImag1QP.P)
			}

			cnt1++
		}

		if cnt1%(QiOverF>>1) != 0 {
			ringQ.ReduceLvl(levelQ, tmpReal0QP.Q, tmpReal0QP.Q)
			ringQ.ReduceLvl(levelQ, tmpReal1QP.Q, tmpReal1QP.Q)
			ringQ.ReduceLvl(levelQ, tmpImag0QP.Q, tmpImag0QP.Q)
			ringQ.ReduceLvl(levelQ, tmpImag1QP.Q, tmpImag1QP.Q)
		}

		if cnt1%(PiOverF>>1) != 0 {
			ringP.ReduceLvl(levelP, tmpReal0QP.P, tmpReal0QP.P)
			ringP.ReduceLvl(levelP, tmpReal1QP.P, tmpReal1QP.P)
			ringP.ReduceLvl(levelP, tmpImag0QP.P, tmpImag0QP.P)
			ringP.ReduceLvl(levelP, tmpImag1QP.P, tmpImag1QP.P)
		}

		// If j != 0, then rotates ((tmp0QP.Q, tmp0QP.P), (tmp1QP.Q, tmp1QP.P)) by N1*j and adds the result on ((c0QP.Q, c0QP.P), (c1QP.Q, c1QP.P))
		if j != 0 {

			// Hoisting of the ModDown of sum(sum(phi(d1) * plaintext))
			ks.Baseconverter.ModDownQPtoQNTT(levelQ, levelP, tmpReal1QP.Q, tmpReal1QP.P, tmpReal1QP.Q) // c1 * plaintext + sum(phi(d1) * plaintext) + phi(c1) * plaintext mod Q
			ks.Baseconverter.ModDownQPtoQNTT(levelQ, levelP, tmpImag1QP.Q, tmpImag1QP.P, tmpImag1QP.Q)

			galEl := eval.params.GaloisElementForColumnRotationBy(N1 * j)

			rtk, generated := eval.rtks.Keys[galEl]
			if !generated {
				panic("switching key not available")
			}

			index := eval.permuteNTTIndex[galEl]

			tmpReal1QP.Q.IsNTT = true
			tmpImag1QP.Q.IsNTT = true
			ks.SwitchKeysInPlaceNoModDown(levelQ, tmpReal1QP.Q, rtk, c0RealQP.Q, c0RealQP.P, c1RealQP.Q, c1RealQP.P) // Switchkey(P*phi(tmpRes_1)) = (d0, d1) in base QP
			ks.SwitchKeysInPlaceNoModDown(levelQ, tmpImag1QP.Q, rtk, c0ImagQP.Q, c0ImagQP.P, c1ImagQP.Q, c1ImagQP.P)

			ringQP.AddLvl(levelQ, levelP, c0RealQP, tmpReal0QP, c0RealQP)
			ringQP.AddLvl(levelQ, levelP, c0ImagQP, tmpImag0QP, c0ImagQP)

			// Outer loop rotations
			if cnt0 == 0 {
				ringQP.PermuteNTTWithIndexLvl(levelQ, levelP, c0RealQP, index, c0OutRealQP)
				ringQP.PermuteNTTWithIndexLvl(levelQ, levelP, c1RealQP, index, c1OutRealQP)
				ringQP.PermuteNTTWithIndexLvl(levelQ, levelP, c0ImagQP, index, c0OutImagQP)
				ringQP.PermuteNTTWithIndexLvl(levelQ, levelP, c1ImagQP, index, c1OutImagQP)
			} else {
				ringQP.PermuteNTTWithIndexAndAddNoModLvl(levelQ, levelP, c0RealQP, index, c0OutRealQP)
				ringQP.PermuteNTTWithIndexAndAddNoModLvl(levelQ, levelP, c1RealQP, index, c1OutRealQP)
				ringQP.PermuteNTTWithIndexAndAddNoModLvl(levelQ, levelP, c0ImagQP, index, c0OutImagQP)
				ringQP.PermuteNTTWithIndexAndAddNoModLvl(levelQ, levelP, c1ImagQP, index, c1OutImagQP)
			}

			// Else directly adds on ((c0QP.Q, c0QP.P), (c1QP.Q, c1QP.P))
		} else {
			if cnt0 == 0 {
				ringQP.CopyValuesLvl(levelQ, levelP, tmpReal0QP, c0OutRealQP)
				ringQP.CopyValuesLvl(levelQ, levelP, tmpReal1QP, c1OutRealQP)
				ringQP.CopyValuesLvl(levelQ, levelP, tmpImag0QP, c0OutImagQP)
				ringQP.CopyValuesLvl(levelQ, levelP, tmpImag1QP, c1OutImagQP)
			} else {
				ringQP.AddNoModLvl(levelQ, levelP, c0OutRealQP, tmpReal0QP, c0OutRealQP)
				ringQP.AddNoModLvl(levelQ, levelP, c1OutRealQP, tmpReal1QP, c1OutRealQP)
				ringQP.AddNoModLvl(levelQ, levelP, c0OutImagQP, tmpImag0QP, c0OutImagQP)
				ringQP.AddNoModLvl(levelQ, levelP, c1OutImagQP, tmpImag1QP, c1OutImagQP)
			}
		}

		if cnt0%QiOverF == QiOverF-1 {
			ringQ.ReduceLvl(levelQ, ctOut.Real.Value[0], ctOut.Real.Value[0])
			ringQ.ReduceLvl(levelQ, ctOut.Real.Value[1], ctOut.Real.Value[1])
			ringQ.ReduceLvl(levelQ, ctOut.Imag.Value[0], ctOut.Imag.Value[0])
			ringQ.ReduceLvl(levelQ, ctOut.Imag.Value[1], ctOut.Imag.Value[1])
		}

		if cnt0%PiOverF == PiOverF-1 {
			ringP.ReduceLvl(levelP, c0OutRealQP.P, c0OutRealQP.P)
			ringP.ReduceLvl(levelP, c1OutRealQP.P, c1OutRealQP.P)
			ringP.ReduceLvl(levelP, c0OutImagQP.P, c0OutImagQP.P)
			ringP.ReduceLvl(levelP, c1OutImagQP.P, c1OutImagQP.P)
		}

		cnt0++
	}

	if cnt0%QiOverF != 0 {
		ringQ.ReduceLvl(levelQ, ctOut.Real.Value[0], ctOut.Real.Value[0])
		ringQ.ReduceLvl(levelQ, ctOut.Real.Value[1], ctOut.Real.Value[1])
		ringQ.ReduceLvl(levelQ, ctOut.Imag.Value[0], ctOut.Imag.Value[0])
		ringQ.ReduceLvl(levelQ, ctOut.Imag.Value[1], ctOut.Imag.Value[1])
	}

	if cnt0%PiOverF != 0 {
		ringP.ReduceLvl(levelP, c0OutRealQP.P, c0OutRealQP.P)
		ringP.ReduceLvl(levelP, c1OutRealQP.P, c1OutRealQP.P)
		ringP.ReduceLvl(levelP, c0OutImagQP.P, c0OutImagQP.P)
		ringP.ReduceLvl(levelP, c1OutImagQP.P, c1OutImagQP.P)
	}

	ks.Baseconverter.ModDownQPtoQNTT(levelQ, levelP, ctOut.Real.Value[0], c0OutRealQP.P, ctOut.Real.Value[0]) // sum(phi(c0 * P + d0_QP))/P
	ks.Baseconverter.ModDownQPtoQNTT(levelQ, levelP, ctOut.Real.Value[1], c1OutRealQP.P, ctOut.Real.Value[1]) // sum(phi(d1_QP))/P
	ks.Baseconverter.ModDownQPtoQNTT(levelQ, levelP, ctOut.Imag.Value[0], c0OutImagQP.P, ctOut.Imag.Value[0]) // sum(phi(c0 * P + d0_QP))/P
	ks.Baseconverter.ModDownQPtoQNTT(levelQ, levelP, ctOut.Imag.Value[1], c1OutImagQP.P, ctOut.Imag.Value[1]) // sum(phi(d1_QP))/P

	ctOut.Real.Scale = matrix.Scale * ctIn.Real.Scale
	ctOut.Imag.Scale = ctOut.Real.Scale

	ctInRotRealQP = nil
	ctInRotImagQP = nil
	runtime.GC()
}
