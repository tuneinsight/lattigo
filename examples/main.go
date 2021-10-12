package main

import (
	"fmt"
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
	"math/bits"
	"runtime"
	"time"
)

// LogN of the ring degree of the used parameters
var LogN = 15

func main() {
	RLWE()
	//SimulatedComplex()
	//Complex()
}

// RLWE -> LWE -> RLWE repacking
func RLWE() {
	Q := []uint64{0x80000000080001}
	P := []uint64{0x1fffffffffe00001}

	RLWEParams := rlwe.ParametersLiteral{
		LogN:     12,
		Q:        Q,
		P:        P,
		Sigma:    rlwe.DefaultSigma,
		RingType: rlwe.RingStandard,
	}
	scale := float64(1 << 40)

	params, _ := rlwe.NewParametersFromLiteral(RLWEParams)
	ks := rlwe.NewKeySwitcher(params)
	kgen := rlwe.NewKeyGenerator(params)

	sk := kgen.GenSecretKey()
	encryptor := rlwe.NewEncryptor(params, sk)
	decryptor := rlwe.NewDecryptor(params, sk)

	rotations := []int{}
	for i := 1; i < params.N(); i <<= 1 {
		rotations = append(rotations, i)
	}

	rtks := kgen.GenRotationKeysForRotations(rotations, true, sk)

	ringQ := params.RingQ()

	plaintext := rlwe.NewPlaintext(params, params.MaxLevel())
	plaintext.Value.IsNTT = true

	for i := 0; i < ringQ.N; i++ {
		plaintext.Value.Coeffs[0][i] = uint64(float64(i+1) / float64(ringQ.N) * scale)
	}
	ringQ.NTT(plaintext.Value, plaintext.Value)

	ciphertext := rlwe.NewCiphertextNTT(params, 1, plaintext.Level())

	encryptor.Encrypt(plaintext, ciphertext)

	//Extract LWE
	ringQ.InvNTT(ciphertext.Value[0], ciphertext.Value[0])
	ringQ.InvNTT(ciphertext.Value[1], ciphertext.Value[1])
	LWE := ExtractLWESamples(ciphertext.Value[0], ciphertext.Value[1], ringQ)
	ringQ.NTT(ciphertext.Value[0], ciphertext.Value[0])
	ringQ.NTT(ciphertext.Value[1], ciphertext.Value[1])

	ciphertexts := make([]*rlwe.Ciphertext, len(LWE))

	X := ringQ.NewPoly()
	Xpow := ringQ.NewPoly()

	X.Coeffs[0][1] = 1
	Xpow.Coeffs[0][0] = 1

	ringQ.NTT(X, X)
	ringQ.NTT(Xpow, Xpow)
	ringQ.MForm(X, X)
	ringQ.MForm(Xpow, Xpow)

	level := params.MaxLevel()

	permuteNTTIndex := make(map[uint64][]uint64, len(rtks.Keys))
	for galEl := range rtks.Keys {
		permuteNTTIndex[galEl] = ringQ.PermuteNTTIndex(galEl)
	}

	galEls := make(map[uint64]uint64)
	for i := 0; i < params.LogN()-1; i++ {
		galEls[1<<i] = params.GaloisElementForColumnRotationBy(1 << i)
	}

	acc := ringQ.NewPoly()
	tmp := rlwe.NewCiphertextNTT(params, 1, plaintext.Level())

	_ = level
	_ = tmp
	_ = ks

	for i := range ciphertexts[:] {

		// Alocates ciphertext
		ciphertexts[i] = rlwe.NewCiphertextNTT(params, 1, plaintext.Level())
		ciphertexts[i].Value[0].Coeffs[0][0] = LWE[i].b

		// Copy coefficients multiplied by X^{N-1} in reverse order:
		// a_{0} -a_{N-1} -a2_{N-2} ... -a_{1}
		tmp0 := acc.Coeffs[0]
		tmp1 := LWE[i].a
		tmp0[0] = tmp1[0]
		for j := 1; j < ringQ.N; j++ {
			tmp0[j] = ringQ.Modulus[0] - tmp1[ringQ.N-j]
		}

		copy(ciphertexts[i].Value[1].Coeffs[0], acc.Coeffs[0])

		// Switches to NTT domain
		ringQ.NTT(ciphertexts[i].Value[0], ciphertexts[i].Value[0])
		ringQ.NTT(ciphertexts[i].Value[1], ciphertexts[i].Value[1])
	}

	now := time.Now()
	ciphertext = PackLWEs(ciphertexts[:4096], ks, rtks, permuteNTTIndex, params, decryptor, plaintext)
	fmt.Printf("Done : %s\n", time.Since(now))

	// Trace

	for j := 11; j < 11; j++ {
		Rotate(ciphertext, galEls[uint64(1<<j)], permuteNTTIndex, params, ks, rtks, tmp)
		ringQ.Add(ciphertext.Value[0], tmp.Value[0], ciphertext.Value[0])
		ringQ.Add(ciphertext.Value[1], tmp.Value[1], ciphertext.Value[1])
	}

	DecryptAndPrint(decryptor, ringQ.N, ringQ, ciphertext, plaintext, scale)
}

// PackLWEs repacks LWE ciphertexts into a RLWE ciphertext
func PackLWEs(ciphertexts []*rlwe.Ciphertext, ks *rlwe.KeySwitcher, rtks *rlwe.RotationKeySet, permuteNTTIndex map[uint64][]uint64, params rlwe.Parameters, decryptor rlwe.Decryptor, plaintext *rlwe.Plaintext) *rlwe.Ciphertext {

	ringQ := params.RingQ()

	L := bits.Len64(uint64(len(ciphertexts))) - 1

	if L == 0 {
		// Multiplies by N^-1
		//ring.MulScalarMontgomeryVec(ciphertexts[0].Value[0].Coeffs[0], ciphertexts[0].Value[0].Coeffs[0], ringQ.NttNInv[0], ringQ.Modulus[0], ringQ.MredParams[0])
		//ring.MulScalarMontgomeryVec(ciphertexts[0].Value[1].Coeffs[0], ciphertexts[0].Value[1].Coeffs[0], ringQ.NttNInv[0], ringQ.Modulus[0], ringQ.MredParams[0])
		return ciphertexts[0]
	}

	odd := make([]*rlwe.Ciphertext, len(ciphertexts)>>1)
	even := make([]*rlwe.Ciphertext, len(ciphertexts)>>1)

	for i := 0; i < len(ciphertexts)>>1; i++ {
		odd[i] = ciphertexts[2*i]
		even[i] = ciphertexts[2*i+1]
	}

	ctEven := PackLWEs(odd, ks, rtks, permuteNTTIndex, params, decryptor, plaintext)
	ctOdd := PackLWEs(even, ks, rtks, permuteNTTIndex, params, decryptor, plaintext)

	//DecryptAndPrint(decryptor, ringQ.N, ringQ, ctEven, plaintext, 1<<40)
	//DecryptAndPrint(decryptor, ringQ.N, ringQ, ctOdd, plaintext, 1<<40)

	tmpEven := ctEven.CopyNew()

	//X^(N/2^L)
	XpowNoverL := ringQ.NewPoly()
	XpowNoverL.Coeffs[0][ringQ.N/(1<<L)] = 1
	ringQ.NTT(XpowNoverL, XpowNoverL)
	ringQ.MForm(XpowNoverL, XpowNoverL)

	// ctOdd * X^(N/2^L)
	ringQ.MulCoeffsMontgomery(ctOdd.Value[0], XpowNoverL, ctOdd.Value[0])
	ringQ.MulCoeffsMontgomery(ctOdd.Value[1], XpowNoverL, ctOdd.Value[1])

	// ctEven + ctOdd * X^(N/2^L)
	ringQ.Add(ctEven.Value[0], ctOdd.Value[0], ctEven.Value[0])
	ringQ.Add(ctEven.Value[1], ctOdd.Value[1], ctEven.Value[1])

	// phi(ctEven - ctOdd * X^(N/2^L), 2^L+1)
	ringQ.Sub(tmpEven.Value[0], ctOdd.Value[0], tmpEven.Value[0])
	ringQ.Sub(tmpEven.Value[1], ctOdd.Value[1], tmpEven.Value[1])

	if L == 1 {
		Rotate(tmpEven, uint64(2*ringQ.N-1), permuteNTTIndex, params, ks, rtks, tmpEven)
	} else {
		Rotate(tmpEven, params.GaloisElementForColumnRotationBy(1<<(L-2)), permuteNTTIndex, params, ks, rtks, tmpEven)
	}

	// ctEven + ctOdd * X^(N/2^L) + phi(ctEven - ctOdd * X^(N/2^L), 2^L+1)
	ringQ.Add(ctEven.Value[0], tmpEven.Value[0], ctEven.Value[0])
	ringQ.Add(ctEven.Value[1], tmpEven.Value[1], ctEven.Value[1])

	return ctEven

}

// DecryptAndPrint decrypts and prints the first N values.
func DecryptAndPrint(decryptor rlwe.Decryptor, N int, ringQ *ring.Ring, ciphertext *rlwe.Ciphertext, plaintext *rlwe.Plaintext, scale float64) {
	decryptor.Decrypt(ciphertext, plaintext)
	ringQ.InvNTT(plaintext.Value, plaintext.Value)

	v := make([]float64, N)

	for j := 0; j < N; j++ {
		if plaintext.Value.Coeffs[0][j] >= ringQ.Modulus[0]>>1 {
			v[j] = -float64(ringQ.Modulus[0] - plaintext.Value.Coeffs[0][j])
		} else {
			v[j] = float64(plaintext.Value.Coeffs[0][j])
		}

		v[j] /= scale
	}

	for j := 0; j < N; j++ {
		fmt.Printf("%12.4f ", v[j])
	}
	fmt.Printf("\n")
}

// Rotate rotates a ciphertext
func Rotate(ctIn *rlwe.Ciphertext, galEl uint64, permuteNTTindex map[uint64][]uint64, params rlwe.Parameters, ks *rlwe.KeySwitcher, rtks *rlwe.RotationKeySet, ctOut *rlwe.Ciphertext) {
	ringQ := params.RingQ()
	rtk, _ := rtks.GetRotationKey(galEl)
	level := utils.MinInt(ctIn.Level(), ctOut.Level())
	index := permuteNTTindex[galEl]
	ks.SwitchKeysInPlace(level, ctIn.Value[1], rtk, ks.Pool[1].Q, ks.Pool[2].Q)
	ringQ.AddLvl(level, ks.Pool[1].Q, ctIn.Value[0], ks.Pool[1].Q)
	ringQ.PermuteNTTWithIndexLvl(level, ks.Pool[1].Q, index, ctOut.Value[0])
	ringQ.PermuteNTTWithIndexLvl(level, ks.Pool[2].Q, index, ctOut.Value[1])
}

// LWESample is a struct for RNS LWE samples
type LWESample struct {
	b uint64
	a []uint64
}

// DecryptLWE decrypts an LWE sample
func DecryptLWE(ringQ *ring.Ring, lwe LWESample, scale float64, skInvNTT *ring.Poly) float64 {

	tmp := ringQ.NewPoly()
	pol := new(ring.Poly)
	pol.Coeffs = [][]uint64{lwe.a}
	ringQ.MulCoeffsMontgomery(pol, skInvNTT, tmp)
	qi := ringQ.Modulus[0]
	tmp0 := tmp.Coeffs[0]
	tmp1 := lwe.b
	for j := 0; j < ringQ.N; j++ {
		tmp1 = ring.CRed(tmp1+tmp0[j], qi)
	}

	if tmp1 >= ringQ.Modulus[0]>>1 {
		tmp1 = ringQ.Modulus[0] - tmp1
		return -float64(tmp1) / scale
	}

	return float64(tmp1) / scale
}

// ExtractLWESamples extracts LWE samples from a R-LWE sample
func ExtractLWESamples(b, a *ring.Poly, ringQ *ring.Ring) (LWE []LWESample) {

	N := ringQ.N

	LWE = make([]LWESample, N)

	// Copy coefficients multiplied by X^{N-1} in reverse order:
	// a_{0} -a_{N-1} -a2_{N-2} ... -a_{1}
	acc := ringQ.NewPoly()
	for i, qi := range ringQ.Modulus {
		tmp0 := acc.Coeffs[i]
		tmp1 := a.Coeffs[i]
		tmp0[0] = tmp1[0]
		for j := 1; j < N; j++ {
			tmp0[j] = qi - tmp1[ringQ.N-j]
		}
	}

	pol := b

	// Real values
	for i := 0; i < N; i++ {

		LWE[i].b = pol.Coeffs[0][i]
		LWE[i].a = make([]uint64, N)
		copy(LWE[i].a, acc.Coeffs[0])

		// Multiplies the accumulator by X
		MulBySmallMonomial(ringQ, acc, 1)
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

// SimulatedComplex implements complex arithmetic using RCKKS
func SimulatedComplex() {
	// Schemes parameters are created from scratch
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:     LogN,
		LogQ:     []int{55, 40, 40, 40},
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
	_ = ctB

	/*
		var ctC Ciphertext
		var tot time.Duration
		for i := 0; i < 1; i++{
			now := time.Now()
			ctC = eval.MulRelinNew(ctA, ctB)
			eval.Rescale(ctC.Real, eval.params.Scale(), ctC.Real)
			eval.Rescale(ctC.Imag, eval.params.Scale(), ctC.Imag)
			tot += time.Since(now)
		}
		fmt.Println("RCKKS - LogN :", params.LogN())
		fmt.Printf("Done : %s\n", tot/100.)
	*/

	ctC := ctA

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
