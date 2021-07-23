package main

import (
	"fmt"
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
	"math"
	"time"
)

func main() {
	var err error
	

	var paramsRLWE ckks.Parameters
	if paramsRLWE, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:     10,
		LogSlots: 7,
		Scale:    1 << 30,
		Sigma:    rlwe.DefaultSigma,
		Q: []uint64{
			0xffff820001,       // 40 Q0
			0x2000000a0001,     // 45 CtS
			0x2000000e0001,     // 45 CtS
			0x10000140001,      // 40
			0xffffe80001,       // 40
			0xffffc40001,       // 40
			0x100003e0001,      // 40
			0xffffb20001,       // 40
			0x10000500001,      // 40
			0xffff940001,       // 40
			0xffff8a0001,       // 40
			0xfffffffff840001,  // 60 Sine (double angle)
			0x1000000000860001, // 60 Sine (double angle)
			0xfffffffff6a0001,  // 60 Sine
			0x1000000000980001, // 60 Sine
			0xfffffffff5a0001,  // 60 Sine
			0x1000000000b00001, // 60 Sine
			0x1000000000ce0001, // 60 Sine
			0xfffffffff2a0001,  // 60 Sine
			0x100000000060001,  // 58 Repack & Change of basis
		},
		P: []uint64{
			0x1fffffffffe00001, // Pi 61
		},
	}); err != nil {
		panic(err)
	}

	var paramsLWE ckks.Parameters
	if paramsLWE, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:     paramsRLWE.LogN(),
		LogSlots: paramsRLWE.LogSlots(),
		Scale:    paramsRLWE.Scale(),
		Sigma:    paramsRLWE.Sigma(),
		Q: []uint64{
			paramsRLWE.Q()[0],       // 40 Q0
		},
		P: paramsRLWE.P(),
	}); err != nil {
		panic(err)
	}

	SlotsToCoeffsParameters := ckks.EncodingMatricesParameters{
		LevelStart:  2,
		BSGSRatio:   16.0,
		BitReversed: false,
		ScalingFactor: [][]float64{
			{paramsRLWE.QiFloat64(1)},
			{paramsRLWE.QiFloat64(2)},
		},
	}

	EvalModParameters := ckks.EvalModParameters{
		Q:             paramsRLWE.Q()[0],
		LevelStart:    paramsRLWE.MaxLevel() - 1,
		SineType:      ckks.Cos1,
		MessageRatio:  256.0,
		K:             16,
		SineDeg:       63,
		DoubleAngle:   2,
		ArcSineDeg:    0,
		ScalingFactor: 1 << 60,
	}

	EvalModPoly := EvalModParameters.GenPoly()

	//ringQRLWE := paramsRLWE.RingQ()
	ringQLWE := paramsLWE.RingQ()
	Q := paramsRLWE.Q()[0]

	// RLWE Parameters
	encoder := ckks.NewEncoder(paramsRLWE)
	kgenRLWE := ckks.NewKeyGenerator(paramsRLWE)
	skRLWE := kgenRLWE.GenSecretKey()
	skRLWE2 := kgenRLWE.GenSecretKey()
	encryptor := ckks.NewEncryptor(paramsRLWE, skRLWE)
	decryptor := ckks.NewDecryptor(paramsRLWE, skRLWE)
	decryptor2 := ckks.NewDecryptor(paramsRLWE, skRLWE2)

	SlotsToCoeffsMatrix := encoder.GenHomomorphicEncodingMatrices(SlotsToCoeffsParameters, 1.0)
	rotations := SlotsToCoeffsParameters.Rotations(paramsRLWE.LogN(), paramsRLWE.LogSlots())
	rotKey := kgenRLWE.GenRotationKeysForRotations(rotations, true, skRLWE)
	rlk := kgenRLWE.GenRelinearizationKey(skRLWE, 2)
	eval := ckks.NewEvaluator(paramsRLWE, rlwe.EvaluationKey{Rlk: rlk, Rtks: rotKey})

	// LWE Parameters
	kgenLWE := ckks.NewKeyGenerator(paramsLWE)
	skLWE := kgenLWE.GenSecretKeySparse(64)
	decryptorLWE := ckks.NewDecryptor(paramsLWE, skLWE)
	_= decryptorLWE

	// RLWE -> LWE SWK
	swkRLWEDimToLWEDim := kgenRLWE.GenSwitchingKey(skRLWE, skRLWE2)
	
	// Random complex plaintext encrypted
	values := make([]complex128, paramsRLWE.Slots())
	for i := range values {
		values[i] = complex(float64(i+1)/float64(paramsRLWE.Slots()), 1+float64(i+1)/float64(paramsRLWE.Slots()))
	}

	plaintext := ckks.NewPlaintext(paramsRLWE, paramsRLWE.MaxLevel(), paramsRLWE.Scale())
	// Must encode with 2*Slots because a real vector is returned
	encoder.Encode(plaintext, values, utils.MinInt(paramsRLWE.LogSlots()+1, paramsRLWE.LogN()-1))
	ct := encryptor.EncryptNew(plaintext)

	// ******** STEP 1 : HOMOMORPHIC DECODING *******
	ct = eval.SlotsToCoeffs(ct, nil, SlotsToCoeffsMatrix)

	// Decrypt and print coefficient domain
	coeffsFloat := encoder.DecodeCoeffsPublic(decryptor.DecryptNew(ct), 0)
	valuesFloat := make([]complex128, paramsRLWE.Slots())
	gap := paramsRLWE.N() / (2 * paramsRLWE.Slots())
	for i, idx := 0, 0; i < paramsRLWE.Slots(); i, idx = i+1, idx+gap {
		valuesFloat[i] = complex(coeffsFloat[idx], coeffsFloat[idx+(paramsRLWE.N()>>1)])
	}
	ckks.SliceBitReverseInPlaceComplex128(valuesFloat, paramsRLWE.Slots())

	// ******** STEP 2 : RLWE -> LWE EXTRACTION *************

	// Scale the message to Delta = Q/MessageRatio
	scale := math.Exp2(math.Round(math.Log2(float64(EvalModPoly.Q) / EvalModPoly.MessageRatio)))
	eval.ScaleUp(ct, math.Round(scale/ct.Scale), ct)

	fmt.Println(encoder.DecodeCoeffsPublic(decryptor.DecryptNew(ct), 0)[:8])


	//Switch to lower dimension
	eval.SwitchKeys(ct, swkRLWEDimToLWEDim, ct)

	fmt.Println(encoder.DecodeCoeffsPublic(decryptor2.DecryptNew(ct), 0)[:8])

	ctLWE := new(ckks.Ciphertext)
	// RLWE -> LWE Extraction
	lweReal, lweImag := ExtractLWESamplesBitReversed(ctLWE, paramsLWE)

	// Encode the secret-key
	skLWEInvNTT := ringQLWE.NewPoly()
	ring.CopyValues(skLWE.Value, skLWEInvNTT)
	ringQLWE.InvNTT(skLWEInvNTT, skLWEInvNTT)

	// Visual of some values
	fmt.Println(valuesFloat[0], valuesFloat[paramsLWE.Slots()-1])
	fmt.Println(complex(DecryptLWE(ringQLWE, lweReal[0], scale, skLWEInvNTT), DecryptLWE(ringQLWE, lweImag[0], scale, skLWEInvNTT)),
		complex(DecryptLWE(ringQLWE, lweReal[paramsLWE.Slots()-1], scale, skLWEInvNTT), DecryptLWE(ringQLWE, lweImag[paramsLWE.Slots()-1], scale, skLWEInvNTT)))


	// ********* STEP 3 : LWE -> RLWE REPACKING

	ptLWE := ckks.NewPlaintext(paramsRLWE, paramsRLWE.MaxLevel(), 1.0)

	// Encode the LWE samples
	lweEncoded := make([]complex128, paramsRLWE.Slots())
	for i := 0; i < paramsRLWE.Slots(); i++ {
		lweEncoded[i] = complex(float64(lweReal[i].b), float64(lweImag[i].b))
		lweEncoded[i] *= complex(math.Pow(1/(float64(EvalModPoly.K)*EvalModPoly.QDiff()), 1.0), 0) // pre-scaling for Cheby
	}

	encoder.EncodeNTT(ptLWE, lweEncoded, paramsRLWE.LogSlots())

	// Encodes and Encrypts skLWE
	ringQLWE.InvMFormLvl(0, skLWEInvNTT, skLWEInvNTT)
	fmt.Println("Encrypt SK")
	skFloat := make([]complex128, paramsRLWE.N())
	for i, s := range skLWEInvNTT.Coeffs[0] {
		if s >= Q>>1 {
			skFloat[i] = -complex(float64(Q-s), 0)
		} else {
			skFloat[i] = complex(float64(s), 0)
		}

		skFloat[i] *= complex(math.Pow(1.0/(float64(EvalModPoly.K)*EvalModPoly.QDiff()), 0.5), 0) // sqrt(pre-scaling for Cheby)
	}

	ringQLWE.InvMFormLvl(0, skLWEInvNTT, skLWEInvNTT)


	ptSk := ckks.NewPlaintext(paramsRLWE, paramsRLWE.MaxLevel(), float64(paramsRLWE.Q()[paramsRLWE.MaxLevel()]))
	encoder.Encode(ptSk, skFloat, paramsRLWE.LogSlots())
	ctSk := encryptor.EncryptNew(ptSk)


	// Encodes A

	// Compute			     skLeft 			     skRight
	//    	   __________  	   _       __________ 	    _
	//   	  |		     |	  | |	  |		     |	   | |
	//   	  |		     |    | |	  |		     |     | |
	//   n    |	 ALeft   |  x | |  +  |  ARight  |  x  | | 	= 	A x sk
	//   	  |		     |	  | |	  |		     |     | |
	//   	  |__________|	  |_|	  |__________|     |_|
	//
	//			  N/n          1           N/n 		    1

	// Constructs matrix

	fmt.Println("Encode A")
	AVectors := make([][]complex128, paramsLWE.N())
	for i := range AVectors {
		tmp := make([]complex128, paramsLWE.N())
		for j := 0; j < paramsLWE.N(); j++ {
			tmp[j] = complex(float64(lweReal[i].a[j]), float64(lweImag[i].a[j]))
			tmp[j] *= complex(math.Pow(1/(float64(EvalModPoly.K)*EvalModPoly.QDiff()), 0.5), 0) // sqrt(pre-scaling for Cheby)
		}

		AVectors[i] = tmp
	}

	// Diagonalize
	AMatDiag := make(map[int][]complex128)
	for i := 0; i < paramsLWE.N(); i++ {
		tmp := make([]complex128, paramsLWE.N())
		for j := 0; j < paramsLWE.N(); j++ {
			tmp[j] = AVectors[j][(j+i)%paramsLWE.N()]
		}

		AMatDiag[i] = tmp
	}

	ptMatDiag := encoder.EncodeDiagMatrixBSGSAtLvl(paramsRLWE.MaxLevel(), AMatDiag, 1.0, 16.0, paramsLWE.LogN())


	fmt.Println("GenRepackKeys")
	rotKeyRepack := kgenRLWE.GenRotationKeysForRotations(paramsRLWE.RotationsForDiagMatrixMult(ptMatDiag), false, skRLWE)

	evalRepack := ckks.NewEvaluator(paramsRLWE, rlwe.EvaluationKey{Rlk: rlk, Rtks: rotKeyRepack})

	fmt.Println("Start Repacking")
	startTotal := time.Now()
	ctAs := evalRepack.LinearTransform(ctSk, ptMatDiag)[0]
	fmt.Printf("Done : %s\n", time.Since(startTotal))

	fmt.Printf("Done : %s\n", time.Since(startTotal))

	eval.Rescale(ctAs, 1.0, ctAs)
	eval.Add(ctAs, ptLWE, ctAs)

	ctAs.Scale = scale

	ctAsConj := eval.ConjugateNew(ctAs)
	ctAsReal := eval.AddNew(ctAs, ctAsConj)
	ctAsImag := eval.SubNew(ctAs, ctAsConj)

	ctAsReal.Scale = ctAsReal.Scale * 2
	ctAsImag.Scale = ctAsImag.Scale * 2

	// Scale the message up to Sine/MessageRatio
	eval.ScaleUp(ctAsReal, math.Round((EvalModPoly.ScalingFactor/EvalModPoly.MessageRatio)/ctAsReal.Scale), ctAsReal)
	eval.ScaleUp(ctAsImag, math.Round((EvalModPoly.ScalingFactor/EvalModPoly.MessageRatio)/ctAsImag.Scale), ctAsImag)

	// EvalMod
	ctAsReal = eval.EvalMod(ctAsReal, EvalModPoly)
	eval.DivByi(ctAsImag, ctAsImag)
	ctAsImag = eval.EvalMod(ctAsImag, EvalModPoly)
	eval.MultByi(ctAsImag, ctAsImag)
	eval.Add(ctAsReal, ctAsImag, ctAsReal)

	v := encoder.DecodePublic(decryptor.DecryptNew(ctAsReal), paramsRLWE.LogSlots(), 0)

	fmt.Println(v[0], v[paramsRLWE.Slots()-1])

}

//DecryptLWE decrypts an LWE sample
func DecryptLWE(ringQ *ring.Ring, lwe RNSLWESample, scale float64, skInvNTT *ring.Poly) float64 {

	tmp := ringQ.NewPolyLvl(0)
	pol := new(ring.Poly)
	pol.Coeffs = [][]uint64{lwe.a}
	ringQ.MulCoeffsMontgomeryLvl(0, pol, skInvNTT, tmp)
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

// RNSLWESample is a struct for RNS LWE samples
type RNSLWESample struct {
	b uint64
	a []uint64
}

// ExtractLWESamplesBitReversed extracts LWE samples from a R-LWE sample
func ExtractLWESamplesBitReversed(ct *ckks.Ciphertext, params ckks.Parameters) (LWEReal, LWEImag []RNSLWESample) {

	ringQ := params.RingQ()

	LWEReal = make([]RNSLWESample, params.Slots())
	LWEImag = make([]RNSLWESample, params.Slots())

	// Copies coefficients multiplied by X^{N-1} in reverse order :
	// a_{0} -a_{N-1} -a2_{N-2} ... -a_{1}
	acc := ringQ.NewPolyLvl(ct.Level())
	for i, qi := range ringQ.Modulus[:ct.Level()+1] {
		tmp0 := acc.Coeffs[i]
		tmp1 := ct.Value[1].Coeffs[i]
		tmp0[0] = tmp1[0]
		for j := 1; j < ringQ.N; j++ {
			tmp0[j] = qi - tmp1[ringQ.N-j]
		}
	}

	pol := ct.Value[0]

	gap := params.N() / (2 * params.Slots()) // Gap between plaintext coefficient if sparse packed

	// Real values
	for i, idx := 0, 0; i < params.Slots(); i, idx = i+1, idx+gap {

		iRev := utils.BitReverse64(uint64(i), uint64(params.LogSlots()))

		LWEReal[iRev].b = pol.Coeffs[0][idx]
		LWEReal[iRev].a = make([]uint64, params.N())
		copy(LWEReal[iRev].a, acc.Coeffs[0])

		// Multiplies the accumulator by X^{N/(2*slots)}
		MulBySmallMonomial(ringQ, acc, gap)
	}

	// Imaginary values
	for i, idx := 0, 0; i < params.Slots(); i, idx = i+1, idx+gap {

		iRev := utils.BitReverse64(uint64(i), uint64(params.LogSlots()))

		LWEImag[iRev].b = pol.Coeffs[0][idx+(params.N()>>1)]
		LWEImag[iRev].a = make([]uint64, params.N())
		copy(LWEImag[iRev].a, acc.Coeffs[0])

		// Multiplies the accumulator by X^{N/(2*slots)}
		MulBySmallMonomial(ringQ, acc, gap)
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
