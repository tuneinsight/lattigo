package main

import (
	"fmt"
	"math"
	"time"

	"github.com/tuneinsight/lattigo/v3/ckks"
	ckksAdvanced "github.com/tuneinsight/lattigo/v3/ckks/advanced"
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// This example is an implementation of the RLWE -> LWE extraction followed by an LWE -> RLWE repacking
// (bridge between CKKS and FHEW ciphertext) based on "Pegasus: Bridging Polynomial and Non-polynomial
// Evaluations in Homomorphic Encryption".
// It showcases advanced tools of the CKKS scheme, such as homomorphic decoding and homomorphic modular reduction.

func main() {

	// Ring Learning With Error parameters
	fmt.Printf("Gen RLWE Parameters... ")
	start := time.Now()
	paramsRLWE := genRLWEParameters()
	fmt.Printf("Done (%s)\n", time.Since(start))

	// Learning With Error parameters
	fmt.Printf("Gen LWE Parameters... ")
	start = time.Now()
	paramsLWE := genLWEParameters(paramsRLWE)
	fmt.Printf("Done (%s)\n", time.Since(start))

	fmt.Printf("RLWE Params : logN=%2d, logQP=%3d\n", paramsRLWE.LogN(), paramsRLWE.LogQP())
	fmt.Printf("LWE  Params : logN=%2d, logQP=%3d\n", paramsLWE.LogN(), paramsLWE.LogQP())

	// Homomorphic decoding parameters
	SlotsToCoeffsParameters := ckksAdvanced.EncodingMatrixLiteral{
		LogN:                paramsRLWE.LogN(),
		LogSlots:            paramsRLWE.LogSlots(),
		Scaling:             1.0,
		LinearTransformType: ckksAdvanced.SlotsToCoeffs,
		LevelStart:          2,     // starting level
		BSGSRatio:           4.0,   // ratio between n1/n2 for n1*n2 = slots
		BitReversed:         false, // bit-reversed input
		ScalingFactor: [][]float64{ // Decomposition level of the encoding matrix
			{paramsRLWE.QiFloat64(1)}, // Scale of the second matrix
			{paramsRLWE.QiFloat64(2)}, // Scale of the first matrix
		},
	}

	// Homomorphic modular reduction parameters
	EvalModParameters := ckksAdvanced.EvalModLiteral{
		Q:             paramsRLWE.Q()[0],         // Modulus
		LevelStart:    paramsRLWE.MaxLevel() - 1, // Starting level of the procedure
		SineType:      ckksAdvanced.Cos1,         // Type of approximation
		MessageRatio:  256.0,                     // Q/|m|
		K:             16,                        // Interval of approximation
		SineDeg:       63,                        // Degree of approximation
		DoubleAngle:   2,                         // Number of double angle evaluation
		ArcSineDeg:    0,                         // Degree of arcsine Taylor polynomial
		ScalingFactor: 1 << 60,                   // Scaling factor during the procedure
	}

	// Generates the homomorphic modular reduction polynomial approximation
	fmt.Printf("Gen EvalMod Poly... ")
	start = time.Now()
	EvalModPoly := ckksAdvanced.NewEvalModPolyFromLiteral(EvalModParameters)
	fmt.Printf("Done (%s)\n", time.Since(start))

	// RLWE Parameters
	start = time.Now()
	encoder := ckks.NewEncoder(paramsRLWE)
	kgenRLWE := ckks.NewKeyGenerator(paramsRLWE)
	skRLWE := kgenRLWE.GenSecretKey()
	encryptor := ckks.NewEncryptor(paramsRLWE, skRLWE)
	decryptor := ckks.NewDecryptor(paramsRLWE, skRLWE)

	fmt.Printf("Gen SlotsToCoeffs Matrices... ")
	start = time.Now()
	SlotsToCoeffsMatrix := ckksAdvanced.NewHomomorphicEncodingMatrixFromLiteral(SlotsToCoeffsParameters, encoder)
	fmt.Printf("Done (%s)\n", time.Since(start))

	fmt.Printf("Gen Evaluation Keys:\n")
	fmt.Printf("	Decoding Keys... ")
	start = time.Now()
	rotKey := kgenRLWE.GenRotationKeysForRotations(SlotsToCoeffsParameters.Rotations(), true, skRLWE)
	fmt.Printf("Done (%s)\n", time.Since(start))
	fmt.Printf("	Relinearization Key... ")
	start = time.Now()
	rlk := kgenRLWE.GenRelinearizationKey(skRLWE, 2)
	fmt.Printf("Done (%s)\n", time.Since(start))

	fmt.Printf("	Repacking Keys... ")
	nonzerodiags := make([]int, paramsRLWE.Slots())
	for i := range nonzerodiags {
		nonzerodiags[i] = i
	}
	rotationsRepack := paramsRLWE.RotationsForLinearTransform(nonzerodiags, paramsRLWE.Slots(), 4.0)
	rotationsRepack = append(rotationsRepack, paramsRLWE.RotationsForTrace(paramsRLWE.LogSlots(), paramsLWE.LogN())...)
	rotKeyRepack := kgenRLWE.GenRotationKeysForRotations(rotationsRepack, false, skRLWE)

	fmt.Printf("Done (%s)\n", time.Since(start))

	eval := ckksAdvanced.NewEvaluator(paramsRLWE, rlwe.EvaluationKey{Rlk: rlk, Rtks: rotKey})

	// LWE Parameters
	kgenLWE := ckks.NewKeyGenerator(paramsLWE)
	skLWE := kgenLWE.GenSecretKeyWithHammingWeight(64)

	// RLWE -> LWE Switching key
	fmt.Printf("	RLWE -> LWE Switching Key... ")
	start = time.Now()
	swkRLWEDimToLWEDim := kgenRLWE.GenSwitchingKey(skRLWE, skLWE)
	fmt.Printf("Done (%s)\n", time.Since(start))

	// Encodes and Encrypts skLWE
	fmt.Printf("Encode & Encrypt SK LWE... ")
	start = time.Now()
	skLWEInvNTT := paramsLWE.RingQ().NewPoly()
	ring.CopyValues(skLWE.Value.Q, skLWEInvNTT)
	paramsLWE.RingQ().InvNTT(skLWEInvNTT, skLWEInvNTT)
	Q := paramsRLWE.Q()[0]
	paramsLWE.RingQ().InvMFormLvl(0, skLWEInvNTT, skLWEInvNTT)
	skFloat := make([]float64, paramsLWE.N())
	for i, s := range skLWEInvNTT.Coeffs[0] {
		if s >= Q>>1 {
			skFloat[i] = -float64(Q - s)
		} else {
			skFloat[i] = float64(s)
		}

		skFloat[i] *= math.Pow(1.0/(EvalModPoly.K()*EvalModPoly.QDiff()), 0.5) // sqrt(pre-scaling for Cheby)
	}

	paramsLWE.RingQ().MFormLvl(0, skLWEInvNTT, skLWEInvNTT)
	ptSk := ckks.NewPlaintext(paramsRLWE, paramsRLWE.MaxLevel(), paramsRLWE.QiFloat64(paramsRLWE.MaxLevel()))
	encoder.Encode(skFloat, ptSk, paramsLWE.LogN())
	ctSk := encryptor.EncryptNew(ptSk)
	fmt.Printf("Done (%s)\n", time.Since(start))

	// ********** PLAINTEXT GENERATION & ENCRYPTION **************

	// Random complex plaintext encrypted
	fmt.Printf("Gen Plaintext & Encrypt... ")
	start = time.Now()
	values := make([]complex128, paramsRLWE.Slots())
	for i := range values {
		values[i] = complex(float64(i+1)/float64(paramsRLWE.Slots()), 1+float64(i+1)/float64(paramsRLWE.Slots()))
	}

	plaintext := ckks.NewPlaintext(paramsRLWE, paramsRLWE.MaxLevel(), paramsRLWE.DefaultScale())
	// Must encode with 2*Slots because a real vector is returned
	encoder.Encode(values, plaintext, paramsRLWE.LogSlots())
	ct := encryptor.EncryptNew(plaintext)
	fmt.Printf("Done (%s)\n", time.Since(start))

	// ******** STEP 1 : HOMOMORPHIC DECODING *******
	fmt.Printf("Homomorphic Decoding... ")
	start = time.Now()
	ct = eval.SlotsToCoeffsNew(ct, nil, SlotsToCoeffsMatrix)
	fmt.Printf("Done (%s)\n", time.Since(start))

	// ******** STEP 2 : RLWE -> LWE EXTRACTION *************

	fmt.Printf("RLWE -> LWE Extraction... ")
	start = time.Now()
	// Scale the message to Delta = Q/MessageRatio
	scale := math.Exp2(math.Round(math.Log2(float64(EvalModParameters.Q) / EvalModParameters.MessageRatio)))
	eval.ScaleUp(ct, math.Round(scale/ct.Scale), ct)

	// Switch from RLWE parameters to LWE parameters
	ctTmp := eval.SwitchKeysNew(ct, swkRLWEDimToLWEDim)
	ctLWE := ckks.NewCiphertext(paramsLWE, 1, 0, ctTmp.Scale)

	// Switch the ciphertext outside of the NTT domain for the LWE extraction
	for i := range ctLWE.Value {
		paramsRLWE.RingQ().InvNTTLvl(0, ctTmp.Value[i], ctTmp.Value[i])
	}

	rlwe.SwitchCiphertextRingDegree(ctTmp.El(), ctLWE.El())

	// RLWE -> LWE Extraction
	lweReal, lweImag := ExtractLWESamplesBitReversed(ctLWE, paramsLWE)
	fmt.Printf("Done (%s)\n", time.Since(start))

	// Visual of some values
	fmt.Println("Visual Comparison :")
	fmt.Printf("Slot %4d : RLWE %f LWE %f\n", 0, values[0], complex(DecryptLWE(paramsLWE.RingQ(), lweReal[0], scale, skLWEInvNTT), DecryptLWE(paramsLWE.RingQ(), lweImag[0], scale, skLWEInvNTT)))
	fmt.Printf("Slot %4d : RLWE %f LWE %f\n", paramsLWE.Slots()-1, values[paramsLWE.Slots()-1], complex(DecryptLWE(paramsLWE.RingQ(), lweReal[paramsLWE.Slots()-1], scale, skLWEInvNTT), DecryptLWE(paramsLWE.RingQ(), lweImag[paramsLWE.Slots()-1], scale, skLWEInvNTT)))

	// ********* STEP 3 : LWE -> RLWE REPACKING
	fmt.Printf("Encode LWE Samples... ")
	start = time.Now()
	ptLWE := ckks.NewPlaintext(paramsRLWE, paramsRLWE.MaxLevel(), 1.0)

	// Encode the LWE samples as a vector
	lweEncoded := make([]complex128, paramsRLWE.Slots())
	for i := 0; i < paramsRLWE.Slots(); i++ {
		lweEncoded[i] = complex(float64(lweReal[i].b), float64(lweImag[i].b))
		lweEncoded[i] *= complex(math.Pow(1/(EvalModPoly.K()*EvalModPoly.QDiff()), 1.0), 0) // pre-scaling for Cheby
	}

	encoder.Encode(lweEncoded, ptLWE, paramsRLWE.LogSlots())
	fmt.Printf("Done (%s)\n", time.Since(start))

	// Encode A

	// Compute               skLeft                  skRight
	//         __________      _       __________       _
	//        |          |    | |     |          |     | |
	//        |          |    | |     |          |     | |
	//   n    |  ALeft   |  x | |  +  |  ARight  |  x  | | 	=   A x sk
	//        |          |    | |     |          |     | |
	//        |__________|    |_|     |__________|     |_|
	//
	//            N/n          1          N/n           1

	fmt.Printf("Encode A... ")
	start = time.Now()
	AVectors := make([][]complex128, paramsLWE.Slots())
	for i := range AVectors {
		tmp := make([]complex128, paramsLWE.N())
		for j := 0; j < paramsLWE.N(); j++ {
			tmp[j] = complex(float64(lweReal[i].a[j]), float64(lweImag[i].a[j]))
			tmp[j] *= complex(math.Pow(1/(EvalModPoly.K()*EvalModPoly.QDiff()), 0.5), 0) // sqrt(pre-scaling for Cheby)
		}

		AVectors[i] = tmp
	}

	// Diagonalize
	AMatDiag := make(map[int][]complex128)
	for i := 0; i < paramsLWE.Slots(); i++ {
		tmp := make([]complex128, paramsLWE.N())
		for j := 0; j < paramsLWE.N(); j++ {
			tmp[j] = AVectors[j%paramsLWE.Slots()][(j+i)%paramsLWE.N()]
		}
		AMatDiag[i] = tmp
	}

	linTransf := ckks.GenLinearTransformBSGS(encoder, AMatDiag, paramsRLWE.MaxLevel(), 1.0, 4.0, paramsLWE.LogN())
	fmt.Printf("Done (%s)\n", time.Since(start))

	evalRepack := ckks.NewEvaluator(paramsRLWE, rlwe.EvaluationKey{Rlk: rlk, Rtks: rotKeyRepack})

	fmt.Printf("Homomorphic Partial Decryption : pt = A x sk + encode(LWE) + I(X)*Q... ")
	start = time.Now()
	ctAs := evalRepack.LinearTransformNew(ctSk, linTransf)[0]                // A_left * sk || A_right * sk
	ctAs = evalRepack.TraceNew(ctAs, paramsLWE.LogSlots(), paramsLWE.LogN()) // A * sk || A * sk
	evalRepack.MultByConst(ctAs, paramsLWE.N()/paramsLWE.Slots(), ctAs)
	eval.Rescale(ctAs, 1.0, ctAs)
	eval.Add(ctAs, ptLWE, ctAs) // A * sk || A * sk + LWE_real || LWE_imag = RLWE + I(X) * Q
	ctAs.Scale = scale
	fmt.Printf("Done (%s)\n", time.Since(start))

	fmt.Printf("Homomorphic Modular Reduction : pt mod Q... ")
	start = time.Now()
	// Extract imaginary part : RLWE_real + I(X)*Q ; RLWE_imag + I(X)*Q
	ctAsConj := eval.ConjugateNew(ctAs)
	ctAsReal := eval.AddNew(ctAs, ctAsConj)
	ctAsImag := eval.SubNew(ctAs, ctAsConj)
	ctAsReal.Scale = ctAsReal.Scale * 2                                                                                   // Divide by 2
	ctAsImag.Scale = ctAsImag.Scale * 2                                                                                   // Divide by 2
	eval.ScaleUp(ctAsReal, math.Round((EvalModPoly.ScalingFactor()/EvalModPoly.MessageRatio())/ctAsReal.Scale), ctAsReal) // Scale the real message up to Sine/MessageRatio
	eval.ScaleUp(ctAsImag, math.Round((EvalModPoly.ScalingFactor()/EvalModPoly.MessageRatio())/ctAsImag.Scale), ctAsImag) // Scale the imag message up to Sine/MessageRatio
	ctAsReal = eval.EvalModNew(ctAsReal, EvalModPoly)                                                                     // Real mod Q
	eval.DivByi(ctAsImag, ctAsImag)
	ctAsImag = eval.EvalModNew(ctAsImag, EvalModPoly) // (-i*imag mod Q)*i
	eval.MultByi(ctAsImag, ctAsImag)
	eval.Add(ctAsReal, ctAsImag, ctAsReal) // Repack both imag and real parts
	fmt.Printf("Done (%s)\n", time.Since(start))

	fmt.Println("Visual Comparison :")
	v := encoder.DecodePublic(decryptor.DecryptNew(ctAsReal), paramsRLWE.LogSlots(), 0)
	fmt.Printf("Slot %4d : Want %f Have %f\n", 0, values[0], v[0])
	fmt.Printf("Slot %4d : Want %f Have %f\n", paramsRLWE.Slots()-1, values[paramsRLWE.Slots()-1], v[paramsRLWE.Slots()-1])

}

func genRLWEParameters() (paramsRLWE ckks.Parameters) {
	var err error
	if paramsRLWE, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:         15,
		LogSlots:     9,
		DefaultScale: 1 << 30,
		Sigma:        rlwe.DefaultSigma,
		Q: []uint64{
			0xffff820001,       // 40 Q0
			0x2000000a0001,     // 45
			0x2000000e0001,     // 45
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
			0x1fffffffffe00001, // 61
			0x1fffffffffc80001, // 61
			0x1fffffffffb40001, // 61
		},
	}); err != nil {
		panic(err)
	}
	return
}

func genLWEParameters(paramsRLWE ckks.Parameters) (paramsLWE ckks.Parameters) {
	var err error
	if paramsLWE, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:         10,
		LogSlots:     paramsRLWE.LogSlots(),
		DefaultScale: paramsRLWE.DefaultScale(),
		Sigma:        paramsRLWE.Sigma(),
		Q:            paramsRLWE.Q()[:1], // 40 Q0
		P:            paramsRLWE.P()[:1], // Pi 61
	}); err != nil {
		panic(err)
	}
	return
}

// DecryptLWE decrypts an LWE sample
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

	// Copy coefficients multiplied by X^{N-1} in reverse order:
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

		// Multiply the accumulator by X^{N/(2*slots)}
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
