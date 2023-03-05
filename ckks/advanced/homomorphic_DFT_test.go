package advanced

import (
	"flag"
	"runtime"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring/distribution"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")

var minPrec float64 = 15

func TestHomomorphicEncoding(t *testing.T) {
	var err error

	if runtime.GOARCH == "wasm" {
		t.Skip("skipping homomorphic encoding tests for GOARCH=wasm")
	}

	ParametersLiteral := ckks.ParametersLiteral{
		LogN: 13,
		Q: []uint64{
			0x10000000006e0001, // 60 Q0
			0x2000000a0001,     // 45
			0x2000000e0001,     // 45
			0x1fffffc20001,     // 45
		},
		P: []uint64{
			0x1fffffffffe00001, // Pi 61
			0x1fffffffffc80001, // Pi 61
		},

		Xs:       &distribution.Ternary{H: 192},
		LogScale: 45,
	}

	testHomomorphicDFTMatrixLiteralMarshalling(t)

	var params ckks.Parameters
	if params, err = ckks.NewParametersFromLiteral(ParametersLiteral); err != nil {
		panic(err)
	}

	for _, logSlots := range []int{params.MaxLogSlots() - 1, params.MaxLogSlots()} {
		for _, testSet := range []func(params ckks.Parameters, logSlots int, t *testing.T){
			testCoeffsToSlots,
			testSlotsToCoeffs,
		} {
			testSet(params, logSlots, t)
			runtime.GC()
		}
	}
}

func testHomomorphicDFTMatrixLiteralMarshalling(t *testing.T) {
	t.Run("Marshalling", func(t *testing.T) {
		m := HomomorphicDFTMatrixLiteral{
			LogSlots:        15,
			Type:            Decode,
			LevelStart:      12,
			LogBSGSRatio:    2,
			Levels:          []int{1, 1, 1},
			BitReversed:     true,
			RepackImag2Real: true,
		}

		data, err := m.MarshalBinary()
		assert.Nil(t, err)

		mNew := new(HomomorphicDFTMatrixLiteral)
		if err := mNew.UnmarshalBinary(data); err != nil {
			assert.Nil(t, err)
		}
		assert.Equal(t, m, *mNew)
	})
}

func testCoeffsToSlots(params ckks.Parameters, LogSlots int, t *testing.T) {

	slots := 1 << LogSlots

	var sparse bool = LogSlots < params.MaxLogSlots()

	packing := "FullPacking"
	if LogSlots < params.MaxLogSlots() {
		packing = "SparsePacking"
	}

	params2NLit := params.ParametersLiteral()
	params2NLit.LogN++

	var params2N ckks.Parameters
	var err error
	if params2N, err = ckks.NewParametersFromLiteral(params2NLit); err != nil {
		t.Fatal(err)
	}

	ecd2N := ckks.NewEncoder(params2N)

	t.Run("CoeffsToSlots/"+packing, func(t *testing.T) {

		// This test tests the homomorphic encoding
		// It first generates a vector of complex values of size slots
		//
		// vReal + i*vImag
		//
		// Then encode coefficient-wise and encrypts the vectors :
		//
		// Enc(bitReverse(vReal)||bitReverse(vImg))
		//
		// And applies the homomorphic Encoding (will merge both vectors if there was two)
		//
		// Enc(iFFT(vReal+ i*vImag))
		//
		// And returns the result in one ciphertext if the ciphertext can store it else in two ciphertexts
		//
		// Enc(Ecd(vReal) || Ecd(vImag)) or Enc(Ecd(vReal)) and Enc(Ecd(vImag))
		//
		// Then checks that Dcd(Dec(Enc(Ecd(vReal)))) = vReal and Dcd(Dec(Enc(Ecd(vImag)))) = vImag

		Levels := make([]int, params.MaxLevel())
		for i := range Levels {
			Levels[i] = 1
		}

		CoeffsToSlotsParametersLiteral := HomomorphicDFTMatrixLiteral{
			LogSlots:        LogSlots,
			Type:            Encode,
			RepackImag2Real: true,
			LevelStart:      params.MaxLevel(),
			Levels:          Levels,
		}

		kgen := ckks.NewKeyGenerator(params)
		sk := kgen.GenSecretKeyNew()
		encoder := ckks.NewEncoder(params)
		encryptor := ckks.NewEncryptor(params, sk)
		decryptor := ckks.NewDecryptor(params, sk)

		// Generates the encoding matrices
		CoeffsToSlotMatrices := NewHomomorphicDFTMatrixFromLiteral(CoeffsToSlotsParametersLiteral, encoder)

		// Gets Galois elements
		galEls := CoeffsToSlotsParametersLiteral.GaloisElements(params)

		// Instantiates the EvaluationKeySet
		evk := rlwe.NewEvaluationKeySet()

		// Generates and adds the keys
		for _, galEl := range galEls {
			evk.GaloisKeys[galEl] = kgen.GenGaloisKeyNew(galEl, sk)
		}

		// Also adds the conjugate key
		evk.GaloisKeys[params.GaloisElementForRowRotation()] = kgen.GenGaloisKeyNew(params.GaloisElementForRowRotation(), sk)

		// Creates an evaluator with the rotation keys
		eval := NewEvaluator(params, evk)

		// Generates the vector of random complex values
		values := make([]complex128, slots)
		for i := range values {
			values[i] = complex(sampling.RandFloat64(-1, 1), sampling.RandFloat64(-1, 1))
		}

		// Splits between real and imaginary
		valuesReal := make([]float64, slots)
		for i := range valuesReal {
			valuesReal[i] = real(values[i])
		}

		valuesImag := make([]float64, slots)
		for i := range valuesImag {
			valuesImag[i] = imag(values[i])
		}

		// Applies bit-reverse on the original complex vector
		utils.BitReverseInPlaceSlice(values, params.Slots())

		// Maps to a float vector
		// Add gaps if sparse packing
		valuesFloat := make([]float64, params.N())
		gap := params.N() / (2 * slots)
		for i, jdx, idx := 0, params.N()>>1, 0; i < slots; i, jdx, idx = i+1, jdx+gap, idx+gap {
			valuesFloat[idx] = real(values[i])
			valuesFloat[jdx] = imag(values[i])
		}

		// Encodes coefficient-wise and encrypts the test vector
		pt := ckks.NewPlaintext(params, params.MaxLevel())
		pt.LogSlots = LogSlots

		pt.EncodingDomain = rlwe.CoefficientsDomain
		encoder.Encode(valuesFloat, pt)
		pt.EncodingDomain = rlwe.SlotsDomain

		ct := encryptor.EncryptNew(pt)

		// Applies the homomorphic DFT
		ct0, ct1 := eval.CoeffsToSlotsNew(ct, CoeffsToSlotMatrices)

		// Checks against the original coefficients
		if sparse {

			ct0.EncodingDomain = rlwe.CoefficientsDomain

			coeffsReal := make([]float64, params.N())

			encoder.Decode(decryptor.DecryptNew(ct0), coeffsReal)

			// Plaintext circuit
			vec := make([]complex128, 2*slots)

			// Embed real vector into the complex vector (trivial)
			for i, j := 0, slots; i < slots; i, j = i+1, j+1 {
				vec[i] = complex(valuesReal[i], 0)
				vec[j] = complex(valuesImag[i], 0)
			}

			// IFFT
			encoder.IFFT(vec, LogSlots+1)

			// Extract complex vector into real vector
			vecReal := make([]float64, params.N())
			for i, idx, jdx := 0, 0, params.N()>>1; i < 2*slots; i, idx, jdx = i+1, idx+gap/2, jdx+gap/2 {
				vecReal[idx] = real(vec[i])
				vecReal[jdx] = imag(vec[i])
			}

			// Compares
			verifyTestVectors(params, ecd2N, nil, vecReal, coeffsReal, t)

		} else {

			ct0.EncodingDomain = rlwe.CoefficientsDomain
			ct1.EncodingDomain = rlwe.CoefficientsDomain

			coeffsReal := make([]float64, params.N())
			coeffsImag := make([]float64, params.N())

			encoder.Decode(decryptor.DecryptNew(ct0), coeffsReal)
			encoder.Decode(decryptor.DecryptNew(ct1), coeffsImag)

			vec0 := make([]complex128, slots)
			vec1 := make([]complex128, slots)

			// Embed real vector into the complex vector (trivial)
			for i := 0; i < slots; i++ {
				vec0[i] = complex(valuesReal[i], 0)
				vec1[i] = complex(valuesImag[i], 0)
			}

			// IFFT
			encoder.IFFT(vec0, LogSlots)
			encoder.IFFT(vec1, LogSlots)

			// Extract complex vectors into real vectors
			vecReal := make([]float64, params.N())
			vecImag := make([]float64, params.N())
			for i, j := 0, slots; i < slots; i, j = i+1, j+1 {
				vecReal[i], vecReal[j] = real(vec0[i]), imag(vec0[i])
				vecImag[i], vecImag[j] = real(vec1[i]), imag(vec1[i])
			}

			verifyTestVectors(params, ecd2N, nil, vecReal, coeffsReal, t)
			verifyTestVectors(params, ecd2N, nil, vecImag, coeffsImag, t)
		}
	})
}

func testSlotsToCoeffs(params ckks.Parameters, LogSlots int, t *testing.T) {

	slots := 1 << LogSlots

	var sparse bool = LogSlots < params.LogN()-1

	packing := "FullPacking"
	if LogSlots < params.LogN()-1 {
		packing = "SparsePacking"
	}

	t.Run("SlotsToCoeffs/"+packing, func(t *testing.T) {

		// This test tests the homomorphic decoding
		// It first generates a complex vector of size 2*slots
		// if 2*slots == N, then two vectors are generated, one for the real part, one for the imaginary part :
		//
		// vReal and vReal (both floating point vectors because the encoding always result in a real vector)
		//
		// Then encode and encrypts the vectors :
		//
		// Enc(Ecd(vReal)) and Enc(Ecd(vImag))
		//
		// And applies the homomorphic decoding (will merge both vectors if there was two)
		//
		// Enc(FFT(Ecd(vReal) + i*Ecd(vImag)))
		//
		// The result should be the decoding of the initial vectors bit-reversed
		//
		// Enc(FFT(Ecd(vReal) + i*Ecd(vImag))) = Enc(BitReverse(Dcd(Ecd(vReal + i*vImag))))
		//
		// The first N/2 slots of the plaintext will be the real part while the last N/2 the imaginary part
		// In case of 2*slots < N, then there is a gap of N/(2*slots) between each values

		Levels := make([]int, params.MaxLevel())
		for i := range Levels {
			Levels[i] = 1
		}

		SlotsToCoeffsParametersLiteral := HomomorphicDFTMatrixLiteral{
			LogSlots:        LogSlots,
			Type:            Decode,
			RepackImag2Real: true,
			LevelStart:      params.MaxLevel(),
			Levels:          Levels,
		}

		kgen := ckks.NewKeyGenerator(params)
		sk := kgen.GenSecretKeyNew()
		encoder := ckks.NewEncoder(params)
		encryptor := ckks.NewEncryptor(params, sk)
		decryptor := ckks.NewDecryptor(params, sk)

		// Generates the encoding matrices
		SlotsToCoeffsMatrix := NewHomomorphicDFTMatrixFromLiteral(SlotsToCoeffsParametersLiteral, encoder)

		// Gets the Galois elements
		galEls := SlotsToCoeffsParametersLiteral.GaloisElements(params)

		// Instantiates the EvaluationKeySet
		evk := rlwe.NewEvaluationKeySet()

		// Generates and adds the keys
		for _, galEl := range galEls {
			evk.GaloisKeys[galEl] = kgen.GenGaloisKeyNew(galEl, sk)
		}

		// Also adds the conjugate key
		evk.GaloisKeys[params.GaloisElementForRowRotation()] = kgen.GenGaloisKeyNew(params.GaloisElementForRowRotation(), sk)

		// Creates an evaluator with the rotation keys
		eval := NewEvaluator(params, evk)

		// Generates the n first slots of the test vector (real part to encode)
		valuesReal := make([]complex128, slots)
		for i := range valuesReal {
			valuesReal[i] = complex(sampling.RandFloat64(-1, 1), 0)
		}

		// Generates the n first slots of the test vector (imaginary part to encode)
		valuesImag := make([]complex128, slots)
		for i := range valuesImag {
			valuesImag[i] = complex(sampling.RandFloat64(-1, 1), 0)
		}

		// If sparse, there there is the space to store both vectors in one
		if sparse {
			for i := range valuesReal {
				valuesReal[i] = complex(real(valuesReal[i]), real(valuesImag[i]))
			}
			LogSlots++
		}

		// Encodes and encrypts the test vectors
		plaintext := ckks.NewPlaintext(params, params.MaxLevel())
		plaintext.LogSlots = LogSlots
		encoder.Encode(valuesReal, plaintext)
		ct0 := encryptor.EncryptNew(plaintext)
		var ct1 *rlwe.Ciphertext
		if !sparse {
			encoder.Encode(valuesImag, plaintext)
			ct1 = encryptor.EncryptNew(plaintext)
		}

		// Applies the homomorphic DFT
		res := eval.SlotsToCoeffsNew(ct0, ct1, SlotsToCoeffsMatrix)

		// Decrypt and decode in the coefficient domain
		coeffsFloat := make([]float64, params.N())
		res.EncodingDomain = rlwe.CoefficientsDomain

		encoder.Decode(decryptor.DecryptNew(res), coeffsFloat)

		// Extracts the coefficients and construct the complex vector
		// This is simply coefficient ordering
		valuesTest := make([]complex128, slots)
		gap := params.N() / (2 * slots)
		for i, idx := 0, 0; i < slots; i, idx = i+1, idx+gap {
			valuesTest[i] = complex(coeffsFloat[idx], coeffsFloat[idx+(params.N()>>1)])
		}

		// The result is always returned as a single complex vector, so if full-packing (2 initial vectors)
		// then repacks both vectors together
		if !sparse {
			for i := range valuesReal {
				valuesReal[i] += complex(0, real(valuesImag[i]))
			}
		}

		// Result is bit-reversed, so applies the bit-reverse permutation on the reference vector
		utils.BitReverseInPlaceSlice(valuesReal, params.Slots())

		verifyTestVectors(params, encoder, decryptor, valuesReal, valuesTest, t)
	})
}

func verifyTestVectors(params ckks.Parameters, encoder *ckks.Encoder, decryptor rlwe.Decryptor, valuesWant, element interface{}, t *testing.T) {

	precStats := ckks.GetPrecisionStats(params, encoder, decryptor, valuesWant, element, nil, false)

	if *printPrecisionStats {
		t.Log(precStats.String())
	}

	rf64, _ := precStats.MeanPrecision.Real.Float64()
	if64, _ := precStats.MeanPrecision.Imag.Float64()

	require.GreaterOrEqual(t, rf64, minPrec)
	require.GreaterOrEqual(t, if64, minPrec)
}
