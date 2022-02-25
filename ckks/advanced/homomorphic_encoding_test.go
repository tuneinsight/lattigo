package advanced

import (
	"flag"
	"runtime"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v3/ckks"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")

var minPrec float64 = 15

func TestHomomorphicEncoding(t *testing.T) {
	var err error

	if runtime.GOARCH == "wasm" {
		t.Skip("skipping homomorphic encoding tests for GOARCH=wasm")
	}

	ParametersLiteral := ckks.ParametersLiteral{
		LogN:         13,
		LogSlots:     12,
		DefaultScale: 1 << 45,
		Sigma:        rlwe.DefaultSigma,
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
	}

	testEncodingMatrixLiteralMarshalling(t)

	var params ckks.Parameters
	if params, err = ckks.NewParametersFromLiteral(ParametersLiteral); err != nil {
		panic(err)
	}

	for _, testSet := range []func(params ckks.Parameters, t *testing.T){
		testCoeffsToSlots,
		testSlotsToCoeffs,
	} {
		testSet(params, t)
		runtime.GC()
	}

	ParametersLiteral.LogSlots--
	if params, err = ckks.NewParametersFromLiteral(ParametersLiteral); err != nil {
		panic(err)
	}

	for _, testSet := range []func(params ckks.Parameters, t *testing.T){
		testCoeffsToSlots,
		testSlotsToCoeffs,
	} {
		testSet(params, t)
		runtime.GC()
	}
}

func testEncodingMatrixLiteralMarshalling(t *testing.T) {
	t.Run("Marshalling", func(t *testing.T) {
		m := EncodingMatrixLiteral{
			LinearTransformType: CoeffsToSlots,
			LevelStart:          12,
			BSGSRatio:           16.0,
			BitReversed:         false,
			ScalingFactor: [][]float64{
				{0x100000000060001},
				{0xfffffffff00001},
				{0xffffffffd80001},
				{0x1000000002a0001},
			},
		}

		data, err := m.MarshalBinary()
		assert.Nil(t, err)

		mNew := new(EncodingMatrixLiteral)
		if err := mNew.UnmarshalBinary(data); err != nil {
			assert.Nil(t, err)
		}
		assert.Equal(t, m, *mNew)
	})
}

func testCoeffsToSlots(params ckks.Parameters, t *testing.T) {

	packing := "FullPacking"
	if params.LogSlots() < params.LogN()-1 {
		packing = "SparsePacking"
	}

	t.Run("CoeffsToSlots/"+packing, func(t *testing.T) {

		// This test tests the homomorphic encoding
		// It first generates a vector of complex values of size params.Slots()
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
		// And returns the result in one ciphextext if the ciphertext can store it else in two ciphertexts
		//
		// Enc(Ecd(vReal) || Ecd(vImag)) or Enc(Ecd(vReal)) and Enc(Ecd(vImag))
		//
		// Then checks that Dcd(Dec(Enc(Ecd(vReal)))) = vReal and Dcd(Dec(Enc(Ecd(vImag)))) = vImag

		CoeffsToSlotsParametersLiteral := EncodingMatrixLiteral{
			LogN:                params.LogN(),
			LogSlots:            params.LogSlots(),
			Scaling:             1.0 / float64(2*params.Slots()),
			LinearTransformType: CoeffsToSlots,
			RepackImag2Real:     true,
			LevelStart:          params.MaxLevel(),
			BSGSRatio:           16.0,
			BitReversed:         false,
			ScalingFactor: [][]float64{
				{params.QiFloat64(params.MaxLevel() - 2)},
				{params.QiFloat64(params.MaxLevel() - 1)},
				{params.QiFloat64(params.MaxLevel() - 0)},
			},
		}

		kgen := ckks.NewKeyGenerator(params)
		sk := kgen.GenSecretKey()
		encoder := ckks.NewEncoder(params)
		encryptor := ckks.NewEncryptor(params, sk)
		decryptor := ckks.NewDecryptor(params, sk)

		// Generates the encoding matrices
		CoeffsToSlotMatrices := NewHomomorphicEncodingMatrixFromLiteral(CoeffsToSlotsParametersLiteral, encoder)

		// Gets the rotations indexes for CoeffsToSlots
		rotations := CoeffsToSlotsParametersLiteral.Rotations()

		// Generates the rotation keys
		rotKey := kgen.GenRotationKeysForRotations(rotations, true, sk)

		// Creates an evaluator with the rotation keys
		eval := NewEvaluator(params, rlwe.EvaluationKey{Rlk: nil, Rtks: rotKey})

		// Generates the vector of random complex values
		values := make([]complex128, params.Slots())
		for i := range values {
			values[i] = complex(utils.RandFloat64(-1, 1), utils.RandFloat64(-1, 1))
		}

		// Splits between real and imaginary
		valuesReal := make([]complex128, params.Slots())
		for i := range valuesReal {
			valuesReal[i] = complex(real(values[i]), 0)
		}

		valuesImag := make([]complex128, params.Slots())
		for i := range valuesImag {
			valuesImag[i] = complex(imag(values[i]), 0)
		}

		// Applies bit-reverse on the original complex vector
		ckks.SliceBitReverseInPlaceComplex128(values, params.Slots())

		// Maps to a float vector
		// Add gaps if sparse packing
		valuesFloat := make([]float64, params.N())
		gap := params.N() / (2 * params.Slots())
		for i, jdx, idx := 0, params.N()>>1, 0; i < params.Slots(); i, jdx, idx = i+1, jdx+gap, idx+gap {
			valuesFloat[idx] = real(values[i])
			valuesFloat[jdx] = imag(values[i])
		}

		// Encodes coefficient-wise and encrypts the test vector
		plaintext := ckks.NewPlaintext(params, params.MaxLevel(), params.DefaultScale())
		encoder.EncodeCoeffs(valuesFloat, plaintext)
		ciphertext := encryptor.EncryptNew(plaintext)

		// Applies the homomorphic DFT
		ct0, ct1 := eval.CoeffsToSlotsNew(ciphertext, CoeffsToSlotMatrices)

		// Checks against the original coefficients
		var coeffsReal, coeffsImag []complex128
		if params.LogSlots() < params.LogN()-1 {
			coeffsRealImag := encoder.DecodePublic(decryptor.DecryptNew(ct0), params.LogSlots()+1, 0)
			coeffsReal = coeffsRealImag[:params.Slots()]
			coeffsImag = coeffsRealImag[params.Slots():]
		} else {
			coeffsReal = encoder.DecodePublic(decryptor.DecryptNew(ct0), params.LogSlots(), 0)
			coeffsImag = encoder.DecodePublic(decryptor.DecryptNew(ct1), params.LogSlots(), 0)
		}

		verifyTestVectors(params, encoder, nil, valuesReal, coeffsReal, params.LogSlots(), 0, t)
		verifyTestVectors(params, encoder, nil, valuesImag, coeffsImag, params.LogSlots(), 0, t)
	})
}

func testSlotsToCoeffs(params ckks.Parameters, t *testing.T) {

	packing := "FullPacking"
	if params.LogSlots() < params.LogN()-1 {
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

		SlotsToCoeffsParametersLiteral := EncodingMatrixLiteral{
			LogN:                params.LogN(),
			LogSlots:            params.LogSlots(),
			Scaling:             1.0,
			LinearTransformType: SlotsToCoeffs,
			RepackImag2Real:     true,
			LevelStart:          params.MaxLevel(),
			BSGSRatio:           16.0,
			BitReversed:         false,
			ScalingFactor: [][]float64{
				{params.QiFloat64(params.MaxLevel() - 2)},
				{params.QiFloat64(params.MaxLevel() - 1)},
				{params.QiFloat64(params.MaxLevel() - 0)},
			},
		}

		kgen := ckks.NewKeyGenerator(params)
		sk := kgen.GenSecretKey()
		encoder := ckks.NewEncoder(params)
		encryptor := ckks.NewEncryptor(params, sk)
		decryptor := ckks.NewDecryptor(params, sk)

		// Generates the encoding matrices
		SlotsToCoeffsMatrix := NewHomomorphicEncodingMatrixFromLiteral(SlotsToCoeffsParametersLiteral, encoder)

		// Gets the rotations indexes for SlotsToCoeffs
		rotations := SlotsToCoeffsParametersLiteral.Rotations()

		// Generates the rotation keys
		rotKey := kgen.GenRotationKeysForRotations(rotations, true, sk)

		// Creates an evaluator with the rotation keys
		eval := NewEvaluator(params, rlwe.EvaluationKey{Rlk: nil, Rtks: rotKey})

		// Generates the n first slots of the test vector (real part to encode)
		valuesReal := make([]complex128, params.Slots())
		for i := range valuesReal {
			valuesReal[i] = complex(float64(i+1)/float64(params.Slots()), 0)
		}

		// Generates the n first slots of the test vector (imaginary part to encode)
		valuesImag := make([]complex128, params.Slots())
		for i := range valuesImag {
			valuesImag[i] = complex(float64(i+1)/float64(params.Slots()), 0)
		}

		// If sparse, there there is the space to store both vectors in one
		if params.LogSlots() < params.LogN()-1 {
			for i := range valuesReal {
				valuesReal[i] += complex(0, real(valuesImag[i]))
			}
		}

		// Encodes and encrypts the test vectors
		logSlots := params.LogSlots()
		if params.LogSlots() < params.LogN()-1 {
			logSlots++
		}

		plaintext := ckks.NewPlaintext(params, params.MaxLevel(), params.DefaultScale())
		encoder.Encode(valuesReal, plaintext, logSlots)
		ct0 := encryptor.EncryptNew(plaintext)
		var ct1 *ckks.Ciphertext
		if params.LogSlots() == params.LogN()-1 {
			encoder.Encode(valuesImag, plaintext, logSlots)
			ct1 = encryptor.EncryptNew(plaintext)
		}

		// Applies the homomorphic DFT
		res := eval.SlotsToCoeffsNew(ct0, ct1, SlotsToCoeffsMatrix)

		// Decrypt and decode in the coefficient domain
		coeffsFloat := encoder.DecodeCoeffsPublic(decryptor.DecryptNew(res), 0)

		// Extracts the coefficients and construct the complex vector
		// This is simply coefficient ordering
		valuesTest := make([]complex128, params.Slots())
		gap := params.N() / (2 * params.Slots())
		for i, idx := 0, 0; i < params.Slots(); i, idx = i+1, idx+gap {
			valuesTest[i] = complex(coeffsFloat[idx], coeffsFloat[idx+(params.N()>>1)])
		}

		// The result is always returned as a single complex vector, so if full-packing (2 initial vectors)
		// then repacks both vectors together
		if params.LogSlots() == params.LogN()-1 {
			for i := range valuesReal {
				valuesReal[i] += complex(0, real(valuesImag[i]))
			}
		}

		// Result is bit-reversed, so applies the bit-reverse permutation on the reference vector
		ckks.SliceBitReverseInPlaceComplex128(valuesReal, params.Slots())

		verifyTestVectors(params, encoder, decryptor, valuesReal, valuesTest, params.LogSlots(), 0, t)
	})
}

func verifyTestVectors(params ckks.Parameters, encoder ckks.Encoder, decryptor ckks.Decryptor, valuesWant []complex128, element interface{}, logSlots int, bound float64, t *testing.T) {

	precStats := ckks.GetPrecisionStats(params, encoder, decryptor, valuesWant, element, logSlots, bound)
	if *printPrecisionStats {
		t.Log(precStats.String())
	}
	require.GreaterOrEqual(t, precStats.MeanPrecision.Real, minPrec)
	require.GreaterOrEqual(t, precStats.MeanPrecision.Imag, minPrec)
}
