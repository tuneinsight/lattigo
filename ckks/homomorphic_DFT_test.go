package ckks

import (
	"math/big"
	"math/rand"
	"runtime"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

func TestHomomorphicDFT(t *testing.T) {
	var err error

	if runtime.GOARCH == "wasm" {
		t.Skip("skipping homomorphic DFT tests for GOARCH=wasm")
	}

	ParametersLiteral := ParametersLiteral{
		LogN:              10,
		LogQ:              []int{60, 45, 45, 45, 45, 45, 45, 45},
		LogP:              []int{61, 61},
		Xs:                ring.Ternary{H: 192},
		LogPlaintextScale: 90,
	}

	testHomomorphicDFTMatrixLiteralMarshalling(t)

	var params Parameters
	if params, err = NewParametersFromLiteral(ParametersLiteral); err != nil {
		t.Fatal(err)
	}

	for _, logSlots := range []int{params.PlaintextLogDimensions().Cols - 1, params.PlaintextLogDimensions().Cols} {
		for _, testSet := range []func(params Parameters, logSlots int, t *testing.T){
			testHomomorphicEncoding,
			testHomomorphicDecoding,
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

func testHomomorphicEncoding(params Parameters, LogSlots int, t *testing.T) {

	slots := 1 << LogSlots

	var sparse bool = LogSlots < params.PlaintextLogDimensions().Cols

	packing := "FullPacking"
	if sparse {
		packing = "SparsePacking"
	}

	var params2N Parameters
	var err error
	if params2N, err = NewParametersFromLiteral(ParametersLiteral{
		LogN:              params.LogN() + 1,
		LogQ:              []int{60},
		LogP:              []int{61},
		LogPlaintextScale: params.LogPlaintextScale(),
	}); err != nil {
		t.Fatal(err)
	}

	ecd2N := NewEncoder(params2N)

	t.Run("Encode/"+packing, func(t *testing.T) {

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
		// And returns the result in one ciphextext if the ciphertext can store it else in two ciphertexts
		//
		// Enc(Ecd(vReal) || Ecd(vImag)) or Enc(Ecd(vReal)) and Enc(Ecd(vImag))
		//
		// Then checks that Dcd(Dec(Enc(Ecd(vReal)))) = vReal and Dcd(Dec(Enc(Ecd(vImag)))) = vImag

		Levels := make([]int, params.MaxDepth())
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

		kgen := NewKeyGenerator(params)
		sk := kgen.GenSecretKeyNew()
		encoder := NewEncoder(params)
		encryptor, err := NewEncryptor(params, sk)
		require.NoError(t, err)
		decryptor, err := NewDecryptor(params, sk)
		require.NoError(t, err)

		// Generates the encoding matrices
		CoeffsToSlotMatrices, err := NewHomomorphicDFTMatrixFromLiteral(CoeffsToSlotsParametersLiteral, encoder)
		require.NoError(t, err)

		// Gets Galois elements
		galEls := append(CoeffsToSlotsParametersLiteral.GaloisElements(params), params.GaloisElementOrderTwoOrthogonalSubgroup())

		// Generates and adds the keys
		gks, err := kgen.GenGaloisKeysNew(galEls, sk)
		require.NoError(t, err)

		// Instantiates the EvaluationKeySet
		evk := rlwe.NewMemEvaluationKeySet(nil, gks...)

		// Creates an evaluator with the rotation keys
		eval := NewEvaluator(params, evk)

		prec := params.PlaintextPrecision()

		// Generates the vector of random complex values
		values := make([]*bignum.Complex, slots)
		r := rand.New(rand.NewSource(0))
		for i := range values {
			values[i] = bignum.NewComplex().SetPrec(prec)
			values[i][0].SetFloat64(2*r.Float64() - 1)
			values[i][1].SetFloat64(2*r.Float64() - 1)
		}

		// Splits between real and imaginary
		valuesReal := make([]*big.Float, slots)
		for i := range valuesReal {
			valuesReal[i] = new(big.Float).Set(values[i][0])
		}

		valuesImag := make([]*big.Float, slots)
		for i := range valuesImag {
			valuesImag[i] = new(big.Float).Set(values[i][1])
		}

		// Applies bit-reverse on the original complex vector
		utils.BitReverseInPlaceSlice(values, slots)

		// Maps to a float vector
		// Add gaps if sparse packing
		valuesFloat := make([]*big.Float, params.N())
		gap := params.N() / (2 * slots)
		for i, jdx, idx := 0, params.N()>>1, 0; i < slots; i, jdx, idx = i+1, jdx+gap, idx+gap {
			valuesFloat[idx] = values[i][0]
			valuesFloat[jdx] = values[i][1]
		}

		// Encodes coefficient-wise and encrypts the test vector
		pt := NewPlaintext(params, params.MaxLevel())
		pt.PlaintextLogDimensions = ring.Dimensions{Rows: 0, Cols: LogSlots}

		pt.EncodingDomain = rlwe.CoeffsDomain
		if err = encoder.Encode(valuesFloat, pt); err != nil {
			t.Fatal(err)
		}
		pt.EncodingDomain = rlwe.SlotsDomain

		ct, err := encryptor.EncryptNew(pt)
		require.NoError(t, err)

		// Applies the homomorphic DFT
		ct0, ct1, err := eval.CoeffsToSlotsNew(ct, CoeffsToSlotMatrices)
		require.NoError(t, err)

		// Checks against the original coefficients
		if sparse {

			ct0.EncodingDomain = rlwe.CoeffsDomain

			have := make([]*big.Float, params.N())

			if err = encoder.Decode(decryptor.DecryptNew(ct0), have); err != nil {
				t.Fatal(err)
			}

			// Plaintext circuit
			vec := make([]*bignum.Complex, 2*slots)

			// Embed real vector into the complex vector (trivial)
			for i, j := 0, slots; i < slots; i, j = i+1, j+1 {
				vec[i] = bignum.NewComplex().SetPrec(prec)
				vec[i][0].Set(valuesReal[i])
				vec[j] = bignum.NewComplex().SetPrec(prec)
				vec[j][0].Set(valuesImag[i])
			}

			// IFFT
			if err = encoder.IFFT(vec, LogSlots+1); err != nil {
				t.Fatal(err)
			}

			// Extract complex vector into real vector
			want := make([]*big.Float, params.N())
			for i, idx, jdx := 0, 0, params.N()>>1; i < 2*slots; i, idx, jdx = i+1, idx+gap/2, jdx+gap/2 {
				want[idx] = vec[i][0]
				want[jdx] = vec[i][1]
			}

			// Compares
			verifyTestVectors(params, ecd2N, nil, want, have, nil, t)

		} else {

			ct0.EncodingDomain = rlwe.CoeffsDomain
			ct1.EncodingDomain = rlwe.CoeffsDomain

			haveReal := make([]*big.Float, params.N())
			if err = encoder.Decode(decryptor.DecryptNew(ct0), haveReal); err != nil {
				t.Fatal(err)
			}

			haveImag := make([]*big.Float, params.N())
			if err = encoder.Decode(decryptor.DecryptNew(ct1), haveImag); err != nil {
				t.Fatal(err)
			}

			vec0 := make([]*bignum.Complex, slots)
			vec1 := make([]*bignum.Complex, slots)

			// Embed real vector into the complex vector (trivial)
			for i := 0; i < slots; i++ {
				vec0[i] = bignum.NewComplex().SetPrec(prec)
				vec0[i][0].Set(valuesReal[i])
				vec1[i] = bignum.NewComplex().SetPrec(prec)
				vec1[i][0].Set(valuesImag[i])
			}

			// IFFT
			if err = encoder.IFFT(vec0, LogSlots); err != nil {
				t.Fatal(err)
			}
			if err = encoder.IFFT(vec1, LogSlots); err != nil {
				t.Fatal(err)
			}

			// Extract complex vectors into real vectors
			wantReal := make([]*big.Float, params.N())
			wantImag := make([]*big.Float, params.N())
			for i, j := 0, slots; i < slots; i, j = i+1, j+1 {
				wantReal[i], wantReal[j] = vec0[i][0], vec0[i][1]
				wantImag[i], wantImag[j] = vec1[i][0], vec1[i][1]
			}

			verifyTestVectors(params, ecd2N, nil, wantReal, haveReal, nil, t)
			verifyTestVectors(params, ecd2N, nil, wantImag, haveImag, nil, t)
		}
	})
}

func testHomomorphicDecoding(params Parameters, LogSlots int, t *testing.T) {

	slots := 1 << LogSlots

	var sparse bool = LogSlots < params.PlaintextLogDimensions().Cols

	packing := "FullPacking"
	if sparse {
		packing = "SparsePacking"
	}

	t.Run("Decode/"+packing, func(t *testing.T) {

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

		Levels := make([]int, params.MaxDepth())
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

		kgen := NewKeyGenerator(params)
		sk := kgen.GenSecretKeyNew()
		encoder := NewEncoder(params)
		encryptor, err := NewEncryptor(params, sk)
		require.NoError(t, err)
		decryptor, err := NewDecryptor(params, sk)
		require.NoError(t, err)

		// Generates the encoding matrices
		SlotsToCoeffsMatrix, err := NewHomomorphicDFTMatrixFromLiteral(SlotsToCoeffsParametersLiteral, encoder)
		require.NoError(t, err)

		// Gets the Galois elements
		galEls := append(SlotsToCoeffsParametersLiteral.GaloisElements(params), params.GaloisElementOrderTwoOrthogonalSubgroup())

		// Generates and adds the keys
		gks, err := kgen.GenGaloisKeysNew(galEls, sk)
		require.NoError(t, err)

		// Instantiates the EvaluationKeySet
		evk := rlwe.NewMemEvaluationKeySet(nil, gks...)

		// Creates an evaluator with the rotation keys
		eval := NewEvaluator(params, evk)

		prec := params.PlaintextPrecision()

		// Generates the n first slots of the test vector (real part to encode)
		valuesReal := make([]*bignum.Complex, slots)
		for i := range valuesReal {
			valuesReal[i] = bignum.NewComplex().SetPrec(prec)
			valuesReal[i][0].SetFloat64(sampling.RandFloat64(-1, 1))
		}

		// Generates the n first slots of the test vector (imaginary part to encode)
		valuesImag := make([]*bignum.Complex, slots)
		for i := range valuesImag {
			valuesImag[i] = bignum.NewComplex().SetPrec(prec)
			valuesImag[i][0].SetFloat64(sampling.RandFloat64(-1, 1))
		}

		// If sparse, there there is the space to store both vectors in one
		if sparse {
			for i := range valuesReal {
				valuesReal[i][1].Add(valuesReal[i][1], valuesImag[i][0])
			}
			LogSlots++
		}

		// Encodes and encrypts the test vectors
		plaintext := NewPlaintext(params, params.MaxLevel())
		plaintext.PlaintextLogDimensions = ring.Dimensions{Rows: 0, Cols: LogSlots}
		if err = encoder.Encode(valuesReal, plaintext); err != nil {
			t.Fatal(err)
		}
		ct0, err := encryptor.EncryptNew(plaintext)
		require.NoError(t, err)

		var ct1 *rlwe.Ciphertext
		if !sparse {
			if err = encoder.Encode(valuesImag, plaintext); err != nil {
				t.Fatal(err)
			}
			var err error
			ct1, err = encryptor.EncryptNew(plaintext)
			require.NoError(t, err)
		}

		// Applies the homomorphic DFT
		res, err := eval.SlotsToCoeffsNew(ct0, ct1, SlotsToCoeffsMatrix)
		require.NoError(t, err)

		// Decrypt and decode in the coefficient domain
		coeffsFloat := make([]*big.Float, params.N())
		res.EncodingDomain = rlwe.CoeffsDomain

		if err = encoder.Decode(decryptor.DecryptNew(res), coeffsFloat); err != nil {
			t.Fatal(err)
		}

		// Extracts the coefficients and construct the complex vector
		// This is simply coefficient ordering
		valuesTest := make([]*bignum.Complex, slots)
		gap := params.N() / (2 * slots)
		for i, idx := 0, 0; i < slots; i, idx = i+1, idx+gap {
			valuesTest[i] = &bignum.Complex{coeffsFloat[idx], coeffsFloat[idx+(params.N()>>1)]}
		}

		// The result is always returned as a single complex vector, so if full-packing (2 initial vectors)
		// then repacks both vectors together
		if !sparse {
			for i := range valuesReal {
				valuesReal[i][1].Add(valuesReal[i][1], valuesImag[i][0])
			}
		}

		// Result is bit-reversed, so applies the bit-reverse permutation on the reference vector
		utils.BitReverseInPlaceSlice(valuesReal, slots)

		verifyTestVectors(params, encoder, decryptor, valuesReal, valuesTest, nil, t)
	})
}
