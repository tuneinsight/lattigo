// Package main implements an example of Blind Rotation (a.k.a. Lookup Table) evaluation.
// These packages can be used to implement all the functionalities of the TFHE scheme.
package main

import (
	"fmt"
	"time"

	"github.com/tuneinsight/lattigo/v6/core/rgsw/blindrot"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils"
)

// Function to evaluate
func sign(x float64) float64 {
	if x >= 0 {
		return 1
	}

	return -1
}

func main() {
	// RLWE parameters of the Blind Rotation
	// N=1024, Q=0x7fff801 -> ~2^128 ROP-security
	paramsBR, err := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		LogN:    10,
		Q:       []uint64{0x7fff801},
		NTTFlag: true,
	})

	if err != nil {
		panic(err)
	}

	// RLWE parameters of the samples
	// N=512, Q=0x3001 -> ~2^128 ROP-security
	paramsLWE, err := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		LogN:    9,
		Q:       []uint64{0x3001},
		NTTFlag: true,
	})

	if err != nil {
		panic(err)
	}

	// Set the parameters for the blind rotation keys
	evkParams := rlwe.EvaluationKeyParameters{BaseTwoDecomposition: utils.Pointy(7)}

	// Scale of the RLWE samples
	scaleLWE := float64(paramsLWE.Q()[0]) / 4.0

	// Scale of the test poly
	scaleBR := float64(paramsBR.Q()[0]) / 4.0

	// Number of values samples stored in the RLWE sample
	slots := 32

	// Test poly
	testPoly := blindrot.InitTestPolynomial(sign, rlwe.NewScale(scaleBR), paramsBR.RingQ(), -1, 1)

	// Index map of which test poly to evaluate on which slot
	testPolyMap := make(map[int]*ring.Poly)
	for i := 0; i < slots; i++ {
		testPolyMap[i] = &testPoly
	}

	// RLWE secret for the samples
	skLWE := rlwe.NewKeyGenerator(paramsLWE).GenSecretKeyNew()

	// RLWE encryptor for the samples
	encryptorLWE := rlwe.NewEncryptor(paramsLWE, skLWE)

	// Values to encrypt in the RLWE sample
	values := make([]float64, slots)
	for i := 0; i < slots; i++ {
		values[i] = (-1.0 + float64(2*i)/float64(slots))
	}

	// Encode multiples values in a single RLWE
	ptLWE := rlwe.NewPlaintext(paramsLWE, paramsLWE.MaxLevel())
	for i := range values {
		if values[i] < 0 {
			ptLWE.Value.Coeffs[0][i] = paramsLWE.Q()[0] - uint64(-values[i]*scaleLWE)
		} else {
			ptLWE.Value.Coeffs[0][i] = uint64(values[i] * scaleLWE)
		}
	}

	paramsLWE.RingQ().NTT(ptLWE.Value, ptLWE.Value)

	// Encrypt the multiples values in a single RLWE
	ctLWE := rlwe.NewCiphertext(paramsLWE, 1, paramsLWE.MaxLevel())
	if err = encryptorLWE.Encrypt(ptLWE, ctLWE); err != nil {
		panic(err)
	}

	// Evaluator for the Blind Rotations
	eval := blindrot.NewEvaluator(paramsBR, paramsLWE)

	// Secret of the RGSW ciphertexts encrypting the bits of skLWE
	skBR := rlwe.NewKeyGenerator(paramsBR).GenSecretKeyNew()

	// Collection of RGSW ciphertexts encrypting the bits of skLWE under skBR
	blindeRotateKey := blindrot.GenEvaluationKeyNew(paramsBR, skBR, paramsLWE, skLWE, evkParams)

	// Evaluation of BlindRotate(ctLWE) = testPoly(X) * X^{dec{ctLWE}}
	// Returns one RLWE sample per slot in ctLWE

	now := time.Now()
	ctsBR, err := eval.Evaluate(ctLWE, testPolyMap, blindeRotateKey)
	if err != nil {
		panic(err)
	}
	fmt.Printf("Done: %s (avg/BlindRotation %3.1f [ms])\n", time.Since(now), float64(time.Since(now).Milliseconds())/float64(slots))

	// Decrypts, decodes and compares
	q := paramsBR.Q()[0]
	qHalf := q >> 1
	decryptorBR := rlwe.NewDecryptor(paramsBR, skBR)
	ptBR := rlwe.NewPlaintext(paramsBR, paramsBR.MaxLevel())
	for i := 0; i < slots; i++ {

		decryptorBR.Decrypt(ctsBR[i], ptBR)

		if ptBR.IsNTT {
			paramsBR.RingQ().INTT(ptBR.Value, ptBR.Value)
		}

		c := ptBR.Value.Coeffs[0][0]

		var a float64
		if c >= qHalf {
			a = -float64(q-c) / scaleBR
		} else {
			a = float64(c) / scaleBR
		}

		fmt.Printf("%7.4f - %7.4f\n", a, values[i])
	}
}
