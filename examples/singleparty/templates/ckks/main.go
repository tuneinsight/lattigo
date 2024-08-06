// Package main is a template encrypted arithmetic with floating point values, with a set of example parameters, key generation, encoding, encryption, decryption and decoding.
package main

import (
	"fmt"
	"math/rand"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func main() {
	var err error
	var params ckks.Parameters

	// 128-bit secure parameters enabling depth-7 circuits.
	// LogN:14, LogQP: 431.
	if params, err = ckks.NewParametersFromLiteral(
		ckks.ParametersLiteral{
			LogN:            14,                                    // log2(ring degree)
			LogQ:            []int{55, 45, 45, 45, 45, 45, 45, 45}, // log2(primes Q) (ciphertext modulus)
			LogP:            []int{61},                             // log2(primes P) (auxiliary modulus)
			LogDefaultScale: 45,                                    // log2(scale)
		}); err != nil {
		panic(err)
	}

	// Key Generator
	kgen := rlwe.NewKeyGenerator(params)

	// Secret Key
	sk := kgen.GenSecretKeyNew()

	// Encoder
	ecd := ckks.NewEncoder(params)

	// Encryptor
	enc := rlwe.NewEncryptor(params, sk)

	// Decryptor
	dec := rlwe.NewDecryptor(params, sk)

	// Vector of plaintext values
	values := make([]float64, params.MaxSlots())

	// Source for sampling random plaintext values (not cryptographically secure)
	/* #nosec G404 */
	r := rand.New(rand.NewSource(0))

	// Populates the vector of plaintext values
	for i := range values {
		values[i] = 2*r.Float64() - 1 // uniform in [-1, 1]
	}

	// Allocates a plaintext at the max level.
	// Default rlwe.MetaData:
	// - IsBatched = true (slots encoding)
	// - Scale = params.DefaultScale()
	pt := ckks.NewPlaintext(params, params.MaxLevel())

	// Encodes the vector of plaintext values
	if err = ecd.Encode(values, pt); err != nil {
		panic(err)
	}

	// Encrypts the vector of plaintext values
	var ct *rlwe.Ciphertext
	if ct, err = enc.EncryptNew(pt); err != nil {
		panic(err)
	}

	// Allocates a vector for the reference values
	want := make([]float64, params.MaxSlots())
	copy(want, values)

	PrintPrecisionStats(params, ct, want, ecd, dec)
}

// PrintPrecisionStats decrypts, decodes and prints the precision stats of a ciphertext.
func PrintPrecisionStats(params ckks.Parameters, ct *rlwe.Ciphertext, want []float64, ecd *ckks.Encoder, dec *rlwe.Decryptor) {

	var err error

	// Decrypts the vector of plaintext values
	pt := dec.DecryptNew(ct)

	// Decodes the plaintext
	have := make([]float64, params.MaxSlots())
	if err = ecd.Decode(pt, have); err != nil {
		panic(err)
	}

	// Pretty prints some values
	fmt.Printf("Have: ")
	for i := 0; i < 4; i++ {
		fmt.Printf("%20.15f ", have[i])
	}
	fmt.Printf("...\n")

	fmt.Printf("Want: ")
	for i := 0; i < 4; i++ {
		fmt.Printf("%20.15f ", want[i])
	}
	fmt.Printf("...\n")

	// Pretty prints the precision stats
	fmt.Println(ckks.GetPrecisionStats(params, ecd, dec, have, want, 0, false).String())
}
