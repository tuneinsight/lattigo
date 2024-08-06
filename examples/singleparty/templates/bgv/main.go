// Package main is a template encrypted modular arithmetic integers, with a set of example parameters, key generation, encoding, encryption, decryption and decoding.
package main

import (
	"fmt"
	"math/rand"
	"slices"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/bgv"
)

func main() {
	var err error
	var params bgv.Parameters

	// 128-bit secure parameters enabling depth-7 circuits.
	// LogN:14, LogQP: 431.
	if params, err = bgv.NewParametersFromLiteral(
		bgv.ParametersLiteral{
			LogN:             14,                                    // log2(ring degree)
			LogQ:             []int{55, 45, 45, 45, 45, 45, 45, 45}, // log2(primes Q) (ciphertext modulus)
			LogP:             []int{61},                             // log2(primes P) (auxiliary modulus)
			PlaintextModulus: 0x10001,                               // log2(scale)
		}); err != nil {
		panic(err)
	}

	// Key Generator
	kgen := rlwe.NewKeyGenerator(params)

	// Secret Key
	sk := kgen.GenSecretKeyNew()

	// Encoder
	ecd := bgv.NewEncoder(params)

	// Encryptor
	enc := rlwe.NewEncryptor(params, sk)

	// Decryptor
	dec := rlwe.NewDecryptor(params, sk)

	// Vector of plaintext values
	values := make([]uint64, params.MaxSlots())

	// Source for sampling random plaintext values (not cryptographically secure)
	/* #nosec G404 */
	r := rand.New(rand.NewSource(0))

	// Populates the vector of plaintext values
	T := params.PlaintextModulus()
	for i := range values {
		values[i] = r.Uint64() % T
	}

	// Allocates a plaintext at the max level.
	// Default rlwe.MetaData:
	// - IsBatched = true (slots encoding)
	// - Scale = params.DefaultScale()
	pt := bgv.NewPlaintext(params, params.MaxLevel())

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
	want := make([]uint64, params.MaxSlots())
	copy(want, values)

	PrintPrecisionStats(params, ct, want, ecd, dec)
}

// PrintPrecisionStats decrypts, decodes and prints the precision stats of a ciphertext.
func PrintPrecisionStats(params bgv.Parameters, ct *rlwe.Ciphertext, want []uint64, ecd *bgv.Encoder, dec *rlwe.Decryptor) {

	var err error

	// Decrypts the vector of plaintext values
	pt := dec.DecryptNew(ct)

	// Decodes the plaintext
	have := make([]uint64, params.MaxSlots())
	if err = ecd.Decode(pt, have); err != nil {
		panic(err)
	}

	// Pretty prints some values
	fmt.Printf("Have: ")
	for i := 0; i < 4; i++ {
		fmt.Printf("%d ", have[i])
	}
	fmt.Printf("...\n")

	fmt.Printf("Want: ")
	for i := 0; i < 4; i++ {
		fmt.Printf("%d ", want[i])
	}
	fmt.Printf("...\n")

	if !slices.Equal(want, have) {
		panic("wrong result: bad decryption or encrypted/plaintext circuits do not match")
	}
}
