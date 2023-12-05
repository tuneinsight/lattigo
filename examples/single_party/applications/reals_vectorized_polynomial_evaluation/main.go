// Package main implements an example of vectorized polynomial evaluation.
package main

import (
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
)

func main() {
	var err error
	var params hefloat.Parameters

	// 128-bit secure parameters enabling depth-7 circuits.
	// LogN:14, LogQP: 431.
	if params, err = hefloat.NewParametersFromLiteral(
		hefloat.ParametersLiteral{
			LogN:            14,                                    // log2(ring degree)
			LogQ:            []int{55, 45, 45, 45, 45, 45, 45, 45}, // log2(primes Q) (ciphertext modulus)
			LogP:            []int{61},                             // log2(primes P) (auxiliary modulus)
			LogDefaultScale: 45,                                    // log2(scale)
			RingType:        ring.ConjugateInvariant,
		}); err != nil {
		panic(err)
	}

	// Key Generator
	kgen := rlwe.NewKeyGenerator(params)

	// Secret Key
	sk := kgen.GenSecretKeyNew()

	// Encoder
	ecd := hefloat.NewEncoder(params)

	// Encryptor
	enc := rlwe.NewEncryptor(params, sk)

	// Decryptor
	dec := rlwe.NewDecryptor(params, sk)

	// Relinearization Key
	rlk := kgen.GenRelinearizationKeyNew(sk)

	// Evaluation Key Set with the Relinearization Key
	evk := rlwe.NewMemEvaluationKeySet(rlk)

	// Evaluator
	eval := hefloat.NewEvaluator(params, evk)

	// Samples values in [-K, K]
	K := 25.0

	// Allocates a plaintext at the max level.
	pt := hefloat.NewPlaintext(params, params.MaxLevel())

	// Vector of plaintext values
	values := make([]float64, pt.Slots())

	// Populates the vector of plaintext values
	for i := range values {
		values[i] = sampling.RandFloat64(-K, K)
	}

	// Encodes the vector of plaintext values
	if err = ecd.Encode(values, pt); err != nil {
		panic(err)
	}

	// Encrypts the vector of plaintext values
	var ct *rlwe.Ciphertext
	if ct, err = enc.EncryptNew(pt); err != nil {
		panic(err)
	}

	// f(x)
	sigmoid := func(x float64) (y float64) {
		return 1 / (math.Exp(-x) + 1)
	}

	// g0(x) = f'(x) * (f(x)-0)
	g0 := func(x float64) (y float64) {
		y = sigmoid(x)
		return y * (1 - y) * (y - 0)
	}

	// g1(x) = f'(x) * (f(x)-1)
	g1 := func(x float64) (y float64) {
		y = sigmoid(x)
		return y * (1 - y) * (y - 1)
	}

	// Defines on which slots g0(x) and g1(x) have to be evaluated
	even := make([]int, ct.Slots()>>1) // List of all even slots
	odd := make([]int, ct.Slots()>>1)  // List of all odd slots
	for i := 0; i < ct.Slots()>>1; i++ {
		even[i] = 2 * i
		odd[i] = 2*i + 1
	}

	mapping := map[int][]int{
		0: even, // g0(x) is evaluated on all even slots
		1: odd,  // g1(x) is evaluated on all odd slots
	}

	// Vectorized Chebyhsev approximation of g0(x) and g1(x) in the domain [-K, K] of degree 63.
	var polys hefloat.PolynomialVector
	if polys, err = hefloat.NewPolynomialVector([]bignum.Polynomial{
		GetChebyshevPoly(K, 63, g0),
		GetChebyshevPoly(K, 63, g1),
	}, mapping); err != nil {
		panic(err)
	}

	// Instantiates the polynomial evaluator
	polyEval := hefloat.NewPolynomialEvaluator(params, eval)

	// Retrieves the vectorized change of basis y = scalar * x + constant
	scalar, constant := polys.ChangeOfBasis(ct.Slots())

	// Performes the vectorized change of basis Standard -> Chebyshev
	if err := eval.Mul(ct, scalar, ct); err != nil {
		panic(err)
	}

	if err := eval.Add(ct, constant, ct); err != nil {
		panic(err)
	}

	if err := eval.Rescale(ct, ct); err != nil {
		panic(err)
	}

	// Evaluates the vectorized polynomial
	if ct, err = polyEval.Evaluate(ct, polys, params.DefaultScale()); err != nil {
		panic(err)
	}

	// Allocates a vector for the reference values
	want := make([]float64, ct.Slots())
	for i := 0; i < ct.Slots()>>1; i++ {
		want[2*i+0], _ = polys.Value[0].Evaluate(values[2*i+0])[0].Float64()
		want[2*i+1], _ = polys.Value[1].Evaluate(values[2*i+1])[0].Float64()
		//want[2*i+0] = sigmoidDerivLabel0(values[2*i+0])
		//want[2*i+1] = sigmoidDerivLabel1(values[2*i+1])
	}

	// Decrypts and print the stats about the precision.
	PrintPrecisionStats(params, ct, want, ecd, dec)
}

// GetChebyshevPoly returns the Chebyshev polynomial approximation of f the
// in the interval [-K, K] for the given degree.
func GetChebyshevPoly(K float64, degree int, f64 func(x float64) (y float64)) bignum.Polynomial {

	FBig := func(x *big.Float) (y *big.Float) {
		xF64, _ := x.Float64()
		return new(big.Float).SetPrec(x.Prec()).SetFloat64(f64(xF64))
	}

	var prec uint = 128

	interval := bignum.Interval{
		A:     *bignum.NewFloat(-K, prec),
		B:     *bignum.NewFloat(K, prec),
		Nodes: degree,
	}

	// Returns the polynomial.
	return bignum.ChebyshevApproximation(FBig, interval)
}

// PrintPrecisionStats decrypts, decodes and prints the precision stats of a ciphertext.
func PrintPrecisionStats(params hefloat.Parameters, ct *rlwe.Ciphertext, want []float64, ecd *hefloat.Encoder, dec *rlwe.Decryptor) {

	var err error

	// Decrypts the vector of plaintext values
	pt := dec.DecryptNew(ct)

	// Decodes the plaintext
	have := make([]float64, ct.Slots())
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
	fmt.Println(hefloat.GetPrecisionStats(params, ecd, dec, have, want, 0, false).String())
}
