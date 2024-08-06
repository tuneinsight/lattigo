// Package main implements an example of smooth function approximation using minimax polynomial interpolation.
package main

import (
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
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
			RingType:        ring.ConjugateInvariant,
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

	// Relinearization Key
	rlk := kgen.GenRelinearizationKeyNew(sk)

	// Evaluation Key Set with the Relinearization Key
	evk := rlwe.NewMemEvaluationKeySet(rlk)

	// Evaluator
	eval := ckks.NewEvaluator(params, evk)

	// Samples values in [-K, K]
	K := 25.0

	// Allocates a plaintext at the max level.
	pt := ckks.NewPlaintext(params, params.MaxLevel())

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

	sigmoid := func(x float64) (y float64) {
		return 1 / (math.Exp(-x) + 1)
	}

	// Minimax approximation of the sigmoid in the domain [-K, K] of degree 63.
	poly := polynomial.NewPolynomial(GetMinimaxPoly(K, 63, sigmoid))

	// Instantiates the polynomial evaluator
	polyEval := polynomial.NewEvaluator(params, eval)

	// Retrieves the change of basis y = scalar * x + constant
	scalar, constant := poly.ChangeOfBasis()

	// Performes the change of basis Standard -> Chebyshev
	if err := eval.Mul(ct, scalar, ct); err != nil {
		panic(err)
	}

	if err := eval.Add(ct, constant, ct); err != nil {
		panic(err)
	}

	if err := eval.Rescale(ct, ct); err != nil {
		panic(err)
	}

	// Evaluates the polynomial
	if ct, err = polyEval.Evaluate(ct, poly, params.DefaultScale()); err != nil {
		panic(err)
	}

	// Allocates a vector for the reference values and
	// evaluates the same circuit on the plaintext values
	want := make([]float64, ct.Slots())
	for i := range want {
		want[i], _ = poly.Evaluate(values[i])[0].Float64()
		//want[i] = sigmoid(values[i])
	}

	// Decrypts and print the stats about the precision.
	PrintPrecisionStats(params, ct, want, ecd, dec)
}

// GetMinimaxPoly returns the minimax polynomial approximation of f the
// in the interval [-K, K] for the given degree.
func GetMinimaxPoly(K float64, degree int, f64 func(x float64) (y float64)) bignum.Polynomial {

	FBig := func(x *big.Float) (y *big.Float) {
		xF64, _ := x.Float64()
		return new(big.Float).SetPrec(x.Prec()).SetFloat64(f64(xF64))
	}

	// Bit-precision of the arbitrary precision arithmetic used by the minimax solver
	var prec uint = 160

	// Minimax (Remez) approximation of sigmoid
	r := bignum.NewRemez(bignum.RemezParameters{
		// Function to Approximate
		Function: FBig,

		// Polynomial basis of the approximation
		Basis: bignum.Chebyshev,

		// Approximation in [A, B] of degree Nodes.
		Intervals: []bignum.Interval{
			{
				A:     *bignum.NewFloat(-K, prec),
				B:     *bignum.NewFloat(K, prec),
				Nodes: degree,
			},
		},

		// Bit-precision of the solver
		Prec: prec,

		// Scan step for root finding
		ScanStep: bignum.NewFloat(1/16.0, prec),
		// Optimizes the scan-step for root finding
		OptimalScanStep: true,
	})

	// Max 10 iters, and normalized min/max error of 1e-15
	fmt.Printf("Minimax Approximation of Degree %d\n", degree)
	r.Approximate(10, 1e-15)
	fmt.Println()

	// Shoes the coeffs with 50 decimals of precision
	fmt.Printf("Minimax Chebyshev Coefficients [%f, %f]\n", -K, K)
	r.ShowCoeffs(16)
	fmt.Println()

	// Shows the min and max error with 50 decimals of precision
	fmt.Println("Minimax Error")
	r.ShowError(16)
	fmt.Println()

	// Returns the polynomial.
	return bignum.NewPolynomial(bignum.Chebyshev, r.Coeffs, [2]float64{-K, K})
}

// PrintPrecisionStats decrypts, decodes and prints the precision stats of a ciphertext.
func PrintPrecisionStats(params ckks.Parameters, ct *rlwe.Ciphertext, want []float64, ecd *ckks.Encoder, dec *rlwe.Decryptor) {

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
	fmt.Println(ckks.GetPrecisionStats(params, ecd, dec, have, want, 0, false).String())
}
