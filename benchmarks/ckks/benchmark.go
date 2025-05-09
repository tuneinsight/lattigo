package main

import (
	"fmt"
	"math"
	"time"

	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

type Benchmark struct {
	OpName  string
	Level   int
	AvgTime time.Duration
	StdDev  time.Duration
}

// Run op runtime times and take the average time
func measureOp(runtime int, level int, opName string, op func()) Benchmark {
	var sum time.Duration
	times := make([]time.Duration, runtime)

	// Collect all measurements
	for i := 0; i < runtime; i++ {
		start := time.Now()
		op()
		elapsed := time.Since(start)
		times[i] = elapsed
		sum += elapsed
	}

	avgTime := sum / time.Duration(runtime)

	// Calculate standard deviation
	var sqDiffSum float64
	avgNano := float64(avgTime.Nanoseconds())
	for _, t := range times {
		diff := float64(t.Nanoseconds()) - avgNano
		sqDiffSum += diff * diff
	}
	stdDev := time.Duration(math.Sqrt(sqDiffSum / float64(runtime)))

	return Benchmark{
		OpName:  opName,
		Level:   level,
		AvgTime: avgTime,
		StdDev:  stdDev,
	}
}

// Run all operations for all levels and return the average time for each operation
func benchmarkAllLevels(params ckks.Parameters, btpParams bootstrapping.Parameters) []Benchmark {
	times := 10
	benchmarks := []Benchmark{}

	kgen := ckks.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	pk := kgen.GenPublicKeyNew(sk)
	rlk := kgen.GenRelinearizationKeyNew(sk)
	evk := rlwe.NewMemEvaluationKeySet(rlk)

	enc := rlwe.NewEncryptor(params, pk)
	ecd := ckks.NewEncoder(params)

	eval := ckks.NewEvaluator(params, evk)

	// Rotation by 5 positions to the left
	galEl := []uint64{
		params.GaloisElement(5),
		params.GaloisElementForComplexConjugation(),
	}
	galKeys := kgen.GenGaloisKeysNew(galEl, sk)
	rotEval := ckks.NewEvaluator(params, rlwe.NewMemEvaluationKeySet(rlk, galKeys...))

	btpEvk, _, _ := btpParams.GenEvaluationKeys(sk)
	btpEval, _ := bootstrapping.NewEvaluator(btpParams, btpEvk)

	values := make([]complex128, params.MaxSlots())
	for i := range values {
		values[i] = sampling.RandComplex128(-1, 1)
	}

	for level := params.MaxLevel(); level > 0; level-- {
		pt1 := ckks.NewPlaintext(params, level)
		pt2 := ckks.NewPlaintext(params, level)
		ecd.Encode(values, pt1)
		ecd.Encode(values, pt2)

		ct1, _ := enc.EncryptNew(pt1)
		ct2, _ := enc.EncryptNew(pt2)
		ctout := ckks.NewCiphertext(params, 1, level)

		benchmarks = append(benchmarks, measureOp(times, level,
			"PlainAdd",
			func() { eval.AddNew(ct1, pt2) }))

		benchmarks = append(benchmarks, measureOp(times, level,
			"CipherAdd",
			func() { eval.AddNew(ct1, ct2) }))

		benchmarks = append(benchmarks, measureOp(times, level,
			"PlainMult",
			func() { eval.MulRelinNew(ct1, pt2) }))

		benchmarks = append(benchmarks, measureOp(times, level,
			"CipherMult",
			func() {
				res, err := eval.MulRelinNew(ct1, ct2)
				if err != nil {
					panic(err)
				}
				eval.Rescale(res, res)
			}))

		benchmarks = append(benchmarks, measureOp(times, level,
			"Rotate",
			func() {
				_, err := rotEval.RotateNew(ct1, 5)
				if err != nil {
					panic(err)
				}
			}))

		benchmarks = append(benchmarks, measureOp(2*times, level,
			"Rescale",
			func() {
				err := eval.Rescale(ct1, ctout)
				if err != nil {
					panic(err)
				}
			}))

		benchmarks = append(benchmarks, measureOp(times, level,
			"Bootstrap",
			func() {
				_, err := btpEval.Bootstrap(ct1)
				if err != nil {
					panic(err)
				}
			}))

		benchmarks = append(benchmarks, measureOp(times, level,
			"ModSwitch",
			func() {
				eval.DropLevel(ct1, 1)
			}))

		fmt.Println(".")
	}

	return benchmarks
}

func main() {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            16,                                                // Log2 of the ring degree
		LogQ:            []int{55, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40}, // Log2 of the ciphertext prime moduli
		LogP:            []int{61, 61, 61},                                 // Log2 of the key-switch auxiliary prime moduli
		LogDefaultScale: 40,                                                // Log2 of the scale
		Xs:              ring.Ternary{H: 192},
	})
	if err != nil {
		panic(err)
	}

	btpParams, err := bootstrapping.NewParametersFromLiteral(params, bootstrapping.ParametersLiteral{
		LogN: utils.Pointy(16),
		LogP: []int{61, 61, 61, 61},
	})
	if err != nil {
		panic(err)
	}

	// Run benchmarks
	results := benchmarkAllLevels(params, btpParams)
	for _, result := range results {
		fmt.Printf("Level %d\n", result.Level)
		fmt.Printf("%s:\n", result.OpName)
		fmt.Printf("  Average Time: %v\n", result.AvgTime)
		fmt.Printf("  Std Dev: %v\n", result.StdDev)
	}
}
