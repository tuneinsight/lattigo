package bootstrapping

import (
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ckks/advanced"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

type FrontendParameters struct {
	LogN                   uint64
	LogScale               uint64
	UsableLevels           uint64
	BootstrappingPrecision uint64
}

func IsPrime(value uint64) bool {
	// Completely accurate for value < 2^64, i.e., all uint64's
	return new(big.Int).SetUint64(value).ProbablyPrime(0)
}

func MakeFrontendParameters(logN uint64, logScale uint64, usableLevels uint64, bootstrappingPrecision uint64) FrontendParameters {
	if max_log_scale := uint64(63); logScale > max_log_scale {
		error_string := fmt.Sprintf("logScale larger than %d", max_log_scale)
		panic(error_string)
	}
	fp := FrontendParameters{LogN: logN, LogScale: logScale, UsableLevels: usableLevels, BootstrappingPrecision: bootstrappingPrecision}

	if fp.GetQ0()/(uint64(1)<<10) < uint64(1)<<logScale {
		panic("Ciphertext scale when bootstrapping must be an exact power of two less than or equal to Q0/2^10")
	}
	return fp
}

func Concatenate(slices [][]uint64) []uint64 {
	var result []uint64
	for _, slice := range slices {
		result = append(result, slice...)
	}
	return result
}

// Moduli used for P, during bootstrapping, and q0 cannot be used for usable levels
var reservedModuli = Concatenate([][]uint64{qSuffixBootPrecision19, qSuffixBootPrecision26, qSuffixBootPrecision32, []uint64{q0BootPrecision19, q0BootPrecision26, q0BootPrecision32}, pBootPrecision19, pBootPrecision26, pBootPrecision32})

var qSuffixBootPrecision19 = []uint64{
	0x1000000000b00001, // 60 StC  (30)
	0x1000000000ce0001, // 60 StC  (30+30)
	0x80000000440001,   // 55 Sine (double angle)
	0x7fffffffba0001,   // 55 Sine (double angle)
	0x80000000500001,   // 55 Sine
	0x7fffffffaa0001,   // 55 Sine
	0x800000005e0001,   // 55 Sine
	0x7fffffff7e0001,   // 55 Sine
	0x7fffffff380001,   // 55 Sine
	0x80000000ca0001,   // 55 Sine
	0x200000000e0001,   // 53 CtS
	0x20000000140001,   // 53 CtS
	0x20000000280001,   // 53 CtS
	0x1fffffffd80001,   // 53 CtS
}

var qSuffixBootPrecision26 = []uint64{
	0x7fffe60001,       // 39 StC
	0x7fffe40001,       // 39 StC
	0x7fffe00001,       // 39 StC
	0xfffffffff840001,  // 60 Sine (double angle)
	0x1000000000860001, // 60 Sine (double angle)
	0xfffffffff6a0001,  // 60 Sine
	0x1000000000980001, // 60 Sine
	0xfffffffff5a0001,  // 60 Sine
	0x1000000000b00001, // 60 Sine
	0x1000000000ce0001, // 60 Sine
	0xfffffffff2a0001,  // 60 Sine
	0x100000000060001,  // 56 CtS
	0xfffffffff00001,   // 56 CtS
	0xffffffffd80001,   // 56 CtS
	0x1000000002a0001,  // 56 CtS
}

var qSuffixBootPrecision32 = []uint64{
	0x3ffffe80001,      // 42 StC
	0x3ffffd20001,      // 42 StC
	0x3ffffca0001,      // 42 StC
	0xffffffffffc0001,  // 60 ArcSine
	0xfffffffff240001,  // 60 ArcSine
	0x1000000000f00001, // 60 ArcSine
	0xfffffffff840001,  // 60 Double angle
	0x1000000000860001, // 60 Double angle
	0xfffffffff6a0001,  // 60 Sine
	0x1000000000980001, // 60 Sine
	0xfffffffff5a0001,  // 60 Sine
	0x1000000000b00001, // 60 Sine
	0x1000000000ce0001, // 60 Sine
	0xfffffffff2a0001,  // 60 Sine
	0x400000000360001,  // 58 CtS
	0x3ffffffffbe0001,  // 58 CtS
	0x400000000660001,  // 58 CtS
	0x4000000008a0001,  // 58 CtS
}

var pBootPrecision19 = []uint64{
	0x1fffffffffe00001, // Pi 61
	0x1fffffffffc80001, // Pi 61
	0x1fffffffffb40001, // Pi 61
	0x1fffffffff500001, // Pi 61
	0x1fffffffff420001, // Pi 61
}

var pBootPrecision26 = []uint64{
	0x1fffffffffe00001, // Pi 61
	0x1fffffffffc80001, // Pi 61
	0x1fffffffffb40001, // Pi 61
	0x1fffffffff500001, // Pi 61
	0x1fffffffff420001, // Pi 61
}

var pBootPrecision32 = []uint64{
	0x1fffffffffe00001, // Pi 61
	0x1fffffffffc80001, // Pi 61
	0x1fffffffffb40001, // Pi 61
	0x1fffffffff500001, // Pi 61
}

var parametersBootPrecision19 = Parameters{
	EphemeralSecretWeight: 32,
	SlotsToCoeffsParameters: advanced.EncodingMatrixLiteral{
		LinearTransformType: advanced.SlotsToCoeffs,
		RepackImag2Real:     true,
		LevelStart:          0,
		BSGSRatio:           2.0,
		BitReversed:         false,
		ScalingFactor: [][]float64{
			{1073741824.0},
			{1073741824.0062866, 1073741824.0062866},
		},
	},
	EvalModParameters: advanced.EvalModLiteral{
		Q:             0x80000000080001,
		LevelStart:    0,
		SineType:      advanced.Cos1,
		MessageRatio:  256.0,
		K:             16,
		SineDeg:       30,
		DoubleAngle:   3,
		ArcSineDeg:    0,
		ScalingFactor: 1 << 55,
	},
	CoeffsToSlotsParameters: advanced.EncodingMatrixLiteral{
		LinearTransformType: advanced.CoeffsToSlots,
		RepackImag2Real:     true,
		LevelStart:          0,
		BSGSRatio:           2.0,
		BitReversed:         false,
		ScalingFactor: [][]float64{
			{0x200000000e0001},
			{0x20000000140001},
			{0x20000000280001},
			{0x1fffffffd80001},
		},
	},
}

var parametersBootPrecision26 = Parameters{
	EphemeralSecretWeight: 32,
	SlotsToCoeffsParameters: advanced.EncodingMatrixLiteral{
		LinearTransformType: advanced.SlotsToCoeffs,
		RepackImag2Real:     true,
		LevelStart:          12,
		BSGSRatio:           2.0,
		BitReversed:         false,
		ScalingFactor: [][]float64{
			{0x7fffe60001},
			{0x7fffe40001},
			{0x7fffe00001},
		},
	},
	EvalModParameters: advanced.EvalModLiteral{
		Q:             0x10000000006e0001,
		LevelStart:    20,
		SineType:      advanced.Cos1,
		MessageRatio:  256.0,
		K:             16,
		SineDeg:       30,
		DoubleAngle:   3,
		ArcSineDeg:    0,
		ScalingFactor: 1 << 60,
	},
	CoeffsToSlotsParameters: advanced.EncodingMatrixLiteral{
		LinearTransformType: advanced.CoeffsToSlots,
		RepackImag2Real:     true,
		LevelStart:          24,
		BSGSRatio:           2.0,
		BitReversed:         false,
		ScalingFactor: [][]float64{
			{0x100000000060001},
			{0xfffffffff00001},
			{0xffffffffd80001},
			{0x1000000002a0001},
		},
	},
}

var parametersBootPrecision32 = Parameters{
	EphemeralSecretWeight: 32,
	SlotsToCoeffsParameters: advanced.EncodingMatrixLiteral{
		LinearTransformType: advanced.SlotsToCoeffs,
		RepackImag2Real:     true,
		LevelStart:          0,
		BSGSRatio:           2.0,
		BitReversed:         false,
		ScalingFactor: [][]float64{
			{0x3ffffe80001},
			{0x3ffffd20001},
			{0x3ffffca0001},
		},
	},
	EvalModParameters: advanced.EvalModLiteral{
		Q:             0x10000000006e0001,
		LevelStart:    0,
		SineType:      advanced.Cos1,
		MessageRatio:  4.0,
		K:             16,
		SineDeg:       30,
		DoubleAngle:   3,
		ArcSineDeg:    7,
		ScalingFactor: 1 << 60,
	},
	CoeffsToSlotsParameters: advanced.EncodingMatrixLiteral{
		LinearTransformType: advanced.CoeffsToSlots,
		RepackImag2Real:     true,
		LevelStart:          0,
		BSGSRatio:           2.0,
		BitReversed:         false,
		ScalingFactor: [][]float64{
			{0x400000000360001},
			{0x3ffffffffbe0001},
			{0x400000000660001},
			{0x4000000008a0001},
		},
	},
}

func Contains(slice []uint64, value uint64) bool {
	for _, elem := range slice {
		if elem == value {
			return true
		}
	}
	return false
}

func Divides(slice []uint64, value uint64) bool {
	for _, elem := range slice {
		if elem%value == 0 {
			return true
		}
	}
	return false
}

var q0BootPrecision19 = uint64(0x80000000080001)
var q0BootPrecision26 = uint64(0x10000000006e0001)
var q0BootPrecision32 = uint64(0x10000000006e0001)

func (fp FrontendParameters) getUsableLevelPrimes() []uint64 {
	var primes []uint64
	starting_qh := uint64(1) << (fp.LogScale - fp.LogN - 1) // == scale / (2 * N)
	ntt_friendly_multiple := uint64(1) << (fp.LogN + 1)     // == 2 * N
	high_candidate_qh := starting_qh
	low_candidate_qh := starting_qh
	candidate := uint64(0)
	for len(primes) < int(fp.UsableLevels) {
		if high_candidate_qh-starting_qh > starting_qh-low_candidate_qh {
			candidate = low_candidate_qh*ntt_friendly_multiple + 1
			low_candidate_qh--
		} else {
			candidate = high_candidate_qh*ntt_friendly_multiple + 1
			high_candidate_qh++
		}

		if IsPrime(candidate) && !Divides(reservedModuli, candidate) {
			primes = append(primes, candidate)
		}
	}
	return primes
}

func (fp FrontendParameters) GetSlotsToCoeffsLevelStart() int {
	switch fp.BootstrappingPrecision {
	case 19:
		return 2 + int(fp.UsableLevels)
	case 26:
		return 3 + int(fp.UsableLevels)
	case 32:
		return 3 + int(fp.UsableLevels)
	default:
		panic("")
	}
}

func (fp FrontendParameters) GetEvalModLevelStart() int {
	switch fp.BootstrappingPrecision {
	case 19:
		return 10 + int(fp.UsableLevels)
	case 26:
		return 11 + int(fp.UsableLevels)
	case 32:
		return 14 + int(fp.UsableLevels)
	default:
		panic("")
	}
}

func (fp FrontendParameters) GetCoeffsToSlotsLevelStart() int {
	switch fp.BootstrappingPrecision {
	case 19:
		return 14 + int(fp.UsableLevels)
	case 26:
		return 15 + int(fp.UsableLevels)
	case 32:
		return 18 + int(fp.UsableLevels)
	default:
		panic("")
	}
}

func (fp FrontendParameters) GetQ0() uint64 {
	switch fp.BootstrappingPrecision {
	case 19:
		return q0BootPrecision19
	case 26:
		return q0BootPrecision26
	case 32:
		return q0BootPrecision32
	default:
		panic("")
	}
}

func (fp FrontendParameters) GetSuffix() []uint64 {
	switch fp.BootstrappingPrecision {
	case 19:
		return qSuffixBootPrecision19
	case 26:
		return qSuffixBootPrecision26
	case 32:
		return qSuffixBootPrecision32
	default:
		panic("")
	}
}

func (fp FrontendParameters) GetQ() []uint64 {
	Q := []uint64{fp.GetQ0()}
	Q = append(Q, fp.getUsableLevelPrimes()...)
	Q = append(Q, fp.GetSuffix()...)
	return Q
}

func (fp FrontendParameters) GetP() []uint64 {
	switch fp.BootstrappingPrecision {
	case 19:
		return pBootPrecision19
	case 26:
		return pBootPrecision26
	case 32:
		return pBootPrecision32
	default:
		panic("")
	}
}

func (fp FrontendParameters) GetParameters() Parameters {
	var parameters Parameters
	switch fp.BootstrappingPrecision {
	case 19:
		parameters = parametersBootPrecision19
	case 26:
		parameters = parametersBootPrecision26
	case 32:
		parameters = parametersBootPrecision32
	default:
		panic("")
	}
	parameters.SlotsToCoeffsParameters.LevelStart = fp.GetSlotsToCoeffsLevelStart()
	parameters.EvalModParameters.LevelStart = fp.GetEvalModLevelStart()
	parameters.CoeffsToSlotsParameters.LevelStart = fp.GetCoeffsToSlotsLevelStart()
	return parameters
}

func (fp FrontendParameters) GetDefaultParametersSparse() defaultParametersLiteral {
	Q := fp.GetQ()

	return defaultParametersLiteral{
		ckks.ParametersLiteral{
			LogN:         int(fp.LogN),
			LogSlots:     int(fp.LogN) - 1,
			DefaultScale: math.Pow(2, float64(fp.LogScale)),
			Sigma:        rlwe.DefaultSigma,
			H:            192,
			Q:            Q,
			P:            fp.GetP(),
		},
		fp.GetParameters(),
	}
}
