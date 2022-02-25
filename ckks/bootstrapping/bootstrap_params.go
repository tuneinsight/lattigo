package bootstrapping

import (
	"github.com/tuneinsight/lattigo/v3/ckks/advanced"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// Parameters is a struct for the default bootstrapping parameters
type Parameters struct {
	SlotsToCoeffsParameters advanced.EncodingMatrixLiteral
	EvalModParameters       advanced.EvalModLiteral
	CoeffsToSlotsParameters advanced.EncodingMatrixLiteral
	EphemeralSecretWeight   int // Hamming weight of the ephemeral secret. If 0, no ephemeral secret is used during the bootstrapping.
}

// MarshalBinary encode the target Parameters on a slice of bytes.
func (p *Parameters) MarshalBinary() (data []byte, err error) {
	data = []byte{}
	tmp := []byte{}

	if tmp, err = p.SlotsToCoeffsParameters.MarshalBinary(); err != nil {
		return nil, err
	}

	data = append(data, uint8(len(tmp)))
	data = append(data, tmp...)

	if tmp, err = p.EvalModParameters.MarshalBinary(); err != nil {
		return nil, err
	}

	data = append(data, uint8(len(tmp)))
	data = append(data, tmp...)

	if tmp, err = p.CoeffsToSlotsParameters.MarshalBinary(); err != nil {
		return nil, err
	}

	data = append(data, uint8(len(tmp)))
	data = append(data, tmp...)

	tmp = make([]byte, 4)
	tmp[0] = uint8(p.EphemeralSecretWeight >> 24)
	tmp[1] = uint8(p.EphemeralSecretWeight >> 16)
	tmp[2] = uint8(p.EphemeralSecretWeight >> 8)
	tmp[3] = uint8(p.EphemeralSecretWeight >> 0)
	data = append(data, tmp...)
	return
}

// UnmarshalBinary decodes a slice of bytes on the target Parameters.
func (p *Parameters) UnmarshalBinary(data []byte) (err error) {

	pt := 0
	dLen := int(data[pt])

	if err := p.SlotsToCoeffsParameters.UnmarshalBinary(data[pt+1 : pt+dLen+1]); err != nil {
		return err
	}

	pt += dLen
	pt++
	dLen = int(data[pt])

	if err := p.EvalModParameters.UnmarshalBinary(data[pt+1 : pt+dLen+1]); err != nil {
		return err
	}

	pt += dLen
	pt++
	dLen = int(data[pt])

	if err := p.CoeffsToSlotsParameters.UnmarshalBinary(data[pt+1 : pt+dLen+1]); err != nil {
		return err
	}

	pt += dLen
	pt++
	dLen = int(data[pt])

	p.EphemeralSecretWeight = int(data[pt])<<24 | int(data[pt+1])<<16 | int(data[pt+2])<<8 | int(data[pt+3])

	return
}

// RotationsForBootstrapping returns the list of rotations performed during the Bootstrapping operation.
func (p *Parameters) RotationsForBootstrapping(params ckks.Parameters) (rotations []int) {

	logN := params.LogN()
	logSlots := params.LogSlots()

	// List of the rotation key values to needed for the bootstrapp
	rotations = []int{}

	//SubSum rotation needed X -> Y^slots rotations
	for i := logSlots; i < logN-1; i++ {
		if !utils.IsInSliceInt(1<<i, rotations) {
			rotations = append(rotations, 1<<i)
		}
	}

	p.CoeffsToSlotsParameters.LogN = logN
	p.SlotsToCoeffsParameters.LogN = logN

	p.CoeffsToSlotsParameters.LogSlots = logSlots
	p.SlotsToCoeffsParameters.LogSlots = logSlots

	rotations = append(rotations, p.CoeffsToSlotsParameters.Rotations()...)
	rotations = append(rotations, p.SlotsToCoeffsParameters.Rotations()...)

	return
}
<<<<<<< btp_eprint
=======
<<<<<<< btp_eprint
=======

// DefaultCKKSParameters are default parameters for the bootstrapping.
// To be used in conjonction with DefaultParameters.
var DefaultCKKSParameters = []ckks.ParametersLiteral{
	{
		LogN:         16,
		LogSlots:     15,
		DefaultScale: 1 << 40,
		H:            192,
		Sigma:        rlwe.DefaultSigma,
		Q: []uint64{
			0x10000000006e0001, // 60 Q0
			0x10000140001,      // 40
			0xffffe80001,       // 40
			0xffffc40001,       // 40
			0x100003e0001,      // 40
			0xffffb20001,       // 40
			0x10000500001,      // 40
			0xffff940001,       // 40
			0xffff8a0001,       // 40
			0xffff820001,       // 40
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
			0x100000000060001,  // 58 CtS
			0xfffffffff00001,   // 58 CtS
			0xffffffffd80001,   // 58 CtS
			0x1000000002a0001,  // 58 CtS
		},
		P: []uint64{
			0x1fffffffffe00001, // Pi 61
			0x1fffffffffc80001, // Pi 61
			0x1fffffffffb40001, // Pi 61
			0x1fffffffff500001, // Pi 61
			0x1fffffffff420001, // Pi 61
		},
	},
	{
		LogN:         16,
		LogSlots:     15,
		DefaultScale: 1 << 45,
		H:            192,
		Sigma:        rlwe.DefaultSigma,
		Q: []uint64{
			0x10000000006e0001, // 60 Q0
			0x2000000a0001,     // 45
			0x2000000e0001,     // 45
			0x1fffffc20001,     // 45
			0x200000440001,     // 45
			0x200000500001,     // 45
			0x3ffffe80001,      //42 StC
			0x3ffffd20001,      //42 StC
			0x3ffffca0001,      //42 StC
			0xffffffffffc0001,  // ArcSine
			0xfffffffff240001,  // ArcSine
			0x1000000000f00001, // ArcSine
			0xfffffffff840001,  // Double angle
			0x1000000000860001, // Double angle
			0xfffffffff6a0001,  // Sine
			0x1000000000980001, // Sine
			0xfffffffff5a0001,  // Sine
			0x1000000000b00001, // Sine
			0x1000000000ce0001, // Sine
			0xfffffffff2a0001,  // Sine
			0x400000000360001,  // 58 CtS
			0x3ffffffffbe0001,  // 58 CtS
			0x400000000660001,  // 58 CtS
			0x4000000008a0001,  // 58 CtS
		},
		P: []uint64{
			0x1fffffffffe00001, // Pi 61
			0x1fffffffffc80001, // Pi 61
			0x1fffffffffb40001, // Pi 61
			0x1fffffffff500001, // Pi 61
		},
	},
	{
		LogN:         16,
		LogSlots:     15,
		DefaultScale: 1 << 30,
		H:            192,
		Sigma:        rlwe.DefaultSigma,
		Q: []uint64{
			0x80000000080001,   // 55 Q0
			0xffffffffffc0001,  // 60
			0x10000000006e0001, // 60
			0xfffffffff840001,  // 60
			0x1000000000860001, // 60
			0xfffffffff6a0001,  // 60
			0x1000000000980001, // 60
			0xfffffffff5a0001,  // 60
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
		},
		P: []uint64{
			0x1fffffffffe00001, // Pi 61
			0x1fffffffffc80001, // Pi 61
			0x1fffffffffb40001, // Pi 61
			0x1fffffffff500001, // Pi 61
			0x1fffffffff420001, // Pi 61
		},
	},
	{
		LogN:         16,
		LogSlots:     15,
		DefaultScale: 1 << 40,
		H:            32768,
		Sigma:        rlwe.DefaultSigma,
		Q: []uint64{
			0x4000000120001, // 60 Q0
			0x10000140001,
			0xffffe80001,
			0xffffc40001,
			0x100003e0001,
			0xffffb20001,
			0x10000500001,
			0xffff940001,
			0xffff8a0001,
			0xffff820001,
			0x100000000060001,  // 56 StC (28 + 28)
			0xffa0001,          // 28 StC
			0xffffffffffc0001,  // 60 Sine (double angle)
			0x10000000006e0001, // 60 Sine (double angle)
			0xfffffffff840001,  // 60 Sine (double angle)
			0x1000000000860001, // 60 Sine (double angle)
			0xfffffffff6a0001,  // 60 Sine
			0x1000000000980001, // 60 Sine
			0xfffffffff5a0001,  // 60 Sine
			0x1000000000b00001, // 60 Sine
			0x1000000000ce0001, // 60 Sine
			0xfffffffff2a0001,  // 60 Sine
			0xfffffffff240001,  // 60 Sine
			0x1000000000f00001, // 60 Sine
			0x200000000e0001,   // 53 CtS
			0x20000000140001,   // 53 CtS
			0x20000000280001,   // 53 CtS
			0x1fffffffd80001,   // 53 CtS
		},
		P: []uint64{
			0x1fffffffffe00001, // Pi 61
			0x1fffffffffc80001, // Pi 61
			0x1fffffffffb40001, // Pi 61
			0x1fffffffff500001, // Pi 61
			0x1fffffffff420001, // Pi 61
			0x1fffffffff380001, // Pi 61
		},
	},
	{
		LogN:         15,
		LogSlots:     14,
		DefaultScale: 1 << 25,
		H:            192,
		Sigma:        rlwe.DefaultSigma,
		Q: []uint64{
			0x1fff90001,       // 32 Q0
			0x4000000420001,   // 50
			0x1fc0001,         // 25
			0xffffffffffc0001, // 60 StC (30+30)
			0x4000000120001,   // 50 Sine
			0x40000001b0001,   // 50 Sine
			0x3ffffffdf0001,   // 50 Sine
			0x4000000270001,   // 50 Sine
			0x3ffffffd20001,   // 50 Sine
			0x3ffffffcd0001,   // 50 Sine
			0x4000000350001,   // 50 Sine
			0x3ffffffc70001,   // 50 Sine
			0x1fffffff50001,   // 49 CtS
			0x1ffffffea0001,   // 49 CtS
		},
		P: []uint64{
			0x7fffffffe0001, // 51
			0x8000000110001, // 51
		},
	},
}

// DefaultParameters are default bootstrapping params for the bootstrapping.
var DefaultParameters = []Parameters{

	// SET I
	// 1546
	{
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
			K:             25,
			SineDeg:       63,
			DoubleAngle:   2,
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
	},

	// SET II
	// 1547
	{
		SlotsToCoeffsParameters: advanced.EncodingMatrixLiteral{
			LinearTransformType: advanced.SlotsToCoeffs,
			RepackImag2Real:     true,
			LevelStart:          8,
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
			LevelStart:    19,
			SineType:      advanced.Cos1,
			MessageRatio:  4.0,
			K:             25,
			SineDeg:       63,
			DoubleAngle:   2,
			ArcSineDeg:    7,
			ScalingFactor: 1 << 60,
		},
		CoeffsToSlotsParameters: advanced.EncodingMatrixLiteral{
			LinearTransformType: advanced.CoeffsToSlots,
			RepackImag2Real:     true,
			LevelStart:          23,
			BSGSRatio:           2.0,
			BitReversed:         false,
			ScalingFactor: [][]float64{
				{0x400000000360001},
				{0x3ffffffffbe0001},
				{0x400000000660001},
				{0x4000000008a0001},
			},
		},
	},

	// SET III
	// 1553
	{
		SlotsToCoeffsParameters: advanced.EncodingMatrixLiteral{
			LinearTransformType: advanced.SlotsToCoeffs,
			RepackImag2Real:     true,
			LevelStart:          9,
			BSGSRatio:           2.0,
			BitReversed:         false,
			ScalingFactor: [][]float64{
				{1073741824.0},
				{1073741824.0062866, 1073741824.0062866},
			},
		},
		EvalModParameters: advanced.EvalModLiteral{
			Q:             0x80000000080001,
			LevelStart:    17,
			SineType:      advanced.Cos1,
			MessageRatio:  256.0,
			K:             25,
			SineDeg:       63,
			DoubleAngle:   2,
			ArcSineDeg:    0,
			ScalingFactor: 1 << 55,
		},
		CoeffsToSlotsParameters: advanced.EncodingMatrixLiteral{
			LinearTransformType: advanced.CoeffsToSlots,
			RepackImag2Real:     true,
			LevelStart:          21,
			BSGSRatio:           2.0,
			BitReversed:         false,
			ScalingFactor: [][]float64{
				{0x200000000e0001},
				{0x20000000140001},
				{0x20000000280001},
				{0x1fffffffd80001},
			},
		},
	},

	// Set IV
	// 1792
	{
		SlotsToCoeffsParameters: advanced.EncodingMatrixLiteral{
			LinearTransformType: advanced.SlotsToCoeffs,
			RepackImag2Real:     true,
			LevelStart:          11,
			BSGSRatio:           2.0,
			BitReversed:         false,
			ScalingFactor: [][]float64{
				{268435456.0007324, 268435456.0007324},
				{0xffa0001},
			},
		},
		EvalModParameters: advanced.EvalModLiteral{
			Q:             0x4000000120001,
			LevelStart:    23,
			SineType:      advanced.Cos2,
			MessageRatio:  256.0,
			K:             325,
			SineDeg:       255,
			DoubleAngle:   4,
			ArcSineDeg:    0,
			ScalingFactor: 1 << 60,
		},
		CoeffsToSlotsParameters: advanced.EncodingMatrixLiteral{
			LinearTransformType: advanced.CoeffsToSlots,
			RepackImag2Real:     true,
			LevelStart:          27,
			BSGSRatio:           2.0,
			BitReversed:         false,
			ScalingFactor: [][]float64{
				{0x200000000e0001},
				{0x20000000140001},
				{0x20000000280001},
				{0x1fffffffd80001},
			},
		},
	},

	// Set V
	// 768
	{
		SlotsToCoeffsParameters: advanced.EncodingMatrixLiteral{
			LinearTransformType: advanced.SlotsToCoeffs,
			RepackImag2Real:     true,
			LevelStart:          3,
			BSGSRatio:           2.0,
			BitReversed:         false,
			ScalingFactor: [][]float64{
				{1073741823.9998779, 1073741823.9998779},
			},
		},
		EvalModParameters: advanced.EvalModLiteral{
			Q:             0x1fff90001,
			LevelStart:    11,
			SineType:      advanced.Cos1,
			MessageRatio:  256.0,
			K:             25,
			SineDeg:       63,
			DoubleAngle:   2,
			ArcSineDeg:    0,
			ScalingFactor: 1 << 50,
		},
		CoeffsToSlotsParameters: advanced.EncodingMatrixLiteral{
			LinearTransformType: advanced.CoeffsToSlots,
			RepackImag2Real:     true,
			LevelStart:          13,
			BSGSRatio:           2.0,
			BitReversed:         false,
			ScalingFactor: [][]float64{
				{0x1fffffff50001},
				{0x1ffffffea0001},
			},
		},
	},
}
>>>>>>> [ckks/advanced]: better StC & CtS
>>>>>>> [ckks/advanced]: better StC & CtS
