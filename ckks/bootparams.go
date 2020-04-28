package ckks

import (
	"math"
)

type BootParams struct {
	Parameters

	SinType   SinType // Chose betwenn [Sin(2*pi*x)] or [cos(2*pi*x/r) with double angle formula]
	SinRange  uint64  // K parameter (interpolation in the range -K to K)
	SinDeg    uint64  // Degree of the interpolation
	SinRescal uint64  // Number of rescale and double angle formula (only applies for cos)
	BabySplit uint64  // L parameter of the Baby-step giant-step algorithm (the smallest the more precision but the more non-scalr multiplications)

	CtSDepth uint64 // Depth of the Coeffs To Slots
	StCDepth uint64 // Depth of the Slots To Coeffs

	CtSRescale bool // Rescaling by 1/(2*SinRange*N*2^{SinRescal}) during Coeffs To Slots
	StCRescale bool // Rescaling by 2^45/(initial scale) during Slots To Coeffs

	SinDepth uint64 // Automatically set
}

type SinType uint64

const (
	Sin = SinType(0)
	Cos = SinType(1)
)

func (b *BootParams) Gen() {

	if b.SinType == SinType(Sin) && b.SinRescal != 0 {
		panic("BootParams: cannot use double angle formul for SinType = Sin -> must use SinType = Cos")
	}

	b.SinDepth = uint64(math.Ceil(math.Log2(float64(b.SinDeg))) + float64(b.SinRescal) + 1)

	if !b.CtSRescale {
		b.SinDepth++
	}

	b.Parameters.Gen()
}

var BootstrappParams = []*BootParams{

	// 1430 Sin
	{Parameters: Parameters{
		LogN:     16,
		LogSlots: 10,
		LogModuli: LogModuli{
			LogQi: []uint64{55, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 50, 25, 55, 55, 55, 55, 55, 55, 55, 55, 55, 45, 45, 45},
			LogPi: []uint64{55, 55, 55, 55},
		},
		Scale: 1 << 45,
		Sigma: 3.2,
	},
		SinType:    Sin,
		SinRange:   15,
		SinDeg:     127,
		SinRescal:  0,
		BabySplit:  3,
		CtSDepth:   3,
		StCDepth:   3,
		CtSRescale: false,
		StCRescale: false},

	// 1435 Sin
	{Parameters: Parameters{
		LogN:     16,
		LogSlots: 10,
		LogModuli: LogModuli{
			LogQi: []uint64{55, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 50, 25, 55, 55, 55, 55, 55, 55, 55, 55, 50, 50, 50},
			LogPi: []uint64{55, 55, 55, 55},
		},
		Scale: 1 << 45,
		Sigma: 3.2,
	},
		SinType:    Sin,
		SinRange:   15,
		SinDeg:     127,
		SinRescal:  0,
		BabySplit:  3,
		CtSDepth:   3,
		StCDepth:   3,
		CtSRescale: true,
		StCRescale: false},

	// 1440 Sin
	{Parameters: Parameters{
		LogN:     16,
		LogSlots: 10,
		LogModuli: LogModuli{
			LogQi: []uint64{55, 60, 60, 60, 60, 60, 60, 60, 60, 30, 55, 55, 55, 55, 55, 55, 55, 55, 55, 45, 45, 45},
			LogPi: []uint64{61, 61, 61, 61},
		},
		Scale: 1 << 30,
		Sigma: 3.2,
	},
		SinType:    Sin,
		SinRange:   15,
		SinDeg:     127,
		SinRescal:  0,
		BabySplit:  3,
		CtSDepth:   3,
		StCDepth:   3,
		CtSRescale: false,
		StCRescale: true},

	// 1430 Sin
	{Parameters: Parameters{
		LogN:     16,
		LogSlots: 10,
		LogModuli: LogModuli{
			LogQi: []uint64{55, 60, 60, 60, 60, 60, 60, 60, 60, 60, 55, 55, 55, 55, 55, 55, 55, 55, 50, 50, 50},
			LogPi: []uint64{61, 61, 61, 61},
		},
		Scale: 1 << 30,
		Sigma: 3.2,
	},
		SinType:    Sin,
		SinRange:   15,
		SinDeg:     127,
		SinRescal:  0,
		BabySplit:  3,
		CtSDepth:   3,
		StCDepth:   3,
		CtSRescale: true,
		StCRescale: true},

	// 1425 cos
	{Parameters: Parameters{
		LogN:     16,
		LogSlots: 10,
		LogModuli: LogModuli{
			LogQi: []uint64{55, 60, 60, 60, 60, 60, 60, 60, 60, 55, 55, 55, 55, 55, 55, 55, 55, 55, 50, 50, 50},
			LogPi: []uint64{61, 61, 61, 61},
		},
		Scale: 1 << 30,
		Sigma: 3.2,
	},
		SinType:    Cos,
		SinRange:   15,
		SinDeg:     38,
		SinRescal:  2,
		BabySplit:  2,
		CtSDepth:   3,
		StCDepth:   3,
		CtSRescale: true,
		StCRescale: true},
}
