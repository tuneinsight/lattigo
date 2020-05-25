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

	CtSLevel []uint64 // Level of the Coeffs To Slots
	StCLevel []uint64 // Level of the Slots To Coeffs

	SinDepth uint64 // Automatically set
}

type SinType uint64

const (
	Sin = SinType(0)
	Cos = SinType(1)
)

func (b *BootParams) Gen() error { // TODO "Generate" ?

	if b.SinType == SinType(Sin) && b.SinRescal != 0 {
		panic("BootParams: cannot use double angle formul for SinType = Sin -> must use SinType = Cos")
	}

	b.SinDepth = uint64(math.Ceil(math.Log2(float64(b.SinDeg))) + float64(b.SinRescal) + 1)

	return b.Parameters.Gen()
}

//TODO : hardcode moduli chain
var BootstrappParams = []*BootParams{

	// 1398 Sin - 550
	// 10 : 20.67
	// 11 : 20.09
	// 12 : 19.32
	// 13 : 18.16
	// 14 : 17.60
	// 15 : 18.31
	{Parameters: Parameters{
		LogN:     16,
		LogSlots: 15,
		LogModuli: LogModuli{
			LogQi: []uint64{55, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 55, 56, 28, 55, 55, 55, 55, 55, 55, 55, 53, 53, 53},
			LogPi: []uint64{55, 55, 55, 55},
		},
		Scale: 1 << 45,
		Sigma: 3.2,
	},
		SinType:   Sin,
		SinRange:  15,
		SinDeg:    127,
		SinRescal: 0,
		CtSLevel:  []uint64{23, 22, 21},
		StCLevel:  []uint64{13, 12, 12},
	},

	// 1384 Sin - 505
	// 10 : 20.57
	// 11 : 20.00
	// 12 : 19.50
	// 13 : 19.01
	// 14 : 18.50
	// 15 : 18.20
	{Parameters: Parameters{
		LogN:     16,
		LogSlots: 15,
		LogModuli: LogModuli{
			LogQi: []uint64{55, 60, 60, 60, 60, 55, 55, 55, 60, 60, 55, 55, 55, 55, 55, 55, 55, 53, 53, 53},
			LogPi: []uint64{61, 61, 61, 61},
		},
		Scale: 1 << 30,
		Sigma: 3.2,
	},
		SinType:   Sin,
		SinRange:  15,
		SinDeg:    127,
		SinRescal: 0,
		CtSLevel:  []uint64{19, 18, 17},
		StCLevel:  []uint64{9, 9, 8},
	},

	// 1408 Cos - 550
	// 10 : 20.91
	// 11 : 19.93
	// 12 : 19.54
	// 13 : 18.27
	// 14 : 17.60 40
	// 15 : 20.10 40
	{Parameters: Parameters{
		LogN:     16,
		LogSlots: 15,
		LogModuli: LogModuli{
			LogQi: []uint64{55, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 56, 28, 55, 55, 55, 55, 55, 55, 55, 55, 53, 53, 53},
			LogPi: []uint64{55, 55, 55, 55},
		},
		Scale: 1 << 45,
		Sigma: 3.2,
	},
		SinType:   Cos,
		SinRange:  16,
		SinDeg:    40,
		SinRescal: 2,
		CtSLevel:  []uint64{23, 22, 21},
		StCLevel:  []uint64{12, 11, 11},
	},

	// 1439 Cos - 505
	// 10 : 21.44
	// 11 : 21.00
	// 12 : 20.46
	// 13 : 19.89
	// 14 : 19.37
	// 15 : 19.51
	{Parameters: Parameters{
		LogN:     16,
		LogSlots: 15,
		LogModuli: LogModuli{
			LogQi: []uint64{55, 60, 60, 60, 60, 60, 60, 60, 60, 60, 55, 55, 55, 55, 55, 55, 55, 55, 53, 53, 53},
			LogPi: []uint64{61, 61, 61, 61},
		},
		Scale: 1 << 30,
		Sigma: 3.2,
	},
		SinType:   Cos,
		SinRange:  16,
		SinDeg:    40,
		SinRescal: 2,
		CtSLevel:  []uint64{20, 19, 18},
		StCLevel:  []uint64{9, 9, 8},
	},
}
