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

func (b *BootParams) Copy() *BootParams {
	paramsCopy := &BootParams{
		Parameters: *b.Parameters.Copy(),
		SinType:    b.SinType,
		SinRange:   b.SinRange,
		SinDeg:     b.SinDeg,
		SinRescal:  b.SinRescal,
		BabySplit:  b.BabySplit,
		CtSLevel:   make([]uint64, len(b.CtSLevel)),
		StCLevel:   make([]uint64, len(b.StCLevel)),
		SinDepth:   b.SinDepth,
	}
	copy(paramsCopy.CtSLevel, b.CtSLevel)
	copy(paramsCopy.StCLevel, b.StCLevel)
	return paramsCopy
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
		Moduli: Moduli{
			Qi: []uint64{
				0x80000000080001,
				0x2000000a0001,
				0x2000000e0001,
				0x1fffffc20001,
				0x200000440001,
				0x200000500001,
				0x200000620001,
				0x1fffff980001,
				0x2000006a0001,
				0x1fffff7e0001,
				0x200000860001,
				0x200000a60001,
				0x100000000060001,
				0xffa0001,
				0x7fffffffba0001,
				0x80000000500001,
				0x7fffffffaa0001,
				0x800000005e0001,
				0x7fffffff7e0001,
				0x7fffffff380001,
				0x80000000ca0001,
				0x200000000e0001,
				0x20000000140001,
				0x20000000280001,
			},
			Pi: []uint64{
				0x80000000e00001,
				0x7ffffffef00001,
				0x800000011c0001,
				0x7ffffffeba0001,
			},
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
	// 15 : 18.28
	{Parameters: Parameters{
		LogN:     16,
		LogSlots: 15,
		Moduli: Moduli{
			Qi: []uint64{
				0x80000000080001,
				0xffffffffffc0001,
				0x10000000006e0001,
				0xfffffffff840001,
				0x1000000000860001,
				0xfffffffff6a0001,
				0x1000000000980001,
				0xfffffffff5a0001,
				0x1000000000b00001,
				0x1000000000ce0001,
				0x7fffffffaa0001,
				0x800000005e0001,
				0x7fffffff7e0001,
				0x7fffffff380001,
				0x80000000ca0001,
				0x80000000e00001,
				0x7ffffffef00001,
				0x200000000e0001,
				0x20000000140001,
				0x20000000280001,
			},
			Pi: []uint64{
				0x1fffffffffe00001,
				0x1fffffffffc80001,
				0x1fffffffffb40001,
				0x1fffffffff500001,
			},
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
