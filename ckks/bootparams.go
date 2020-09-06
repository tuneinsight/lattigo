package ckks

import (
	"math"
)

type BootParams struct {
	Parameters
	LogModuli

	H uint64 // Hamming weight of the secret key

	SinType      SinType  // Chose betwenn [Sin(2*pi*x)] or [cos(2*pi*x/r) with double angle formula]
	SinRange     uint64   // K parameter (interpolation in the range -K to K)
	SinDeg       uint64   // Degree of the interpolation
	SinRescal    uint64   // Number of rescale and double angle formula (only applies for cos)
	CtSLevel     []uint64 // Level of the Coeffs To Slots
	StCLevel     []uint64 // Level of the Slots To Coeffs
	MaxN1N2Ratio float64  // n1/n2 ratio for the bsgs algo for matrix x vector eval

	SinDepth uint64 // Automatically set
}

type SinType uint64

const (
	Sin = SinType(0)
	Cos = SinType(1)
)

func (b *BootParams) Gen() (err error) { // TODO "Generate" ?

	if b.SinType == SinType(Sin) && b.SinRescal != 0 {
		panic("BootParams: cannot use double angle formul for SinType = Sin -> must use SinType = Cos")
	}

	b.SinDepth = uint64(math.Ceil(math.Log2(float64(b.SinDeg))) + float64(b.SinRescal))

	scale := b.Parameters.scale
	logSlots := b.Parameters.logSlots

	var p *Parameters

	if b.qi != nil {
		if p, err = NewParametersFromModuli(b.logN, b.Moduli); err != nil {
			return err
		}
	} else {
		if p, err = NewParametersFromLogModuli(b.logN, b.LogModuli); err != nil {
			return err
		}
	}

	b.Parameters = *p

	b.SetScale(scale)
	b.SetLogSlots(logSlots)

	return nil
}

func (b *BootParams) Copy() *BootParams {
	paramsCopy := &BootParams{
		Parameters:   *b.Parameters.Copy(),
		LogModuli:    b.LogModuli.Copy(),
		SinType:      b.SinType,
		SinRange:     b.SinRange,
		SinDeg:       b.SinDeg,
		SinRescal:    b.SinRescal,
		CtSLevel:     make([]uint64, len(b.CtSLevel)),
		StCLevel:     make([]uint64, len(b.StCLevel)),
		MaxN1N2Ratio: b.MaxN1N2Ratio,
		SinDepth:     b.SinDepth,
	}
	copy(paramsCopy.CtSLevel, b.CtSLevel)
	copy(paramsCopy.StCLevel, b.StCLevel)
	return paramsCopy
}

var BootstrappParams = []*BootParams{

	// SET I
	// h = 128
	// 1398 Sin - 550
	{Parameters: Parameters{
		logN:     16,
		logSlots: 15,
		Moduli: Moduli{
			qi: []uint64{
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
			pi: []uint64{
				0x80000000e00001,
				0x7ffffffef00001,
				0x800000011c0001,
				0x7ffffffeba0001,
			},
		},
		scale: 1 << 45,
		sigma: DefaultSigma,
	},
		H:            128,
		SinType:      Sin,
		SinRange:     15,
		SinDeg:       127,
		SinRescal:    0,
		CtSLevel:     []uint64{23, 22, 21},
		StCLevel:     []uint64{13, 12, 12},
		MaxN1N2Ratio: 16.0,
	},

	// SET II
	// 1525 Cos - 550
	{Parameters: Parameters{
		logN:     16,
		logSlots: 15,
		scale:    1 << 45,
		sigma:    DefaultSigma,
	},
		LogModuli: LogModuli{
			LogQi: []uint64{55, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 56, 28, 55, 55, 55, 55, 55, 55, 55, 55, 53, 53, 53, 53},
			LogPi: []uint64{56, 56, 56, 56, 56},
		},
		H:            196,
		SinType:      Cos,
		SinRange:     21,
		SinDeg:       52,
		SinRescal:    2,
		CtSLevel:     []uint64{24, 23, 22, 21},
		StCLevel:     []uint64{12, 11, 11},
		MaxN1N2Ratio: 16.0,
	},

	// SET III
	// 1384 Sin - 505
	{Parameters: Parameters{
		logN:     16,
		logSlots: 15,
		Moduli: Moduli{
			qi: []uint64{
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
			pi: []uint64{
				0x1fffffffffe00001,
				0x1fffffffffc80001,
				0x1fffffffffb40001,
				0x1fffffffff500001,
			},
		},
		scale: 1 << 30,
		sigma: DefaultSigma,
	},
		H:            128,
		SinType:      Sin,
		SinRange:     15,
		SinDeg:       127,
		SinRescal:    0,
		CtSLevel:     []uint64{19, 18, 17},
		StCLevel:     []uint64{9, 9, 8},
		MaxN1N2Ratio: 16.0,
	},

	// SET IV
	// 1439 Cos - 505
	{Parameters: Parameters{
		logN:     16,
		logSlots: 10,
		scale:    1 << 30,
		sigma:    DefaultSigma,
	},
		LogModuli: LogModuli{
			LogQi: []uint64{55, 60, 60, 60, 60, 60, 60, 60, 60, 60, 55, 55, 55, 55, 55, 55, 55, 55, 53, 53, 53},
			LogPi: []uint64{61, 61, 61, 61},
		},
		H:            128,
		SinType:      Cos,
		SinRange:     17,
		SinDeg:       42,
		SinRescal:    2,
		CtSLevel:     []uint64{20, 19, 18},
		StCLevel:     []uint64{9, 9, 8},
		MaxN1N2Ratio: 16.0,
	},

	// SET V
	// 1553 Cos - 505
	{Parameters: Parameters{
		logN:     16,
		logSlots: 15,
		scale:    1 << 30,
		sigma:    DefaultSigma,
	},
		LogModuli: LogModuli{
			LogQi: []uint64{55, 60, 60, 60, 60, 60, 60, 60, 60, 60, 55, 55, 55, 55, 55, 55, 55, 55, 53, 53, 53, 53},
			LogPi: []uint64{61, 61, 61, 61, 61},
		},
		H:            192,
		SinType:      Cos,
		SinRange:     21,
		SinDeg:       52,
		SinRescal:    2,
		CtSLevel:     []uint64{21, 20, 19, 18},
		StCLevel:     []uint64{9, 9, 8},
		MaxN1N2Ratio: 16.0,
	},
}
