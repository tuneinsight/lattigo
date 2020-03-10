package ckks

import (
	"math"
)

type BootParams struct {
	Parameters

	sinType   SinType
	sinDeg    uint64
	sinRescal uint64
	sinDepth  uint64

	ctsDepth uint64
	stcDepth uint64

	ctsRescale bool
	stcRescale bool
}

type SinType uint64

const (
	Sin = uint64(0)
	Cos = uint64(1)
)

func (b *BootParams) GenFromLogModuli() {

	if b.sinType == SinType(Sin) && b.sinRescal != 0 {
		panic("BootParams: cannot use double angle formul for sinType = Sin -> must use sinType = Cos")
	}

	b.sinDepth = uint64(math.Ceil(math.Log2(float64(b.sinDeg))) + float64(b.sinRescal))

	if !b.ctsRescale {
		b.sinDepth++
	}

	b.Parameters.GenFromLogModuli()
}

var BootstrappParams = []*BootParams{

	// 1430 sin
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
		sinType:    SinType(Sin),
		sinDeg:     127,
		sinRescal:  0,
		ctsDepth:   3,
		stcDepth:   3,
		ctsRescale: false,
		stcRescale: false},

	// 1435 sin
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
		sinType:    SinType(Sin),
		sinDeg:     127,
		sinRescal:  0,
		ctsDepth:   3,
		stcDepth:   3,
		ctsRescale: true,
		stcRescale: false},

	// 1440 sin
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
		sinType:    SinType(Sin),
		sinDeg:     127,
		sinRescal:  0,
		ctsDepth:   3,
		stcDepth:   3,
		ctsRescale: false,
		stcRescale: true},

	// 1430 sin
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
		sinType:    SinType(Sin),
		sinDeg:     127,
		sinRescal:  0,
		ctsDepth:   3,
		stcDepth:   3,
		ctsRescale: true,
		stcRescale: true},

	// 1425 cos
	{Parameters: Parameters{
		LogN:     16,
		LogSlots: 10,
		LogModuli: LogModuli{
			LogQi: []uint64{55, 60, 60, 60, 60, 60, 60, 60, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 50, 50, 50},
			LogPi: []uint64{61, 61, 61, 61},
		},
		Scale: 1 << 30,
		Sigma: 3.2,
	},
		sinType:    SinType(Cos),
		sinDeg:     30,
		sinRescal:  3,
		ctsDepth:   3,
		stcDepth:   3,
		ctsRescale: true,
		stcRescale: true},
}
