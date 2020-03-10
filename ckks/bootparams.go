package ckks

type BootParams struct {
	Parameters
	sinDepth   uint64
	ctsDepth   uint64
	stcDepth   uint64
	ctsRescale bool
	stcRescale bool
}

func (b *BootParams) GenFromLogModuli() {
	b.Parameters.GenFromLogModuli()
}

var BootstrappParams = []*BootParams{

	// 1430
	{Parameters: Parameters{LogN: 16,
		LogSlots: 10,
		LogModuli: LogModuli{
			LogQi: []uint64{55, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 50, 25, 55, 55, 55, 55, 55, 55, 55, 55, 55, 45, 45, 45},
			LogPi: []uint64{55, 55, 55, 55},
		},
		Scale: 1 << 45,
		Sigma: 3.2,
	},
		sinDepth:   9,
		ctsDepth:   3,
		stcDepth:   3,
		ctsRescale: false,
		stcRescale: false},

	// 1435
	{Parameters: Parameters{LogN: 16,
		LogSlots: 10,
		LogModuli: LogModuli{
			LogQi: []uint64{55, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 50, 25, 55, 55, 55, 55, 55, 55, 55, 55, 50, 50, 50},
			LogPi: []uint64{55, 55, 55, 55},
		},
		Scale: 1 << 45,
		Sigma: 3.2,
	},
		sinDepth:   8,
		ctsDepth:   3,
		stcDepth:   3,
		ctsRescale: true,
		stcRescale: false},

	// 1440
	{Parameters: Parameters{LogN: 16,
		LogSlots: 10,
		LogModuli: LogModuli{
			LogQi: []uint64{55, 60, 60, 60, 60, 60, 60, 60, 60, 30, 55, 55, 55, 55, 55, 55, 55, 55, 55, 45, 45, 45},
			LogPi: []uint64{61, 61, 61, 61},
		},
		Scale: 1 << 30,
		Sigma: 3.2,
	},
		sinDepth:   9,
		ctsDepth:   3,
		stcDepth:   3,
		ctsRescale: false,
		stcRescale: true},

	// 1430
	{Parameters: Parameters{LogN: 16,
		LogSlots: 10,
		LogModuli: LogModuli{
			LogQi: []uint64{55, 60, 60, 60, 60, 60, 60, 60, 60, 60, 55, 55, 55, 55, 55, 55, 55, 55, 50, 50, 50},
			LogPi: []uint64{61, 61, 61, 61},
		},
		Scale: 1 << 30,
		Sigma: 3.2,
	},
		sinDepth:   8,
		ctsDepth:   3,
		stcDepth:   3,
		ctsRescale: true,
		stcRescale: true},

	// 1425 cos
	{Parameters: Parameters{LogN: 16,
		LogSlots: 10,
		LogModuli: LogModuli{
			LogQi: []uint64{55, 60, 60, 60, 60, 60, 60, 60, 60, 55, 55, 55, 55, 55, 55, 55, 55, 55, 50, 50, 50},
			LogPi: []uint64{61, 61, 61, 61},
		},
		Scale: 1 << 30,
		Sigma: 3.2,
	},
		sinDepth:   9,
		ctsDepth:   3,
		stcDepth:   3,
		ctsRescale: true,
		stcRescale: true},
}
