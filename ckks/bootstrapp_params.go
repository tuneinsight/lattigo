package ckks

// BootstrappParams is a struct for the default bootstrapping parameters
type BootstrappParams struct {
	H            uint64   // Hamming weight of the secret key
	SinType      SinType  // Chose betwenn [Sin(2*pi*x)] or [cos(2*pi*x/r) with double angle formula]
	SinRange     uint64   // K parameter (interpolation in the range -K to K)
	SinDeg       uint64   // Degree of the interpolation
	SinRescal    uint64   // Number of rescale and double angle formula (only applies for cos)
	CtSLevel     []uint64 // Level of the Coeffs To Slots
	StCLevel     []uint64 // Level of the Slots To Coeffs
	MaxN1N2Ratio float64  // n1/n2 ratio for the bsgs algo for matrix x vector eval
}

// CtSDepth returns the number of levels allocated to CoeffsToSlots
func (b *BootstrappParams) CtSDepth() uint64 {
	return uint64(len(b.CtSLevel))
}

// StCDepth returns the number of levels allocated to SlotToCoeffs
func (b *BootstrappParams) StCDepth() uint64 {
	return uint64(len(b.StCLevel))
}

// SinType is the type of function used during the bootstrapping
// for the homomorphic modular reduction
type SinType uint64

// Sin and Cos are the two proposed functions for SinType
const (
	Sin = SinType(0)
	Cos = SinType(1)
)

// Copy return a new BootstrappParams which is a copy of the target
func (b *BootstrappParams) Copy() *BootstrappParams {
	paramsCopy := &BootstrappParams{
		H:            b.H,
		SinType:      b.SinType,
		SinRange:     b.SinRange,
		SinDeg:       b.SinDeg,
		SinRescal:    b.SinRescal,
		CtSLevel:     make([]uint64, len(b.CtSLevel)),
		StCLevel:     make([]uint64, len(b.StCLevel)),
		MaxN1N2Ratio: b.MaxN1N2Ratio,
	}
	copy(paramsCopy.CtSLevel, b.CtSLevel)
	copy(paramsCopy.StCLevel, b.StCLevel)
	return paramsCopy
}

// DefaultBootstrappSchemeParams are default scheme params for the bootstrapping
var DefaultBootstrappSchemeParams = []*Parameters{
	{
		logN:     16,
		logSlots: 15,
		Moduli: Moduli{
			qi: []uint64{
				0x80000000080001,  // 55 Q0
				0x2000000a0001,    // 45
				0x2000000e0001,    // 45
				0x1fffffc20001,    // 45
				0x200000440001,    // 45
				0x200000500001,    // 45
				0x200000620001,    // 45
				0x1fffff980001,    // 45
				0x2000006a0001,    // 45
				0x1fffff7e0001,    // 45
				0x200000860001,    // 45
				0x200000a60001,    // 45
				0x100000000060001, // 56 StC (28 + 28)
				0xffa0001,         // 28 StC
				0x7fffffffba0001,  // 55 Sine
				0x80000000500001,  // 55 Sine
				0x7fffffffaa0001,  // 55 Sine
				0x800000005e0001,  // 55 Sine
				0x7fffffff7e0001,  // 55 Sine
				0x7fffffff380001,  // 55 Sine
				0x80000000ca0001,  // 55 Sine
				0x200000000e0001,  // 53 CtS
				0x20000000140001,  // 53 CtS
				0x20000000280001,  // 53 CtS
			},
			pi: []uint64{
				0x80000000e00001, // 55
				0x7ffffffef00001, // 55
				0x800000011c0001, // 55
				0x7ffffffeba0001, // 55
			},
		},
		scale: 1 << 45,
		sigma: DefaultSigma,
	},

	{
		logN:     16,
		logSlots: 15,
		Moduli: Moduli{
			qi: []uint64{
				0x80000000080001,  // 55 Q0
				0x2000000a0001,    // 45
				0x2000000e0001,    // 45
				0x1fffffc20001,    // 45
				0x200000440001,    // 45
				0x200000500001,    // 45
				0x200000620001,    // 45
				0x1fffff980001,    // 45
				0x2000006a0001,    // 45
				0x1fffff7e0001,    // 45
				0x200000860001,    // 45
				0x100000000060001, // 56 StC (28 + 28)
				0xffa0001,         // 28 StC
				0x80000000440001,  // 55 Sine (double angle)
				0x7fffffffba0001,  // 55 Sine (double angle)
				0x80000000500001,  // 55 Sine
				0x7fffffffaa0001,  // 55 Sine
				0x800000005e0001,  // 55 Sine
				0x7fffffff7e0001,  // 55 Sine
				0x7fffffff380001,  // 55 Sine
				0x80000000ca0001,  // 55 Sine
				0x200000000e0001,  // 53 CtS
				0x20000000140001,  // 53 CtS
				0x20000000280001,  // 53 CtS
				0x1fffffffd80001,  // 53 CtS
			},
			pi: []uint64{
				0xfffffffff00001,  // 56
				0xffffffffd80001,  // 56
				0x1000000002a0001, // 56
				0xffffffffd20001,  // 56
				0x100000000480001, // 56
			},
		},
		scale: 1 << 45,
		sigma: DefaultSigma,
	},

	{
		logN:     16,
		logSlots: 15,
		Moduli: Moduli{
			qi: []uint64{
				0x80000000080001,   // 55 Q0
				0xffffffffffc0001,  // 60
				0x10000000006e0001, // 60
				0xfffffffff840001,  // 60
				0x1000000000860001, // 60
				0xfffffffff6a0001,  // 60
				0x1000000000980001, // 60
				0xfffffffff5a0001,  // 60
				0x1000000000b00001, // 60 StC (30)
				0x1000000000ce0001, // 60 StC (30 + 30)
				0x7fffffffaa0001,   // 55 Sine
				0x800000005e0001,   // 55 Sine
				0x7fffffff7e0001,   // 55 Sine
				0x7fffffff380001,   // 55 Sine
				0x80000000ca0001,   // 55 Sine
				0x80000000e00001,   // 55 Sine
				0x7ffffffef00001,   // 55 Sine
				0x200000000e0001,   // 53 CtS
				0x20000000140001,   // 53 CtS
				0x20000000280001,   // 53 CtS
			},
			pi: []uint64{
				0x1fffffffffe00001, // 61
				0x1fffffffffc80001, // 61
				0x1fffffffffb40001, // 61
				0x1fffffffff500001, // 61
			},
		},
		scale: 1 << 30,
		sigma: DefaultSigma,
	},

	{
		logN:     16,
		logSlots: 15,
		Moduli: Moduli{
			qi: []uint64{
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
			},
			pi: []uint64{
				0x1fffffffffe00001, // Pi 61
				0x1fffffffffc80001, // Pi 61
				0x1fffffffffb40001, // Pi 61
				0x1fffffffff500001, // Pi 61
			},
		},
		scale: 1 << 30,
		sigma: DefaultSigma,
	},

	{
		logN:     16,
		logSlots: 15,
		Moduli: Moduli{
			qi: []uint64{
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
			pi: []uint64{
				0x1fffffffffe00001, // Pi 61
				0x1fffffffffc80001, // Pi 61
				0x1fffffffffb40001, // Pi 61
				0x1fffffffff500001, // Pi 61
				0x1fffffffff420001, // Pi 61
			},
		},
		scale: 1 << 30,
		sigma: DefaultSigma,
	},
}

// DefaultBootstrappParams are default bootstrapping params for the bootstrapping
var DefaultBootstrappParams = []*BootstrappParams{

	// SET I
	// h = 128
	// 1398 Sin - 550
	{
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
	{
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
	{
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
	{
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
	{
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

// DefaultBootstrappSchemeParamsShort are insecure params
// for quick correctness testing of the bootstrapping
var DefaultBootstrappSchemeParamsShort = []*Parameters{

	{
		logN:     14,
		logSlots: 13,
		Moduli: Moduli{
			qi: []uint64{
				0x80000000080001,  // 55 Q0
				0x2000000a0001,    // 45
				0x2000000e0001,    // 45
				0x1fffffc20001,    // 45
				0x200000440001,    // 45
				0x200000500001,    // 45
				0x200000620001,    // 45
				0x1fffff980001,    // 45
				0x2000006a0001,    // 45
				0x1fffff7e0001,    // 45
				0x200000860001,    // 45
				0x100000000060001, // 56 StC (28 + 28)
				0xffa0001,         // 28 StC
				0x80000000440001,  // 55 Sine (double angle)
				0x7fffffffba0001,  // 55 Sine (double angle)
				0x80000000500001,  // 55 Sine
				0x7fffffffaa0001,  // 55 Sine
				0x800000005e0001,  // 55 Sine
				0x7fffffff7e0001,  // 55 Sine
				0x7fffffff380001,  // 55 Sine
				0x80000000ca0001,  // 55 Sine
				0x200000000e0001,  // 53 CtS
				0x20000000140001,  // 53 CtS
				0x20000000280001,  // 53 CtS
				0x1fffffffd80001,  // 53 CtS
			},
			pi: []uint64{
				0xfffffffff00001,  // 56
				0xffffffffd80001,  // 56
				0x1000000002a0001, // 56
				0xffffffffd20001,  // 56
				0x100000000480001, // 56
			},
		},
		scale: 1 << 45,
		sigma: DefaultSigma,
	},

	{
		logN:     14,
		logSlots: 12,
		Moduli: Moduli{
			qi: []uint64{
				0x80000000080001,  // 55 Q0
				0x2000000a0001,    // 45
				0x2000000e0001,    // 45
				0x1fffffc20001,    // 45
				0x200000440001,    // 45
				0x200000500001,    // 45
				0x200000620001,    // 45
				0x1fffff980001,    // 45
				0x2000006a0001,    // 45
				0x1fffff7e0001,    // 45
				0x200000860001,    // 45
				0x100000000060001, // 56 StC (28 + 28)
				0xffa0001,         // 28 StC
				0x80000000440001,  // 55 Sine (double angle)
				0x7fffffffba0001,  // 55 Sine (double angle)
				0x80000000500001,  // 55 Sine
				0x7fffffffaa0001,  // 55 Sine
				0x800000005e0001,  // 55 Sine
				0x7fffffff7e0001,  // 55 Sine
				0x7fffffff380001,  // 55 Sine
				0x80000000ca0001,  // 55 Sine
				0x200000000e0001,  // 53 CtS
				0x20000000140001,  // 53 CtS
				0x20000000280001,  // 53 CtS
				0x1fffffffd80001,  // 53 CtS
			},
			pi: []uint64{
				0xfffffffff00001,  // 56
				0xffffffffd80001,  // 56
				0x1000000002a0001, // 56
				0xffffffffd20001,  // 56
				0x100000000480001, // 56
			},
		},
		scale: 1 << 45,
		sigma: DefaultSigma,
	},
}

// DefaultBootstrappParamsShort are default bootstrapping params for the
// DefaultBootstrappSchemeParamsShort scheme params
var DefaultBootstrappParamsShort = []*BootstrappParams{

	// SET II
	// 1525 Cos - 550
	{
		H:            196,
		SinType:      Cos,
		SinRange:     21,
		SinDeg:       52,
		SinRescal:    2,
		CtSLevel:     []uint64{24, 23, 22, 21},
		StCLevel:     []uint64{12, 11, 11},
		MaxN1N2Ratio: 16.0,
	},

	// SET II
	// 1525 Cos - 550
	{
		H:            196,
		SinType:      Cos,
		SinRange:     21,
		SinDeg:       52,
		SinRescal:    2,
		CtSLevel:     []uint64{24, 23, 22, 21},
		StCLevel:     []uint64{12, 11, 11},
		MaxN1N2Ratio: 16.0,
	},
}
