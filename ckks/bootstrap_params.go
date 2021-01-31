package ckks

// BootstrappingParameters is a struct for the default bootstrapping parameters
type BootstrappingParameters struct {
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
func (b *BootstrappingParameters) CtSDepth() uint64 {
	return uint64(len(b.CtSLevel))
}

// StCDepth returns the number of levels allocated to SlotToCoeffs
func (b *BootstrappingParameters) StCDepth() uint64 {
	return uint64(len(b.StCLevel))
}

// SinType is the type of function used during the bootstrapping
// for the homomorphic modular reduction
type SinType uint64

// Sin and Cos are the two proposed functions for SinType
const (
	Sin  = SinType(0) // Standard Chebyshev approximation of (1/2pi) * sin(2pix)
	Cos1 = SinType(1) // Special approximation (Han and Ki) of pow((1/2pi), 1/2^r) * cos(2pi(x-0.25)/2^r)
	Cos2 = SinType(2) // Standard Chebyshev approximation of pow((1/2pi), 1/2^r) * cos(2pi(x-0.25)/2^r)
)

// Copy return a new BootstrapParams which is a copy of the target
func (b *BootstrappingParameters) Copy() *BootstrappingParameters {
	paramsCopy := &BootstrappingParameters{
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

// DefaultBootstrapSchemeParams are default scheme params for the bootstrapping
var DefaultBootstrapSchemeParams = []*Parameters{

	{
		logN:     14,
		logSlots: 13,
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
			0x80000000e00001, 
			0x7ffffffef00001, 
			0x800000011c0001,
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
		scale: 1 << 45,
		sigma: DefaultSigma,
	},

	{
		logN:     16,
		logSlots: 15,
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
		scale: 1 << 30,
		sigma: DefaultSigma,
	},

	{
		logN:     16,
		logSlots: 15,
		qi: []uint64{
			0x80000000080001,   // 55 Q0
			0x2000000a0001,     // 45
			0x2000000e0001,     // 45
			0x1fffffc20001,     // 45
			0x200000440001,     // 45
			0x200000500001,     // 45
			0x200000620001,     // 45
			0x1fffff980001,     // 45
			0x2000006a0001,     // 45
			0x1fffff7e0001,     // 45
			0x100000000060001,  // 56 StC (28 + 28)
			0xffa0001,          // 28 StC
			0xffffffffffc0001,  // 60 Sine (double angle)
			0x10000000006e0001, // 60 Sine (double angle)
			0xfffffffff840001,  // 60 Sine (double angle)
			0x1000000000860001, // 60 Sine
			0xfffffffff6a0001,  // 60 Sine
			0x1000000000980001, // 60 Sine
			0xfffffffff5a0001,  // 60 Sine
			0x1000000000b00001, // 60 Sine
			0x1000000000ce0001, // 60 Sine
			0xfffffffff2a0001,  // 60 Sine
			0xfffffffff240001,  // 60 Sine
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
			0x1fffffffff380001, // Pi 61
		},
		scale: 1 << 45,
		sigma: DefaultSigma,
	},

	{
		logN:     15,
		logSlots: 14,
		qi: []uint64{
			0x7fffb0001,       // 35 Q0
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
		pi: []uint64{
			0x7e40000000001, // 50
			0x7c80000000001, // 50
		},
		scale: 1 << 25,
		sigma: DefaultSigma,
	},
}

// DefaultBootstrapParams are default bootstrapping params for the bootstrapping
var DefaultBootstrapParams = []*BootstrappingParameters{

	// SET II
	// 1525 - 550
	{
		H:            96,
		SinType:      Cos1,
		SinRange:     21,
		SinDeg:       52,
		SinRescal:    2,
		CtSLevel:     []uint64{27, 26, 25, 24},
		StCLevel:     []uint64{12, 11, 11},
		MaxN1N2Ratio: 16.0,
	},

	// SET V
	// 1553 - 505
	{
		H:            192,
		SinType:      Cos1,
		SinRange:     21,
		SinDeg:       52,
		SinRescal:    2,
		CtSLevel:     []uint64{21, 20, 19, 18},
		StCLevel:     []uint64{9, 9, 8},
		MaxN1N2Ratio: 16.0,
	},

	// Set VII
	// 1773 - 460
	{
		H:            32768,
		SinType:      Cos2,
		SinRange:     257,
		SinDeg:       250,
		SinRescal:    3,
		CtSLevel:     []uint64{26, 25, 24, 23},
		StCLevel:     []uint64{11, 10, 10},
		MaxN1N2Ratio: 16.0,
	},

	// Set IV
	// 768 - 110
	{
		H:            192,
		SinType:      Cos1,
		SinRange:     21,
		SinDeg:       52,
		SinRescal:    2,
		CtSLevel:     []uint64{13, 12},
		StCLevel:     []uint64{3, 3},
		MaxN1N2Ratio: 16.0,
	},
}
