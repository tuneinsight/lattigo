package ckks

import(
	"math"
)





// SinType is the type of function used during the bootstrapping
// for the homomorphic modular reduction
type SinType uint64

// Sin and Cos are the two proposed functions for SinType
const (
	Sin  = SinType(0) // Standard Chebyshev approximation of (1/2pi) * sin(2pix)
	Cos1 = SinType(1) // Special approximation (Han and Ki) of pow((1/2pi), 1/2^r) * cos(2pi(x-0.25)/2^r)
	Cos2 = SinType(2) // Standard Chebyshev approximation of pow((1/2pi), 1/2^r) * cos(2pi(x-0.25)/2^r)
)




// DefaultBootstrapSchemeParams are default scheme params for the bootstrapping
var DefaultBootstrapSchemeParams = []*Parameters{
	/*
	{
		logN:     14,
		logSlots: 13,
		qi: []uint64{
			0x4000000120001,  // 55 Q0
			0x2000000a0001,    // 45
			0x2000000e0001,    // 45
			0x1fffffc20001,    // 45
			0x200000440001,    // 45
			0x200000500001,    // 45
			0x200000620001,    // 45
			0x1fffff980001,    // 45
			0x10004a0001, 
			0x1000500001,
			0x1000960001,
			0x4000000f20001, 
			0x40000010a0001, 
			0x4000001260001,
			0x3ffffffd20001,  // 55 Sine (double angle)
			0x4000000420001,  // 55 Sine (double angle)
			0x3ffffffb80001,  // 55 Sine
			0x4000000660001,  // 55 Sine
			0x40000007e0001,  // 55 Sine
			0x4000000800001,  // 55 Sine
			0x40000008a0001,  // 55 Sine
			0x4000000de0001,  // 55 Sine
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
	*/
}


// BootstrappingParameters is a struct for the default bootstrapping parameters
type BootstrappingParameters struct {
	ResidualModuli
	KeySwitchModuli
	SlotToCoeffsModuli
	SineEvalModuli
	CoeffsToSlotsModuli
	LogN uint64
	LogSlots uint64
	Scale float64
	Sigma float64
	H            uint64   // Hamming weight of the secret key
	SinType      SinType  // Chose betwenn [Sin(2*pi*x)] or [cos(2*pi*x/r) with double angle formula]
	MessageRatio float64  // Ratio between Q0 and m, i.e. Q[0]/|m|
	SinRange     uint64   // K parameter (interpolation in the range -K to K)
	SinDeg       uint64   // Degree of the interpolation
	SinRescal    uint64   // Number of rescale and double angle formula (only applies for cos)
	ArcSineDeg   uint64   // Degree of the Taylor arcsine composed with f(2*pi*x) (if zero then not used)
	MaxN1N2Ratio float64  // n1/n2 ratio for the bsgs algo for matrix x vector eval
}

func (b *BootstrappingParameters) Params() (p *Parameters, err error){
	Qi := append(b.ResidualModuli, b.SlotToCoeffsModuli.Qi...)
	Qi = append(Qi, b.SineEvalModuli.Qi...)
	Qi = append(Qi, b.CoeffsToSlotsModuli.Qi...)
	
	if p, err = NewParametersFromModuli(b.LogN, &Moduli{Qi, b.KeySwitchModuli}); err != nil{
		return nil, err
	}

	p.SetScale(b.Scale)
	p.SetLogSlots(b.LogSlots)
	p.SetSigma(b.Sigma)
	return 
}

// Copy return a new BootstrapParams which is a copy of the target
func (b *BootstrappingParameters) Copy() *BootstrappingParameters {
	paramsCopy := &BootstrappingParameters{
		H:            b.H,
		SinType:      b.SinType,
		MessageRatio: b.MessageRatio,
		SinRange:     b.SinRange,
		SinDeg:       b.SinDeg,
		SinRescal:    b.SinRescal,
		ArcSineDeg:   b.ArcSineDeg,
		MaxN1N2Ratio: b.MaxN1N2Ratio,
	}

	paramsCopy.ResidualModuli = make([]uint64, len(b.ResidualModuli))
	copy(paramsCopy.ResidualModuli, b.ResidualModuli)

	paramsCopy.CoeffsToSlotsModuli.Qi = make([]uint64, b.CtSDepth())
	copy(paramsCopy.CoeffsToSlotsModuli.Qi, b.CoeffsToSlotsModuli.Qi)

	paramsCopy.CoeffsToSlotsModuli.ScalingFactor = make([]float64, b.CtSDepth())
	copy(paramsCopy.CoeffsToSlotsModuli.ScalingFactor, b.CoeffsToSlotsModuli.ScalingFactor)

	paramsCopy.SineEvalModuli.Qi = make([]uint64, b.CtSDepth())
	copy(paramsCopy.SineEvalModuli.Qi, b.SineEvalModuli.Qi)

	paramsCopy.SineEvalModuli.ScalingFactor = b.SineEvalModuli.ScalingFactor

	paramsCopy.SlotToCoeffsModuli.Qi = make([]uint64, b.StCDepth())
	copy(paramsCopy.SlotToCoeffsModuli.Qi, b.SlotToCoeffsModuli.Qi)

	paramsCopy.SlotToCoeffsModuli.ScalingFactor = make([]float64, b.StCDepth())
	copy(paramsCopy.SlotToCoeffsModuli.ScalingFactor, b.SlotToCoeffsModuli.ScalingFactor)

	return paramsCopy
}

type ResidualModuli []uint64
type KeySwitchModuli []uint64

type CoeffsToSlotsModuli struct{
	Qi []uint64
	ScalingFactor []float64
}

type SineEvalModuli struct{
	Qi []uint64
	ScalingFactor float64
}

type SlotToCoeffsModuli struct{
	Qi []uint64
	ScalingFactor []float64
}

// SineDepth returns the depth of the SineEval. If true, then also 
// counts the double angle formula.
func (b *BootstrappingParameters) SineEvalDepth(withRescale bool) uint64{
	depth := uint64(math.Ceil(math.Log2(float64(b.SinDeg+1))))

	if withRescale{
		depth += b.SinRescal
	}

	return depth
}

// ArcSineDepth returns the depth of the arcsine polynomial.
func (b *BootstrappingParameters) ArcSineDepth() uint64{
	return uint64(math.Ceil(math.Log2(float64(b.ArcSineDeg+1))))
} 


// CtSDepth returns the number of levels allocated to CoeffsToSlots.
func (b *BootstrappingParameters) CtSDepth() uint64 {
	return uint64(len(b.CoeffsToSlotsModuli.Qi))
}

// StCDepth returns the number of levels allocated to SlotToCoeffs.
func (b *BootstrappingParameters) StCDepth() uint64 {
	return uint64(len(b.SlotToCoeffsModuli.Qi))
}

// DefaultBootstrapParams are default bootstrapping params for the bootstrapping.
var DefaultBootstrapParams = []*BootstrappingParameters{

	// SET II
	// 1525 - 550
	{	
		LogN:     14,
		LogSlots: 13,
		Scale: 1 << 45,
		Sigma: DefaultSigma,
		ResidualModuli: []uint64{
			0x4000000120001,  // 55 Q0
			0x2000000a0001,    // 45
			0x2000000e0001,    // 45
			0x1fffffc20001,    // 45
			0x200000440001,    // 45
			0x200000500001,    // 45
			0x200000620001,    // 45
			0x1fffff980001,    // 45
		},
		KeySwitchModuli: []uint64{
			0xfffffffff00001,  // 56
			0xffffffffd80001,  // 56
			0x1000000002a0001, // 56
			0xffffffffd20001,  // 56
			0x100000000480001, // 56
		},
		SlotToCoeffsModuli:SlotToCoeffsModuli{
			Qi :[]uint64{
				0x10004a0001, 
				0x1000500001,
				0x1000960001,
			},
			ScalingFactor: []float64{
				0x10004a0001,
				0x1000500001,
				0x1000960001,
			},
		},
		SineEvalModuli:SineEvalModuli{
			Qi :[]uint64{
				0x4000000f20001,  // 50 Arcsine 
				0x40000010a0001,  // 50 Arcsine 
				0x4000001260001,  // 50 Arcsine
				0x3ffffffd20001,  // 50 Sine (double angle)
				0x4000000420001,  // 50 Sine (double angle)
				0x3ffffffb80001,  // 50 Sine
				0x4000000660001,  // 50 Sine
				0x40000007e0001,  // 50 Sine
				0x4000000800001,  // 50 Sine
				0x40000008a0001,  // 50 Sine
				0x4000000de0001,  // 50 Sine
			},
			ScalingFactor: 1<<50,
		},
		CoeffsToSlotsModuli:CoeffsToSlotsModuli{
			Qi :[]uint64{
				0x200000000e0001,  // 53 CtS
				0x20000000140001,  // 53 CtS
				0x20000000280001,  // 53 CtS
				0x1fffffffd80001,  // 53 CtS
			},
			ScalingFactor: []float64{
				0x200000000e0001,
				0x20000000140001, 
				0x20000000280001, 
				0x1fffffffd80001,
			},
		},
		H:            192,
		SinType:      Cos1,
		MessageRatio: 4.0,
		SinRange:     25,
		SinDeg:       63,
		SinRescal:    2,
		ArcSineDeg:   7,
		MaxN1N2Ratio: 16.0,
	},

	/*
	// SET V
	// 1553 - 505
	{
		H:            192,
		SinType:      Cos1,
		MessageRatio:        1024.0,
		SinRange:     25,
		SinDeg:       63,
		SinRescal:    2,
		ArcSineDeg:   0,
		CtSLevel:     []uint64{21, 20, 19, 18},
		StCLevel:     []uint64{9, 9, 8},
		MaxN1N2Ratio: 16.0,
	},

	// Set VII
	// 1773 - 460
	{
		H:            32768,
		SinType:      Cos2,
		MessageRatio:        1024.0,
		SinRange:     325,
		SinDeg:       255,
		SinRescal:    4,
		ArcSineDeg:   0,
		CtSLevel:     []uint64{26, 25, 24, 23},
		StCLevel:     []uint64{11, 10, 10},
		MaxN1N2Ratio: 16.0,
	},

	// Set IV
	// 768 - 110
	{
		H:            192,
		SinType:      Cos1,
		MessageRatio:        1024.0,
		SinRange:     25,
		SinDeg:       63,
		SinRescal:    2,
		ArcSineDeg:   0,
		CtSLevel:     []uint64{13, 12},
		StCLevel:     []uint64{3, 3},
		MaxN1N2Ratio: 16.0,
	},
	*/
}
