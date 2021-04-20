package ckks

import (
	"errors"
	"fmt"
	"math"
	"math/big"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// Name of the different default parameter sets
const (
	// PN12QP109 is the index in DefaultParams for logQP = 109
	PN12QP109 = iota
	// PN13QP218 is the index in DefaultParams for logQP = 218
	PN13QP218
	// PN14QP438 is the index in DefaultParams for logQP = 438
	PN14QP438
	// PN15QP880 is the index in DefaultParams for logQP = 880
	PN15QP880
	// PN16QP1761 is the index in DefaultParams for logQP = 1761
	PN16QP1761

	// PN12QP101pq is the index in DefaultParams for logQP = 101 (post quantum)
	PN12QP101pq
	// PN13QP202pq is the index in DefaultParams for logQP = 202 (post quantum)
	PN13QP202pq
	// PN14QP411pq is the index in DefaultParams for logQP = 411 (post quantum)
	PN14QP411pq
	// PN15QP827pq is the index in DefaultParams for logQP = 827 (post quantum)
	PN15QP827pq
	// PN16QP1654pq is the index in DefaultParams for logQP = 1654 (post quantum)
	PN16QP1654pq
)

// DefaultParams is a set of default CKKS parameters ensuring 128 bit security.
var DefaultParams = []*ParametersDef{

	//LogQi = 109
	{LogN: 12,
		LogSlots: 11,
		Q: []uint64{0x200000e001, // 37 + 32
			0x100006001},
		P:     []uint64{0x3ffffea001}, // 38
		Scale: 1 << 32,
		Sigma: rlwe.DefaultSigma,
	},

	//LogQi = 218
	{LogN: 13,
		LogSlots: 12,
		Q: []uint64{0x1fffec001, // 33 + 5 x 30
			0x3fff4001,
			0x3ffe8001,
			0x40020001,
			0x40038001,
			0x3ffc0001},
		P:     []uint64{0x800004001}, // 35
		Scale: 1 << 30,
		Sigma: rlwe.DefaultSigma,
	},

	//LogQiP = 438
	{LogN: 14,
		LogSlots: 13,
		Q: []uint64{0x200000008001, 0x400018001, // 45 + 9 x 34
			0x3fffd0001, 0x400060001,
			0x400068001, 0x3fff90001,
			0x400080001, 0x4000a8001,
			0x400108001, 0x3ffeb8001},
		P:     []uint64{0x7fffffd8001, 0x7fffffc8001}, // 43, 43
		Scale: 1 << 34,
		Sigma: rlwe.DefaultSigma,
	},

	//LogQi = 880
	{LogN: 15,
		LogSlots: 14,
		Q: []uint64{0x4000000120001, 0x10000140001, 0xffffe80001, // 50 + 17 x 40
			0x10000290001, 0xffffc40001, 0x100003e0001,
			0x10000470001, 0x100004b0001, 0xffffb20001,
			0x10000500001, 0x10000650001, 0xffff940001,
			0xffff8a0001, 0xffff820001, 0xffff780001,
			0x10000890001, 0xffff750001, 0x10000960001},
		P:     []uint64{0x40000001b0001, 0x3ffffffdf0001, 0x4000000270001}, // 50, 50, 50
		Scale: 1 << 40,
		Sigma: rlwe.DefaultSigma,
	},

	//LogQi = 1761
	{LogN: 16,
		LogSlots: 15,
		Q: []uint64{0x80000000080001, 0x2000000a0001, 0x2000000e0001, 0x1fffffc20001, // 55 + 33 x 45
			0x200000440001, 0x200000500001, 0x200000620001, 0x1fffff980001,
			0x2000006a0001, 0x1fffff7e0001, 0x200000860001, 0x200000a60001,
			0x200000aa0001, 0x200000b20001, 0x200000c80001, 0x1fffff360001,
			0x200000e20001, 0x1fffff060001, 0x200000fe0001, 0x1ffffede0001,
			0x1ffffeca0001, 0x1ffffeb40001, 0x200001520001, 0x1ffffe760001,
			0x2000019a0001, 0x1ffffe640001, 0x200001a00001, 0x1ffffe520001,
			0x200001e80001, 0x1ffffe0c0001, 0x1ffffdee0001, 0x200002480001,
			0x1ffffdb60001, 0x200002560001},
		P:     []uint64{0x80000000440001, 0x7fffffffba0001, 0x80000000500001, 0x7fffffffaa0001}, // 4 x 55
		Scale: 1 << 45,
		Sigma: rlwe.DefaultSigma,
	},

	//LogQi = 101.00001186816735
	{LogN: 12,
		LogSlots: 11,
		Q:        []uint64{0x800004001, 0x40002001}, // 35 + 30
		P:        []uint64{0x1000002001},            // 36
		Scale:    1 << 30,
		Sigma:    rlwe.DefaultSigma,
	},

	//LogQi = 201.9936341352857
	{LogN: 13,
		LogSlots: 12,
		Q:        []uint64{0x1fffec001, 0x8008001, 0x8020001, 0x802c001, 0x7fa8001, 0x7f74001}, // 33 + 5 x 27
		P:        []uint64{0x400018001},                                                        // 34
		Scale:    1 << 27,
		Sigma:    rlwe.DefaultSigma,
	},

	//LogQiP = 411.0000787495673
	{LogN: 14,
		LogSlots: 13,
		Q: []uint64{0x10000048001, 0x200038001, 0x1fff90001, 0x200080001, 0x1fff60001,
			0x2000b8001, 0x200100001, 0x1fff00001, 0x1ffef0001, 0x200128001}, // 40 + 9 x 33

		P:     []uint64{0x1ffffe0001, 0x1ffffc0001}, // 37, 37
		Scale: 1 << 33,
		Sigma: rlwe.DefaultSigma,
	},

	//LogQi = 827.0000771955918
	{LogN: 15,
		LogSlots: 14,
		Q: []uint64{0x400000060001, 0x4000170001, 0x3fffe80001, 0x40002f0001, 0x4000300001,
			0x3fffcf0001, 0x40003f0001, 0x3fffc10001, 0x4000450001, 0x3fffb80001,
			0x3fffb70001, 0x40004a0001, 0x3fffb20001, 0x4000510001, 0x3fffaf0001,
			0x4000540001, 0x4000560001, 0x4000590001}, // 46 + 17 x 38
		P:     []uint64{0x2000000a0001, 0x2000000e0001, 0x2000001d0001}, // 3 x 45
		Scale: 1 << 38,
		Sigma: rlwe.DefaultSigma,
	},

	//LogQi = 1653.999999
	{LogN: 16,
		LogSlots: 15,
		Q: []uint64{0x80000000080001, 0x2000000a0001, 0x2000000e0001, 0x1fffffc20001, 0x200000440001,
			0x200000500001, 0x200000620001, 0x1fffff980001, 0x2000006a0001, 0x1fffff7e0001,
			0x200000860001, 0x200000a60001, 0x200000aa0001, 0x200000b20001, 0x200000c80001,
			0x1fffff360001, 0x200000e20001, 0x1fffff060001, 0x200000fe0001, 0x1ffffede0001,
			0x1ffffeca0001, 0x1ffffeb40001, 0x200001520001, 0x1ffffe760001, 0x2000019a0001,
			0x1ffffe640001, 0x200001a00001, 0x1ffffe520001, 0x200001e80001, 0x1ffffe0c0001,
			0x1ffffdee0001, 0x200002480001}, // 55 + 31 x 45
		P:     []uint64{0x7fffffffe0001, 0x80000001c0001, 0x80000002c0001, 0x7ffffffd20001}, // 4 x 51
		Scale: 1 << 45,
		Sigma: rlwe.DefaultSigma,
	},
}

type Parameters interface {
	rlwe.Parameters

	LogQLvl(level uint64) uint64
	LogSlots() uint64
	MaxLevel() uint64
	MaxLogSlots() uint64
	QLvl(level uint64) *big.Int
	Scale() float64
	Slots() uint64
}

type ParametersDef struct {
	Q        []uint64
	P        []uint64
	LogN     uint64 // Ring degree (power of 2)
	LogSlots uint64
	Scale    float64
	Sigma    float64 // Gaussian sampling variance
}

// parameters represents a given parameter set for the CKKS cryptosystem.
type parameters struct {
	rlwe.Parameters

	logSlots uint64
	scale    float64
}

func NewParametersFromParamDef(paramDef *ParametersDef) (*parameters, error) {
	m := new(rlwe.Moduli)
	m.Qi = make([]uint64, len(paramDef.Q))
	copy(m.Qi, paramDef.Q)
	m.Pi = make([]uint64, len(paramDef.P))
	copy(m.Pi, paramDef.P)
	return NewParametersFromModuli(paramDef.LogN, m, paramDef.Sigma, paramDef.LogSlots, paramDef.Scale)
}

// NewParametersFromModuli creates a new Parameters struct and returns a pointer to it.
func NewParametersFromModuli(logN uint64, m *rlwe.Moduli, sigma float64, logSlots uint64, scale float64) (p *parameters, err error) {

	var rlweParams *rlwe.ParametersStruct
	if rlweParams, err = rlwe.NewRLWEParameters(logN, m.Qi, m.Pi, sigma); err != nil {
		return nil, err
	}

	return &parameters{rlweParams, logSlots, scale}, nil

}

// NewParametersFromLogModuli creates a new Parameters struct and returns a pointer to it.
func NewParametersFromLogModuli(logN uint64, lm *rlwe.LogModuli, sigma float64, logSlots uint64, scale float64) (p *parameters, err error) {

	if err = rlwe.CheckLogModuli(lm); err != nil {
		return nil, err
	}

	// If LogModuli is valid and then generates the moduli
	return NewParametersFromModuli(logN, rlwe.GenModuli(lm, logN), sigma, logSlots, scale)
}

func GetDefaultParameters(paramsId int) *parameters {
	if paramsId >= len(DefaultParams) {
		panic(fmt.Errorf("paramsId %d does not exist", paramsId))
	}
	params, err := NewParametersFromParamDef(DefaultParams[paramsId])
	if err != nil {
		panic(err)
	}
	return params
}

// LogSlots returns the log of the number of slots
func (p *parameters) LogSlots() uint64 {
	return p.logSlots
}

// MaxLevel returns the maximum ciphertext level
func (p *parameters) MaxLevel() uint64 {
	return p.QCount() - 1
}

// Slots returns number of available plaintext slots
func (p *parameters) Slots() uint64 {
	return 1 << p.logSlots
}

// MaxSlots returns the theoretical maximum of plaintext slots allowed by the ring degree
func (p *parameters) MaxSlots() uint64 {
	return p.N() >> 1
}

// MaxLogSlots returns the log of the maximum number of slots enabled by the parameters
func (p *parameters) MaxLogSlots() uint64 {
	return p.LogN() - 1
}

// Scale returns the default plaintext/ciphertext scale
func (p *parameters) Scale() float64 {
	return p.scale
}

// LogQLvl returns the size of the modulus Q in bits at a specific level
func (p *parameters) LogQLvl(level uint64) uint64 {
	tmp := p.QLvl(level)
	return uint64(tmp.BitLen())
}

// QLvl returns the product of the moduli at the given level as a big.Int
func (p *parameters) QLvl(level uint64) *big.Int {
	tmp := ring.NewUint(1)
	for _, qi := range p.Q()[:level+1] {
		tmp.Mul(tmp, ring.NewUint(qi))
	}
	return tmp
}

// // Copy creates a copy of the target parameters.
// func (p *parameters) Copy() (paramsCopy *parameters) {

// 	paramsCopy = new(parameters)
// 	paramsCopy.Parameters = p.Parameters.Copy()
// 	paramsCopy.logSlots = p.logSlots
// 	paramsCopy.scale = p.scale
// 	paramsCopy.sigma = p.sigma
// 	paramsCopy.qi = make([]uint64, len(p.qi))
// 	copy(paramsCopy.qi, p.qi)
// 	paramsCopy.pi = make([]uint64, len(p.pi))
// 	copy(paramsCopy.pi, p.pi)
// 	return
// }

// Equals compares two sets of parameters for equality.
func (p *parameters) Equals(other Parameters) (res bool) {

	if p == other {
		return true
	}

	res = p.LogN() == other.LogN()
	res = res && (p.logSlots == other.LogSlots())
	res = res && (p.scale == other.Scale())
	res = res && (p.Sigma() == other.Sigma())
	res = res && utils.EqualSliceUint64(p.Q(), other.Q())
	res = res && utils.EqualSliceUint64(p.P(), other.P())
	return
}

// MarshalBinary returns a []byte representation of the parameter set.
func (p *parameters) MarshalBinary() ([]byte, error) {
	if p.LogN() == 0 { // if N is 0, then p is the zero value
		return []byte{}, nil
	}

	// Data 21 byte + QPiCount * 8 byte:
	// 1 byte : logN
	// 1 byte : logSlots
	// 8 byte : scale
	// 8 byte : sigma
	// 1 byte : #qi
	// 1 byte : #pi
	b := utils.NewBuffer(make([]byte, 0, 20+(p.QPCount())<<3))

	b.WriteUint8(uint8(p.LogN()))
	b.WriteUint8(uint8(p.logSlots))
	b.WriteUint64(math.Float64bits(p.scale))
	b.WriteUint64(math.Float64bits(p.Sigma()))
	b.WriteUint8(uint8(p.QCount()))
	b.WriteUint8(uint8(p.PCount()))
	b.WriteUint64Slice(p.Q())
	b.WriteUint64Slice(p.P())

	return b.Bytes(), nil
}

// UnmarshalBinary decodes a []byte into a parameter set struct
func (p *parameters) UnmarshalBinary(data []byte) (err error) {

	if len(data) < 20 {
		return errors.New("invalid parameters encoding")
	}

	b := utils.NewBuffer(data)

	logN := uint64(b.ReadUint8())
	p.logSlots = uint64(b.ReadUint8())

	if p.logSlots > logN-1 {
		return fmt.Errorf("LogSlots larger than %d", rlwe.MaxLogN-1)
	}

	p.scale = math.Float64frombits(b.ReadUint64())
	sigma := math.Float64frombits(b.ReadUint64())

	lenQi := b.ReadUint8()
	lenPi := b.ReadUint8()

	qi := make([]uint64, lenQi)
	pi := make([]uint64, lenPi)

	b.ReadUint64Slice(qi)
	b.ReadUint64Slice(pi)

	rlweParams, err := rlwe.NewRLWEParameters(logN, qi, pi, sigma)
	if err != nil {
		return err
	}
	p.Parameters = rlweParams

	return nil
}
