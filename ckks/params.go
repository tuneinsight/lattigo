package ckks

import (
	"errors"
	"fmt"
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
	"math"
	"math/bits"
)

type Moduli struct {
	Q []uint64
	P []uint64
}

type LogModuli struct {
	LogQi []uint64
	LogPi []uint64
}

type Parameters struct {
	Moduli
	LogModuli
	LogN     uint64
	LogSlots uint64
	LogQP    uint64
	MaxLevel uint64
	Scale    float64
	Sigma    float64
	Alpha    uint64
	Beta     uint64
}

func init() {
	for _, params := range DefaultParams {
		params.Gen()
	}
}

func (p *Parameters) Gen() {
	if p.Q == nil && p.P == nil {

		p.Q, p.P = GenModuli(p)

		tmp := ring.NewUint(p.Q[0])

		for i := 1; i < len(p.Q); i++ {
			tmp.Mul(tmp, ring.NewUint(p.Q[i]))
		}

		for i := 0; i < len(p.P); i++ {
			tmp.Mul(tmp, ring.NewUint(p.P[i]))
		}

		p.LogQP = uint64(tmp.BitLen())

	} else if p.LogQi == nil && p.LogPi == nil {

		p.LogQi = make([]uint64, len(p.Q), len(p.Q))
		for i := range p.Q {
			p.LogQi[i] = uint64(bits.Len64(p.Q[i]) - 1)
		}

		p.LogPi = make([]uint64, len(p.P), len(p.P))
		for i := range p.P {
			p.LogPi[i] = uint64(bits.Len64(p.P[i]) - 1)
		}

	} else {
		panic("Invalid parameters : must set correctly either Moduli or LogModuli")
	}

	p.MaxLevel = uint64(len(p.Q) - 1)

	p.Alpha = uint64(len(p.P))
	p.Beta = uint64(math.Ceil(float64(len(p.Q)) / float64(len(p.P))))
}

// Copy creates a copy of the target parameters.
func (p *Parameters) Copy() (paramsCopy *Parameters) {

	paramsCopy = new(Parameters)
	paramsCopy.LogN = p.LogN
	paramsCopy.LogSlots = p.LogSlots
	paramsCopy.Scale = p.Scale
	paramsCopy.Sigma = p.Sigma

	paramsCopy.LogQi = make([]uint64, len(p.LogQi))
	copy(paramsCopy.LogQi, p.LogQi)

	paramsCopy.LogPi = make([]uint64, len(p.LogPi))
	copy(paramsCopy.LogPi, p.LogPi)

	paramsCopy.Q = make([]uint64, len(p.Q))
	copy(paramsCopy.Q, p.Q)

	paramsCopy.P = make([]uint64, len(p.P))
	copy(paramsCopy.P, p.P)

	paramsCopy.LogQP = p.LogQP
	paramsCopy.MaxLevel = p.MaxLevel

	paramsCopy.Alpha = p.Alpha
	paramsCopy.Beta = p.Beta

	return
}

// DefaultParams is a set of default CKKS parameters ensuring 128 bit security.
var DefaultParams = map[uint64]*Parameters{

	//LogQi = 109
	12: {LogN: 12,
		LogSlots: 11,
		LogModuli: LogModuli{
			LogQi: []uint64{37, 32},
			LogPi: []uint64{38},
		},
		Scale: 1 << 32,
		Sigma: 3.2},

	//LogQi = 218
	13: {LogN: 13,
		LogSlots: 12,
		LogModuli: LogModuli{
			LogQi: []uint64{33, 30, 30, 30, 30, 30},
			LogPi: []uint64{35},
		},
		Scale: 1 << 30,
		Sigma: 3.2},

	//LogQiP = 438
	14: {LogN: 14,
		LogSlots: 13,
		LogModuli: LogModuli{
			LogQi: []uint64{45, 34, 34, 34, 34, 34, 34, 34, 34, 34},
			LogPi: []uint64{43, 43},
		},
		Scale: 1 << 34,
		Sigma: 3.2},

	//LogQi = 880
	15: {LogN: 15,
		LogSlots: 14,
		LogModuli: LogModuli{
			LogQi: []uint64{50, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40},
			LogPi: []uint64{50, 50, 50},
		},
		Scale: 1 << 40,
		Sigma: 3.2},

	//LogQi = 1761
	16: {LogN: 16,
		LogSlots: 15,
		LogModuli: LogModuli{
			LogQi: []uint64{55, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
			LogPi: []uint64{55, 55, 55, 55},
		},
		Scale: 1 << 45,
		Sigma: 3.2},
}

// MaxLogN is the largest supported polynomial modulus degree
const MaxLogN = 16

// MaxModuliCount is the largest supported number of moduli in the RNS representation
const MaxModuliCount = 34

const MaxModuliSize = 60

// Equals compares two sets of parameters for equality
func (p *Parameters) Equals(other *Parameters) bool {
	if p == other {
		return true
	}
	return p.LogN == other.LogN && p.LogSlots == other.LogSlots && utils.EqualSliceUint64(p.LogQi, other.LogQi) && utils.EqualSliceUint64(p.LogPi, other.LogPi) && p.Scale == other.Scale && p.Sigma == other.Sigma
}

// MarshalBinary returns a []byte representation of the parameter set
func (p *Parameters) MarshalBinary() ([]byte, error) {

	if p.LogN == 0 { // if N is 0, then p is the zero value
		return []byte{}, nil
	}

	if p.LogN > MaxLogN {
		return nil, errors.New("polynomial degree is too large")
	}

	if p.LogSlots >= p.LogN {
		return nil, errors.New("slot number is too large")
	}

	tmp := make([]uint8, len(p.LogQi)+len(p.LogPi), len(p.LogQi)+len(p.LogPi))

	for i, qi := range p.LogQi {
		if qi > MaxModuliSize {
			return nil, errors.New("invalid LogQi, values must be smaller or equal to 60 bits")
		}
		tmp[i] = uint8(qi)
	}

	for i, pi := range p.LogPi {
		if pi > MaxModuliSize {
			return nil, errors.New("invalid LogPi, values must be smaller or equal to 60 bits")
		}
		tmp[i+len(p.LogQi)] = uint8(pi)
	}

	b := utils.NewBuffer(make([]byte, 0, 20+len(p.LogQi)+len(p.LogPi)))
	b.WriteUint8(uint8(p.LogN))
	b.WriteUint8(uint8(p.LogSlots))
	b.WriteUint8(uint8(len(p.LogQi)))
	b.WriteUint8(uint8(len(p.LogPi)))
	b.WriteUint64(uint64(p.Scale))
	b.WriteUint64(uint64(p.Sigma * (1 << 32)))
	b.WriteUint8Slice(tmp)

	return b.Bytes(), nil
}

// UnMarshalBinary decodes a []byte into a parameter set struct
func (p *Parameters) UnMarshalBinary(data []byte) error {
	if len(data) < 3 {
		return errors.New("invalid parameters encoding")
	}
	b := utils.NewBuffer(data)

	p.LogN = uint64(b.ReadUint8())
	if p.LogN > MaxLogN {
		return errors.New("polynomial degree is too large")
	}

	p.LogSlots = uint64(b.ReadUint8())
	if p.LogSlots >= p.LogN {
		return errors.New("slot number is too large")
	}

	lenLogQi := uint64(b.ReadUint8())
	if lenLogQi > MaxModuliCount {
		return fmt.Errorf("len(LogQi) is larger than %d", MaxModuliCount)
	}

	lenLogPi := uint64(b.ReadUint8())
	if lenLogPi > MaxModuliCount {
		return fmt.Errorf("len(LogPi) is larger than %d", MaxModuliCount)
	}

	p.Scale = float64(b.ReadUint64())

	p.Sigma = math.Round((float64(b.ReadUint64())/float64(1<<32))*100) / 100

	tmp := make([]uint8, lenLogQi+lenLogPi, lenLogQi+lenLogPi)

	b.ReadUint8Slice(tmp[:lenLogQi])
	b.ReadUint8Slice(tmp[lenLogQi:])

	p.LogQi = make([]uint64, lenLogQi, lenLogQi)
	p.LogPi = make([]uint64, lenLogPi, lenLogPi)

	for i := range p.LogQi {
		p.LogQi[i] = uint64(tmp[i])
	}

	for i := range p.LogPi {
		p.LogPi[i] = uint64(tmp[i+len(p.LogQi)])
	}

	p.Gen()

	return nil
}
