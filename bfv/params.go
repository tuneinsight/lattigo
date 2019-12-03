package bfv

import (
	"errors"
	"fmt"
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
	"math"
	"math/bits"
)

// MaxLogN is the log2 of the largest supported polynomial modulus degree
const MaxLogN = 16

// MaxModuliCount is the largest supported number of 60 moduli in the RNS representation
const MaxModuliCount = 34

const MaxModuliSize = 60

// Plaintext modulus allowing batching for the corresponding N in ascending bitsize.
var tBatching = map[uint64][]uint64{
	4096: {40961, 114689, 188417, 417793, 1032193, 2056193, 4169729, 8380417, 16760833, 33538049, 67084289, 134176769,
		268369921, 536813569, 1073692673, 2147377153, 4294828033},
	8192: {65537, 114689, 163841, 1032193, 1785857, 4079617, 8273921, 16760833, 33538049, 67043329, 133857281,
		268369921, 536690689, 1073692673, 2147352577, 4294475777},
	16384: {65537, 163841, 786433, 1769473, 3735553, 8257537, 16580609, 33292289, 67043329, 133857281, 268369921,
		536641537, 1073643521, 2147352577, 4294475777},
	32768: {65537, 786433, 1769473, 3735553, 8257537, 16580609, 33292289, 67043329, 132710401, 268369921, 536608769,
		1073479681, 2147352577, 4293918721},
}

type Moduli struct {
	Qi    []uint64 // Ciphertext prime moduli
	Pi    []uint64 // Keys additional prime moduli
	QiMul []uint64 // Ciphertext secondary prime moduli
}

type LogModuli struct {
	LogQi    []uint64 // Ciphertext prime moduli bit-size
	LogPi    []uint64 // Keys additional prime moduli bit-size
	LogQiMul []uint64 // Ciphertext secondary prime moduli bit-size
}

// Parameters represents a given parameter set for the BFV cryptosystem.
type Parameters struct {
	Moduli
	LogModuli
	LogN  uint64 // Ring degree (power of 2)
	T     uint64 // Plaintext modulus
	LogQP uint64
	Sigma float64 // Gaussian sampling variance
	Alpha uint64
	Beta  uint64
}

func init() {
	for _, params := range DefaultParams {
		params.Gen()
	}
}

func (p *Parameters) Gen() {
	if p.Qi == nil && p.Pi == nil && p.QiMul == nil {

		p.Qi, p.Pi, p.QiMul = GenModuli(p)

		tmp := ring.NewUint(1)

		for _, qi := range p.Qi {
			tmp.Mul(tmp, ring.NewUint(qi))
		}

		for _, pi := range p.Pi {
			tmp.Mul(tmp, ring.NewUint(pi))
		}

		p.LogQP = uint64(tmp.BitLen())

	} else if p.LogQi == nil && p.LogPi == nil && p.LogQiMul == nil {

		p.LogQi = make([]uint64, len(p.Qi), len(p.Qi))
		for i := range p.Qi {
			p.LogQi[i] = uint64(bits.Len64(p.Qi[i]) - 1)
		}

		p.LogPi = make([]uint64, len(p.Pi), len(p.Pi))
		for i := range p.Pi {
			p.LogPi[i] = uint64(bits.Len64(p.Pi[i]) - 1)
		}

		p.LogQiMul = make([]uint64, len(p.QiMul), len(p.QiMul))
		for i := range p.QiMul {
			p.LogQiMul[i] = uint64(bits.Len64(p.QiMul[i]) - 1)
		}

	} else {
		panic("invalid parameters : must set correctly either Moduli or LogModuli")
	}

	p.Alpha = uint64(len(p.Pi))
	p.Beta = uint64(math.Ceil(float64(len(p.Qi)) / float64(len(p.Pi))))
}

// Copy creates a copy of the target parameters.
func (p *Parameters) Copy() (paramsCopy *Parameters) {

	paramsCopy = new(Parameters)
	paramsCopy.LogN = p.LogN
	paramsCopy.T = p.T
	paramsCopy.Sigma = p.Sigma

	paramsCopy.Qi = make([]uint64, len(p.Qi))
	copy(paramsCopy.Qi, p.Qi)

	paramsCopy.Pi = make([]uint64, len(p.Pi))
	copy(paramsCopy.Pi, p.Pi)

	paramsCopy.QiMul = make([]uint64, len(p.QiMul))
	copy(paramsCopy.QiMul, p.QiMul)

	paramsCopy.LogQi = make([]uint64, len(p.LogQi))
	copy(paramsCopy.LogQi, p.LogQi)

	paramsCopy.LogPi = make([]uint64, len(p.LogPi))
	copy(paramsCopy.LogPi, p.LogPi)

	paramsCopy.LogQiMul = make([]uint64, len(p.LogQiMul))
	copy(paramsCopy.LogQiMul, p.LogQiMul)

	paramsCopy.LogQP = p.LogQP
	paramsCopy.Alpha = p.Alpha
	paramsCopy.Beta = p.Beta

	return
}

// DefaultParams is a set of default BFV parameters ensuring 128 bit security.
var DefaultParams = map[uint64]*Parameters{

	//logQ1+P = 109
	12: {LogN: 12,
		T: 65537,
		LogModuli: LogModuli{
			LogQi:    []uint64{39, 39},
			LogPi:    []uint64{30},
			LogQiMul: []uint64{60, 60},
		},
		Sigma: 3.2},

	//logQ1+P = 218
	13: {LogN: 13,
		T: 65537,
		LogModuli: LogModuli{
			LogQi:    []uint64{54, 54, 54},
			LogPi:    []uint64{55},
			LogQiMul: []uint64{60, 60, 60},
		},
		Sigma: 3.2},

	//logQ1+P = 438
	14: {LogN: 14,
		T: 65537,
		LogModuli: LogModuli{
			LogQi:    []uint64{56, 55, 55, 54, 54, 54},
			LogPi:    []uint64{55, 55},
			LogQiMul: []uint64{60, 60, 60, 60, 60, 60},
		},
		Sigma: 3.2},

	//logQ1+P = 880
	15: {LogN: 15,
		T: 65537,
		LogModuli: LogModuli{
			LogQi:    []uint64{59, 59, 59, 58, 58, 58, 58, 58, 58, 58, 58, 58},
			LogPi:    []uint64{60, 60, 60},
			LogQiMul: []uint64{60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60},
		},
		Sigma: 3.2},
}

// Equals compares two sets of parameters for equality
func (p *Parameters) Equals(other *Parameters) bool {
	if p == other {
		return true
	}

	if p.Qi == nil && p.Pi == nil && p.QiMul == nil {
		return p.LogN == other.LogN && p.T == other.T && utils.EqualSliceUint64(p.LogQi, other.LogQi) && utils.EqualSliceUint64(p.LogPi, other.LogPi) && utils.EqualSliceUint64(p.LogQiMul, other.LogQiMul) && p.Sigma == other.Sigma
	}

	if p.LogQi == nil && p.LogPi == nil && p.LogQiMul == nil {
		return p.LogN == other.LogN && p.T == other.T && utils.EqualSliceUint64(p.Qi, other.Qi) && utils.EqualSliceUint64(p.Pi, other.Pi) && utils.EqualSliceUint64(p.QiMul, other.QiMul) && p.Sigma == other.Sigma
	}

	return false

}

// MarshalBinary returns a []byte representation of the parameter set
func (p *Parameters) MarshalBinary() ([]byte, error) {
	if p.LogN == 0 { // if N is 0, then p is the zero value
		return []byte{}, nil
	}
	b := utils.NewBuffer(make([]byte, 0, 4+8+8+len(p.LogQi)+len(p.LogPi)+len(p.LogQiMul)))

	if p.LogN > MaxLogN {
		return nil, fmt.Errorf("LogN is larger than %d", MaxLogN)
	}

	if len(p.LogQi) > MaxModuliCount {
		return nil, fmt.Errorf("len(LogQi) is larger than %d", MaxModuliCount)
	}

	if len(p.LogPi) > MaxModuliCount {
		return nil, fmt.Errorf("len(LogPi) is larger than %d", MaxModuliCount)
	}

	if len(p.LogQiMul) > MaxModuliCount {
		return nil, fmt.Errorf("len(LogQiMul) is larger than %d", MaxModuliCount)
	}

	b.WriteUint8(uint8(p.LogN))
	b.WriteUint8(uint8(len(p.LogQi)))
	b.WriteUint8(uint8(len(p.LogPi)))
	b.WriteUint8(uint8(len(p.LogQiMul)))
	b.WriteUint64(p.T)
	b.WriteUint64(uint64(p.Sigma * (1 << 32)))

	tmp := make([]uint8, len(p.LogQi)+len(p.LogPi)+len(p.LogQiMul), len(p.LogQi)+len(p.LogPi)+len(p.LogQiMul))

	for i, qi := range p.LogQi {
		if qi > MaxModuliSize {
			return nil, fmt.Errorf("LogQi n°%d is larger than %d", i, MaxModuliSize)
		}
		tmp[i] = uint8(qi)
	}

	for i, pi := range p.LogPi {
		if pi > MaxModuliSize {
			return nil, fmt.Errorf("LogQi n°%d is larger than %d", i, MaxModuliSize)
		}
		tmp[i+len(p.LogQi)] = uint8(pi)
	}

	for i, qi := range p.LogQiMul {
		if qi > MaxModuliSize {
			return nil, fmt.Errorf("LogQiMul n°%d is larger than %d", i, MaxModuliSize)
		}
		tmp[i+len(p.LogQi)+len(p.LogPi)] = uint8(qi)
	}

	b.WriteUint8Slice(tmp)
	return b.Bytes(), nil
}

// UnmarshalBinary decodes a []byte into a parameter set struct
func (p *Parameters) UnmarshalBinary(data []byte) error {
	if len(data) < 3 {
		return errors.New("invalid parameters encoding")
	}
	b := utils.NewBuffer(data)

	p.LogN = uint64(b.ReadUint8())
	if p.LogN > MaxLogN {
		return fmt.Errorf("LogN is larger than %d", MaxLogN)
	}

	lenLogQi := b.ReadUint8()
	if lenLogQi > MaxModuliCount {
		return fmt.Errorf("len(LogQi) is larger than %d", MaxModuliCount)
	}

	lenLogPi := b.ReadUint8()
	if lenLogPi > MaxModuliCount {
		return fmt.Errorf("len(LogPi) is larger than %d", MaxModuliCount)
	}

	lenLogQiMul := b.ReadUint8()
	if lenLogQiMul > MaxModuliCount {
		return fmt.Errorf("len(LogQiMul) is larger than %d", MaxModuliCount)
	}

	p.T = b.ReadUint64()
	p.Sigma = math.Round((float64(b.ReadUint64())/float64(1<<32))*100) / 100
	p.LogQi = make([]uint64, lenLogQi, lenLogQi)
	p.LogPi = make([]uint64, lenLogPi, lenLogPi)
	p.LogQiMul = make([]uint64, lenLogQiMul, lenLogQiMul)

	tmp := make([]uint8, lenLogQi+lenLogPi+lenLogQiMul, lenLogQi+lenLogPi+lenLogQiMul)

	b.ReadUint8Slice(tmp[:lenLogQi])
	b.ReadUint8Slice(tmp[lenLogQi : lenLogQi+lenLogPi])
	b.ReadUint8Slice(tmp[lenLogQi+lenLogPi:])

	var qi uint64

	for i := range p.LogQi {

		qi = uint64(tmp[i])

		if qi > MaxModuliSize {
			return fmt.Errorf("LogQi n°%d is larger than %d", i, MaxModuliSize)
		}

		p.LogQi[i] = qi
	}

	for i := range p.LogPi {

		qi = uint64(tmp[i+len(p.LogQi)])

		if qi > MaxModuliSize {
			return fmt.Errorf("LogPi n°%d is larger than %d", i, MaxModuliSize)
		}

		p.LogPi[i] = qi
	}

	for i := range p.LogQiMul {

		qi = uint64(tmp[i+len(p.LogQi)+len(p.LogPi)])

		if qi > MaxModuliSize {
			return fmt.Errorf("LogQiMul n°%d is larger than %d", i, MaxModuliSize)
		}

		p.LogQiMul[i] = qi
	}

	p.Gen()

	return nil
}
