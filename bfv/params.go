package bfv

import (
	"errors"
	"fmt"
	"github.com/ldsec/lattigo/utils"
	"math"
)

// MaxN is the largest supported polynomial modulus degree
const MaxLogN = 16

// MaxModuliCount is the largest supported number of 60 moduli in the RNS representation
const MaxModuliCount = 34

// Modulies for 128 security according to http://homomorphicencryption.org/white_papers/security_homomorphic_encryption_white_paper.pdf

// Power of 2 plaintext modulus
var tPow2 = []uint64{2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144,
	524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912,
	1073741824, 2147483648, 4294967296}

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

// Parameters represents a given parameter set for the BFV cryptosystem.
type Parameters struct {
	LogN  uint8
	T     uint64
	Q1    []uint8
	P     []uint8
	Q2    []uint8
	Sigma float64
}

// Copy creates a copy of the target parameters.
func (p *Parameters) Copy() (paramsCopy *Parameters) {

	paramsCopy = new(Parameters)
	paramsCopy.LogN = p.LogN
	paramsCopy.T = p.T

	paramsCopy.Q1 = make([]uint8, len(p.Q1))
	for i := range p.Q1 {
		paramsCopy.Q1[i] = p.Q1[i]
	}

	paramsCopy.P = make([]uint8, len(p.P))
	for i := range p.P {
		paramsCopy.P[i] = p.P[i]
	}

	paramsCopy.Q2 = make([]uint8, len(p.Q2))
	for i := range p.Q2 {
		paramsCopy.Q2[i] = p.Q2[i]
	}

	paramsCopy.Sigma = p.Sigma

	return
}

// DefaultParams is a set of default BFV parameters ensuring 128 bit security.
var DefaultParams = map[uint64]*Parameters{

	//logQ = 109
	12: {LogN: 12,
		T:     65537,
		Q1:    []uint8{39, 39},
		P:     []uint8{30},
		Q2:    []uint8{60, 60},
		Sigma: 3.2},

	//logQ = 218
	13: {LogN: 13,
		T:     65537,
		Q1:    []uint8{54, 54, 54},
		P:     []uint8{55},
		Q2:    []uint8{60, 60, 60},
		Sigma: 3.2},

	//logQ = 438
	14: {LogN: 14,
		T:     65537,
		Q1:    []uint8{56, 55, 55, 54, 54, 54},
		P:     []uint8{55, 55},
		Q2:    []uint8{60, 60, 60, 60, 60, 60},
		Sigma: 3.2},

	//logQ = 880
	15: {LogN: 15,
		T:     65537,
		Q1:    []uint8{59, 59, 59, 58, 58, 58, 58, 58, 58, 58, 58, 58},
		P:     []uint8{60, 60, 60},
		Q2:    []uint8{60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60},
		Sigma: 3.2},
}

// Equals compares two sets of parameters for equality
func (p *Parameters) Equals(other *Parameters) bool {
	if p == other {
		return true
	}
	return p.LogN == other.LogN && p.T == other.T && utils.EqualSliceUint8(p.Q1, other.Q1) && utils.EqualSliceUint8(p.P, other.P) && utils.EqualSliceUint8(p.Q2, other.Q2) && p.Sigma == other.Sigma
}

// MarshalBinary returns a []byte representation of the parameter set
func (p *Parameters) MarshalBinary() ([]byte, error) {
	if p.LogN == 0 { // if N is 0, then p is the zero value
		return []byte{}, nil
	}
	b := utils.NewBuffer(make([]byte, 0, 4+((3+len(p.Q1)+len(p.Q2)+len(p.P))<<3)))
	b.WriteUint8(uint8(p.LogN))
	b.WriteUint8(uint8(len(p.Q1)))
	b.WriteUint8(uint8(len(p.P)))
	b.WriteUint8(uint8(len(p.Q2)))
	b.WriteUint64(p.T)
	b.WriteUint64(uint64(p.Sigma * (1 << 32)))
	b.WriteUint8Slice(p.Q1)
	b.WriteUint8Slice(p.P)
	b.WriteUint8Slice(p.Q2)
	return b.Bytes(), nil
}

// UnmarshalBinary decodes a []byte into a parameter set struct
func (p *Parameters) UnmarshalBinary(data []byte) error {
	if len(data) < 3 {
		return errors.New("invalid parameters encoding")
	}
	b := utils.NewBuffer(data)

	p.LogN = b.ReadUint8()
	if p.LogN > MaxLogN {
		return errors.New("polynomial degree is too large")
	}

	lenQ1 := b.ReadUint8()
	if lenQ1 > MaxModuliCount {
		return fmt.Errorf("len(Q1) is larger than %d", MaxModuliCount)
	}

	lenP := b.ReadUint8()
	if lenP > MaxModuliCount {
		return fmt.Errorf("len(lenP) is larger than %d", MaxModuliCount)
	}

	lenQ2 := b.ReadUint8()
	if lenQ2 > MaxModuliCount {
		return fmt.Errorf("len(Q2) is larger than %d", MaxModuliCount)
	}

	p.T = b.ReadUint64()
	p.Sigma = math.Round((float64(b.ReadUint64())/float64(1<<32))*100) / 100
	p.Q1 = make([]uint8, lenQ1, lenQ1)
	p.P = make([]uint8, lenP, lenP)
	p.Q2 = make([]uint8, lenQ2, lenQ2)

	b.ReadUint8Slice(p.Q1)
	b.ReadUint8Slice(p.P)
	b.ReadUint8Slice(p.Q2)
	return nil
}
