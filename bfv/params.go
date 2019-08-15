package bfv

import (
	"errors"
	"fmt"
	"github.com/lca1/lattigo/ring"
	"github.com/lca1/lattigo/utils"
	"math"
	"math/bits"
)

// MaxN is the largest supported polynomial modulus degree
const MaxN = 1 << 15

// MaxModuliCount is the largest supported number of 60 moduli in the RNS representation
const MaxModuliCount = 256

// Power of 2 plaintext modulus
var Tpow2 = []uint64{2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144,
	524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912,
	1073741824, 2147483648, 4294967296}

// Plaintext modulus allowing batching for the corresponding N in ascending bitsize.
var TBatching = map[uint64][]uint64{
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
	N     uint64
	T     uint64
	Qi    []uint64
	Pi    []uint64
	Sigma float64
}

// DefaultParams is an array default parameters with increasing homomorphic capacity.
// These parameters correspond to 128 bit security level for secret keys in the ternary distribution
// (see //https://projects.csail.mit.edu/HEWorkshop/HomomorphicEncryptionStandard2018.pdf).
var DefaultParams = []Parameters{
	{4096, 65537, ring.Qi60[len(ring.Qi60)-2:], ring.Pi60[len(ring.Pi60)-3:], 3.19},
	{8192, 65537, ring.Qi60[len(ring.Qi60)-4:], ring.Pi60[len(ring.Pi60)-5:], 3.19},
	{16384, 65537, ring.Qi60[len(ring.Qi60)-8:], ring.Pi60[len(ring.Pi60)-9:], 3.19},
	{32768, 65537, ring.Qi60[len(ring.Qi60)-16:], ring.Pi60[len(ring.Pi60)-17:], 3.19},
}

// Equals compares two sets of parameters for equality
func (p *Parameters) Equals(other *Parameters) bool {
	if p == other {
		return true
	}
	return p.N == other.N && EqualSlice(p.Qi, other.Qi) && EqualSlice(p.Pi, other.Pi) && p.Sigma == other.Sigma
}

// MarshalBinary returns a []byte representation of the parameter set
func (p *Parameters) MarshalBinary() ([]byte, error) {
	if p.N == 0 { // if N is 0, then p is the zero value
		return []byte{}, nil
	}
	b := utils.NewBuffer(make([]byte, 0, 3+((2+len(p.Qi)+len(p.Pi))<<3)))
	b.WriteUint8(uint8(bits.Len64(p.N) - 1))
	b.WriteUint8(uint8(len(p.Qi)))
	b.WriteUint8(uint8(len(p.Pi)))
	b.WriteUint64(p.T)
	b.WriteUint64(uint64(p.Sigma * (1 << 32)))
	b.WriteUint64Slice(p.Qi)
	b.WriteUint64Slice(p.Pi)
	return b.Bytes(), nil
}

// UnMarshalBinary decodes a []byte into a parameter set struct
func (p *Parameters) UnMarshalBinary(data []byte) error {
	if len(data) < 3 {
		return errors.New("invalid parameters encoding")
	}
	b := utils.NewBuffer(data)
	p.N = 1 << uint64(b.ReadUint8())
	if p.N > MaxN {
		return errors.New("polynomial degree is too large")
	}
	lenQi := uint64(b.ReadUint8())
	if lenQi > MaxModuliCount {
		return fmt.Errorf("len(Qi) is larger than %d", MaxModuliCount)
	}
	lenPi := uint64(b.ReadUint8())
	if lenPi > MaxModuliCount {
		return fmt.Errorf("len(Pi) is larger than %d", MaxModuliCount)
	}
	p.T = b.ReadUint64()
	p.Sigma = math.Round((float64(b.ReadUint64())/float64(1<<32))*100) / 100
	p.Qi = make([]uint64, lenQi, lenQi)
	p.Pi = make([]uint64, lenPi, lenPi)
	b.ReadUint64Slice(p.Qi)
	b.ReadUint64Slice(p.Pi)
	return nil
}
