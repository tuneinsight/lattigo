package bfv

import (
	"errors"
	"fmt"
	"github.com/ldsec/lattigo/utils"
	"math"
	"math/bits"
)

// MaxN is the largest supported polynomial modulus degree
const MaxN = 1 << 16

// MaxModuliCount is the largest supported number of 60 moduli in the RNS representation
const MaxModuliCount = 34

// Modulies for 128 security according to http://homomorphicencryption.org/white_papers/security_homomorphic_encryption_white_paper.pdf

var logN13Q218 = []uint64{0x7ffffffffb4001, 0x3fffffffef8001, 0x3fffffffeb8001}

var logN14Q438 = []uint64{0x7fffffffe90001, 0x7fffffffd58001, 0x7fffffffbf0001, 0x7fffffffbd0001,
	0x7fffffffba0001, 0x7fffffffb58001}

var logN15Q881 = []uint64{0x7ffffffffe70001, 0x7ffffffffe10001, 0x7ffffffffcc0001, 0x7ffffffffba0001,
	0x7ffffffffb00001, 0x7ffffffff630001, 0x7ffffffff510001, 0x7ffffffff3f0001,
	0x7ffffffff350001, 0x7ffffffff320001, 0x7ffffffff2c0001, 0x7fffffffe90001}

var logN16Q1770 = []uint64{0x7ffffffffcc0001, 0x7ffffffffba0001, 0x7ffffffffb00001, 0x7ffffffff320001,
	0x7ffffffff2c0001, 0x7ffffffff240001, 0x7fffffffefa0001, 0x7fffffffede0001,
	0x7fffffffe900001, 0x7fffffffe3c0001, 0x7fffffffe240001, 0x7fffffffddc0001,
	0x7fffffffdbe0001, 0x7fffffffd740001, 0x7fffffffd640001, 0x7fffffffd1a0001,
	0x7fffffffd0a0001, 0x7fffffffd080001, 0x7fffffffcda0001, 0x7fffffffccc0001,
	0x7fffffffcbc0001, 0x7fffffffcae0001, 0x7fffffffc980001, 0x7fffffffc480001,
	0x7fffffffc020001, 0x7fffffffbcc0001}

var Pi60 = []uint64{0xffffffffe400001, 0xffffffffd000001, 0xffffffffa200001, 0xffffffff9600001,
	0xfffffffeb200001, 0xfffffffea400001, 0xfffffffe8000001, 0xfffffffe3e00001,
	0xfffffffe2200001, 0xfffffffe0800001, 0xfffffffdd400001, 0xfffffffd9000001,
	0xfffffffcea00001, 0xfffffffcdc00001, 0xfffffffc7200001, 0xfffffffc5a00001,
	0xfffffffc5400001, 0xfffffffc4200001, 0xfffffffc2e00001, 0xfffffffbfa00001,
	0xfffffffbf200001, 0xfffffffbce00001, 0xfffffffba400001, 0xfffffffba000001,
	0xfffffffb8c00001, 0xfffffffb1400001, 0xfffffffafc00001, 0xfffffffaf800001,
	0xfffffffa6600001, 0xfffffffa5000001, 0xfffffff9ee00001, 0xfffffff9d600001,
	0xfffffff9ba00001, 0xfffffff99a00001, 0xfffffff94800001, 0xfffffff91000001,
	0xfffffff90600001, 0xfffffff8e600001, 0xfffffff8a400001, 0xfffffff88200001}

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
	N               uint64
	T               uint64
	Qi              []uint64
	KeySwitchPrimes []uint64
	Pi              []uint64
	Sigma           float64
}

// DefaultParams is an array default parameters with increasing homomorphic capacity.
// These parameters correspond to 128 bit security level for secret keys in the ternary distribution
// (see //https://projects.csail.mit.edu/HEWorkshop/HomomorphicEncryptionStandard2018.pdf).
var DefaultParams = []Parameters{
	{8192, 65537, logN13Q218, []uint64{0x7fffffffeac001}, Pi60[len(Pi60)-len(logN13Q218):], 3.19},
	{16384, 65537, logN14Q438, []uint64{0x3fffffffeb8001, 0x3fffffffef8001}, Pi60[len(Pi60)-len(logN14Q438):], 3.19},
	{32768, 65537, logN15Q881, []uint64{0x7ffffffff240001, 0x7ffffffff230001, 0x7fffffffefa0001}, Pi60[len(Pi60)-len(logN15Q881):], 3.19},
}

// Equals compares two sets of parameters for equality
func (p *Parameters) Equals(other *Parameters) bool {
	if p == other {
		return true
	}
	return p.N == other.N && utils.EqualSliceUint64(p.Qi, other.Qi) && utils.EqualSliceUint64(p.Pi, other.Pi) && p.Sigma == other.Sigma
}

// MarshalBinary returns a []byte representation of the parameter set
func (p *Parameters) MarshalBinary() ([]byte, error) {
	if p.N == 0 { // if N is 0, then p is the zero value
		return []byte{}, nil
	}
	b := utils.NewBuffer(make([]byte, 0, 4+((3+len(p.Qi)+len(p.Pi)+len(p.KeySwitchPrimes))<<3)))
	b.WriteUint8(uint8(bits.Len64(p.N) - 1))
	b.WriteUint8(uint8(len(p.Qi)))
	b.WriteUint8(uint8(len(p.KeySwitchPrimes)))
	b.WriteUint8(uint8(len(p.Pi)))
	b.WriteUint64(p.T)
	b.WriteUint64(uint64(p.Sigma * (1 << 32)))
	b.WriteUint64Slice(p.Qi)
	b.WriteUint64Slice(p.KeySwitchPrimes)
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

	lenKeySwitchPrimes := uint64(b.ReadUint8())
	if lenKeySwitchPrimes > MaxModuliCount {
		return fmt.Errorf("len(lenKeySwitchPrimes) is larger than %d", MaxModuliCount)
	}

	lenPi := uint64(b.ReadUint8())
	if lenPi > MaxModuliCount {
		return fmt.Errorf("len(Pi) is larger than %d", MaxModuliCount)
	}

	p.T = b.ReadUint64()
	p.Sigma = math.Round((float64(b.ReadUint64())/float64(1<<32))*100) / 100
	p.Qi = make([]uint64, lenQi, lenQi)
	p.KeySwitchPrimes = make([]uint64, lenKeySwitchPrimes, lenKeySwitchPrimes)
	p.Pi = make([]uint64, lenPi, lenPi)

	b.ReadUint64Slice(p.Qi)
	b.ReadUint64Slice(p.KeySwitchPrimes)
	b.ReadUint64Slice(p.Pi)
	return nil
}
