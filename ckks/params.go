package ckks

import (
	"errors"
	"fmt"
	"github.com/ldsec/lattigo/utils"
	"math"
)

// Parameters is a struct storing the necessary parameters to instantiate a new ckkscontext.
//
// logN : the ring degree such that N = 2^logN
//
// Q : a list of moduli in bit size, such that prod(moduli) <= logQ, where logQ is the maximum bit size of the product of all the moduli
// for a given security parameter. This moduli chain allows to customize the value by which the ciphertexts will be rescaled allong the homomorphic computations,
// and can enhance performances for optimized applications.
//
// logScale : the default scale of the scheme such that scale = 2^logScale
//
// sigma : the variance used by the ckkscontext to sample gaussian polynomials.
type Parameters struct {
	LogN  uint8
	Q     []uint8
	P     []uint8
	Scale float64
	Sigma float64
}

// Copy creates a copy of the target parameters.
func (p *Parameters) Copy() (paramsCopy *Parameters) {

	paramsCopy = new(Parameters)
	paramsCopy.LogN = p.LogN

	paramsCopy.Q = make([]uint8, len(p.Q))
	for i := range p.Q {
		paramsCopy.Q[i] = p.Q[i]
	}

	paramsCopy.P = make([]uint8, len(p.P))
	for i := range p.P {
		paramsCopy.P[i] = p.P[i]
	}

	paramsCopy.Scale = p.Scale
	paramsCopy.Sigma = p.Sigma

	return
}

// DefaultParams is a set of default CKKS parameters ensuring 128 bit security.
var DefaultParams = map[uint64]*Parameters{

	//logQ = 109
	12: {LogN: 12,
		Q:     []uint8{37, 32},
		P:     []uint8{38},
		Scale: 1 << 32,
		Sigma: 3.2},

	//logQ = 218
	13: {LogN: 13,
		Q:     []uint8{33, 30, 30, 30, 30, 30},
		P:     []uint8{35},
		Scale: 1 << 30,
		Sigma: 3.2},

	//logQP = 438
	14: {LogN: 14,
		Q:     []uint8{45, 34, 34, 34, 34, 34, 34, 34, 34, 34},
		P:     []uint8{43, 43},
		Scale: 1 << 34,
		Sigma: 3.2},

	//logQ = 880
	15: {LogN: 15,
		Q:     []uint8{50, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40},
		P:     []uint8{50, 50, 50},
		Scale: 1 << 40,
		Sigma: 3.2},

	//logQ = 1761
	16: {LogN: 16,
		Q:     []uint8{55, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
		P:     []uint8{55, 55, 55, 55},
		Scale: 1 << 45,
		Sigma: 3.2},
}

// MaxLogN is the largest supported polynomial modulus degree
const MaxLogN = 16

// MaxModuliCount is the largest supported number of moduli in the RNS representation
const MaxModuliCount = 34

// Equals compares two sets of parameters for equality
func (p *Parameters) Equals(other *Parameters) bool {
	if p == other {
		return true
	}
	return p.LogN == other.LogN && utils.EqualSliceUint8(p.Q, other.Q) && p.Scale == other.Scale && p.Sigma == other.Sigma
}

// MarshalBinary returns a []byte representation of the parameter set
func (p *Parameters) MarshalBinary() ([]byte, error) {
	if p.LogN == 0 { // if N is 0, then p is the zero value
		return []byte{}, nil
	}
	b := utils.NewBuffer(make([]byte, 0, 7+((2+len(p.Q))<<3)))
	b.WriteUint8(uint8(p.LogN))
	b.WriteUint8(uint8(len(p.Q)))
	b.WriteUint8(uint8(len(p.P)))
	b.WriteUint64(uint64(p.Scale))
	b.WriteUint64(uint64(p.Sigma * (1 << 32)))
	b.WriteUint8Slice(p.Q)
	b.WriteUint8Slice(p.P)
	return b.Bytes(), nil
}

// UnMarshalBinary decodes a []byte into a parameter set struct
func (p *Parameters) UnMarshalBinary(data []byte) error {
	if len(data) < 3 {
		return errors.New("invalid parameters encoding")
	}
	b := utils.NewBuffer(data)
	p.LogN = b.ReadUint8()
	if p.LogN > MaxLogN {
		return errors.New("polynomial degree is too large")
	}
	lenQ := uint64(b.ReadUint8())
	if lenQ > MaxModuliCount {
		return fmt.Errorf("len(Q) is larger than %d", MaxModuliCount)
	}
	lenP := uint64(b.ReadUint8())
	if lenP > MaxModuliCount {
		return fmt.Errorf("len(P) is larger than %d", MaxModuliCount)
	}

	p.Scale = float64(b.ReadUint64())

	p.Sigma = math.Round((float64(b.ReadUint64())/float64(1<<32))*100) / 100
	p.Q = make([]uint8, lenQ, lenQ)
	b.ReadUint8Slice(p.Q)
	p.P = make([]uint8, lenP, lenP)
	b.ReadUint8Slice(p.P)
	return nil
}
