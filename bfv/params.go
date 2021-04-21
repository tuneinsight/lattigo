package bfv

import (
	"errors"
	"fmt"
	"math"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

const (
	// PN12QP109 is a set of parameters with N = 2^12 and log(QP) = 109
	PN12QP109 = iota
	// PN13QP218 is a set of parameters with N = 2^13 and log(QP) = 218
	PN13QP218
	// PN14QP438 is a set of parameters with N = 2^14 and log(QP) = 438
	PN14QP438
	// PN15QP880 is a set of parameters with N = 2^15 and log(QP) = 880
	PN15QP880

	// PN12QP101pq is the index in DefaultParams for logQP = 101 (post quantum)
	PN12QP101pq
	// PN13QP202pq is the index in DefaultParams for logQP = 202 (post quantum)
	PN13QP202pq
	// PN14QP411pq is the index in DefaultParams for logQP = 411 (post quantum)
	PN14QP411pq
	// PN15QP827pq is the index in DefaultParams for logQP = 827 (post quantum)
	PN15QP827pq
)

// DefaultParams is a set of default BFV parameters ensuring 128 bit security.
var DefaultParams = []ParametersLiteral{

	{
		LogN:  12,
		T:     65537,
		Q:     []uint64{0x7ffffec001, 0x8000016001}, // 39 + 39 bits
		P:     []uint64{0x40002001},                 // 30 bits
		Sigma: rlwe.DefaultSigma,
	},

	{
		LogN:  13,
		T:     65537,
		Q:     []uint64{0x3fffffffef8001, 0x4000000011c001, 0x40000000120001}, // 54 + 54 + 54 bits
		P:     []uint64{0x7ffffffffb4001},                                     // 55 bits
		Sigma: rlwe.DefaultSigma,
	},

	{
		LogN: 14,
		T:    65537,
		Q: []uint64{0x100000000060001, 0x80000000068001, 0x80000000080001,
			0x3fffffffef8001, 0x40000000120001, 0x3fffffffeb8001}, // 56 + 55 + 55 + 54 + 54 + 54 bits
		P:     []uint64{0x80000000130001, 0x7fffffffe90001}, // 55 + 55 bits
		Sigma: rlwe.DefaultSigma,
	},

	{
		LogN: 15,
		T:    65537,
		Q: []uint64{0x7ffffffffe70001, 0x7ffffffffe10001, 0x7ffffffffcc0001, // 59 + 59 + 59 bits
			0x400000000270001, 0x400000000350001, 0x400000000360001, // 58 + 58 + 58 bits
			0x3ffffffffc10001, 0x3ffffffffbe0001, 0x3ffffffffbd0001, // 58 + 58 + 58 bits
			0x4000000004d0001, 0x400000000570001, 0x400000000660001}, // 58 + 58 + 58 bits
		P:     []uint64{0xffffffffffc0001, 0x10000000001d0001, 0x10000000006e0001}, // 60 + 60 + 60 bits
		Sigma: rlwe.DefaultSigma,
	},

	{ // LogQP = 101.00005709794536
		LogN:  12,
		T:     65537,
		Q:     []uint64{0x800004001, 0x800008001}, // 2*35
		P:     []uint64{0x80014001},               // 1*31
		Sigma: rlwe.DefaultSigma,
	},

	{ // LogQP = 201.99999999994753
		LogN:  13,
		T:     65537,
		Q:     []uint64{0x7fffffffe0001, 0x7fffffffcc001, 0x3ffffffffc001}, // 2*51 + 50
		P:     []uint64{0x4000000024001},                                   // 50
		Sigma: rlwe.DefaultSigma,
	},

	{ // LogQP = 410.9999999999886
		LogN:  14,
		T:     65537,
		Q:     []uint64{0x7fffffffff18001, 0x8000000000f8001, 0x7ffffffffeb8001, 0x800000000158001, 0x7ffffffffe70001}, // 5*59
		P:     []uint64{0x7ffffffffe10001, 0x400000000068001},                                                          // 59+58
		Sigma: rlwe.DefaultSigma,
	},

	{ // LogQP = 826.9999999999509
		LogN: 15,
		T:    65537,
		Q: []uint64{0x7ffffffffe70001, 0x7ffffffffe10001, 0x7ffffffffcc0001, 0x7ffffffffba0001, 0x8000000004a0001,
			0x7ffffffffb00001, 0x800000000890001, 0x8000000009d0001, 0x7ffffffff630001, 0x800000000a70001,
			0x7ffffffff510001}, // 11*59
		P:     []uint64{0x800000000b80001, 0x800000000bb0001, 0xffffffffffc0001}, // 2*59+60
		Sigma: rlwe.DefaultSigma,
	},
}

type Parameters struct {
	rlwe.Parameters
	t uint64
}

// ParametersLiteral represents a given parameter set for the BFV cryptosystem.
type ParametersLiteral struct {
	LogN  int // Log Ring degree (power of 2)
	Q     []uint64
	P     []uint64
	T     uint64  // Plaintext modulus
	Sigma float64 // Gaussian sampling standard deviation
}

func NewParametersFromParamDef(paramDef ParametersLiteral) (Parameters, error) {
	m := new(rlwe.Moduli)
	m.Qi = make([]uint64, len(paramDef.Q))
	copy(m.Qi, paramDef.Q)
	m.Pi = make([]uint64, len(paramDef.P))
	copy(m.Pi, paramDef.P)
	return NewParametersFromModuli(paramDef.LogN, m, paramDef.T)
}

// NewParametersFromModuli creates a new Parameters struct and returns a pointer to it.
func NewParametersFromModuli(logN int, m *rlwe.Moduli, t uint64) (p Parameters, err error) {

	var rlweParams rlwe.Parameters
	if rlweParams, err = rlwe.NewRLWEParameters(logN, m.Qi, m.Pi, rlwe.DefaultSigma); err != nil {
		return Parameters{}, err
	}

	return Parameters{rlweParams, t}, nil
}

// NewParametersFromLogModuli creates a new Parameters struct and returns a pointer to it.
func NewParametersFromLogModuli(logN int, lm *rlwe.LogModuli, t uint64) (p Parameters, err error) {

	if err = rlwe.CheckLogModuli(lm); err != nil {
		return Parameters{}, err
	}

	// If LogModuli is valid and then generates the moduli
	return NewParametersFromModuli(logN, rlwe.GenModuli(lm, logN), t)
}

// T returns the plaintext coefficient modulus t
func (p Parameters) T() uint64 {
	return p.t
}

func (p Parameters) RingT() *ring.Ring {
	ringQP, err := ring.NewRing(p.N(), []uint64{p.t})
	if err != nil {
		panic(err) // Parameter type invariant
	}
	return ringQP
}

// Equals compares two sets of parameters for equality.
func (p *Parameters) Equals(other Parameters) (res bool) {

	res = p.LogN() == other.LogN()
	res = res && (p.t == other.T())
	res = res && (p.Sigma() == other.Sigma())

	res = res && utils.EqualSliceUint64(p.Q(), other.Q())
	res = res && utils.EqualSliceUint64(p.P(), other.P())

	return
}

// MarshalBinary returns a []byte representation of the parameter set.
func (p *Parameters) MarshalBinary() ([]byte, error) {
	if p.LogN() == 0 { // if N is 0, then p is the zero value
		return []byte{}, nil
	}

	// data : 19 byte + len(QPi) * 8 byte
	// 1 byte : logN
	// 1 byte : #pi
	// 1 byte : #pi
	// 8 byte : t
	// 8 byte : sigma
	b := utils.NewBuffer(make([]byte, 0, 19+(len(p.Q())+len(p.P()))<<3))

	b.WriteUint8(uint8(p.LogN()))
	b.WriteUint8(uint8(len(p.Q())))
	b.WriteUint8(uint8(len(p.P())))
	b.WriteUint64(p.T())
	b.WriteUint64(math.Float64bits(p.Sigma()))
	b.WriteUint64Slice(p.Q())
	b.WriteUint64Slice(p.P())

	return b.Bytes(), nil
}

// UnmarshalBinary decodes a []byte into a parameter set struct.
func (p *Parameters) UnmarshalBinary(data []byte) error {
	if len(data) < 19 {
		return errors.New("invalid parameters encoding")
	}
	b := utils.NewBuffer(data)

	logN := int(b.ReadUint8())

	if logN > rlwe.MaxLogN {
		return fmt.Errorf("logN larger than %d", rlwe.MaxLogN)
	}

	lenQi := b.ReadUint8()
	lenPi := b.ReadUint8()

	p.t = b.ReadUint64()
	sigma := math.Float64frombits(b.ReadUint64())
	qi := make([]uint64, lenQi)
	pi := make([]uint64, lenPi)

	b.ReadUint64Slice(qi)
	b.ReadUint64Slice(pi)

	var err error
	if p.Parameters, err = rlwe.NewRLWEParameters(logN, qi, pi, sigma); err != nil {
		return err
	}

	if err != nil {
		return err
	}
	return nil
}
