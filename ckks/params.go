package ckks

import (
	"errors"
	"fmt"
	"math"
	"math/bits"

	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
)

// MaxLogN is the log2 of the largest supported polynomial modulus degree.
const MaxLogN = 16

// MaxModuliCount is the largest supported number of moduli in the RNS representation.
const MaxModuliCount = 34

// MaxModuliSize is the largest bit-length supported for the moduli in the RNS representation.
const MaxModuliSize = 60

// DefaultSigma is the default error distribution standard deviation
const DefaultSigma = 3.2

const (
	PN12QP109 = iota
	PN13QP218
	PN14QP438
	PN15QP880
	PN16QP1761
	PN16BootCheby
)

type DefaultParam struct {
	LogModuli
	LogN     uint64
	LogSlots uint64
	Scale    float64
}

// DefaultParams is a set of default CKKS parameters ensuring 128 bit security.
var DefaultParams = []*DefaultParam{

	//LogQi = 109
	{LogN: 12,
		LogSlots: 11,
		LogModuli: LogModuli{
			LogQi: []uint64{37, 32},
			LogPi: []uint64{38},
		},
		Scale: 1 << 32,
	},

	//LogQi = 218
	{LogN: 13,
		LogSlots: 12,
		LogModuli: LogModuli{
			LogQi: []uint64{33, 30, 30, 30, 30, 30},
			LogPi: []uint64{35},
		},
		Scale: 1 << 30,
	},

	//LogQiP = 438
	{LogN: 14,
		LogSlots: 13,
		LogModuli: LogModuli{
			LogQi: []uint64{45, 34, 34, 34, 34, 34, 34, 34, 34, 34},
			LogPi: []uint64{43, 43},
		},
		Scale: 1 << 34,
	},

	//LogQi = 880
	{LogN: 15,
		LogSlots: 14,
		LogModuli: LogModuli{
			LogQi: []uint64{50, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40},
			LogPi: []uint64{50, 50, 50},
		},
		Scale: 1 << 40,
	},

	//LogQi = 1761
	{LogN: 16,
		LogSlots: 15,
		LogModuli: LogModuli{
			LogQi: []uint64{55, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
			LogPi: []uint64{55, 55, 55, 55},
		},
		Scale: 1 << 45,
	},

	//LogQi = 1761
	{LogN: 16,
		LogSlots: 15,
		LogModuli: LogModuli{
			LogQi: []uint64{55, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
			LogPi: []uint64{61},
		},
		Scale: 1 << 55,
	},
}

// Moduli stores the NTT primes of the RNS representation.
type Moduli struct {
	qi []uint64 // Ciphertext prime moduli
	pi []uint64 // Keys additional prime moduli
}

// Qi returns a new slice with the factors of the ciphertext modulus q
func (m *Moduli) Qi() []uint64 {
	qi := make([]uint64, len(m.qi))
	copy(qi, m.qi)
	return qi
}

// QiCount returns the number of factors of the ciphertext modulus q
func (m *Moduli) QiCount() uint64 {
	return uint64(len(m.qi))
}

// Pi returns a new slice with the factors of the ciphertext modulus extention p
func (m *Moduli) Pi() []uint64 {
	pi := make([]uint64, len(m.pi))
	copy(pi, m.pi)
	return pi
}

// PiCount returns the number of factors of the ciphertext modulus extention p
func (m *Moduli) PiCount() uint64 {
	return uint64(len(m.pi))
}

// PiCount returns the number of factors of the ciphertext modulus extention p
func (m *Moduli) QPiCount() uint64 {
	return uint64(len(m.qi) + len(m.pi))
}

// LogQP returns the size of the extended modulus QP in bits
func (m *Moduli) LogQP() uint64 {
	tmp := ring.NewUint(1)
	for _, qi := range m.qi {
		tmp.Mul(tmp, ring.NewUint(qi))
	}
	for _, pi := range m.pi {
		tmp.Mul(tmp, ring.NewUint(pi))
	}
	return uint64(tmp.BitLen())
}

// Copy creates a copy of the target Moduli.
func (m *Moduli) Copy() Moduli {

	qi := make([]uint64, len(m.qi))
	copy(qi, m.qi)

	pi := make([]uint64, len(m.pi))
	copy(pi, m.pi)

	return Moduli{qi, pi}
}

// LogModuli generates a LogModuli struct from the parameters' Moduli struct and returns it.
func (p *Parameters) LogModuli() (lm LogModuli) {

	lm.LogQi = make([]uint64, len(p.qi), len(p.qi))
	for i := range p.qi {
		lm.LogQi[i] = uint64(math.Round(math.Log2(float64(p.qi[i]))))
	}

	lm.LogPi = make([]uint64, len(p.pi), len(p.pi))
	for i := range p.pi {
		lm.LogPi[i] = uint64(math.Round(math.Log2(float64(p.pi[i]))))
	}

	return lm
}

// LogModuli stores the bit-length of the NTT primes of the RNS representation.
type LogModuli struct {
	LogQi []uint64 // Ciphertext prime moduli bit-size
	LogPi []uint64 // Keys additional prime moduli bit-size
}

// Copy creates a copy of the target Moduli.
func (m *LogModuli) Copy() LogModuli {

	LogQi := make([]uint64, len(m.LogQi))
	copy(LogQi, m.LogQi)

	LogPi := make([]uint64, len(m.LogPi))
	copy(LogPi, m.LogPi)

	return LogModuli{LogQi, LogPi}
}

// Parameters represents a given parameter set for the BFV cryptosystem.
type Parameters struct {
	Moduli
	logN     uint64 // Ring degree (power of 2)
	logSlots uint64
	scale    float64
	sigma    float64 // Gaussian sampling variance

	logQP uint64
	alpha uint64
	beta  uint64
}

// NewParametersFromModuli creates a new Parameters struct and returns a pointer to it.
func NewParametersFromModuli(logN uint64, m Moduli) (p *Parameters, err error) {
	p = new(Parameters)

	if (logN < 3) || (logN > MaxLogN) {
		return nil, fmt.Errorf("invalid polynomial ring log degree: %d", logN)
	}

	p.logN = logN

	if err = checkModuli(m, logN); err != nil {
		return nil, err
	}

	p.Moduli = m.Copy()
	p.logQP = p.Moduli.LogQP()

	if p.PiCount() != 0 {
		p.alpha = p.PiCount()
		p.beta = uint64(math.Ceil(float64(p.QiCount()) / float64(p.PiCount())))
	}

	p.sigma = DefaultSigma

	return p, nil

}

// NewParametersFromLogModuli creates a new Parameters struct and returns a pointer to it.
func NewParametersFromLogModuli(logN uint64, lm LogModuli) (p *Parameters, err error) {

	if err = checkLogModuli(lm); err != nil {
		return nil, err
	}

	// If LogModuli is valid and then generates the moduli
	return NewParametersFromModuli(logN, genModuli(lm, logN))
}

// NewPolyQ returns a new empty polynomial of degree 2^LogN in basis Qi.
func (p *Parameters) NewPolyQ() *ring.Poly {
	return ring.NewPoly(p.N(), p.QiCount())
}

// NewPolyP returns a new empty polynomial of degree 2^LogN in basis Pi.
func (p *Parameters) NewPolyP() *ring.Poly {
	return ring.NewPoly(p.N(), p.PiCount())
}

// NewPolyQP returns a new empty polynomial of degree 2^LogN in basis Qi + Pi.
func (p *Parameters) NewPolyQP() *ring.Poly {
	return ring.NewPoly(p.N(), p.QPiCount())
}

// N returns the ring degree
func (p *Parameters) N() uint64 {
	return 1 << p.logN
}

// LogN returns the log of the degree of the polynomial ring
func (p *Parameters) LogN() uint64 {
	return p.logN
}

// LogSlots returns the log of the number of slots
func (p *Parameters) LogSlots() uint64 {
	return p.logSlots
}

// MaxLevel returns the maximum ciphertext level
func (p *Parameters) MaxLevel() uint64 {
	return p.QiCount() - 1
}

// Levels returns then number of total levels enabled by the parameters
func (p *Parameters) Levels() uint64 {
	return p.QiCount()
}

// LogN returns the log of the degree of the polynomial ring
func (p *Parameters) Slots() uint64 {
	return 1 << p.logSlots
}

// LogN returns the log of the degree of the polynomial ring
func (p *Parameters) MaxSlots() uint64 {
	return p.N() >> 1
}

func (p *Parameters) LogMaxSlots() uint64 {
	return p.logN - 1
}

// MaxLogSlots returns the log of the maximum number of slots enabled by the parameters
func (p *Parameters) MaxLogSlots() uint64 {
	return p.logN - 1
}

// Sigma returns standard deviation of the noise distribution
func (p *Parameters) Sigma() float64 {
	return p.sigma
}

// Alpha returns the number of moduli in in P
func (p *Parameters) Alpha() uint64 {
	return p.alpha
}

// Beta returns the number of element in the RNS decomposition basis: Ceil(lenQi / lenPi)
func (p *Parameters) Beta() uint64 {
	return p.beta
}

// LogQP returns the size of the extanded modulus QP in bits
func (p *Parameters) LogQP() uint64 {
	return p.logQP
}

// WithT returns a copy of this parameter struct with the plaintext modulus set to t
func (p *Parameters) Scale() float64 {
	return p.scale
}

// WithT returns a copy of this parameter struct with the plaintext modulus set to t
func (p *Parameters) SetScale(scale float64) {
	p.scale = scale
}

// SetLogSlots sets the value logSlots of the parameters.
func (p *Parameters) SetLogSlots(logSlots uint64) (err error) {
	if (logSlots == 0) || (logSlots > p.MaxLogSlots()) {
		return fmt.Errorf("slots cannot be greater than LogN-1")
	}

	p.logSlots = logSlots

	return nil
}

// SetSigma sets the value sigma of the parameters
func (p *Parameters) SetSigma(sigma float64) {
	p.sigma = sigma
}

// Copy creates a copy of the target parameters.
func (p *Parameters) Copy() (paramsCopy *Parameters) {

	paramsCopy = new(Parameters)
	paramsCopy.logN = p.logN
	paramsCopy.logSlots = p.logSlots
	paramsCopy.scale = p.scale
	paramsCopy.sigma = p.sigma
	paramsCopy.Moduli = p.Moduli.Copy()
	paramsCopy.logQP = p.logQP
	paramsCopy.alpha = p.alpha
	paramsCopy.beta = p.beta

	return
}

// Equals compares two sets of parameters for equality.
func (p *Parameters) Equals(other *Parameters) (res bool) {

	if p == other {
		return true
	}

	res = p.logN == other.logN
	res = res && (p.logSlots == other.logSlots)
	res = res && (p.scale == other.scale)
	res = res && (p.sigma == other.sigma)

	res = res && utils.EqualSliceUint64(p.qi, other.qi)
	res = res && utils.EqualSliceUint64(p.pi, other.pi)

	res = res && (p.alpha == other.alpha)
	res = res && (p.beta == other.beta)
	res = res && (p.logQP == other.logQP)

	return
}

// MarshalBinary returns a []byte representation of the parameter set.
func (p *Parameters) MarshalBinary() ([]byte, error) {
	if p.logN == 0 { // if N is 0, then p is the zero value
		return []byte{}, nil
	}

	b := utils.NewBuffer(make([]byte, 0, 21+(p.QPiCount())<<3))

	b.WriteUint8(uint8(p.logN))
	b.WriteUint8(uint8(p.logSlots))
	b.WriteUint64(math.Float64bits(p.scale))
	b.WriteUint64(math.Float64bits(p.sigma))
	b.WriteUint8(uint8(p.QiCount()))
	b.WriteUint8(uint8(p.PiCount()))
	b.WriteUint64Slice(p.qi)
	b.WriteUint64Slice(p.pi)

	return b.Bytes(), nil
}

// UnmarshalBinary decodes a []byte into a parameter set struct
func (p *Parameters) UnmarshalBinary(data []byte) (err error) {

	if len(data) < 3 {
		return errors.New("invalid parameters encoding")
	}

	b := utils.NewBuffer(data)

	p.logN = uint64(b.ReadUint8())

	if p.logN > MaxLogN {
		return fmt.Errorf("LogN larger than %d", MaxLogN)
	}

	p.logSlots = uint64(b.ReadUint8())

	if p.logSlots > p.logN-1 {
		return fmt.Errorf("LogSlots larger than %d", MaxLogN-1)
	}

	p.scale = math.Float64frombits(b.ReadUint64())
	p.sigma = math.Float64frombits(b.ReadUint64())

	lenQi := b.ReadUint8()
	lenPi := b.ReadUint8()

	p.qi = make([]uint64, lenQi, lenQi)
	p.pi = make([]uint64, lenPi, lenPi)

	b.ReadUint64Slice(p.qi)
	b.ReadUint64Slice(p.pi)

	if lenPi != 0 {
		p.alpha = uint64(lenPi)
		p.beta = uint64(math.Ceil(float64(lenQi) / float64(lenPi)))
	}

	if err = checkModuli(p.Moduli, p.logN); err != nil {
		return err
	}

	p.logQP = p.Moduli.LogQP()

	return nil
}

func checkModuli(m Moduli, logN uint64) error {

	if len(m.qi) > MaxModuliCount {
		return fmt.Errorf("#Qi is larger than %d", MaxModuliCount)
	}

	if len(m.pi) > MaxModuliCount {
		return fmt.Errorf("#Pi is larger than %d", MaxModuliCount)
	}

	for i, qi := range m.qi {
		if uint64(bits.Len64(qi)-1) > MaxModuliSize+1 {
			return fmt.Errorf("Qi bit-size (i=%d) is larger than %d", i, MaxModuliSize)
		}
	}

	for i, pi := range m.pi {
		if uint64(bits.Len64(pi)-1) > MaxModuliSize+2 {
			return fmt.Errorf("Pi bit-size (i=%d) is larger than %d", i, MaxModuliSize)
		}
	}

	N := uint64(1 << logN)

	for i, qi := range m.qi {
		if !ring.IsPrime(qi) || qi&((N<<1)-1) != 1 {
			return fmt.Errorf("Qi (i=%d) is not an NTT prime", i)
		}
	}

	for i, pi := range m.pi {
		if !ring.IsPrime(pi) || pi&((N<<1)-1) != 1 {
			return fmt.Errorf("Pi (i=%d) is not an NTT prime", i)
		}
	}

	return nil
}

func checkLogModuli(m LogModuli) error {

	if len(m.LogQi) > MaxModuliCount {
		return fmt.Errorf("#LogQi is larger than %d", MaxModuliCount)
	}

	if len(m.LogPi) > MaxModuliCount {
		return fmt.Errorf("#LogPi is larger than %d", MaxModuliCount)
	}

	for i, qi := range m.LogQi {
		if qi > MaxModuliSize {
			return fmt.Errorf("LogQi (i=%d) is larger than %d", i, MaxModuliSize)
		}
	}

	for i, pi := range m.LogPi {
		if pi > MaxModuliSize+1 {
			return fmt.Errorf("LogPi (i=%d) is larger than %d", i, MaxModuliSize)
		}
	}

	return nil
}

func genModuli(lm LogModuli, logN uint64) (m Moduli) {

	// Extracts all the different primes bit size and maps their number
	primesbitlen := make(map[uint64]uint64)
	for _, qi := range lm.LogQi {
		primesbitlen[qi]++
	}

	for _, pj := range lm.LogPi {
		primesbitlen[pj]++
	}

	// For each bit-size, finds that many primes
	primes := make(map[uint64][]uint64)
	for key, value := range primesbitlen {
		primes[key] = ring.GenerateNTTPrimes(key, logN, value)
	}

	// Assigns the primes to the ckks moduli chain
	m.qi = make([]uint64, len(lm.LogQi))
	for i, qi := range lm.LogQi {
		m.qi[i] = primes[qi][0]
		primes[qi] = primes[qi][1:]
	}

	// Assigns the primes to the special primes list for the extended ring
	m.pi = make([]uint64, len(lm.LogPi))
	for i, pj := range lm.LogPi {
		m.pi[i] = primes[pj][0]
		primes[pj] = primes[pj][1:]
	}

	return m
}
