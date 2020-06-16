package ckks

import (
	"errors"
	"fmt"
	"github.com/ldsec/lattigo/ring"
	"github.com/ldsec/lattigo/utils"
	"math"
	"math/bits"
)

// MaxLogN is the log2 of the largest supported polynomial modulus degree.
const MaxLogN = 16

// MaxModuliCount is the largest supported number of moduli in the RNS representation.
const MaxModuliCount = 34

// MaxModuliSize is the largest bit-length supported for the moduli in the RNS representation.
const MaxModuliSize = 60

func init() {
	for _, params := range DefaultParams {
		params.GenFromLogModuli()
	}
}

const (
	// PN12QP109 is the index in DefaultParams for logQ1 = 109
	PN12QP109 = iota
	// PN13QP218 is the index in DefaultParams for logQ1 = 218
	PN13QP218
	// PN14QP438 is the index in DefaultParams for logQ1 = 438
	PN14QP438
	// PN15QP880 is the index in DefaultParams for logQ1 = 880
	PN15QP880
	// PN16QP1761 is the index in DefaultParams for logQ1 = 1761
	PN16QP1761
)

// DefaultParams is a set of default CKKS parameters ensuring 128 bit security.
var DefaultParams = []*Parameters{

	//LogQi = 109
	{LogN: 12,
		LogSlots: 11,
		LogModuli: LogModuli{
			LogQi: []uint64{37, 32},
			LogPi: []uint64{38},
		},
		Scale: 1 << 32,
		Sigma: 3.2},

	//LogQi = 218
	{LogN: 13,
		LogSlots: 12,
		LogModuli: LogModuli{
			LogQi: []uint64{33, 30, 30, 30, 30, 30},
			LogPi: []uint64{35},
		},
		Scale: 1 << 30,
		Sigma: 3.2},

	//LogQiP = 438
	{LogN: 14,
		LogSlots: 13,
		LogModuli: LogModuli{
			LogQi: []uint64{45, 34, 34, 34, 34, 34, 34, 34, 34, 34},
			LogPi: []uint64{43, 43},
		},
		Scale: 1 << 34,
		Sigma: 3.2},

	//LogQi = 880
	{LogN: 15,
		LogSlots: 14,
		LogModuli: LogModuli{
			LogQi: []uint64{50, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40},
			LogPi: []uint64{50, 50, 50},
		},
		Scale: 1 << 40,
		Sigma: 3.2},

	//LogQi = 1761
	{LogN: 16,
		LogSlots: 15,
		LogModuli: LogModuli{
			LogQi: []uint64{55, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
			LogPi: []uint64{55, 55, 55, 55},
		},
		Scale: 1 << 45,
		Sigma: 3.2},
}

// Moduli stores the NTT primes of the RNS representation.
type Moduli struct {
	Qi []uint64 // Ciphertext prime moduli
	Pi []uint64 // Keys additional prime moduli
}

// Copy creates a copy of the target Moduli.
func (m *Moduli) Copy() Moduli {

	Qi := make([]uint64, len(m.Qi))
	copy(Qi, m.Qi)

	Pi := make([]uint64, len(m.Pi))
	copy(Pi, m.Pi)

	return Moduli{Qi, Pi}
}

// LogModuli stores the bit-length of the NTT primes of the RNS representation.
type LogModuli struct {
	LogQi []uint64 // Ciphertext prime moduli bit-size
	LogPi []uint64 // Keys additional prime moduli bit-size
}

// Copy creates a copy of the target LogModuli.
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
	LogModuli
	LogN     uint64 // Ring degree (power of 2)
	LogSlots uint64
	Scale    float64
	Sigma    float64 // Gaussian sampling variance

	logQP uint64
	alpha uint64
	beta  uint64

	isValid bool
}

// NewParametersFromModuli generates a new set or bfv parameters from the input parameters.
func NewParametersFromModuli(LogN, LogSlots uint64, Scale float64, moduli Moduli, Sigma float64) (params *Parameters) {

	if LogN > MaxLogN {
		panic(fmt.Errorf("LogN is larger than %d", MaxLogN))
	}

	if LogSlots >= LogN {
		panic(fmt.Errorf("LogSlots is larger than LogN-1"))
	}

	params = new(Parameters)
	params.LogN = LogN
	params.LogSlots = LogSlots
	params.Scale = Scale
	params.Sigma = Sigma
	params.Moduli = moduli.Copy()
	params.GenFromModuli()
	return
}

// NewParametersFromLogModuli generates a new set or bfv parameters from the input parameters.
func NewParametersFromLogModuli(LogN, LogSlots uint64, Scale float64, logModuli LogModuli, Sigma float64) (params *Parameters) {

	if LogN > MaxLogN {
		panic(fmt.Errorf("LogN is larger than %d", MaxLogN))
	}

	if LogSlots >= LogN {
		panic(fmt.Errorf("LogSlots is larger than LogN-1"))
	}

	params = new(Parameters)
	params.LogN = LogN
	params.LogSlots = LogSlots
	params.Scale = Scale
	params.Sigma = Sigma
	params.LogModuli = logModuli.Copy()
	params.GenFromLogModuli()
	return
}

// Alpha returns #Pi.
func (p *Parameters) Alpha() uint64 {
	return p.alpha
}

// Beta returns ceil(#Qi/#Pi).
func (p *Parameters) Beta() uint64 {
	return p.beta
}

// LogQP returns the bitlength of prod(Qi) * prod(Pi)
func (p *Parameters) LogQP() uint64 {
	return p.logQP
}

// MaxLevel returns #Qi -1
func (p *Parameters) MaxLevel() uint64 {
	return uint64(len(p.Qi) - 1)
}

// IsValid returns a true if the parameters are complete and valid, else false.
func (p *Parameters) IsValid() bool {
	return p.isValid
}

// NewPolyQ returns a new empty polynomial of degree 2^LogN in basis Qi.
func (p *Parameters) NewPolyQ() *ring.Poly {
	return ring.NewPoly(1<<p.LogN, uint64(len(p.Qi)))
}

// NewPolyP returns a new empty polynomial of degree 2^LogN in basis Pi.
func (p *Parameters) NewPolyP() *ring.Poly {
	return ring.NewPoly(1<<p.LogN, uint64(len(p.Pi)))
}

// NewPolyQP returns a new empty polynomial of degree 2^LogN in basis Qi + Pi.
func (p *Parameters) NewPolyQP() *ring.Poly {
	return ring.NewPoly(1<<p.LogN, uint64(len(p.Qi)+len(p.Pi)))
}

// Copy creates a copy of the target parameters.
func (p *Parameters) Copy() (paramsCopy *Parameters) {

	paramsCopy = new(Parameters)
	paramsCopy.LogN = p.LogN
	paramsCopy.LogSlots = p.LogSlots
	paramsCopy.Scale = p.Scale
	paramsCopy.Sigma = p.Sigma
	paramsCopy.Moduli = p.Moduli.Copy()
	paramsCopy.LogModuli = p.LogModuli.Copy()
	paramsCopy.logQP = p.logQP
	paramsCopy.alpha = p.alpha
	paramsCopy.beta = p.beta
	paramsCopy.isValid = p.isValid

	return
}

// Equals compares two sets of parameters for equality.
func (p *Parameters) Equals(other *Parameters) (res bool) {

	if p == other {
		return true
	}

	res = res && (p.LogN == other.LogN)
	res = res && (p.LogSlots == other.LogSlots)
	res = res && (p.Scale == other.Scale)
	res = res && (p.Sigma == other.Sigma)

	res = res && utils.EqualSliceUint64(p.Qi, other.Qi)
	res = res && utils.EqualSliceUint64(p.Pi, other.Pi)
	res = res && utils.EqualSliceUint64(p.LogQi, other.LogQi)
	res = res && utils.EqualSliceUint64(p.LogPi, other.LogPi)

	res = res && (p.alpha == other.alpha)
	res = res && (p.beta == other.beta)
	res = res && (p.logQP == other.logQP)

	res = res && (p.isValid == other.isValid)

	fmt.Println(res)

	return
}

// MarshalBinary returns a []byte representation of the parameter set.
func (p *Parameters) MarshalBinary() ([]byte, error) {
	if p.LogN == 0 { // if N is 0, then p is the zero value
		return []byte{}, nil
	}

	if !p.IsValid() {
		return nil, errors.New("cannot MarshalBinary: parameters not generated or invalid")
	}

	b := utils.NewBuffer(make([]byte, 0, 21+(len(p.LogQi)+len(p.LogPi))<<3))

	b.WriteUint8(uint8(p.LogN))
	b.WriteUint8(uint8(p.LogSlots))
	b.WriteUint64(math.Float64bits(p.Scale))
	b.WriteUint64(math.Float64bits(p.Sigma))
	b.WriteUint8(uint8(len(p.Qi)))
	b.WriteUint8(uint8(len(p.Pi)))
	b.WriteUint64Slice(p.Qi)
	b.WriteUint64Slice(p.Pi)

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
		return fmt.Errorf("LogN larger than %d", MaxLogN)
	}

	p.LogSlots = uint64(b.ReadUint8())

	if p.LogSlots > p.LogN-1 {
		return fmt.Errorf("LogSlots larger than %d", MaxLogN-1)
	}

	p.Scale = math.Float64frombits(b.ReadUint64())
	p.Sigma = math.Float64frombits(b.ReadUint64())

	lenLogQi := b.ReadUint8()
	lenLogPi := b.ReadUint8()

	p.Qi = make([]uint64, lenLogQi, lenLogQi)
	p.Pi = make([]uint64, lenLogPi, lenLogPi)

	b.ReadUint64Slice(p.Qi)
	b.ReadUint64Slice(p.Pi)

	if err := p.checkModuli(); err != nil {
		return err
	}

	p.GenFromModuli()

	return nil
}

// GenFromModuli generates the parameters using the provided moduli.
func (p *Parameters) GenFromModuli() {

	if err := p.checkModuli(); err != nil {
		panic(err)
	}

	tmp := ring.NewUint(1)

	for _, qi := range p.Qi {
		tmp.Mul(tmp, ring.NewUint(qi))
	}

	for _, pi := range p.Pi {
		tmp.Mul(tmp, ring.NewUint(pi))
	}

	p.logQP = uint64(tmp.BitLen())

	p.LogQi = make([]uint64, len(p.Qi), len(p.Qi))
	for i := range p.Qi {
		p.LogQi[i] = uint64(bits.Len64(p.Qi[i]) - 1)
	}

	p.LogPi = make([]uint64, len(p.Pi), len(p.Pi))
	for i := range p.Pi {
		p.LogPi[i] = uint64(bits.Len64(p.Pi[i]) - 1)
	}

	p.alpha = uint64(len(p.Pi))
	p.beta = uint64(math.Ceil(float64(len(p.Qi)) / float64(len(p.Pi))))

	p.isValid = true
}

// GenFromLogModuli generates the parameters using the given bit-size for the moduli.
func (p *Parameters) GenFromLogModuli() {

	if err := p.checkLogModuli(); err != nil {
		panic(err)
	}

	p.Qi, p.Pi = GenModuli(p)

	p.GenFromModuli()
}

func (p *Parameters) checkModuli() error {

	if len(p.Qi) > MaxModuliCount {
		return fmt.Errorf("#LogQi is larger than %d", MaxModuliCount)
	}

	if len(p.Pi) > MaxModuliCount {
		return fmt.Errorf("#Pi is larger than %d", MaxModuliCount)
	}

	for i, qi := range p.Qi {
		if uint64(bits.Len64(qi)-1) > MaxModuliSize {
			return fmt.Errorf("Qi bit-size (i=%d) is larger than %d", i, MaxModuliSize)
		}
	}

	for i, pi := range p.Pi {
		if uint64(bits.Len64(pi)-1) > MaxModuliSize {
			return fmt.Errorf("Pi bit-size (i=%d) is larger than %d", i, MaxModuliSize)
		}
	}

	N := uint64(1 << p.LogN)

	for i, qi := range p.Qi {
		if !ring.IsPrime(qi) || qi&((N<<1)-1) != 1 {
			return fmt.Errorf("Qi (i=%d) is not an NTT prime", i)
		}
	}

	for i, pi := range p.Pi {
		if !ring.IsPrime(pi) || pi&((N<<1)-1) != 1 {
			return fmt.Errorf("Pi (i=%d) is not an NTT prime", i)
		}
	}

	return nil
}

func (p *Parameters) checkLogModuli() error {

	if len(p.LogQi) > MaxModuliCount {
		return fmt.Errorf("#LogQi is larger than %d", MaxModuliCount)
	}

	if len(p.LogPi) > MaxModuliCount {
		return fmt.Errorf("#LogPi is larger than %d", MaxModuliCount)
	}

	for i, qi := range p.LogQi {
		if qi > MaxModuliSize {
			return fmt.Errorf("LogQi (i=%d) is larger than %d", i, MaxModuliSize)
		}
	}

	for i, pi := range p.LogPi {
		if pi > MaxModuliSize {
			return fmt.Errorf("LogPi (i=%d) is larger than %d", i, MaxModuliSize)
		}
	}

	return nil
}
