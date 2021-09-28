package bfv

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

// Plaintext moduli allowing batching for the corresponding N in ascending bit-size.
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

const (
	PN12QP109 = iota
	PN13QP218
	PN14QP438
	PN15QP880
)

// DefaultParams is a set of default BFV parameters ensuring 128 bit security.
var DefaultParams = []*Parameters{

	//logQ1+P = 109
	{LogN: 12,
		T: 65537,
		LogModuli: LogModuli{
			LogQi:    []uint64{39, 39},
			LogPi:    []uint64{30},
			LogQiMul: []uint64{60, 60},
		},
		Sigma: 3.2},

	//logQ1+P = 218
	{LogN: 13,
		T: 65537,
		LogModuli: LogModuli{
			LogQi:    []uint64{54, 54, 54},
			LogPi:    []uint64{55},
			LogQiMul: []uint64{60, 60, 60},
		},
		Sigma: 3.2},

	//logQ1+P = 438
	{LogN: 14,
		T: 65537,
		LogModuli: LogModuli{
			LogQi:    []uint64{56, 55, 55, 54, 54, 54},
			LogPi:    []uint64{55, 55},
			LogQiMul: []uint64{60, 60, 60, 60, 60, 60},
		},
		Sigma: 3.2},

	//logQ1+P = 880
	{LogN: 15,
		T: 65537,
		LogModuli: LogModuli{
			LogQi:    []uint64{59, 59, 59, 58, 58, 58, 58, 58, 58, 58, 58, 58},
			LogPi:    []uint64{60, 60, 60},
			LogQiMul: []uint64{60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60},
		},
		Sigma: 3.2},
}

// Moduli stores the NTT primes of the RNS representation.
type Moduli struct {
	Qi    []uint64 // Ciphertext prime moduli
	Pi    []uint64 // Keys additional prime moduli
	QiMul []uint64 // Ciphertext secondary prime moduli
}

// Copy creates a copy of the target Moduli.
func (m *Moduli) Copy() Moduli {

	Qi := make([]uint64, len(m.Qi))
	copy(Qi, m.Qi)

	Pi := make([]uint64, len(m.Pi))
	copy(Pi, m.Pi)

	QiMul := make([]uint64, len(m.QiMul))
	copy(QiMul, m.QiMul)

	return Moduli{Qi, Pi, QiMul}
}

// LogModuli stores the bit-length of the NTT primes of the RNS representation.
type LogModuli struct {
	LogQi    []uint64 // Ciphertext prime moduli bit-size
	LogPi    []uint64 // Keys additional prime moduli bit-size
	LogQiMul []uint64 // Ciphertext secondary prime moduli bit-size
}

// Copy creates a copy of the target LogModuli.
func (m *LogModuli) Copy() LogModuli {

	LogQi := make([]uint64, len(m.LogQi))
	copy(LogQi, m.LogQi)

	LogPi := make([]uint64, len(m.LogPi))
	copy(LogPi, m.LogPi)

	LogQiMul := make([]uint64, len(m.LogQiMul))
	copy(LogQiMul, m.LogQiMul)

	return LogModuli{LogQi, LogPi, LogQiMul}
}

// Parameters represents a given parameter set for the BFV cryptosystem.
type Parameters struct {
	Moduli
	LogModuli
	LogN  uint64  // Ring degree (power of 2)
	T     uint64  // Plaintext modulus
	Sigma float64 // Gaussian sampling standard deviation

	logQP uint64
	alpha uint64
	beta  uint64

	isValid bool
}

// NewParametersFromModuli generates a new set or BFV parameters from the input parameters.
func NewParametersFromModuli(LogN, T uint64, moduli Moduli, sigma float64) (params *Parameters) {

	if LogN > MaxLogN {
		panic(fmt.Errorf("cannot NewParametersFromModuli: LogN is larger than %d", MaxLogN))
	}

	params = new(Parameters)
	params.LogN = LogN
	params.T = T
	params.Sigma = sigma
	params.Moduli = moduli.Copy()
	params.GenFromModuli()
	return
}

// NewParametersFromLogModuli generates a new set or BFV parameters from the input parameters.
func NewParametersFromLogModuli(LogN, T uint64, logModuli LogModuli, sigma float64) (params *Parameters) {

	if LogN > MaxLogN {
		panic(fmt.Errorf("cannot NewParametersFromLogModuli: LogN is larger than %d", MaxLogN))
	}

	params = new(Parameters)
	params.LogN = LogN
	params.T = T
	params.Sigma = sigma
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

// LogQP returns the bit-length of prod(Qi) * prod(Pi)
func (p *Parameters) LogQP() uint64 {
	return p.logQP
}

// IsValid returns a true if the parameters are complete and valid, and false otherwise.
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

// Copy creates a copy of the target Parameters.
func (p *Parameters) Copy() (paramsCopy *Parameters) {

	paramsCopy = new(Parameters)
	paramsCopy.LogN = p.LogN
	paramsCopy.T = p.T
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
	res = res && (p.T == other.T)
	res = res && (p.Sigma == other.Sigma)

	res = res && utils.EqualSliceUint64(p.Qi, other.Qi)
	res = res && utils.EqualSliceUint64(p.Pi, other.Pi)
	res = res && utils.EqualSliceUint64(p.QiMul, other.QiMul)
	res = res && utils.EqualSliceUint64(p.LogQi, other.LogQi)
	res = res && utils.EqualSliceUint64(p.LogPi, other.LogPi)
	res = res && utils.EqualSliceUint64(p.LogQiMul, other.LogQiMul)

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

	b := utils.NewBuffer(make([]byte, 0, 20+(len(p.LogQi)+len(p.LogPi)+len(p.LogQiMul))<<3))

	b.WriteUint8(uint8(p.LogN))
	b.WriteUint8(uint8(len(p.Qi)))
	b.WriteUint8(uint8(len(p.Pi)))
	b.WriteUint8(uint8(len(p.QiMul)))
	b.WriteUint64(p.T)
	b.WriteUint64(uint64(p.Sigma * (1 << 32)))
	b.WriteUint64Slice(p.Qi)
	b.WriteUint64Slice(p.Pi)
	b.WriteUint64Slice(p.QiMul)

	return b.Bytes(), nil
}

// UnmarshalBinary decodes a []byte into a parameter set struct.
func (p *Parameters) UnmarshalBinary(data []byte) error {
	if len(data) < 3 {
		return errors.New("invalid parameters encoding")
	}
	b := utils.NewBuffer(data)

	p.LogN = uint64(b.ReadUint8())

	if p.LogN > MaxLogN {
		return fmt.Errorf("LogN larger than %d", MaxLogN)
	}

	lenLogQi := b.ReadUint8()
	lenLogPi := b.ReadUint8()
	lenLogQiMul := b.ReadUint8()

	p.T = b.ReadUint64()
	p.Sigma = math.Round((float64(b.ReadUint64())/float64(1<<32))*100) / 100
	p.Qi = make([]uint64, lenLogQi, lenLogQi)
	p.Pi = make([]uint64, lenLogPi, lenLogPi)
	p.QiMul = make([]uint64, lenLogQiMul, lenLogQiMul)

	b.ReadUint64Slice(p.Qi)
	b.ReadUint64Slice(p.Pi)
	b.ReadUint64Slice(p.QiMul)

	if err := p.checkModuli(); err != nil {
		return err
	}

	p.GenFromModuli()

	return nil
}

// GenFromModuli generates a set of parameters from the moduli chain.
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

	p.LogQiMul = make([]uint64, len(p.QiMul), len(p.QiMul))
	for i := range p.QiMul {
		p.LogQiMul[i] = uint64(bits.Len64(p.QiMul[i]) - 1)
	}

	p.alpha = uint64(len(p.Pi))
	p.beta = uint64(math.Ceil(float64(len(p.Qi)) / float64(len(p.Pi))))

	p.isValid = true
}

// GenFromLogModuli generates a set of parameters, including the actual moduli, from the target bit-sizes of the moduli chain.
func (p *Parameters) GenFromLogModuli() {

	if err := p.checkLogModuli(); err != nil {
		panic(err)
	}

	p.Qi, p.Pi, p.QiMul = GenModuli(p)

	p.GenFromModuli()
}

func (p *Parameters) checkModuli() error {

	if len(p.Qi) > MaxModuliCount {
		return fmt.Errorf("#LogQi is larger than %d", MaxModuliCount)
	}

	if len(p.Pi) > MaxModuliCount {
		return fmt.Errorf("#Pi is larger than %d", MaxModuliCount)
	}

	if len(p.QiMul) > MaxModuliCount {
		return fmt.Errorf("#QiMul is larger than %d", MaxModuliCount)
	}

	for i, qi := range p.Qi {
		if uint64(bits.Len64(qi)-1) > MaxModuliSize {
			return fmt.Errorf("Qi bit-size for i=%d is larger than %d", i, MaxModuliSize)
		}
	}

	for i, pi := range p.Pi {
		if uint64(bits.Len64(pi)-1) > MaxModuliSize {
			return fmt.Errorf("Pi bit-size for i=%d is larger than %d", i, MaxModuliSize)
		}
	}

	for i, qi := range p.QiMul {
		if uint64(bits.Len64(qi)-1) > MaxModuliSize {
			return fmt.Errorf("QiMul bitsize n°%d is larger than %d", i, MaxModuliSize)
		}
	}

	N := uint64(1 << p.LogN)

	for i, qi := range p.Qi {
		if !ring.IsPrime(qi) || qi&((N<<1)-1) != 1 {
			return fmt.Errorf("Qi n°%d is not an NTT prime", i)
		}
	}

	for i, pi := range p.Pi {
		if !ring.IsPrime(pi) || pi&((N<<1)-1) != 1 {
			return fmt.Errorf("Pi n°%d is not an NTT prime", i)
		}
	}

	for i, qi := range p.QiMul {
		if !ring.IsPrime(qi) || qi&((N<<1)-1) != 1 {
			return fmt.Errorf("QiMul n°%d is not an NTT prime", i)
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

	if len(p.LogQiMul) > MaxModuliCount {
		return fmt.Errorf("#LogQiMul is larger than %d", MaxModuliCount)
	}

	for i, qi := range p.LogQi {
		if qi > MaxModuliSize {
			return fmt.Errorf("LogQi for i=%d is larger than %d", i, MaxModuliSize)
		}
	}

	for i, pi := range p.LogPi {
		if pi > MaxModuliSize {
			return fmt.Errorf("LogPi for i=%d is larger than %d", i, MaxModuliSize)
		}
	}

	for i, qi := range p.LogQiMul {
		if qi > MaxModuliSize {
			return fmt.Errorf("LogQiMul for i=%d is larger than %d", i, MaxModuliSize)
		}
	}

	return nil
}
