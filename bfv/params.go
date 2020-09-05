package bfv

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
	// PN12QP109 is a set of parameters with N = 2^12 and log(QP) = 109
	PN12QP109 = iota
	// PN13QP218 is a set of parameters with N = 2^13 and log(QP) = 218
	PN13QP218
	// PN14QP438 is a set of parameters with N = 2^14 and log(QP) = 438
	PN14QP438
	// PN15QP880 is a set of parameters with N = 2^15 and log(QP) = 880
	PN15QP880
)

// DefaultSigma is the default error distribution standard deviation
const DefaultSigma = 3.2

// DefaultParams is a set of default BFV parameters ensuring 128 bit security.
var DefaultParams = []*Parameters{

	{
		logN:  12,
		n:     4096,
		t:     65537,
		logQP: 109,
		Moduli: Moduli{
			qi:    []uint64{0x7ffffec001, 0x8000016001},             // 39 + 39 bits
			pi:    []uint64{0x40002001},                             // 30 bits
			qiMul: []uint64{0xfffffffffffc001, 0x100000000000e001}}, // 60 + 60 bits
		sigma: DefaultSigma,
		alpha: 1,
		beta:  2,
	},

	{
		logN:  13,
		n:     8192,
		t:     65537,
		logQP: 218,
		Moduli: Moduli{
			qi:    []uint64{0x3fffffffef8001, 0x4000000011c001, 0x40000000120001},      // 54 + 54 + 54 bits
			pi:    []uint64{0x7ffffffffb4001},                                          // 55 bits
			qiMul: []uint64{0xfffffffffffc001, 0xffffffffffe8001, 0x1000000000024001}}, // 60 + 60 + 60 bits
		sigma: DefaultSigma,
		alpha: 1,
		beta:  3,
	},

	{
		logN:  14,
		n:     16384,
		t:     65537,
		logQP: 438,
		Moduli: Moduli{
			qi: []uint64{0x100000000060001, 0x80000000068001, 0x80000000080001,
				0x3fffffffef8001, 0x40000000120001, 0x3fffffffeb8001}, // 56 + 55 + 55 + 54 + 54 + 54 bits
			pi: []uint64{0x80000000130001, 0x7fffffffe90001}, // 55 + 55 bits
			qiMul: []uint64{0xffffffffffe8001, 0xffffffffffd8001, 0xffffffffffc0001,
				0x1000000000078001, 0xffffffffff28001, 0xfffffffffe38001}}, // 60 + 60 + 60 + 60 + 60 + 60 bits
		sigma: DefaultSigma,
		alpha: 2,
		beta:  3,
	},

	{
		logN:  15,
		n:     32768,
		t:     65537,
		logQP: 880,
		Moduli: Moduli{
			qi: []uint64{0x7ffffffffe70001, 0x7ffffffffe10001, 0x7ffffffffcc0001, // 59 + 59 + 59 bits
				0x400000000270001, 0x400000000350001, 0x400000000360001, // 58 + 58 + 58 bits
				0x3ffffffffc10001, 0x3ffffffffbe0001, 0x3ffffffffbd0001, // 58 + 58 + 58 bits
				0x4000000004d0001, 0x400000000570001, 0x400000000660001}, // 58 + 58 + 58 bits
			pi: []uint64{0xffffffffffc0001, 0x10000000001d0001, 0x10000000006e0001}, // 60 + 60 + 60 bits
			qiMul: []uint64{0xfffffffff840001, 0x1000000000860001, 0x1000000000870001, // 60 + 60 + 60 bits
				0x1000000000930001, 0xfffffffff6a0001, 0x1000000000980001, // 60 + 60 + 60 bits
				0xfffffffff5a0001, 0xfffffffff550001, 0x1000000000b00001, // 60 + 60 + 60 bits
				0xfffffffff330001, 0x1000000000ce0001, 0xfffffffff2a0001}}, // 60 + 60 + 60 bits
		sigma: DefaultSigma,
		alpha: 3,
		beta:  4,
	},
}

// Moduli stores the NTT primes of the RNS representation.
type Moduli struct {
	qi    []uint64 // Ciphertext prime moduli
	pi    []uint64 // Keys additional prime moduli
	qiMul []uint64 // Ciphertext secondary prime moduli
}

// Qi returns a new slice with the factors of the ciphertext modulus q
func (m *Moduli) Qi() []uint64 {
	qi := make([]uint64, len(m.qi))
	copy(qi, m.qi)
	return qi
}

// QiCount returns the number of factors of the ciphertext modulus q
func (m *Moduli) QiCount() int {
	return len(m.qi)
}

// Pi returns a new slice with the factors of the ciphertext modulus extention p
func (m *Moduli) Pi() []uint64 {
	pi := make([]uint64, len(m.pi))
	copy(pi, m.pi)
	return pi
}

// PiCount returns the number of factors of the ciphertext modulus extention p
func (m *Moduli) PiCount() int {
	return len(m.pi)
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

	qiMul := make([]uint64, len(m.qiMul))
	copy(qiMul, m.qiMul)

	return Moduli{qi, pi, qiMul}
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
	logN  uint64  // Log Ring degree (power of 2)
	n     uint64  // Ring degree
	t     uint64  // Plaintext modulus
	sigma float64 // Gaussian sampling standard deviation

	logQP uint64
	alpha uint64
	beta  uint64
}

// NewParametersFromModuli creates a new Parameters struct and returns a pointer to it.
func NewParametersFromModuli(logN uint64, m Moduli, t uint64) (p *Parameters, err error) {

	p = new(Parameters)

	if logN < 0 || logN > MaxLogN {
		return nil, fmt.Errorf("invalid polynomial ring log degree: %d", logN)
	}

	p.logN = logN

	// Checks if Moduli is valid
	if err = checkModuli(m, logN); err != nil {
		return nil, err
	}

	p.Moduli = m.Copy()

	p.logQP = p.Moduli.LogQP()

	p.n = 1 << p.logN

	if len(p.pi) != 0 {
		p.alpha = uint64(len(p.pi))
		p.beta = uint64(math.Ceil(float64(len(p.qi)) / float64(len(p.pi))))
	}

	p.sigma = DefaultSigma

	p.t = t

	return p, nil
}

// NewParametersFromLogModuli creates a new Parameters struct and returns a pointer to it.
func NewParametersFromLogModuli(logN uint64, lm LogModuli, t uint64) (p *Parameters, err error) {

	if err = checkLogModuli(lm); err != nil {
		return nil, err
	}

	// If LogModuli is valid and then generates the moduli
	return NewParametersFromModuli(logN, genModuli(lm, logN), t)
}

// LogN returns the log of the degree of the polynomial ring
func (p *Parameters) LogN() uint64 {
	return p.logN
}

// T returns the plaintext coefficient modulus t
func (p *Parameters) T() uint64 {
	return p.t
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

// LogQP returns the size of the extanded modulus QP in bits
func (p *Parameters) SetT(T uint64) {
	p.t = T
}

func (p *Parameters) WithT(T uint64) (pCopy *Parameters) {
	pCopy = p.Copy()
	pCopy.SetT(T)
	return
}

// LogModuli generates a LogModuli struct from the parameters' Moduli struct and returns it.
func (p *Parameters) LogModuli() LogModuli {
	var lm LogModuli
	lm.LogQi = make([]uint64, len(p.qi), len(p.qi))
	for i := range p.qi {
		lm.LogQi[i] = uint64(math.Round(math.Log2(float64(p.qi[i]))))
	}

	lm.LogPi = make([]uint64, len(p.pi), len(p.pi))
	for i := range p.pi {
		lm.LogPi[i] = uint64(math.Round(math.Log2(float64(p.pi[i]))))
	}

	lm.LogQiMul = make([]uint64, len(p.qiMul), len(p.qiMul))
	for i := range p.qiMul {
		lm.LogQiMul[i] = uint64(math.Round(math.Log2(float64(p.qiMul[i]))))
	}
	return lm
}

// NewPolyQ returns a new empty polynomial of degree 2^logN in basis qi.
func (p *Parameters) NewPolyQ() *ring.Poly {
	return ring.NewPoly(p.n, uint64(len(p.qi)))
}

// NewPolyP returns a new empty polynomial of degree 2^logN in basis Pi.
func (p *Parameters) NewPolyP() *ring.Poly {
	return ring.NewPoly(p.n, uint64(len(p.pi)))
}

// NewPolyQP returns a new empty polynomial of degree 2^logN in basis qi + Pi.
func (p *Parameters) NewPolyQP() *ring.Poly {
	return ring.NewPoly(p.n, uint64(len(p.qi)+len(p.pi)))
}

// Copy creates a copy of the target Parameters.
func (p *Parameters) Copy() (paramsCopy *Parameters) {

	paramsCopy = new(Parameters)
	paramsCopy.logN = p.logN
	paramsCopy.n = p.n
	paramsCopy.t = p.t
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
	res = res && (p.n == other.n)
	res = res && (p.t == other.t)
	res = res && (p.sigma == other.sigma)

	res = res && utils.EqualSliceUint64(p.qi, other.qi)
	res = res && utils.EqualSliceUint64(p.pi, other.pi)
	res = res && utils.EqualSliceUint64(p.qiMul, other.qiMul)

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

	b := utils.NewBuffer(make([]byte, 0, 20+(len(p.qi)+len(p.pi)+len(p.qiMul))<<3))

	b.WriteUint8(uint8(p.logN))
	b.WriteUint8(uint8(len(p.qi)))
	b.WriteUint8(uint8(len(p.pi)))
	b.WriteUint8(uint8(len(p.qiMul)))
	b.WriteUint64(p.t)
	b.WriteUint64(uint64(p.sigma * (1 << 32)))
	b.WriteUint64Slice(p.qi)
	b.WriteUint64Slice(p.pi)
	b.WriteUint64Slice(p.qiMul)

	return b.Bytes(), nil
}

// UnmarshalBinary decodes a []byte into a parameter set struct.
func (p *Parameters) UnmarshalBinary(data []byte) error {
	if len(data) < 3 {
		return errors.New("invalid parameters encoding")
	}
	b := utils.NewBuffer(data)

	p.logN = uint64(b.ReadUint8())
	p.n = 1 << p.logN

	if p.logN > MaxLogN {
		return fmt.Errorf("logN larger than %d", MaxLogN)
	}

	lenQi := b.ReadUint8()
	lenPi := b.ReadUint8()
	lenQiMul := b.ReadUint8()

	if lenPi != 0 {
		p.alpha = uint64(lenPi)
		p.beta = uint64(math.Ceil(float64(lenQi) / float64(lenPi)))
	}

	p.t = b.ReadUint64()
	p.sigma = math.Round((float64(b.ReadUint64())/float64(1<<32))*100) / 100
	p.qi = make([]uint64, lenQi, lenQi)
	p.pi = make([]uint64, lenPi, lenPi)
	p.qiMul = make([]uint64, lenQiMul, lenQiMul)

	b.ReadUint64Slice(p.qi)
	b.ReadUint64Slice(p.pi)
	b.ReadUint64Slice(p.qiMul)

	err := checkModuli(p.Moduli, p.logN) // TODO: check more than moduli.
	if err != nil {
		return err
	}

	p.logQP = p.Moduli.LogQP()

	return nil
}

func checkModuli(m Moduli, logN uint64) (err error) {

	if len(m.qi) > MaxModuliCount {
		return fmt.Errorf("#qi is larger than %d", MaxModuliCount)
	}

	if len(m.pi) > MaxModuliCount {
		return fmt.Errorf("#Pi is larger than %d", MaxModuliCount)
	}

	if len(m.qiMul) > MaxModuliCount {
		return fmt.Errorf("#qiMul is larger than %d", MaxModuliCount)
	}

	for i, qi := range m.qi {
		if uint64(bits.Len64(qi)-1) > MaxModuliSize+1 {
			return fmt.Errorf("qi bit-size for i=%d is larger than %d", i, MaxModuliSize)
		}
	}

	for i, pi := range m.pi {
		if uint64(bits.Len64(pi)-1) > MaxModuliSize+1 {
			return fmt.Errorf("Pi bit-size for i=%d is larger than %d", i, MaxModuliSize)
		}
	}

	for i, qi := range m.qiMul {
		if uint64(bits.Len64(qi)-1) > MaxModuliSize+1 {
			return fmt.Errorf("qiMul bitsize n°%d is larger than %d", i, MaxModuliSize)
		}
	}

	N := uint64(1 << logN)

	for i, qi := range m.qi {
		if !ring.IsPrime(qi) || qi&((N<<1)-1) != 1 {
			return fmt.Errorf("qi n°%d is not an NTT prime", i)
		}
	}

	for i, pi := range m.pi {
		if !ring.IsPrime(pi) || pi&((N<<1)-1) != 1 {
			return fmt.Errorf("Pi n°%d is not an NTT prime", i)
		}
	}

	for i, qi := range m.qiMul {
		if !ring.IsPrime(qi) || qi&((N<<1)-1) != 1 {
			return fmt.Errorf("qiMul n°%d is not an NTT prime", i)
		}
	}

	return nil
}

func checkLogModuli(lm LogModuli) (err error) {

	// Checks if the parameters are empty
	if lm.LogQi == nil || len(lm.LogQi) == 0 {
		return fmt.Errorf("nil or empty slice provided as LogModuli.LogQi") // TODO: are our algorithm working with empty mult basis ?
	}

	if len(lm.LogQi) > MaxModuliCount {
		return fmt.Errorf("#LogQi is larger than %d", MaxModuliCount)
	}

	if len(lm.LogPi) > MaxModuliCount {
		return fmt.Errorf("#LogPi is larger than %d", MaxModuliCount)
	}

	if len(lm.LogQiMul) > MaxModuliCount {
		return fmt.Errorf("#LogQiMul is larger than %d", MaxModuliCount)
	}

	for i, qi := range lm.LogQi {
		if qi > MaxModuliSize {
			return fmt.Errorf("LogQi for i=%d is larger than %d", i, MaxModuliSize)
		}
	}

	for i, pi := range lm.LogPi {
		if pi > MaxModuliSize {
			return fmt.Errorf("LogPi for i=%d is larger than %d", i, MaxModuliSize)
		}
	}

	for i, qi := range lm.LogQiMul {
		if qi > MaxModuliSize {
			return fmt.Errorf("LogQiMul for i=%d is larger than %d", i, MaxModuliSize)
		}
	}

	return nil
}

// GenModuli generates the appropriate primes from the parameters using generateNTTPrimes such that all primes are different.
func genModuli(lm LogModuli, logN uint64) (m Moduli) {

	// Extracts all the different primes bit-size and maps their number
	primesbitlen := make(map[uint64]uint64)

	for _, qi := range lm.LogQi {
		primesbitlen[qi]++
	}

	for _, pj := range lm.LogPi {
		primesbitlen[pj]++
	}

	for _, qi := range lm.LogQiMul {
		primesbitlen[qi]++
	}

	// For each bit-size, it finds that many primes
	primes := make(map[uint64][]uint64)
	for key, value := range primesbitlen {
		primes[key] = ring.GenerateNTTPrimes(key, logN, value)
	}

	// Assigns the primes to the CKKS moduli chain
	m.qi = make([]uint64, len(lm.LogQi))
	for i, qi := range lm.LogQi {
		m.qi[i] = primes[qi][0]
		primes[qi] = primes[qi][1:]
	}

	// Assigns the primes to the special primes list for the the keys context
	m.pi = make([]uint64, len(lm.LogPi))
	for i, pj := range lm.LogPi {
		m.pi[i] = primes[pj][0]
		primes[pj] = primes[pj][1:]
	}

	m.qiMul = make([]uint64, len(lm.LogQiMul))
	for i, qi := range lm.LogQiMul {
		m.qiMul[i] = primes[qi][0]
		primes[qi] = primes[qi][1:]
	}

	return m
}
