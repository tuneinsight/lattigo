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

// DefaultParams is a set of default BFV parameters ensuring 128 bit security.
var DefaultParams = []*Parameters{

	//logQ1+P = 109
	{
		LogN: 12,
		N:    4096,
		T:    65537,
		Moduli: Moduli{
			Qi:    []uint64{0x7ffffec001, 0x8000016001},             // 39 + 39 bits
			Pi:    []uint64{0x40002001},                             // 30 bits
			QiMul: []uint64{0xfffffffffffc001, 0x100000000000e001}}, // 60 + 60 bits
		Sigma: 3.2,
		Alpha: 1,
		Beta:  2,
	},

	//logQ1+P = 218
	{
		LogN: 13,
		N:    8192,
		T:    65537,
		Moduli: Moduli{
			Qi:    []uint64{0x3fffffffef8001, 0x4000000011c001, 0x40000000120001},      // 54 + 54 + 54 bits
			Pi:    []uint64{0x7ffffffffb4001},                                          // 55 bits
			QiMul: []uint64{0xfffffffffffc001, 0xffffffffffe8001, 0x1000000000024001}}, // 60 + 60 + 60 bits
		Sigma: 3.2,
		Alpha: 1,
		Beta:  3,
	},

	//logQ1+P = 438
	{
		LogN: 14,
		N:    16384,
		T:    65537,
		Moduli: Moduli{
			Qi: []uint64{0x100000000060001, 0x80000000068001, 0x80000000080001,
				0x3fffffffef8001, 0x40000000120001, 0x3fffffffeb8001}, // 56 + 55 + 55 + 54 + 54 + 54 bits
			Pi: []uint64{0x80000000130001, 0x7fffffffe90001}, // 55 + 55 bits
			QiMul: []uint64{0xffffffffffe8001, 0xffffffffffd8001, 0xffffffffffc0001,
				0x1000000000078001, 0xffffffffff28001, 0xfffffffffe38001}}, // 60 + 60 + 60 + 60 + 60 + 60 bits
		Sigma: 3.2,
		Alpha: 2,
		Beta:  3,
	},

	//logQ1+P = 880
	{
		LogN: 15,
		N:    32768,
		T:    65537,
		Moduli: Moduli{
			Qi: []uint64{0x7ffffffffe70001, 0x7ffffffffe10001, 0x7ffffffffcc0001, // 59 + 59 + 59 bits
				0x400000000270001, 0x400000000350001, 0x400000000360001, // 58 + 58 + 58 bits
				0x3ffffffffc10001, 0x3ffffffffbe0001, 0x3ffffffffbd0001, // 58 + 58 + 58 bits
				0x4000000004d0001, 0x400000000570001, 0x400000000660001}, // 58 + 58 + 58 bits
			Pi: []uint64{0xffffffffffc0001, 0x10000000001d0001, 0x10000000006e0001}, // 60 + 60 + 60 bits
			QiMul: []uint64{0xfffffffff840001, 0x1000000000860001, 0x1000000000870001, // 60 + 60 + 60 bits
				0x1000000000930001, 0xfffffffff6a0001, 0x1000000000980001, // 60 + 60 + 60 bits
				0xfffffffff5a0001, 0xfffffffff550001, 0x1000000000b00001, // 60 + 60 + 60 bits
				0xfffffffff330001, 0x1000000000ce0001, 0xfffffffff2a0001}}, // 60 + 60 + 60 bits
		Sigma: 3.2,
		Alpha: 3,
		Beta:  4,
	},
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
	LogN  uint64  // Log Ring degree (power of 2)
	N     uint64  // Ring degree
	T     uint64  // Plaintext modulus
	Sigma float64 // Gaussian sampling standard deviation

	LogQP uint64
	Alpha uint64
	Beta  uint64
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
	paramsCopy.N = p.N
	paramsCopy.T = p.T
	paramsCopy.Sigma = p.Sigma
	paramsCopy.Moduli = p.Moduli.Copy()
	paramsCopy.LogQP = p.LogQP
	paramsCopy.Alpha = p.Alpha
	paramsCopy.Beta = p.Beta

	return
}

// Equals compares two sets of parameters for equality.
func (p *Parameters) Equals(other *Parameters) (res bool) {

	if p == other {
		return true
	}

	res = p.LogN == other.LogN
	res = res && (p.N == other.N)
	res = res && (p.T == other.T)
	res = res && (p.Sigma == other.Sigma)

	res = res && utils.EqualSliceUint64(p.Qi, other.Qi)
	res = res && utils.EqualSliceUint64(p.Pi, other.Pi)
	res = res && utils.EqualSliceUint64(p.QiMul, other.QiMul)

	res = res && (p.Alpha == other.Alpha)
	res = res && (p.Beta == other.Beta)
	res = res && (p.LogQP == other.LogQP)

	return
}

// MarshalBinary returns a []byte representation of the parameter set.
func (p *Parameters) MarshalBinary() ([]byte, error) {
	if p.LogN == 0 { // if N is 0, then p is the zero value
		return []byte{}, nil
	}

	b := utils.NewBuffer(make([]byte, 0, 20+(len(p.Qi)+len(p.Pi)+len(p.QiMul))<<3))

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
	p.N = 1 << p.LogN

	if p.LogN > MaxLogN {
		return fmt.Errorf("LogN larger than %d", MaxLogN)
	}

	lenQi := b.ReadUint8()
	lenPi := b.ReadUint8()
	lenQiMul := b.ReadUint8()

	if lenPi != 0 {
		p.Alpha = uint64(lenPi)
		p.Beta = uint64(math.Ceil(float64(lenQi) / float64(lenPi)))
	}

	p.T = b.ReadUint64()
	p.Sigma = math.Round((float64(b.ReadUint64())/float64(1<<32))*100) / 100
	p.Qi = make([]uint64, lenQi, lenQi)
	p.Pi = make([]uint64, lenPi, lenPi)
	p.QiMul = make([]uint64, lenQiMul, lenQiMul)

	b.ReadUint64Slice(p.Qi)
	b.ReadUint64Slice(p.Pi)
	b.ReadUint64Slice(p.QiMul)

	return checkModuli(&p.Moduli, p.LogN) // TODO: check more than moduli.
}

// GenFromLogModuli creates a new Parameters struct and returns a pointer to it.
func GenFromLogModuli(logN uint64, lm *LogModuli) (p *Parameters, err error) {

	if logN < 0 || logN > MaxLogN {
		return nil, fmt.Errorf("invalid polynomial ring log degree: %d", logN)
	}

	if err = checkLogModuli(lm); err != nil {
		return nil, err
	}

	// If LogModuli is valid and then generates the moduli
	p.Moduli = genModuli(lm, logN)

	// Checks if Moduli is valid
	if err = checkModuli(&p.Moduli, logN); err != nil {
		return nil, err
	}

	tmp := ring.NewUint(1)

	for _, qi := range p.Qi {
		tmp.Mul(tmp, ring.NewUint(qi))
	}

	for _, pi := range p.Pi {
		tmp.Mul(tmp, ring.NewUint(pi))
	}

	p.LogQP = uint64(tmp.BitLen())

	p.N = 1 << p.LogN
	if len(p.Pi) != 0 {
		p.Alpha = uint64(len(p.Pi))
		p.Beta = uint64(math.Ceil(float64(len(p.Qi)) / float64(len(p.Pi))))
	}

	return p, nil
}

// GetLogModuli generates a LogModuli struct from the parameters' Moduli struct and returns it.
func (p *Parameters) GetLogModuli() LogModuli {
	var lm LogModuli
	lm.LogQi = make([]uint64, len(p.Qi), len(p.Qi))
	for i := range p.Qi {
		lm.LogQi[i] = uint64(bits.Len64(p.Qi[i]) - 1)
	}

	lm.LogPi = make([]uint64, len(p.Pi), len(p.Pi))
	for i := range p.Pi {
		lm.LogPi[i] = uint64(bits.Len64(p.Pi[i]) - 1)
	}

	lm.LogQiMul = make([]uint64, len(p.QiMul), len(p.QiMul))
	for i := range p.QiMul {
		lm.LogQiMul[i] = uint64(bits.Len64(p.QiMul[i]) - 1)
	}
	return lm
}

func checkModuli(m *Moduli, logN uint64) (err error) {

	if len(m.Qi) > MaxModuliCount {
		return fmt.Errorf("#LogQi is larger than %d", MaxModuliCount)
	}

	if len(m.Pi) > MaxModuliCount {
		return fmt.Errorf("#Pi is larger than %d", MaxModuliCount)
	}

	if len(m.QiMul) > MaxModuliCount {
		return fmt.Errorf("#QiMul is larger than %d", MaxModuliCount)
	}

	for i, qi := range m.Qi {
		if uint64(bits.Len64(qi)-1) > MaxModuliSize+1 {
			return fmt.Errorf("Qi bit-size for i=%d is larger than %d", i, MaxModuliSize)
		}
	}

	for i, pi := range m.Pi {
		if uint64(bits.Len64(pi)-1) > MaxModuliSize+1 {
			return fmt.Errorf("Pi bit-size for i=%d is larger than %d", i, MaxModuliSize)
		}
	}

	for i, qi := range m.QiMul {
		if uint64(bits.Len64(qi)-1) > MaxModuliSize+1 {
			return fmt.Errorf("QiMul bitsize n°%d is larger than %d", i, MaxModuliSize)
		}
	}

	N := uint64(1 << logN)

	for i, qi := range m.Qi {
		if !ring.IsPrime(qi) || qi&((N<<1)-1) != 1 {
			return fmt.Errorf("Qi n°%d is not an NTT prime", i)
		}
	}

	for i, pi := range m.Pi {
		if !ring.IsPrime(pi) || pi&((N<<1)-1) != 1 {
			return fmt.Errorf("Pi n°%d is not an NTT prime", i)
		}
	}

	for i, qi := range m.QiMul {
		if !ring.IsPrime(qi) || qi&((N<<1)-1) != 1 {
			return fmt.Errorf("QiMul n°%d is not an NTT prime", i)
		}
	}

	return nil
}

func checkLogModuli(lm *LogModuli) (err error) {

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

// GenModuli generates the appropriate primes from the parameters using generateCKKSPrimes such that all primes are different.
func genModuli(lm *LogModuli, logN uint64) (m Moduli) {

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
	m.Qi = make([]uint64, len(lm.LogQi))
	for i, qi := range lm.LogQi {
		m.Qi[i] = primes[qi][0]
		primes[qi] = primes[qi][1:]
	}

	// Assigns the primes to the special primes list for the the keys context
	m.Pi = make([]uint64, len(lm.LogPi))
	for i, pj := range lm.LogPi {
		m.Pi[i] = primes[pj][0]
		primes[pj] = primes[pj][1:]
	}

	m.QiMul = make([]uint64, len(lm.LogQiMul))
	for i, qi := range lm.LogQiMul {
		m.QiMul[i] = primes[qi][0]
		primes[qi] = primes[qi][1:]
	}

	return m
}
