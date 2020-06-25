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

		if err := params.Gen(); err != nil {
			panic(err)
		}

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
	LogN  uint64  // Log Ring degree (power of 2)
	N     uint64  // Ring degree
	T     uint64  // Plaintext modulus
	Sigma float64 // Gaussian sampling standard deviation

	LogQP uint64
	Alpha uint64
	Beta  uint64

	isValid bool
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
	paramsCopy.N = p.N
	paramsCopy.T = p.T
	paramsCopy.Sigma = p.Sigma
	paramsCopy.Moduli = p.Moduli.Copy()
	paramsCopy.LogModuli = p.LogModuli.Copy()
	paramsCopy.LogQP = p.LogQP
	paramsCopy.Alpha = p.Alpha
	paramsCopy.Beta = p.Beta
	paramsCopy.isValid = p.isValid

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
	res = res && utils.EqualSliceUint64(p.LogQi, other.LogQi)
	res = res && utils.EqualSliceUint64(p.LogPi, other.LogPi)
	res = res && utils.EqualSliceUint64(p.LogQiMul, other.LogQiMul)

	res = res && (p.Alpha == other.Alpha)
	res = res && (p.Beta == other.Beta)
	res = res && (p.LogQP == other.LogQP)

	res = res && (p.isValid == other.isValid)

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
	p.N = 1 << p.LogN

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

	return p.Gen()
}

// Gen populates the parameter structs fromt the provided parameters.
func (p *Parameters) Gen() (err error) {

	// Checks if the parameters are empty
	if (len(p.Qi) + len(p.Pi) + len(p.QiMul) + len(p.LogQi) + len(p.LogPi) + len(p.LogQiMul)) == 0 {
		return fmt.Errorf("cannot p.Gen() -> Moduli & LogModuli are both empty -> must set one of them")
	}

	// Checks if both Moduli and LogModuli are set
	if (len(p.Qi)+len(p.Pi)+len(p.QiMul) != 0) && (len(p.LogQi)+len(p.LogPi)+len(p.LogQiMul) != 0) {
		return fmt.Errorf("warning Moduli & LogModuli are both set -> LogModuli will be overwritten")
	}

	// If Moduli is not set, then checks if LogModuli is valid and then generates the moduli
	if len(p.Qi)+len(p.Pi)+len(p.QiMul) == 0 {

		if err = p.checkLogModuli(); err != nil {
			return err
		}

		p.Qi, p.Pi, p.QiMul = GenModuli(p)
	}

	// Checks if Moduli is valid
	if err = p.checkModuli(); err != nil {
		return err
	}

	tmp := ring.NewUint(1)

	for _, qi := range p.Qi {
		tmp.Mul(tmp, ring.NewUint(qi))
	}

	for _, pi := range p.Pi {
		tmp.Mul(tmp, ring.NewUint(pi))
	}

	p.LogQP = uint64(tmp.BitLen())

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

	p.N = 1 << p.LogN
	if len(p.LogPi) != 0 {
		p.Alpha = uint64(len(p.Pi))
		p.Beta = uint64(math.Ceil(float64(len(p.Qi)) / float64(len(p.Pi))))
	}

	p.isValid = true

	return nil
}

func (p *Parameters) checkModuli() (err error) {

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
		if uint64(bits.Len64(qi)-1) > MaxModuliSize+1 {
			return fmt.Errorf("Qi bit-size for i=%d is larger than %d", i, MaxModuliSize)
		}
	}

	for i, pi := range p.Pi {
		if uint64(bits.Len64(pi)-1) > MaxModuliSize+1 {
			return fmt.Errorf("Pi bit-size for i=%d is larger than %d", i, MaxModuliSize)
		}
	}

	for i, qi := range p.QiMul {
		if uint64(bits.Len64(qi)-1) > MaxModuliSize+1 {
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

func (p *Parameters) checkLogModuli() (err error) {

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

// GenModuli generates the appropriate primes from the parameters using generateCKKSPrimes such that all primes are different.
func GenModuli(params *Parameters) (Q []uint64, P []uint64, QMul []uint64) {

	// Extracts all the different primes bit-size and maps their number
	primesbitlen := make(map[uint64]uint64)

	for _, qi := range params.LogQi {
		primesbitlen[qi]++
	}

	for _, pj := range params.LogPi {
		primesbitlen[pj]++
	}

	for _, qi := range params.LogQiMul {
		primesbitlen[qi]++
	}

	// For each bit-size, it finds that many primes
	primes := make(map[uint64][]uint64)
	for key, value := range primesbitlen {
		primes[key] = ring.GenerateNTTPrimes(key, params.LogN, value)
	}

	// Assigns the primes to the CKKS moduli chain
	Q = make([]uint64, len(params.LogQi))
	for i, qi := range params.LogQi {
		Q[i] = primes[qi][0]
		primes[qi] = primes[qi][1:]
	}

	// Assigns the primes to the special primes list for the the keys context
	P = make([]uint64, len(params.LogPi))
	for i, pj := range params.LogPi {
		P[i] = primes[pj][0]
		primes[pj] = primes[pj][1:]
	}

	QMul = make([]uint64, len(params.LogQiMul))
	for i, qi := range params.LogQiMul {
		QMul[i] = primes[qi][0]
		primes[qi] = primes[qi][1:]
	}

	return Q, P, QMul
}
