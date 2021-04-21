package rlwe

import (
	"fmt"
	"math"
	"math/big"
	"math/bits"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// MaxLogN is the log2 of the largest supported polynomial modulus degree.
const MaxLogN = 16

// MinLogN is the log2 of the smallest supported polynomial modulus degree (needed to ensure the NTT correctness).
const MinLogN = 4

// MaxModuliCount is the largest supported number of moduli in the RNS representation.
const MaxModuliCount = 34

// MaxModuliSize is the largest bit-length supported for the moduli in the RNS representation.
const MaxModuliSize = 60

// DefaultSigma is the default error distribution standard deviation
const DefaultSigma = 3.2

// GaloisGen is an integer of order N=2^d modulo M=2N and that spans Z_M with the integer -1.
// The j-th ring automorphism takes the root zeta to zeta^(5j).
const GaloisGen uint64 = 5

type Parameters struct {
	logN  uint64
	qi    []uint64
	pi    []uint64
	sigma float64
}

// Moduli stores the NTT primes of the RNS representation.
type Moduli struct {
	Qi []uint64 // Ciphertext prime moduli
	Pi []uint64 // Keys additional prime moduli
}

// Copy creates a copy of the target Moduli.
func (m *Moduli) Copy() Moduli {

	qi := make([]uint64, len(m.Qi))
	copy(qi, m.Qi)

	pi := make([]uint64, len(m.Pi))
	copy(pi, m.Pi)

	return Moduli{qi, pi}
}

// Print prints the moduli in hexadecimal
func (m *Moduli) Print() {
	for _, qi := range m.Qi {
		fmt.Printf("0x%x,\n", qi)
	}
	fmt.Println()

	for _, pj := range m.Pi {
		fmt.Printf("0x%x,\n", pj)
	}
	fmt.Println()
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

func NewRLWEParameters(logn uint64, q, p []uint64, sigma float64) (Parameters, error) { // TEMPORARY constructor

	if (logn < MinLogN) || (logn > MaxLogN) {
		return Parameters{}, fmt.Errorf("invalid polynomial ring log degree: %d", logn)
	}

	// Checks if Moduli is valid
	if err := CheckModuli(&Moduli{q, p}, logn); err != nil {
		return Parameters{}, err
	}

	params := Parameters{
		logN:  logn,
		pi:    make([]uint64, len(p)),
		qi:    make([]uint64, len(q)),
		sigma: sigma,
	}
	copy(params.qi, q)
	copy(params.pi, p)
	return params, nil
}

// N returns the ring degree
func (p *Parameters) N() uint64 {
	return 1 << p.logN
}

// LogN returns the log of the degree of the polynomial ring
func (p *Parameters) LogN() uint64 {
	return p.logN
}

// Sigma returns standard deviation of the noise distribution
func (p *Parameters) Sigma() float64 {
	return p.sigma
}

// Qi returns a new slice with the factors of the ciphertext modulus q
func (p *Parameters) Q() []uint64 {
	qi := make([]uint64, len(p.qi))
	copy(qi, p.qi)
	return qi
}

// QCount returns the number of factors of the ciphertext modulus Q
func (p *Parameters) QCount() uint64 {
	return uint64(len(p.qi))
}

func (p *Parameters) QBigInt() *big.Int {
	q := big.NewInt(1)
	for _, qi := range p.qi {
		q.Mul(q, new(big.Int).SetUint64(qi))
	}
	return q
}

// Pi returns a new slice with the factors of the ciphertext modulus extension P
func (p *Parameters) P() []uint64 {
	pi := make([]uint64, len(p.pi))
	copy(pi, p.pi)
	return pi
}

// PCount returns the number of factors of the ciphertext modulus extension P
func (p *Parameters) PCount() uint64 {
	return uint64(len(p.pi))
}

func (p *Parameters) PBigInt() *big.Int {
	pInt := big.NewInt(1)
	for _, pi := range p.pi {
		pInt.Mul(pInt, new(big.Int).SetUint64(pi))
	}
	return pInt
}

func (p *Parameters) QP() []uint64 {
	qp := make([]uint64, len(p.qi)+len(p.pi))
	copy(qp, p.qi)
	copy(qp[len(p.qi):], p.pi)
	return qp
}

// QPiCount returns the number of factors of the ciphertext modulus + the modulus extension P
func (p *Parameters) QPCount() uint64 {
	return uint64(len(p.qi) + len(p.pi))
}

func (p *Parameters) QPBigInt() *big.Int {
	pqInt := p.QBigInt()
	pqInt.Mul(pqInt, p.PBigInt())
	return pqInt
}

// LogQP returns the size of the extended modulus QP in bits
func (p *Parameters) LogQP() uint64 {
	tmp := ring.NewUint(1)
	for _, qi := range p.qi {
		tmp.Mul(tmp, ring.NewUint(qi))
	}
	for _, pi := range p.pi {
		tmp.Mul(tmp, ring.NewUint(pi))
	}
	return uint64(tmp.BitLen())
}

// Alpha returns the number of moduli in in P
func (p *Parameters) Alpha() uint64 {
	return p.PCount()
}

// Beta returns the number of element in the RNS decomposition basis: Ceil(lenQi / lenPi)
func (p *Parameters) Beta() uint64 {
	if p.Alpha() != 0 {
		return uint64(math.Ceil(float64(p.QCount()) / float64(p.Alpha())))
	}

	return 0
}

// LogModuli generates a LogModuli struct from the parameters' Moduli struct and returns it.
func (p *Parameters) LogModuli() (lm *LogModuli) {
	lm = new(LogModuli)
	lm.LogQi = make([]uint64, len(p.Q()))
	for i := range p.Q() {
		lm.LogQi[i] = uint64(math.Round(math.Log2(float64(p.Q()[i]))))
	}
	lm.LogPi = make([]uint64, len(p.P()))
	for i := range p.P() {
		lm.LogPi[i] = uint64(math.Round(math.Log2(float64(p.P()[i]))))
	}
	return
}

// Moduli returns a struct Moduli with the moduli of the parameters
func (p *Parameters) Moduli() (m *Moduli) {
	m = new(Moduli)
	m.Qi = make([]uint64, len(p.Q()))
	copy(m.Qi, p.Q())
	m.Pi = make([]uint64, len(p.P()))
	copy(m.Pi, p.P())
	return
}

func (p *Parameters) RingQ() *ring.Ring {
	ringQ, err := ring.NewRing(p.N(), p.qi)
	if err != nil {
		panic(err) // Parameter type invariant
	}
	return ringQ
}

func (p *Parameters) RingP() *ring.Ring {
	if len(p.pi) == 0 {
		return nil
	}
	ringP, err := ring.NewRing(p.N(), p.pi)
	if err != nil {
		panic(err) // Parameter type invariant
	}
	return ringP
}

func (p *Parameters) RingQP() *ring.Ring {
	ringQP, err := ring.NewRing(p.N(), append(p.qi, p.pi...))
	if err != nil {
		panic(err) // Parameter type invariant
	}
	return ringQP
}

// GaloisElementForColumnRotationBy returns the galois element for plaintext
// column rotations by k position to the left. Providing a negative k is
// equivalent to a right rotation.
func (p *Parameters) GaloisElementForColumnRotationBy(k int) uint64 {
	twoN := 1 << (p.logN + 1)
	mask := twoN - 1
	kRed := uint64(k & mask)
	return ring.ModExp(GaloisGen, kRed, uint64(twoN))
}

// GaloisElementForRowRotation returns the galois element for generating the row
// rotation automorphism
func (p *Parameters) GaloisElementForRowRotation() uint64 {
	return (1 << (p.logN + 1)) - 1
}

// GaloisElementsForRowInnerSum returns a list of all galois elements required to
// perform an InnerSum operation. This corresponds to all the left rotations by
// k-positions where k is a power of two and the row-rotation element.
func (p *Parameters) GaloisElementsForRowInnerSum() (galEls []uint64) {
	galEls = make([]uint64, p.logN+1, p.logN+1)
	galEls[0] = p.GaloisElementForRowRotation()
	for i := 0; i < int(p.logN)-1; i++ {
		galEls[i+1] = p.GaloisElementForColumnRotationBy(1 << i)
	}
	return galEls
}

// InverseGaloisElement takes a galois element and returns the galois element
//  corresponding to the inverse automorphism
func (p *Parameters) InverseGaloisElement(galEl uint64) uint64 {
	twoN := uint64(1 << (p.logN + 1))
	return ring.ModExp(galEl, twoN-1, twoN)
}

func CheckModuli(m *Moduli, logN uint64) error {

	if len(m.Qi) > MaxModuliCount {
		return fmt.Errorf("#Qi is larger than %d", MaxModuliCount)
	}

	if len(m.Pi) > MaxModuliCount {
		return fmt.Errorf("#Pi is larger than %d", MaxModuliCount)
	}

	for i, qi := range m.Qi {
		if uint64(bits.Len64(qi)-1) > MaxModuliSize+1 {
			return fmt.Errorf("Qi bit-size (i=%d) is larger than %d", i, MaxModuliSize)
		}
	}

	for i, pi := range m.Pi {
		if uint64(bits.Len64(pi)-1) > MaxModuliSize+2 {
			return fmt.Errorf("Pi bit-size (i=%d) is larger than %d", i, MaxModuliSize)
		}
	}

	N := uint64(1 << logN)

	for i, qi := range m.Qi {
		if !ring.IsPrime(qi) || qi&((N<<1)-1) != 1 {
			return fmt.Errorf("Qi (i=%d) is not an NTT prime", i)
		}
	}

	for i, pi := range m.Pi {
		if !ring.IsPrime(pi) || pi&((N<<1)-1) != 1 {
			return fmt.Errorf("Pi (i=%d) is not an NTT prime", i)
		}
	}

	return nil
}

func (p *Parameters) MarshalBinary() ([]byte, error) {
	if p.LogN() == 0 { // if N is 0, then p is the zero value
		return []byte{}, nil
	}

	// 1 byte : logN
	// 1 byte : #Q
	// 1 byte : #P
	// 8 byte : sigma
	// 8 * (#Q) : Q
	// 8 * (#P) : P
	b := utils.NewBuffer(make([]byte, 0, 11+(len(p.qi)+len(p.pi))<<3))
	b.WriteUint8(uint8(p.logN))
	b.WriteUint8(uint8(len(p.qi)))
	b.WriteUint8(uint8(len(p.pi)))
	b.WriteUint64(math.Float64bits(p.sigma))
	b.WriteUint64Slice(p.qi)
	b.WriteUint64Slice(p.pi)
	return b.Bytes(), nil
}

func (params *Parameters) UnmarshalBinary(data []byte) error {
	if len(data) < 11 {
		return fmt.Errorf("invalid rlwe.Parameter serialization")
	}
	b := utils.NewBuffer(data)
	logN := uint64(b.ReadUint8())
	lenQ := uint64(b.ReadUint8())
	lenP := uint64(b.ReadUint8())
	sigma := math.Float64frombits(b.ReadUint64())

	if err := CheckSizeParams(logN, lenQ, lenP); err != nil {
		return err
	}

	q := make([]uint64, lenQ)
	p := make([]uint64, lenP)
	b.ReadUint64Slice(q)
	b.ReadUint64Slice(p)

	var err error
	*params, err = NewRLWEParameters(logN, q, p, sigma)
	return err
}

func (p *Parameters) Equals(other Parameters) bool {
	res := p.logN == other.logN
	res = res && utils.EqualSliceUint64(p.qi, other.qi)
	res = res && utils.EqualSliceUint64(p.pi, other.pi)
	res = res && (p.sigma == other.sigma)
	return res
}

func CheckSizeParams(logN, lenQ, lenP uint64) error {
	if logN > MaxLogN {
		return fmt.Errorf("logN=%d is larger than MaxLogN=%d", logN, MaxLogN)
	}
	if logN < MinLogN {
		return fmt.Errorf("logN=%d is smaller than MinLogN=%d", logN, MinLogN)
	}
	if lenQ > MaxModuliCount {
		return fmt.Errorf("lenQ=%d is larger than MaxModuliCount=%d", lenQ, MaxModuliCount)
	}
	if lenP > MaxModuliCount {
		return fmt.Errorf("lenP=%d is larger than MaxModuliCount=%d", lenP, MaxModuliCount)
	}
	return nil
}

func CheckLogModuli(m *LogModuli) error {

	// Checks if the parameters are empty
	if m.LogQi == nil || len(m.LogQi) == 0 {
		return fmt.Errorf("nil or empty slice provided as LogModuli.LogQi")
	}

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

func GenModuli(lm *LogModuli, logN uint64) (m *Moduli) {

	m = new(Moduli)

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
		primes[key] = ring.GenerateNTTPrimes(key, 2<<logN, value)
	}

	// Assigns the primes to the CKKS moduli chain
	m.Qi = make([]uint64, len(lm.LogQi))
	for i, qi := range lm.LogQi {
		m.Qi[i] = primes[qi][0]
		primes[qi] = primes[qi][1:]
	}

	// Assigns the primes to the special primes list for the extended ring
	m.Pi = make([]uint64, len(lm.LogPi))
	for i, pj := range lm.LogPi {
		m.Pi[i] = primes[pj][0]
		primes[pj] = primes[pj][1:]
	}

	return m
}
