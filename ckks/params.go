package ckks

import (
	"errors"
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

// Name of the different default parameter sets
const (
	// PN12QP109 is the index in DefaultParams for logQP = 109
	PN12QP109 = iota
	// PN13QP218 is the index in DefaultParams for logQP = 218
	PN13QP218
	// PN14QP438 is the index in DefaultParams for logQP = 438
	PN14QP438
	// PN15QP880 is the index in DefaultParams for logQP = 880
	PN15QP880
	// PN16QP1761 is the index in DefaultParams for logQP = 1761
	PN16QP1761

	// PN12QP101pq is the index in DefaultParams for logQP = 101 (post quantum)
	PN12QP101pq
	// PN13QP202pq is the index in DefaultParams for logQP = 202 (post quantum)
	PN13QP202pq
	// PN14QP411pq is the index in DefaultParams for logQP = 411 (post quantum)
	PN14QP411pq
	// PN15QP827pq is the index in DefaultParams for logQP = 827 (post quantum)
	PN15QP827pq
	// PN16QP1654pq is the index in DefaultParams for logQP = 1654 (post quantum)
	PN16QP1654pq
)

// DefaultParams is a set of default CKKS parameters ensuring 128 bit security.
var DefaultParams = []*Parameters{

	//LogQi = 109
	{logN: 12,
		logSlots: 11,
		qi: []uint64{0x200000e001, // 37 + 32
			0x100006001},
		pi:    []uint64{0x3ffffea001}, // 38
		scale: 1 << 32,
		sigma: DefaultSigma,
	},

	//LogQi = 218
	{logN: 13,
		logSlots: 12,
		qi: []uint64{0x1fffec001, // 33 + 5 x 30
			0x3fff4001,
			0x3ffe8001,
			0x40020001,
			0x40038001,
			0x3ffc0001},
		pi:    []uint64{0x800004001}, // 35
		scale: 1 << 30,
		sigma: DefaultSigma,
	},

	//LogQiP = 438
	{logN: 14,
		logSlots: 13,
		qi: []uint64{0x200000008001, 0x400018001, // 45 + 9 x 34
			0x3fffd0001, 0x400060001,
			0x400068001, 0x3fff90001,
			0x400080001, 0x4000a8001,
			0x400108001, 0x3ffeb8001},
		pi:    []uint64{0x7fffffd8001, 0x7fffffc8001}, // 43, 43
		scale: 1 << 34,
		sigma: DefaultSigma,
	},

	//LogQi = 880
	{logN: 15,
		logSlots: 14,
		qi: []uint64{0x4000000120001, 0x10000140001, 0xffffe80001, // 50 + 17 x 40
			0x10000290001, 0xffffc40001, 0x100003e0001,
			0x10000470001, 0x100004b0001, 0xffffb20001,
			0x10000500001, 0x10000650001, 0xffff940001,
			0xffff8a0001, 0xffff820001, 0xffff780001,
			0x10000890001, 0xffff750001, 0x10000960001},
		pi:    []uint64{0x40000001b0001, 0x3ffffffdf0001, 0x4000000270001}, // 50, 50, 50
		scale: 1 << 40,
		sigma: DefaultSigma,
	},

	//LogQi = 1761
	{logN: 16,
		logSlots: 15,
		qi: []uint64{0x80000000080001, 0x2000000a0001, 0x2000000e0001, 0x1fffffc20001, // 55 + 33 x 45
			0x200000440001, 0x200000500001, 0x200000620001, 0x1fffff980001,
			0x2000006a0001, 0x1fffff7e0001, 0x200000860001, 0x200000a60001,
			0x200000aa0001, 0x200000b20001, 0x200000c80001, 0x1fffff360001,
			0x200000e20001, 0x1fffff060001, 0x200000fe0001, 0x1ffffede0001,
			0x1ffffeca0001, 0x1ffffeb40001, 0x200001520001, 0x1ffffe760001,
			0x2000019a0001, 0x1ffffe640001, 0x200001a00001, 0x1ffffe520001,
			0x200001e80001, 0x1ffffe0c0001, 0x1ffffdee0001, 0x200002480001,
			0x1ffffdb60001, 0x200002560001},
		pi:    []uint64{0x80000000440001, 0x7fffffffba0001, 0x80000000500001, 0x7fffffffaa0001}, // 4 x 55
		scale: 1 << 45,
		sigma: DefaultSigma,
	},

	//LogQi = 101.00001186816735
	{logN: 12,
		logSlots: 11,
		qi:       []uint64{0x800004001, 0x40002001}, // 35 + 30
		pi:       []uint64{0x1000002001},            // 36
		scale:    1 << 30,
		sigma:    DefaultSigma,
	},

	//LogQi = 201.9936341352857
	{logN: 13,
		logSlots: 12,
		qi:       []uint64{0x1fffec001, 0x8008001, 0x8020001, 0x802c001, 0x7fa8001, 0x7f74001}, // 33 + 5 x 27
		pi:       []uint64{0x400018001},                                                        // 34
		scale:    1 << 27,
		sigma:    DefaultSigma,
	},

	//LogQiP = 411.0000787495673
	{logN: 14,
		logSlots: 13,
		qi: []uint64{0x10000048001, 0x200038001, 0x1fff90001, 0x200080001, 0x1fff60001,
			0x2000b8001, 0x200100001, 0x1fff00001, 0x1ffef0001, 0x200128001}, // 40 + 9 x 33

		pi:    []uint64{0x1ffffe0001, 0x1ffffc0001}, // 37, 37
		scale: 1 << 33,
		sigma: DefaultSigma,
	},

	//LogQi = 827.0000771955918
	{logN: 15,
		logSlots: 14,
		qi: []uint64{0x400000060001, 0x4000170001, 0x3fffe80001, 0x40002f0001, 0x4000300001,
			0x3fffcf0001, 0x40003f0001, 0x3fffc10001, 0x4000450001, 0x3fffb80001,
			0x3fffb70001, 0x40004a0001, 0x3fffb20001, 0x4000510001, 0x3fffaf0001,
			0x4000540001, 0x4000560001, 0x4000590001}, // 46 + 17 x 38
		pi:    []uint64{0x2000000a0001, 0x2000000e0001, 0x2000001d0001}, // 3 x 45
		scale: 1 << 38,
		sigma: DefaultSigma,
	},

	//LogQi = 1653.999999
	{logN: 16,
		logSlots: 15,
		qi: []uint64{0x80000000080001, 0x2000000a0001, 0x2000000e0001, 0x1fffffc20001, 0x200000440001,
			0x200000500001, 0x200000620001, 0x1fffff980001, 0x2000006a0001, 0x1fffff7e0001,
			0x200000860001, 0x200000a60001, 0x200000aa0001, 0x200000b20001, 0x200000c80001,
			0x1fffff360001, 0x200000e20001, 0x1fffff060001, 0x200000fe0001, 0x1ffffede0001,
			0x1ffffeca0001, 0x1ffffeb40001, 0x200001520001, 0x1ffffe760001, 0x2000019a0001,
			0x1ffffe640001, 0x200001a00001, 0x1ffffe520001, 0x200001e80001, 0x1ffffe0c0001,
			0x1ffffdee0001, 0x200002480001}, // 55 + 31 x 45
		pi:    []uint64{0x7fffffffe0001, 0x80000001c0001, 0x80000002c0001, 0x7ffffffd20001}, // 4 x 51
		scale: 1 << 45,
		sigma: DefaultSigma,
	},
}

// Moduli stores the NTT primes of the RNS representation.
type Moduli struct {
	Qi []uint64 // Ciphertext prime moduli
	Pi []uint64 // Keys additional prime moduli
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

// Copy creates a copy of the target Moduli.
func (m *Moduli) Copy() Moduli {

	qi := make([]uint64, len(m.Qi))
	copy(qi, m.Qi)

	pi := make([]uint64, len(m.Pi))
	copy(pi, m.Pi)

	return Moduli{qi, pi}
}

// LogModuli stores the bit-length of the NTT primes of the RNS representation.
type LogModuli struct {
	LogQi []int // Ciphertext prime moduli bit-size
	LogPi []int // Keys additional prime moduli bit-size
}

// Copy creates a copy of the target Moduli.
func (m *LogModuli) Copy() LogModuli {

	LogQi := make([]int, len(m.LogQi))
	copy(LogQi, m.LogQi)

	LogPi := make([]int, len(m.LogPi))
	copy(LogPi, m.LogPi)

	return LogModuli{LogQi, LogPi}
}

// Parameters represents a given parameter set for the CKKS cryptosystem.
type Parameters struct {
	qi       []uint64
	pi       []uint64
	logN     int // Ring degree (power of 2)
	logSlots int
	scale    float64
	sigma    float64 // Gaussian sampling variance
}

// NewParametersFromModuli creates a new Parameters struct and returns a pointer to it.
func NewParametersFromModuli(logN int, m *Moduli) (p *Parameters, err error) {
	p = new(Parameters)

	if (logN < MinLogN) || (logN > MaxLogN) {
		return nil, fmt.Errorf("invalid polynomial ring log degree: %d", logN)
	}

	p.logN = logN

	if err = checkModuli(m, logN); err != nil {
		return nil, err
	}

	p.qi = make([]uint64, len(m.Qi))
	copy(p.qi, m.Qi)
	p.pi = make([]uint64, len(m.Pi))
	copy(p.pi, m.Pi)

	p.sigma = DefaultSigma

	return p, nil

}

// NewParametersFromLogModuli creates a new Parameters struct and returns a pointer to it.
func NewParametersFromLogModuli(logN int, lm *LogModuli) (p *Parameters, err error) {

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
func (p *Parameters) N() int {
	return 1 << p.logN
}

// LogN returns the log of the degree of the polynomial ring
func (p *Parameters) LogN() int {
	return p.logN
}

// LogSlots returns the log of the number of slots
func (p *Parameters) LogSlots() int {
	return p.logSlots
}

// MaxLevel returns the maximum ciphertext level
func (p *Parameters) MaxLevel() int {
	return p.QiCount() - 1
}

// Levels returns then number of total levels enabled by the parameters
func (p *Parameters) Levels() int {
	return p.QiCount()
}

// Slots returns number of available plaintext slots
func (p *Parameters) Slots() int {
	return 1 << p.logSlots
}

// MaxSlots returns the theoretical maximum of plaintext slots allowed by the ring degree
func (p *Parameters) MaxSlots() int {
	return p.N() >> 1
}

// MaxLogSlots returns the log of the maximum number of slots enabled by the parameters
func (p *Parameters) MaxLogSlots() int {
	return p.logN - 1
}

// Sigma returns standard deviation of the noise distribution
func (p *Parameters) Sigma() float64 {
	return p.sigma
}

// Scale returns the default plaintext/ciphertext scale
func (p *Parameters) Scale() float64 {
	return p.scale
}

// SetScale sets the default plaintext/ciphertext scale
func (p *Parameters) SetScale(scale float64) {
	p.scale = scale
}

// SetLogSlots sets the value logSlots of the parameters.
func (p *Parameters) SetLogSlots(logSlots int) {
	if (logSlots == 0) || (logSlots > p.MaxLogSlots()) {
		panic(fmt.Errorf("slots cannot be greater than LogN-1"))
	}

	p.logSlots = logSlots
}

// SetSigma sets the value sigma of the parameters
func (p *Parameters) SetSigma(sigma float64) {
	p.sigma = sigma
}

// LogModuli generates a LogModuli struct from the parameters' Moduli struct and returns it.
func (p *Parameters) LogModuli() (lm *LogModuli) {

	lm = new(LogModuli)

	lm.LogQi = make([]int, len(p.qi))
	for i := range p.qi {
		lm.LogQi[i] = int(math.Round(math.Log2(float64(p.qi[i]))))
	}

	lm.LogPi = make([]int, len(p.pi))
	for i := range p.pi {
		lm.LogPi[i] = int(math.Round(math.Log2(float64(p.pi[i]))))
	}

	return
}

// QiOverflowMargin returns floor(2^64 / max(Qi)), i.e. the number of times elements of Z_max{Qi} can
// be added together before overflowing 2^64.
func (p *Parameters) QiOverflowMargin(level int) int {
	return int(math.Exp2(64) / float64(utils.MaxSliceUint64(p.qi[:level+1])))
}

// PiOverflowMargin returns floor(2^64 / max(Pi)), i.e. the number of times elements of Z_max{Pi} can
// be added together before overflowing 2^64.
func (p *Parameters) PiOverflowMargin() int {
	return int(math.Exp2(64) / float64(utils.MaxSliceUint64(p.pi)))
}

// Moduli returns a struct Moduli with the moduli of the parameters
func (p *Parameters) Moduli() (m *Moduli) {
	m = new(Moduli)
	m.Qi = make([]uint64, p.QiCount())
	m.Pi = make([]uint64, p.PiCount())
	copy(m.Qi, p.qi)
	copy(m.Pi, p.pi)
	return
}

// Qi returns a new slice with the factors of the ciphertext modulus q
func (p *Parameters) Qi() []uint64 {
	qi := make([]uint64, len(p.qi))
	copy(qi, p.qi)
	return qi
}

// QiCount returns the number of factors of the ciphertext modulus Q
func (p *Parameters) QiCount() int {
	return len(p.qi)
}

// Pi returns a new slice with the factors of the ciphertext modulus extension P
func (p *Parameters) Pi() []uint64 {
	pi := make([]uint64, len(p.pi))
	copy(pi, p.pi)
	return pi
}

// PiCount returns the number of factors of the ciphertext modulus extension P
func (p *Parameters) PiCount() int {
	return len(p.pi)
}

// QPiCount returns the number of factors of the ciphertext modulus + the modulus extension P
func (p *Parameters) QPiCount() int {
	return len(p.qi) + len(p.pi)
}

// LogQP returns the size of the extended modulus QP in bits
func (p *Parameters) LogQP() int {
	tmp := ring.NewUint(1)
	for _, qi := range p.qi {
		tmp.Mul(tmp, ring.NewUint(qi))
	}
	for _, pi := range p.pi {
		tmp.Mul(tmp, ring.NewUint(pi))
	}
	return tmp.BitLen()
}

// LogQLvl returns the size of the modulus Q in bits at a specific level
func (p *Parameters) LogQLvl(level int) int {
	tmp := p.QLvl(level)
	return tmp.BitLen()
}

// QLvl returns the product of the moduli at the given level as a big.Int
func (p *Parameters) QLvl(level int) *big.Int {
	tmp := ring.NewUint(1)
	for _, qi := range p.qi[:level+1] {
		tmp.Mul(tmp, ring.NewUint(qi))
	}
	return tmp
}

// LogQ returns the size of the modulus Q in bits
func (p *Parameters) LogQ() int {
	return p.LogQLvl(p.QiCount() - 1)
}

// Q returns the product of all the moduli as a big.Int
func (p *Parameters) Q() *big.Int {
	return p.QLvl(p.QiCount() - 1)
}

// LogP returns the size of the modulus P in bits
func (p *Parameters) LogP() int {
	tmp := ring.NewUint(1)
	for _, pi := range p.pi {
		tmp.Mul(tmp, ring.NewUint(pi))
	}
	return tmp.BitLen()
}

// LogQAlpha returns the size in bits of the sum of the norm of
// each element of the special RNS decomposition basis for the
// key-switching.
// LogQAlpha is the size of the element that is multiplied by the
// error during the keyswitching and then divided by P.
// LogQAlpha should be smaller than P or the error added during
// the key-switching wont be negligible.
func (p *Parameters) LogQAlpha() int {

	alpha := p.PiCount()

	if alpha == 0 {
		return 0
	}

	res := ring.NewUint(0)
	var j int
	for i := 0; i < p.QiCount(); i = i + alpha {

		j = i + alpha
		if j > p.QiCount() {
			j = p.QiCount()
		}

		tmp := ring.NewUint(1)
		for _, qi := range p.pi[i:j] {
			tmp.Mul(tmp, ring.NewUint(qi))
		}

		res.Add(res, tmp)
	}

	return res.BitLen()
}

// Alpha returns the number of moduli in in P
func (p *Parameters) Alpha() int {
	return p.PiCount()
}

// Beta returns the number of element in the RNS decomposition basis: Ceil(lenQi / lenPi)
func (p *Parameters) Beta() int {
	if p.Alpha() != 0 {
		return int(math.Ceil(float64(p.QiCount()) / float64(p.Alpha())))
	}

	return 0
}

// GaloisElementForColumnRotationBy returns the galois element for plaintext
// column rotations by k position to the left. Providing a negative k is
// equivalent to a right rotation.
func (p *Parameters) GaloisElementForColumnRotationBy(k int) uint64 {
	twoN := 1 << (p.logN + 1)
	mask := twoN - 1
	kRed := k & mask
	return ring.ModExp(uint64(GaloisGen), kRed, uint64(twoN))
}

// GaloisElementForRowRotation returns the galois element corresponding to a row rotation (conjugate) automorphism
func (p *Parameters) GaloisElementForRowRotation() uint64 {
	return (1 << (p.logN + 1)) - 1
}

// GaloisElementsForRowInnerSum returns a list of galois element corresponding to
// all the left rotations by a k-position where k is a power of two.
func (p *Parameters) GaloisElementsForRowInnerSum() (galEls []uint64) {
	galEls = make([]uint64, p.logN+1, p.logN+1)
	galEls[0] = p.GaloisElementForRowRotation()
	for i := 0; i < int(p.logN)-1; i++ {
		galEls[i+1] = p.GaloisElementForColumnRotationBy(1 << i)
	}
	return galEls
}

// InverseGaloisElement returns the galois element for the inverse automorphism of galEl
func (p *Parameters) InverseGaloisElement(galEl uint64) uint64 {
	twoN := 1 << (p.logN + 1)
	return ring.ModExp(galEl, twoN-1, uint64(twoN))
}

// Copy creates a copy of the target parameters.
func (p *Parameters) Copy() (paramsCopy *Parameters) {

	paramsCopy = new(Parameters)
	paramsCopy.logN = p.logN
	paramsCopy.logSlots = p.logSlots
	paramsCopy.scale = p.scale
	paramsCopy.sigma = p.sigma
	paramsCopy.qi = make([]uint64, len(p.qi))
	copy(paramsCopy.qi, p.qi)
	paramsCopy.pi = make([]uint64, len(p.pi))
	copy(paramsCopy.pi, p.pi)
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
	return
}

// MarshalBinary returns a []byte representation of the parameter set.
func (p *Parameters) MarshalBinary() ([]byte, error) {
	if p.logN == 0 { // if N is 0, then p is the zero value
		return []byte{}, nil
	}

	// Data 21 byte + QPiCount * 8 byte:
	// 1 byte : logN
	// 1 byte : logSlots
	// 8 byte : scale
	// 8 byte : sigma
	// 1 byte : #qi
	// 1 byte : #pi
	b := utils.NewBuffer(make([]byte, 0, 20+(p.QPiCount())<<3))

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

	if len(data) < 20 {
		return errors.New("invalid parameters encoding")
	}

	b := utils.NewBuffer(data)

	p.logN = int(b.ReadUint8())

	if p.logN > MaxLogN {
		return fmt.Errorf("LogN larger than %d", MaxLogN)
	}

	p.logSlots = int(b.ReadUint8())

	if p.logSlots > p.logN-1 {
		return fmt.Errorf("LogSlots larger than %d", MaxLogN-1)
	}

	p.scale = math.Float64frombits(b.ReadUint64())
	p.sigma = math.Float64frombits(b.ReadUint64())

	lenQi := b.ReadUint8()
	lenPi := b.ReadUint8()

	p.qi = make([]uint64, lenQi)
	p.pi = make([]uint64, lenPi)

	b.ReadUint64Slice(p.qi)
	b.ReadUint64Slice(p.pi)

	if err = checkModuli(p.Moduli(), p.logN); err != nil {
		return err
	}

	return nil
}

func checkModuli(m *Moduli, logN int) error {

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

	N := 1 << logN

	for i, qi := range m.Qi {
		if !ring.IsPrime(qi) || qi&uint64((N<<1)-1) != 1 {
			return fmt.Errorf("Qi (i=%d) is not an NTT prime", i)
		}
	}

	for i, pi := range m.Pi {
		if !ring.IsPrime(pi) || pi&uint64((N<<1)-1) != 1 {
			return fmt.Errorf("Pi (i=%d) is not an NTT prime", i)
		}
	}

	return nil
}

func checkLogModuli(m *LogModuli) error {

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

func genModuli(lm *LogModuli, logN int) (m *Moduli) {

	m = new(Moduli)

	// Extracts all the different primes bit size and maps their number
	primesbitlen := make(map[int]int)
	for _, qi := range lm.LogQi {
		primesbitlen[qi]++
	}

	for _, pj := range lm.LogPi {
		primesbitlen[pj]++
	}

	// For each bit-size, finds that many primes
	primes := make(map[int][]uint64)
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
