package ckks

import (
	"errors"
	"fmt"
	"math"
	"math/big"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

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
var DefaultParams = []*ParametersDef{

	//LogQi = 109
	{logN: 12,
		logSlots: 11,
		qi: []uint64{0x200000e001, // 37 + 32
			0x100006001},
		pi:    []uint64{0x3ffffea001}, // 38
		scale: 1 << 32,
		sigma: rlwe.DefaultSigma,
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
		sigma: rlwe.DefaultSigma,
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
		sigma: rlwe.DefaultSigma,
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
		sigma: rlwe.DefaultSigma,
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
		sigma: rlwe.DefaultSigma,
	},

	//LogQi = 101.00001186816735
	{logN: 12,
		logSlots: 11,
		qi:       []uint64{0x800004001, 0x40002001}, // 35 + 30
		pi:       []uint64{0x1000002001},            // 36
		scale:    1 << 30,
		sigma:    rlwe.DefaultSigma,
	},

	//LogQi = 201.9936341352857
	{logN: 13,
		logSlots: 12,
		qi:       []uint64{0x1fffec001, 0x8008001, 0x8020001, 0x802c001, 0x7fa8001, 0x7f74001}, // 33 + 5 x 27
		pi:       []uint64{0x400018001},                                                        // 34
		scale:    1 << 27,
		sigma:    rlwe.DefaultSigma,
	},

	//LogQiP = 411.0000787495673
	{logN: 14,
		logSlots: 13,
		qi: []uint64{0x10000048001, 0x200038001, 0x1fff90001, 0x200080001, 0x1fff60001,
			0x2000b8001, 0x200100001, 0x1fff00001, 0x1ffef0001, 0x200128001}, // 40 + 9 x 33

		pi:    []uint64{0x1ffffe0001, 0x1ffffc0001}, // 37, 37
		scale: 1 << 33,
		sigma: rlwe.DefaultSigma,
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
		sigma: rlwe.DefaultSigma,
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
		sigma: rlwe.DefaultSigma,
	},
}

type Parameters interface {
	rlwe.Parameters

	LogQLvl(level int) int
	LogSlots() int
	MaxLevel() int
	MaxLogSlots() int
	QLvl(level int) *big.Int
	Scale() float64
	Slots() int
}

type ParametersDef struct {
	qi       []uint64
	pi       []uint64
	logN     int // Ring degree (power of 2)
	logSlots int
	scale    float64
	sigma    float64 // Gaussian sampling variance
}

// ParametersStruct represents a given parameter set for the CKKS cryptosystem.
type ParametersStruct struct {
	rlwe.ParametersStruct

	logSlots int
	scale    float64
}

func NewParametersFromParamDef(paramDef *ParametersDef) (*ParametersStruct, error) {
	m := new(rlwe.Moduli)
	m.Qi = make([]uint64, len(paramDef.qi))
	copy(m.Qi, paramDef.qi)
	m.Pi = make([]uint64, len(paramDef.pi))
	copy(m.Pi, paramDef.pi)
	return NewParametersFromModuli(paramDef.logN, m, paramDef.sigma, paramDef.logSlots, paramDef.scale)
}

// NewParametersFromModuli creates a new Parameters struct and returns a pointer to it.
func NewParametersFromModuli(logN int, m *rlwe.Moduli, sigma float64, logSlots int, scale float64) (p *ParametersStruct, err error) {

	var rlweParams *rlwe.ParametersStruct
	if rlweParams, err = rlwe.NewRLWEParameters(logN, m.Qi, m.Pi, sigma); err != nil {
		return nil, err
	}

	return &ParametersStruct{*rlweParams, logSlots, scale}, nil

}

// NewParametersFromLogModuli creates a new Parameters struct and returns a pointer to it.
func NewParametersFromLogModuli(logN int, lm *rlwe.LogModuli, sigma float64, logSlots int, scale float64) (p *ParametersStruct, err error) {

	if err = rlwe.CheckLogModuli(lm); err != nil {
		return nil, err
	}

	// If LogModuli is valid and then generates the moduli
	return NewParametersFromModuli(logN, rlwe.GenModuli(lm, logN), sigma, logSlots, scale)
}

func GetDefaultParameters(paramsId int) *ParametersStruct {
	if paramsId >= len(DefaultParams) {
		panic(fmt.Errorf("paramsId %d does not exist", paramsId))
	}
	params, err := NewParametersFromParamDef(DefaultParams[paramsId])
	if err != nil {
		panic(err)
	}
	return params
}

// // NewPolyQ returns a new empty polynomial of degree 2^LogN in basis Qi.
// func (p *ParametersStruct) NewPolyQ() *ring.Poly {
// 	return ring.NewPoly(p.N(), p.QCount())
// }

// // NewPolyP returns a new empty polynomial of degree 2^LogN in basis Pi.
// func (p *ParametersStruct) NewPolyP() *ring.Poly {
// 	return ring.NewPoly(p.N(), p.PCount())
// }

// // NewPolyQP returns a new empty polynomial of degree 2^LogN in basis Qi + Pi.
// func (p *ParametersStruct) NewPolyQP() *ring.Poly {
// 	return ring.NewPoly(p.N(), p.QPCount())
// }

// // N returns the ring degree
// func (p *ParametersStruct) N() uint64 {
// 	return 1 << p.logN
// }

// // LogN returns the log of the degree of the polynomial ring
// func (p *ParametersStruct) LogN() uint64 {
// 	return p.logN
// }

// LogSlots returns the log of the number of slots
func (p *ParametersStruct) LogSlots() int {
	return p.logSlots
}

// MaxLevel returns the maximum ciphertext level
func (p *ParametersStruct) MaxLevel() int {
	return p.QCount() - 1
}

// // Levels returns then number of total levels enabled by the parameters
// func (p *ParametersStruct) Levels() uint64 {
// 	return p.QCount()
// }

// Slots returns number of available plaintext slots
func (p *ParametersStruct) Slots() int {
	return 1 << p.logSlots
}

// MaxSlots returns the theoretical maximum of plaintext slots allowed by the ring degree
func (p *ParametersStruct) MaxSlots() int {
	return p.N() >> 1
}

// MaxLogSlots returns the log of the maximum number of slots enabled by the parameters
func (p *ParametersStruct) MaxLogSlots() int {
	return p.LogN() - 1
}

// // Sigma returns standard deviation of the noise distribution
// func (p *ParametersStruct) Sigma() float64 {
// 	return p.sigma
// }

// Scale returns the default plaintext/ciphertext scale
func (p *ParametersStruct) Scale() float64 {
	return p.scale
}

// // SetScale sets the default plaintext/ciphertext scale
// func (p *ParametersStruct) SetScale(scale float64) {
// 	p.scale = scale
// }

// // SetLogSlots sets the value logSlots of the parameters.
// func (p *ParametersStruct) SetLogSlots(logSlots uint64) {
// 	if (logSlots == 0) || (logSlots > p.MaxLogSlots()) {
// 		panic(fmt.Errorf("slots cannot be greater than LogN-1"))
// 	}

// 	p.logSlots = logSlots
// }

// // SetSigma sets the value sigma of the parameters
// func (p *ParametersStruct) SetSigma(sigma float64) {
// 	p.sigma = sigma
// }

// // LogModuli generates a LogModuli struct from the parameters' Moduli struct and returns it.
// func (p *ParametersStruct) LogModuli() (lm *rlwe.LogModuli) {

// 	lm = new(rlwe.LogModuli)

// 	lm.LogQi = make([]uint64, len(p.qi))
// 	for i := range p.qi {
// 		lm.LogQi[i] = uint64(math.Round(math.Log2(float64(p.qi[i]))))
// 	}

// 	lm.LogPi = make([]uint64, len(p.pi))
// 	for i := range p.pi {
// 		lm.LogPi[i] = uint64(math.Round(math.Log2(float64(p.pi[i]))))
// 	}

// 	return
// }

// // Moduli returns a struct Moduli with the moduli of the parameters
// func (p *ParametersStruct) Moduli() (m *rlwe.Moduli) {
// 	m = new(rlwe.Moduli)
// 	m.Qi = make([]uint64, p.QCount())
// 	m.Pi = make([]uint64, p.PCount())
// 	copy(m.Qi, p.qi)
// 	copy(m.Pi, p.pi)
// 	return
// }

// // Q returns a new slice with the factors of the ciphertext modulus q
// func (p *ParametersStruct) Q() []uint64 {
// 	qi := make([]uint64, len(p.qi))
// 	copy(qi, p.qi)
// 	return qi
// }

// // QCount returns the number of factors of the ciphertext modulus Q
// func (p *ParametersStruct) QCount() uint64 {
// 	return uint64(len(p.qi))
// }

// // P returns a new slice with the factors of the ciphertext modulus extension P
// func (p *ParametersStruct) P() []uint64 {
// 	pi := make([]uint64, len(p.pi))
// 	copy(pi, p.pi)
// 	return pi
// }

// func (p *ParametersStruct) QP() []uint64 {
// 	qp := make([]uint64, len(p.qi)+len(p.pi))
// 	copy(qp, p.qi)
// 	copy(qp[len(p.qi):], p.pi)
// 	return qp
// }

// // PCount returns the number of factors of the ciphertext modulus extension P
// func (p *ParametersStruct) PCount() uint64 {
// 	return uint64(len(p.pi))
// }

// // QPCount returns the number of factors of the ciphertext modulus + the modulus extension P
// func (p *ParametersStruct) QPCount() uint64 {
// 	return uint64(len(p.qi) + len(p.pi))
// }

// // LogQP returns the size of the extended modulus QP in bits
// func (p *ParametersStruct) LogQP() uint64 {
// 	tmp := ring.NewUint(1)
// 	for _, qi := range p.qi {
// 		tmp.Mul(tmp, ring.NewUint(qi))
// 	}
// 	for _, pi := range p.pi {
// 		tmp.Mul(tmp, ring.NewUint(pi))
// 	}
// 	return uint64(tmp.BitLen())
// }

// LogQLvl returns the size of the modulus Q in bits at a specific level
func (p *ParametersStruct) LogQLvl(level int) int {
	tmp := p.QLvl(level)
	return tmp.BitLen()
}

// QLvl returns the product of the moduli at the given level as a big.Int
func (p *ParametersStruct) QLvl(level int) *big.Int {
	tmp := ring.NewUint(1)
	for _, qi := range p.Q()[:level+1] {
		tmp.Mul(tmp, ring.NewUint(qi))
	}
	return tmp
}

// // LogQ returns the size of the modulus Q in bits
// func (p *ParametersStruct) LogQ() uint64 {
// 	return p.LogQLvl(p.QCount() - 1)
// }

// // QBigInt returns the product of all the moduli as a big.Int
// func (p *ParametersStruct) QBigInt() *big.Int {
// 	return p.QLvl(p.QCount() - 1)
// }

// func (p *ParametersStruct) PBigInt() *big.Int {
// 	pInt := big.NewInt(1)
// 	for _, pi := range p.pi {
// 		pInt.Mul(pInt, new(big.Int).SetUint64(pi))
// 	}
// 	return pInt
// }

// func (p *ParametersStruct) QPBigInt() *big.Int {
// 	pqInt := p.QBigInt()
// 	pqInt.Mul(pqInt, p.PBigInt())
// 	return pqInt
// }

// // LogP returns the size of the modulus P in bits
// func (p *ParametersStruct) LogP() uint64 {
// 	tmp := ring.NewUint(1)
// 	for _, pi := range p.pi {
// 		tmp.Mul(tmp, ring.NewUint(pi))
// 	}
// 	return uint64(tmp.BitLen())
// }

// // LogQAlpha returns the size in bits of the sum of the norm of
// // each element of the special RNS decomposition basis for the
// // key-switching.
// // LogQAlpha is the size of the element that is multiplied by the
// // error during the keyswitching and then divided by P.
// // LogQAlpha should be smaller than P or the error added during
// // the key-switching wont be negligible.
// func (p *ParametersStruct) LogQAlpha() uint64 {

// 	alpha := p.PCount()

// 	if alpha == 0 {
// 		return 0
// 	}

// 	res := ring.NewUint(0)
// 	var j uint64
// 	for i := uint64(0); i < p.QCount(); i = i + alpha {

// 		j = i + alpha
// 		if j > p.QCount() {
// 			j = p.QCount()
// 		}

// 		tmp := ring.NewUint(1)
// 		for _, qi := range p.pi[i:j] {
// 			tmp.Mul(tmp, ring.NewUint(qi))
// 		}

// 		res.Add(res, tmp)
// 	}

// 	return uint64(res.BitLen())
// }

// // Alpha returns the number of moduli in in P
// func (p *ParametersStruct) Alpha() uint64 {
// 	return p.PCount()
// }

// // Beta returns the number of element in the RNS decomposition basis: Ceil(lenQi / lenPi)
// func (p *ParametersStruct) Beta() uint64 {
// 	if p.Alpha() != 0 {
// 		return uint64(math.Ceil(float64(p.QCount()) / float64(p.Alpha())))
// 	}

// 	return 0
// }

// // GaloisElementForColumnRotationBy returns the galois element for plaintext
// // column rotations by k position to the left. Providing a negative k is
// // equivalent to a right rotation.
// func (p *ParametersStruct) GaloisElementForColumnRotationBy(k int) uint64 {
// 	twoN := 1 << (p.logN + 1)
// 	mask := twoN - 1
// 	kRed := uint64(k & mask)
// 	return ring.ModExp(GaloisGen, kRed, uint64(twoN))
// }

// // GaloisElementForRowRotation returns the galois element corresponding to a row rotation (conjugate) automorphism
// func (p *ParametersStruct) GaloisElementForRowRotation() uint64 {
// 	return (1 << (p.logN + 1)) - 1
// }

// // GaloisElementsForRowInnerSum returns a list of all galois elements required to
// // perform an InnerSum operation. This corresponds to all the left rotations by
// // k-positions where k is a power of two and the conjugate element.
// func (p *ParametersStruct) GaloisElementsForRowInnerSum() (galEls []uint64) {
// 	galEls = make([]uint64, p.logN+1, p.logN+1)
// 	galEls[0] = p.GaloisElementForRowRotation()
// 	for i := 0; i < int(p.logN)-1; i++ {
// 		galEls[i+1] = p.GaloisElementForColumnRotationBy(1 << i)
// 	}
// 	return galEls
// }

// // InverseGaloisElement returns the galois element for the inverse automorphism of galEl
// func (p *ParametersStruct) InverseGaloisElement(galEl uint64) uint64 {
// 	twoN := uint64(1 << (p.logN + 1))
// 	return ring.ModExp(galEl, twoN-1, twoN)
// }

// // Copy creates a copy of the target parameters.
// func (p *ParametersStruct) Copy() (paramsCopy *ParametersStruct) {

// 	paramsCopy = new(ParametersStruct)
// 	paramsCopy.logN = p.logN
// 	paramsCopy.logSlots = p.logSlots
// 	paramsCopy.scale = p.scale
// 	paramsCopy.sigma = p.sigma
// 	paramsCopy.qi = make([]uint64, len(p.qi))
// 	copy(paramsCopy.qi, p.qi)
// 	paramsCopy.pi = make([]uint64, len(p.pi))
// 	copy(paramsCopy.pi, p.pi)
// 	return
// }

// Equals compares two sets of parameters for equality.
func (p *ParametersStruct) Equals(other Parameters) (res bool) {
	res = p.LogN() == other.LogN()
	res = res && (p.logSlots == other.LogSlots())
	res = res && (p.scale == other.Scale())
	res = res && (p.Sigma() == other.Sigma())
	res = res && utils.EqualSliceUint64(p.Q(), other.Q())
	res = res && utils.EqualSliceUint64(p.P(), other.P())
	return
}

// MarshalBinary returns a []byte representation of the parameter set.
func (p *ParametersStruct) MarshalBinary() ([]byte, error) {
	if p.LogN() == 0 { // if N is 0, then p is the zero value
		return []byte{}, nil
	}

	// Data 21 byte + QPiCount * 8 byte:
	// 1 byte : logN
	// 1 byte : logSlots
	// 8 byte : scale
	// 8 byte : sigma
	// 1 byte : #qi
	// 1 byte : #pi
	b := utils.NewBuffer(make([]byte, 0, 20+(p.QPCount())<<3))

	b.WriteUint8(uint8(p.LogN()))
	b.WriteUint8(uint8(p.logSlots))
	b.WriteUint64(math.Float64bits(p.scale))
	b.WriteUint64(math.Float64bits(p.Sigma()))
	b.WriteUint8(uint8(p.QCount()))
	b.WriteUint8(uint8(p.PCount()))
	b.WriteUint64Slice(p.Q())
	b.WriteUint64Slice(p.P())

	return b.Bytes(), nil
}

// UnmarshalBinary decodes a []byte into a parameter set struct
func (p *ParametersStruct) UnmarshalBinary(data []byte) (err error) {

	if len(data) < 20 {
		return errors.New("invalid parameters encoding")
	}

	b := utils.NewBuffer(data)

	logN := int(b.ReadUint8())
	p.logSlots = int(b.ReadUint8())

	if p.logSlots > logN-1 {
		return fmt.Errorf("LogSlots larger than %d", rlwe.MaxLogN-1)
	}

	p.scale = math.Float64frombits(b.ReadUint64())
	sigma := math.Float64frombits(b.ReadUint64())

	lenQi := b.ReadUint8()
	lenPi := b.ReadUint8()

	qi := make([]uint64, lenQi)
	pi := make([]uint64, lenPi)

	b.ReadUint64Slice(qi)
	b.ReadUint64Slice(pi)

	rlweParams, err := rlwe.NewRLWEParameters(logN, qi, pi, sigma)
	if err != nil {
		return err
	}
	p.ParametersStruct = *rlweParams

	return nil
}

// func (p *ParametersStruct) RingQ() *ring.Ring {
// 	ringQ, err := ring.NewRing(p.N(), p.qi)
// 	if err != nil {
// 		panic(err) // Parameter type invariant
// 	}
// 	return ringQ
// }

// func (p *ParametersStruct) RingP() *ring.Ring {
// 	if len(p.pi) == 0 {
// 		return nil
// 	}
// 	ringP, err := ring.NewRing(p.N(), p.pi)
// 	if err != nil {
// 		panic(err) // Parameter type invariant
// 	}
// 	return ringP
// }

// func (p *ParametersStruct) RingQP() *ring.Ring {
// 	ringQP, err := ring.NewRing(p.N(), append(p.qi, p.pi...))
// 	if err != nil {
// 		panic(err) // Parameter type invariant
// 	}
// 	return ringQP
// }
