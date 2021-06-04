package ckks

import (
	"encoding/binary"
	"encoding/json"
	"fmt"
	"math"
	"math/big"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// Name of the different default parameter sets
var (
	// PN12QP109 is a default parameter set for logN=12 and logQP=109
	PN12QP109 = ParametersLiteral{
		LogN:     12,
		LogSlots: 11,
		Q: []uint64{0x200000e001, // 37 + 32
			0x100006001},
		P:     []uint64{0x3ffffea001}, // 38
		Scale: 1 << 32,
		Sigma: rlwe.DefaultSigma,
	}

	// PN13QP218 is a default parameter set for logN=13 and logQP=218
	PN13QP218 = ParametersLiteral{
		LogN:     13,
		LogSlots: 12,
		Q: []uint64{0x1fffec001, // 33 + 5 x 30
			0x3fff4001,
			0x3ffe8001,
			0x40020001,
			0x40038001,
			0x3ffc0001},
		P:     []uint64{0x800004001}, // 35
		Scale: 1 << 30,
		Sigma: rlwe.DefaultSigma,
	}
	// PN14QP438 is a default parameter set for logN=14 and logQP=438
	PN14QP438 = ParametersLiteral{
		LogN:     14,
		LogSlots: 13,
		Q: []uint64{0x200000008001, 0x400018001, // 45 + 9 x 34
			0x3fffd0001, 0x400060001,
			0x400068001, 0x3fff90001,
			0x400080001, 0x4000a8001,
			0x400108001, 0x3ffeb8001},
		P:     []uint64{0x7fffffd8001, 0x7fffffc8001}, // 43, 43
		Scale: 1 << 34,
		Sigma: rlwe.DefaultSigma,
	}

	// PN15QP880 is a default parameter set for logN=15 and logQP=880
	PN15QP880 = ParametersLiteral{
		LogN:     15,
		LogSlots: 14,
		Q: []uint64{0x4000000120001, 0x10000140001, 0xffffe80001, // 50 + 17 x 40
			0x10000290001, 0xffffc40001, 0x100003e0001,
			0x10000470001, 0x100004b0001, 0xffffb20001,
			0x10000500001, 0x10000650001, 0xffff940001,
			0xffff8a0001, 0xffff820001, 0xffff780001,
			0x10000890001, 0xffff750001, 0x10000960001},
		P:     []uint64{0x40000001b0001, 0x3ffffffdf0001, 0x4000000270001}, // 50, 50, 50
		Scale: 1 << 40,
		Sigma: rlwe.DefaultSigma,
	}
	// PN16QP1761 is a default parameter set for logN=16 and logQP = 1761
	PN16QP1761 = ParametersLiteral{
		LogN:     16,
		LogSlots: 15,
		Q: []uint64{0x80000000080001, 0x2000000a0001, 0x2000000e0001, 0x1fffffc20001, // 55 + 33 x 45
			0x200000440001, 0x200000500001, 0x200000620001, 0x1fffff980001,
			0x2000006a0001, 0x1fffff7e0001, 0x200000860001, 0x200000a60001,
			0x200000aa0001, 0x200000b20001, 0x200000c80001, 0x1fffff360001,
			0x200000e20001, 0x1fffff060001, 0x200000fe0001, 0x1ffffede0001,
			0x1ffffeca0001, 0x1ffffeb40001, 0x200001520001, 0x1ffffe760001,
			0x2000019a0001, 0x1ffffe640001, 0x200001a00001, 0x1ffffe520001,
			0x200001e80001, 0x1ffffe0c0001, 0x1ffffdee0001, 0x200002480001,
			0x1ffffdb60001, 0x200002560001},
		P:     []uint64{0x80000000440001, 0x7fffffffba0001, 0x80000000500001, 0x7fffffffaa0001}, // 4 x 55
		Scale: 1 << 45,
		Sigma: rlwe.DefaultSigma,
	}

	// PN12QP101pq is a default (post quantum) parameter set for logN=12 and logQP=101
	PN12QP101pq = ParametersLiteral{
		LogN:     12,
		LogSlots: 11,
		Q:        []uint64{0x800004001, 0x40002001}, // 35 + 30
		P:        []uint64{0x1000002001},            // 36
		Scale:    1 << 30,
		Sigma:    rlwe.DefaultSigma,
	}
	// PN13QP202pq is a default (post quantum) parameter set for logN=13 and logQP=202
	PN13QP202pq = ParametersLiteral{
		LogN:     13,
		LogSlots: 12,
		Q:        []uint64{0x1fffec001, 0x8008001, 0x8020001, 0x802c001, 0x7fa8001, 0x7f74001}, // 33 + 5 x 27
		P:        []uint64{0x400018001},                                                        // 34
		Scale:    1 << 27,
		Sigma:    rlwe.DefaultSigma,
	}

	// PN14QP411pq is a default (post quantum) parameter set for logN=14 and logQP=411
	PN14QP411pq = ParametersLiteral{
		LogN:     14,
		LogSlots: 13,
		Q: []uint64{0x10000048001, 0x200038001, 0x1fff90001, 0x200080001, 0x1fff60001,
			0x2000b8001, 0x200100001, 0x1fff00001, 0x1ffef0001, 0x200128001}, // 40 + 9 x 33

		P:     []uint64{0x1ffffe0001, 0x1ffffc0001}, // 37, 37
		Scale: 1 << 33,
		Sigma: rlwe.DefaultSigma,
	}

	// PN15QP827pq is a default (post quantum) parameter set for logN=15 and logQP=827
	PN15QP827pq = ParametersLiteral{
		LogN:     15,
		LogSlots: 14,
		Q: []uint64{0x400000060001, 0x4000170001, 0x3fffe80001, 0x40002f0001, 0x4000300001,
			0x3fffcf0001, 0x40003f0001, 0x3fffc10001, 0x4000450001, 0x3fffb80001,
			0x3fffb70001, 0x40004a0001, 0x3fffb20001, 0x4000510001, 0x3fffaf0001,
			0x4000540001, 0x4000560001, 0x4000590001}, // 46 + 17 x 38
		P:     []uint64{0x2000000a0001, 0x2000000e0001, 0x2000001d0001}, // 3 x 45
		Scale: 1 << 38,
		Sigma: rlwe.DefaultSigma,
	}
	// PN16QP1654pq is a default (post quantum) parameter set for logN=16 and logQP=1654
	PN16QP1654pq = ParametersLiteral{LogN: 16,
		LogSlots: 15,
		Q: []uint64{0x80000000080001, 0x2000000a0001, 0x2000000e0001, 0x1fffffc20001, 0x200000440001,
			0x200000500001, 0x200000620001, 0x1fffff980001, 0x2000006a0001, 0x1fffff7e0001,
			0x200000860001, 0x200000a60001, 0x200000aa0001, 0x200000b20001, 0x200000c80001,
			0x1fffff360001, 0x200000e20001, 0x1fffff060001, 0x200000fe0001, 0x1ffffede0001,
			0x1ffffeca0001, 0x1ffffeb40001, 0x200001520001, 0x1ffffe760001, 0x2000019a0001,
			0x1ffffe640001, 0x200001a00001, 0x1ffffe520001, 0x200001e80001, 0x1ffffe0c0001,
			0x1ffffdee0001, 0x200002480001}, // 55 + 31 x 45
		P:     []uint64{0x7fffffffe0001, 0x80000001c0001, 0x80000002c0001, 0x7ffffffd20001}, // 4 x 51
		Scale: 1 << 45,
		Sigma: rlwe.DefaultSigma,
	}
)

// ParametersLiteral is a literal representation of BFV parameters.  It has public
// fields and is used to express unchecked user-defined parameters literally into
// Go programs. The NewParametersFromLiteral function is used to generate the actual
// checked parameters from the literal representation.
type ParametersLiteral struct {
	LogN     int // Ring degree (power of 2)
	Q        []uint64
	P        []uint64
	LogQ     []int   `json:",omitempty"`
	LogP     []int   `json:",omitempty"`
	Sigma    float64 // Gaussian sampling variance
	LogSlots int
	Scale    float64
}

// DefaultParams is a set of default CKKS parameters ensuring 128 bit security in a classic setting.
var DefaultParams = []ParametersLiteral{PN12QP109, PN13QP218, PN14QP438, PN15QP880, PN16QP1761}

// DefaultPostQuantumParams is a set of default CKKS parameters ensuring 128 bit security in a post-quantum setting.
var DefaultPostQuantumParams = []ParametersLiteral{PN12QP101pq, PN13QP202pq, PN14QP411pq, PN15QP827pq, PN16QP1654pq}

// Parameters represents a parameter set for the CKKS cryptosystem. Its fields are private and
// immutable. See ParametersLiteral for user-specified parameters.
type Parameters struct {
	rlwe.Parameters

	logSlots int
	scale    float64
}

// NewParameters instantiate a set of CKKS parameters from the generic RLWE parameters and the CKKS-specific ones.
// It returns the empty parameters Parameters{} and a non-nil error if the specified parameters are invalid.
func NewParameters(rlweParams rlwe.Parameters, logSlot int, scale float64) (p Parameters, err error) {
	if rlweParams.Equals(rlwe.Parameters{}) {
		return Parameters{}, fmt.Errorf("provided RLWE parameters are invalid")
	}
	if logSlot > rlweParams.LogN()-1 {
		return Parameters{}, fmt.Errorf("logSlot=%d is larger than the logN-1=%d", logSlot, rlweParams.LogN()-1)
	}
	return Parameters{rlweParams, logSlot, scale}, nil
}

// NewParametersFromLiteral instantiate a set of CKKS parameters from a ParametersLiteral specification.
// It returns the empty parameters Parameters{} and a non-nil error if the specified parameters are invalid.
func NewParametersFromLiteral(pl ParametersLiteral) (Parameters, error) {
	rlweParams, err := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{LogN: pl.LogN, Q: pl.Q, P: pl.P, LogQ: pl.LogQ, LogP: pl.LogP, Sigma: pl.Sigma})
	if err != nil {
		return Parameters{}, err
	}
	return NewParameters(rlweParams, pl.LogSlots, pl.Scale)
}

// LogSlots returns the log of the number of slots
func (p Parameters) LogSlots() int {
	return p.logSlots
}

// MaxLevel returns the maximum ciphertext level
func (p Parameters) MaxLevel() int {
	return p.QCount() - 1
}

// Slots returns number of available plaintext slots
func (p Parameters) Slots() int {
	return 1 << p.logSlots
}

// MaxSlots returns the theoretical maximum of plaintext slots allowed by the ring degree
func (p Parameters) MaxSlots() int {
	return p.N() >> 1
}

// MaxLogSlots returns the log of the maximum number of slots enabled by the parameters
func (p Parameters) MaxLogSlots() int {
	return p.LogN() - 1
}

// Scale returns the default plaintext/ciphertext scale
func (p Parameters) Scale() float64 {
	return p.scale
}

// LogQLvl returns the size of the modulus Q in bits at a specific level
func (p Parameters) LogQLvl(level int) int {
	tmp := p.QLvl(level)
	return tmp.BitLen()
}

// QLvl returns the product of the moduli at the given level as a big.Int
func (p Parameters) QLvl(level int) *big.Int {
	tmp := ring.NewUint(1)
	for _, qi := range p.Q()[:level+1] {
		tmp.Mul(tmp, ring.NewUint(qi))
	}
	return tmp
}

// RotationsForInnerSum generates the rotations that will be performed by the
// `Evaluator.InnerSum` operation when performed with parameters `batch` and `n`.
func (p Parameters) RotationsForInnerSum(batch, n int) (rotations []int) {
	rotations = []int{}
	for i := 1; i < n; i++ {
		rotations = append(rotations, i*batch)
	}
	return
}

// RotationsForInnerSumLog generates the rotations that will be performed by the
// `Evaluator.InnerSumLog` operation when performed with parameters `batch` and `n`.
func (p Parameters) RotationsForInnerSumLog(batch, n int) (rotations []int) {

	rotations = []int{}
	var k int
	for i := 1; i < n; i <<= 1 {

		k = i
		k *= batch

		if !utils.IsInSliceInt(k, rotations) && k != 0 {
			rotations = append(rotations, k)
		}

		k = n - (n & ((i << 1) - 1))
		k *= batch

		if !utils.IsInSliceInt(k, rotations) && k != 0 {
			rotations = append(rotations, k)
		}
	}

	return
}

// RotationsForReplicate generates the rotations that will be performed by the
// `Evaluator.Replicate` operation when performed with parameters `batch` and `n`.
func (p Parameters) RotationsForReplicate(batch, n int) (rotations []int) {
	return p.RotationsForInnerSum(-batch, n)
}

// RotationsForReplicateLog generates the rotations that will be performed by the
// `Evaluator.ReplicateLog` operation when performed with parameters `batch` and `n`.
func (p Parameters) RotationsForReplicateLog(batch, n int) (rotations []int) {
	return p.RotationsForInnerSumLog(-batch, n)
}

// RotationsForSubSum generates the rotations that will be performed by the
// `Evaluator.SubSum` operation.
func (p Parameters) RotationsForSubSum(logSlots int) (rotations []int) {
	rotations = []int{}

	logN := p.LogN()

	//SubSum rotation needed X -> Y^slots rotations
	for i := logSlots; i < logN-1; i++ {
		if !utils.IsInSliceInt(1<<i, rotations) {
			rotations = append(rotations, 1<<i)
		}
	}

	return
}

// RotationsForDiagMatrixMult generates of all the rotations needed for a the multiplication
// with the provided diagonal plaintext matrix.
func (p Parameters) RotationsForDiagMatrixMult(matrix *PtDiagMatrix) []int {
	slots := 1 << matrix.LogSlots

	rotKeyIndex := []int{}

	var index int

	N1 := matrix.N1

	if len(matrix.Vec) < 3 {

		for j := range matrix.Vec {

			if !utils.IsInSliceInt(j, rotKeyIndex) {
				rotKeyIndex = append(rotKeyIndex, j)
			}
		}

	} else {

		for j := range matrix.Vec {

			index = ((j / N1) * N1) & (slots - 1)

			if index != 0 && !utils.IsInSliceInt(index, rotKeyIndex) {
				rotKeyIndex = append(rotKeyIndex, index)
			}

			index = j & (N1 - 1)

			if index != 0 && !utils.IsInSliceInt(index, rotKeyIndex) {
				rotKeyIndex = append(rotKeyIndex, index)
			}
		}
	}

	return rotKeyIndex
}

// Equals compares two sets of parameters for equality.
func (p Parameters) Equals(other Parameters) bool {
	res := p.Parameters.Equals(other.Parameters)
	res = res && (p.logSlots == other.LogSlots())
	res = res && (p.scale == other.Scale())
	return res
}

// CopyNew makes a deep copy of the receiver and returns it.
func (p Parameters) CopyNew() Parameters {
	p.Parameters = p.Parameters.CopyNew()
	return p
}

// MarshalBinary returns a []byte representation of the parameter set.
func (p Parameters) MarshalBinary() ([]byte, error) {
	if p.LogN() == 0 { // if N is 0, then p is the zero value
		return []byte{}, nil
	}

	rlweBytes, err := p.Parameters.MarshalBinary()
	if err != nil {
		return nil, err
	}

	// len(rlweBytes) : RLWE parameters
	// 1 byte : logSlots
	// 8 byte : scale
	b := utils.NewBuffer(make([]byte, 0, p.MarshalBinarySize()))
	b.WriteUint8Slice(rlweBytes)
	b.WriteUint8(uint8(p.logSlots))
	b.WriteUint64(math.Float64bits(p.scale))
	return b.Bytes(), nil
}

// UnmarshalBinary decodes a []byte into a parameter set struct
func (p *Parameters) UnmarshalBinary(data []byte) (err error) {
	var rlweParams rlwe.Parameters
	if err := rlweParams.UnmarshalBinary(data); err != nil {
		return err
	}
	logSlots := int(data[len(data)-9])
	scale := math.Float64frombits(binary.BigEndian.Uint64(data[len(data)-8:]))
	*p, err = NewParameters(rlweParams, logSlots, scale)
	return err
}

// MarshalBinarySize returns the length of the []byte encoding of the reciever.
func (p Parameters) MarshalBinarySize() int {
	return p.Parameters.MarshalBinarySize() + 9
}

// MarshalJSON returns a JSON representation of this parameter set. See `Marshal` from the `encoding/json` package.
func (p Parameters) MarshalJSON() ([]byte, error) {
	return json.Marshal(ParametersLiteral{LogN: p.LogN(), Q: p.Q(), P: p.P(), Sigma: p.Sigma(), LogSlots: p.logSlots, Scale: p.scale})
}

// UnmarshalJSON reads a JSON representation of a parameter set into the receiver Parameter. See `Unmarshal` from the `encoding/json` package.
func (p *Parameters) UnmarshalJSON(data []byte) (err error) {
	var params ParametersLiteral
	json.Unmarshal(data, &params)
	*p, err = NewParametersFromLiteral(params)
	return
}
