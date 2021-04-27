package bfv

import (
	"encoding/binary"
	"encoding/json"
	"fmt"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// DefaultParametersID is an enumerated type representing a set of default paramters
// hardcoded in the library.
type DefaultParametersID int

const (
	// PN12QP109 is a set of parameters with N = 2^12 and log(QP) = 109
	PN12QP109 = DefaultParametersID(iota)
	// PN13QP218 is a set of parameters with N = 2^13 and log(QP) = 218
	PN13QP218
	// PN14QP438 is a set of parameters with N = 2^14 and log(QP) = 438
	PN14QP438
	// PN15QP880 is a set of parameters with N = 2^15 and log(QP) = 880
	PN15QP880

	// PN12QP101pq is the index in DefaultParams for logQP = 101 (post quantum)
	PN12QP101pq
	// PN13QP202pq is the index in DefaultParams for logQP = 202 (post quantum)
	PN13QP202pq
	// PN14QP411pq is the index in DefaultParams for logQP = 411 (post quantum)
	PN14QP411pq
	// PN15QP827pq is the index in DefaultParams for logQP = 827 (post quantum)
	PN15QP827pq
)

// ParametersLiteral is a literal representation of BFV parameters.  It has public
// fields and is used to express unchecked user-defined parameters literally into
// Go programs. The NewParametersFromLiteral function is used to generate the actual
// checked parameters from the literal representation.
type ParametersLiteral struct {
	LogN  uint64 // Log Ring degree (power of 2)
	Q     []uint64
	P     []uint64
	LogQ  []uint64 `json:",omitempty"`
	LogP  []uint64 `json:",omitempty"`
	Sigma float64  // Gaussian sampling standard deviation
	T     uint64   // Plaintext modulus
}

// DefaultParams is a set of default BFV parameters ensuring 128 bit security.
var DefaultParams = []ParametersLiteral{

	{
		LogN:  12,
		T:     65537,
		Q:     []uint64{0x7ffffec001, 0x8000016001}, // 39 + 39 bits
		P:     []uint64{0x40002001},                 // 30 bits
		Sigma: rlwe.DefaultSigma,
	},

	{
		LogN:  13,
		T:     65537,
		Q:     []uint64{0x3fffffffef8001, 0x4000000011c001, 0x40000000120001}, // 54 + 54 + 54 bits
		P:     []uint64{0x7ffffffffb4001},                                     // 55 bits
		Sigma: rlwe.DefaultSigma,
	},

	{
		LogN: 14,
		T:    65537,
		Q: []uint64{0x100000000060001, 0x80000000068001, 0x80000000080001,
			0x3fffffffef8001, 0x40000000120001, 0x3fffffffeb8001}, // 56 + 55 + 55 + 54 + 54 + 54 bits
		P:     []uint64{0x80000000130001, 0x7fffffffe90001}, // 55 + 55 bits
		Sigma: rlwe.DefaultSigma,
	},

	{
		LogN: 15,
		T:    65537,
		Q: []uint64{0x7ffffffffe70001, 0x7ffffffffe10001, 0x7ffffffffcc0001, // 59 + 59 + 59 bits
			0x400000000270001, 0x400000000350001, 0x400000000360001, // 58 + 58 + 58 bits
			0x3ffffffffc10001, 0x3ffffffffbe0001, 0x3ffffffffbd0001, // 58 + 58 + 58 bits
			0x4000000004d0001, 0x400000000570001, 0x400000000660001}, // 58 + 58 + 58 bits
		P:     []uint64{0xffffffffffc0001, 0x10000000001d0001, 0x10000000006e0001}, // 60 + 60 + 60 bits
		Sigma: rlwe.DefaultSigma,
	},

	{ // LogQP = 101.00005709794536
		LogN:  12,
		T:     65537,
		Q:     []uint64{0x800004001, 0x800008001}, // 2*35
		P:     []uint64{0x80014001},               // 1*31
		Sigma: rlwe.DefaultSigma,
	},

	{ // LogQP = 201.99999999994753
		LogN:  13,
		T:     65537,
		Q:     []uint64{0x7fffffffe0001, 0x7fffffffcc001, 0x3ffffffffc001}, // 2*51 + 50
		P:     []uint64{0x4000000024001},                                   // 50
		Sigma: rlwe.DefaultSigma,
	},

	{ // LogQP = 410.9999999999886
		LogN:  14,
		T:     65537,
		Q:     []uint64{0x7fffffffff18001, 0x8000000000f8001, 0x7ffffffffeb8001, 0x800000000158001, 0x7ffffffffe70001}, // 5*59
		P:     []uint64{0x7ffffffffe10001, 0x400000000068001},                                                          // 59+58
		Sigma: rlwe.DefaultSigma,
	},

	{ // LogQP = 826.9999999999509
		LogN: 15,
		T:    65537,
		Q: []uint64{0x7ffffffffe70001, 0x7ffffffffe10001, 0x7ffffffffcc0001, 0x7ffffffffba0001, 0x8000000004a0001,
			0x7ffffffffb00001, 0x800000000890001, 0x8000000009d0001, 0x7ffffffff630001, 0x800000000a70001,
			0x7ffffffff510001}, // 11*59
		P:     []uint64{0x800000000b80001, 0x800000000bb0001, 0xffffffffffc0001}, // 2*59+60
		Sigma: rlwe.DefaultSigma,
	},
}

// GetDefaultParameters instantiate and returns the default parameter set identified by DefaultParametersID
func GetDefaultParameters(paramsID DefaultParametersID) Parameters {
	if int(paramsID) >= len(DefaultParams) {
		panic(fmt.Errorf("paramsId %d does not exist", paramsID))
	}
	params, err := NewParametersFromLiteral(DefaultParams[paramsID])
	if err != nil {
		panic(err)
	}
	return params
}

// Parameters represents a parameter set for the BFV cryptosystem. Its fields are private and
// immutable. See ParametersLiteral for user-specified parameters.
type Parameters struct {
	rlwe.Parameters
	t uint64
}

// NewParameters instantiate a set of BFV parameters from the generic RLWE parameters and the BFV-specific ones.
// It returns the empty parameters Parameters{} and a non-nil error if the specified parameters are invalid.
func NewParameters(rlweParams rlwe.Parameters, t uint64) (p Parameters, err error) {
	if rlweParams.Equals(rlwe.Parameters{}) {
		return Parameters{}, fmt.Errorf("provided RLWE parameters are invalid")
	}
	if t > rlweParams.Q()[0] {
		return Parameters{}, fmt.Errorf("t=%d is larger than Q[0]=%d", t, rlweParams.Q()[0])
	}
	return Parameters{rlweParams, t}, nil
}

// NewParametersFromLiteral instantiate a set of BFV parameters from a ParametersLiteral specification.
// It returns the empty parameters Parameters{} and a non-nil error if the specified parameters are invalid.
func NewParametersFromLiteral(pl ParametersLiteral) (Parameters, error) {
	rlweParams, err := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{LogN: pl.LogN, Q: pl.Q, P: pl.P, LogQ: pl.LogQ, LogP: pl.LogP, Sigma: pl.Sigma})
	if err != nil {
		return Parameters{}, err
	}
	return NewParameters(rlweParams, pl.T)
}

// T returns the plaintext coefficient modulus t
func (p Parameters) T() uint64 {
	return p.t
}

// RingT instantiates a new ring.Ring corresponding to the plaintext space ring R_t.
func (p Parameters) RingT() *ring.Ring {
	ringQP, err := ring.NewRing(p.N(), []uint64{p.t})
	if err != nil {
		panic(err) // Parameter type invariant
	}
	return ringQP
}

// Equals compares two sets of parameters for equality.
func (p Parameters) Equals(other Parameters) bool {
	res := p.Parameters.Equals(other.Parameters)
	res = res && (p.t == other.T())
	return res
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
	// 8 byte : T
	var tBytes [8]byte
	binary.BigEndian.PutUint64(tBytes[:], p.t)
	data := append(rlweBytes, tBytes[:]...)
	return data, nil
}

// UnmarshalBinary decodes a []byte into a parameter set struct.
func (p *Parameters) UnmarshalBinary(data []byte) error {
	if err := p.Parameters.UnmarshalBinary(data); err != nil {
		return err
	}
	dataBfv := data[len(data)-8:]
	p.t = binary.BigEndian.Uint64(dataBfv)
	return nil
}

// MarshalJSON returns a JSON representation of this parameter set. See `Marshal` from the `encoding/json` package.
func (p Parameters) MarshalJSON() ([]byte, error) {
	return json.Marshal(ParametersLiteral{LogN: p.LogN(), Q: p.Q(), P: p.P(), Sigma: p.Sigma(), T: p.t})
}

// UnmarshalJSON reads a JSON representation of a parameter set into the receiver Parameter. See `Unmarshal` from the `encoding/json` package.
func (p *Parameters) UnmarshalJSON(data []byte) (err error) {
	var params ParametersLiteral
	json.Unmarshal(data, &params)
	*p, err = NewParametersFromLiteral(params)
	return
}
