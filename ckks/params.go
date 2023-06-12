package ckks

import (
	"encoding/json"
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ring"

	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

const (
	NTTFlag = true
)

// Name of the different default parameter sets
var ()

// ParametersLiteral is a literal representation of CKKS parameters.  It has public
// fields and is used to express unchecked user-defined parameters literally into
// Go programs. The NewParametersFromLiteral function is used to generate the actual
// checked parameters from the literal representation.
//
// Users must set the polynomial degree (in log_2, LogN) and the coefficient modulus, by either setting
// the Q and P fields to the desired moduli chain, or by setting the LogQ and LogP fields to
// the desired moduli sizes (in log_2). Users must also specify a default initial scale for the plaintexts.
//
// Optionally, users may specify the error variance (Sigma), the secrets' density (H), the ring
// type (RingType) and the number of slots (in log_2, LogSlots). If left unset, standard default values for
// these field are substituted at parameter creation (see NewParametersFromLiteral).
type ParametersLiteral struct {
	LogN              int
	Q                 []uint64
	P                 []uint64
	LogQ              []int `json:",omitempty"`
	LogP              []int `json:",omitempty"`
	Pow2Base          int
	Xe                ring.DistributionParameters
	Xs                ring.DistributionParameters
	RingType          ring.Type
	LogPlaintextScale int
}

// RLWEParametersLiteral returns the rlwe.ParametersLiteral from the target ckks.ParameterLiteral.
func (p ParametersLiteral) RLWEParametersLiteral() rlwe.ParametersLiteral {
	return rlwe.ParametersLiteral{
		LogN:           p.LogN,
		Q:              p.Q,
		P:              p.P,
		LogQ:           p.LogQ,
		LogP:           p.LogP,
		Pow2Base:       p.Pow2Base,
		Xe:             p.Xe,
		Xs:             p.Xs,
		RingType:       p.RingType,
		NTTFlag:        NTTFlag,
		PlaintextScale: rlwe.NewScale(math.Exp2(float64(p.LogPlaintextScale))),
	}
}

// Parameters represents a parameter set for the CKKS cryptosystem. Its fields are private and
// immutable. See ParametersLiteral for user-specified parameters.
type Parameters struct {
	rlwe.Parameters
}

// NewParameters instantiate a set of CKKS parameters from the generic RLWE parameters and the CKKS-specific ones.
// It returns the empty parameters Parameters{} and a non-nil error if the specified parameters are invalid.
func NewParameters(rlweParams rlwe.Parameters) (p Parameters, err error) {

	if !rlweParams.NTTFlag() {
		return Parameters{}, fmt.Errorf("provided RLWE parameters are invalid for CKKS scheme (NTTFlag must be true)")
	}

	if rlweParams.Equal(rlwe.Parameters{}) {
		return Parameters{}, fmt.Errorf("provided RLWE parameters are invalid")
	}

	return Parameters{rlweParams}, nil
}

// NewParametersFromLiteral instantiate a set of CKKS parameters from a ParametersLiteral specification.
// It returns the empty parameters Parameters{} and a non-nil error if the specified parameters are invalid.
//
// If the `LogSlots` field is left unset, its value is set to `LogN-1` for the Standard ring and to `LogN` for
// the conjugate-invariant ring.
//
// See `rlwe.NewParametersFromLiteral` for default values of the other optional fields.
func NewParametersFromLiteral(pl ParametersLiteral) (Parameters, error) {
	rlweParams, err := rlwe.NewParametersFromLiteral(pl.RLWEParametersLiteral())
	if err != nil {
		return Parameters{}, err
	}

	return NewParameters(rlweParams)
}

// StandardParameters returns the CKKS parameters corresponding to the receiver
// parameter set. If the receiver is already a standard parameter set
// (i.e., RingType==Standard), then the method returns the receiver.
func (p Parameters) StandardParameters() (pckks Parameters, err error) {
	if p.RingType() == ring.Standard {
		return p, nil
	}
	pckks = p
	pckks.Parameters, err = pckks.Parameters.StandardParameters()
	return
}

// ParametersLiteral returns the ParametersLiteral of the target Parameters.
func (p Parameters) ParametersLiteral() (pLit ParametersLiteral) {
	return ParametersLiteral{
		LogN:              p.LogN(),
		Q:                 p.Q(),
		P:                 p.P(),
		Pow2Base:          p.Pow2Base(),
		Xe:                p.Xe(),
		Xs:                p.Xs(),
		RingType:          p.RingType(),
		LogPlaintextScale: p.LogPlaintextScale(),
	}
}

// MaxLevel returns the maximum ciphertext level
func (p Parameters) MaxLevel() int {
	return p.QCount() - 1
}

// PlaintextDimensions returns the maximum dimension of the matrix that can be SIMD packed in a single plaintext polynomial.
func (p Parameters) PlaintextDimensions() [2]int {
	switch p.RingType() {
	case ring.Standard:
		return [2]int{1, p.N() >> 1}
	case ring.ConjugateInvariant:
		return [2]int{1, p.N()}
	default:
		panic("cannot PlaintextDimensions: invalid ring type")
	}
}

// PlaintextLogDimensions returns the log2 of maximum dimension of the matrix that can be SIMD packed in a single plaintext polynomial.
func (p Parameters) PlaintextLogDimensions() [2]int {
	switch p.RingType() {
	case ring.Standard:
		return [2]int{0, p.LogN() - 1}
	case ring.ConjugateInvariant:
		return [2]int{0, p.LogN()}
	default:
		panic("cannot PlaintextLogDimensions: invalid ring type")
	}
}

// PlaintextSlots returns the total number of entries (`slots`) that a plaintext can store.
// This value is obtained by multiplying all dimensions from PlaintextDimensions.
func (p Parameters) PlaintextSlots() int {
	dims := p.PlaintextDimensions()
	return dims[0] * dims[1]
}

// PlaintextLogSlots returns the total number of entries (`slots`) that a plaintext can store.
// This value is obtained by summing all log dimensions from PlaintextLogDimensions.
func (p Parameters) PlaintextLogSlots() int {
	dims := p.PlaintextLogDimensions()
	return dims[0] + dims[1]
}

// LogPlaintextScale returns the log2 of the default plaintext scaling factor.
func (p Parameters) LogPlaintextScale() int {
	return int(math.Round(math.Log2(p.PlaintextScale().Float64())))
}

// LogQLvl returns the size of the modulus Q in bits at a specific level
func (p Parameters) LogQLvl(level int) int {
	tmp := p.QLvl(level)
	return tmp.BitLen()
}

// QLvl returns the product of the moduli at the given level as a big.Int
func (p Parameters) QLvl(level int) *big.Int {
	tmp := bignum.NewInt(1)
	for _, qi := range p.Q()[:level+1] {
		tmp.Mul(tmp, bignum.NewInt(qi))
	}
	return tmp
}

// Equal compares two sets of parameters for equality.
func (p Parameters) Equal(other rlwe.ParametersInterface) bool {
	switch other := other.(type) {
	case Parameters:
		return p.Parameters.Equal(other.Parameters)
	}

	return false
}

// MarshalBinary returns a []byte representation of the parameter set.
// This representation corresponds to the one returned by MarshalJSON.
func (p Parameters) MarshalBinary() ([]byte, error) {
	return p.MarshalJSON()
}

// UnmarshalBinary decodes a []byte into a parameter set struct
func (p *Parameters) UnmarshalBinary(data []byte) (err error) {
	return p.UnmarshalJSON(data)
}

// MarshalJSON returns a JSON representation of this parameter set. See `Marshal` from the `encoding/json` package.
func (p Parameters) MarshalJSON() ([]byte, error) {
	return json.Marshal(p.ParametersLiteral())
}

// UnmarshalJSON reads a JSON representation of a parameter set into the receiver Parameter. See `Unmarshal` from the `encoding/json` package.
func (p *Parameters) UnmarshalJSON(data []byte) (err error) {
	var params ParametersLiteral
	if err = json.Unmarshal(data, &params); err != nil {
		return
	}
	*p, err = NewParametersFromLiteral(params)
	return
}

func (p *ParametersLiteral) UnmarshalJSON(b []byte) (err error) {
	var pl struct {
		LogN              int
		Q                 []uint64
		P                 []uint64
		LogQ              []int
		LogP              []int
		Pow2Base          int
		Xe                map[string]interface{}
		Xs                map[string]interface{}
		RingType          ring.Type
		LogPlaintextScale int
	}

	err = json.Unmarshal(b, &pl)
	if err != nil {
		return err
	}

	p.LogN = pl.LogN
	p.Q, p.P, p.LogQ, p.LogP = pl.Q, pl.P, pl.LogQ, pl.LogP
	p.Pow2Base = pl.Pow2Base
	if pl.Xs != nil {
		p.Xs, err = ring.ParametersFromMap(pl.Xs)
		if err != nil {
			return err
		}
	}
	if pl.Xe != nil {
		p.Xe, err = ring.ParametersFromMap(pl.Xe)
		if err != nil {
			return err
		}
	}
	p.RingType = pl.RingType
	p.LogPlaintextScale = pl.LogPlaintextScale
	return err
}
