package bgv

import (
	"encoding/json"
	"fmt"
	"math"
	"math/bits"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/ring/distribution"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
)

const (
	NTTFlag = true
)

// ParametersLiteral is a literal representation of BGV parameters.  It has public
// fields and is used to express unchecked user-defined parameters literally into
// Go programs. The NewParametersFromLiteral function is used to generate the actual
// checked parameters from the literal representation.
//
// Users must set the polynomial degree (LogN) and the coefficient modulus, by either setting
// the Q and P fields to the desired moduli chain, or by setting the LogQ and LogP fields to
// the desired moduli sizes. Users must also specify the coefficient modulus in plaintext-space
// (T).
//
// Optionally, users may specify the error variance (Sigma) and secrets' density (H). If left
// unset, standard default values for these field are substituted at parameter creation (see
// NewParametersFromLiteral).
type ParametersLiteral struct {
	LogN     int
	Q        []uint64
	P        []uint64
	LogQ     []int `json:",omitempty"`
	LogP     []int `json:",omitempty"`
	Pow2Base int
	Xe       distribution.Distribution
	Xs       distribution.Distribution
	RingType ring.Type
	T        uint64 // Plaintext modulus
}

// RLWEParametersLiteral returns the rlwe.ParametersLiteral from the target bgv.ParametersLiteral.
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
		RingType:       ring.Standard,
		PlaintextScale: rlwe.NewScaleModT(1, p.T),
		NTTFlag:        NTTFlag,
	}
}

// Parameters represents a parameter set for the BGV cryptosystem. Its fields are private and
// immutable. See ParametersLiteral for user-specified parameters.
type Parameters struct {
	rlwe.Parameters
	ringQMul *ring.Ring
	ringT    *ring.Ring
}

// NewParameters instantiate a set of BGV parameters from the generic RLWE parameters and the BGV-specific ones.
// It returns the empty parameters Parameters{} and a non-nil error if the specified parameters are invalid.
func NewParameters(rlweParams rlwe.Parameters, t uint64) (p Parameters, err error) {

	if !rlweParams.NTTFlag() {
		return Parameters{}, fmt.Errorf("provided RLWE parameters are invalid for BGV scheme (NTTFlag must be true)")
	}

	if t == 0 {
		return Parameters{}, fmt.Errorf("invalid parameters: t = 0")
	}

	if utils.IsInSlice(t, rlweParams.Q()) {
		return Parameters{}, fmt.Errorf("insecure parameters: t|Q")
	}

	if rlweParams.Equal(rlwe.Parameters{}) {
		return Parameters{}, fmt.Errorf("provided RLWE parameters are invalid")
	}

	if t > rlweParams.Q()[0] {
		return Parameters{}, fmt.Errorf("t=%d is larger than Q[0]=%d", t, rlweParams.Q()[0])
	}

	var ringQMul *ring.Ring
	nbQiMul := int(math.Ceil(float64(rlweParams.RingQ().ModulusAtLevel[rlweParams.MaxLevel()].BitLen()+rlweParams.LogN()) / 61.0))
	if ringQMul, err = ring.NewRing(rlweParams.N(), ring.GenerateNTTPrimesP(61, 2*rlweParams.N(), nbQiMul)); err != nil {
		return Parameters{}, err
	}

	// Find the largest cyclotomic order enabled by T
	order := uint64(1 << bits.Len64(t))
	for t&(order-1) != 1 {
		order >>= 1
	}

	if order < 16 {
		return Parameters{}, fmt.Errorf("provided plaintext modulus t has cyclotomic order < 16 (ring degree of minimum 8 is required by the backend)")
	}

	var ringT *ring.Ring
	if ringT, err = ring.NewRing(utils.Min(rlweParams.N(), int(order>>1)), []uint64{t}); err != nil {
		return Parameters{}, fmt.Errorf("provided plaintext modulus t is invalid: %w", err)
	}

	return Parameters{
		Parameters: rlweParams,
		ringQMul:   ringQMul,
		ringT:      ringT,
	}, nil
}

// NewParametersFromLiteral instantiate a set of BGV parameters from a ParametersLiteral specification.
// It returns the empty parameters Parameters{} and a non-nil error if the specified parameters are invalid.
//
// See `rlwe.NewParametersFromLiteral` for default values of the optional fields.
func NewParametersFromLiteral(pl ParametersLiteral) (Parameters, error) {
	rlweParams, err := rlwe.NewParametersFromLiteral(pl.RLWEParametersLiteral())
	if err != nil {
		return Parameters{}, err
	}
	return NewParameters(rlweParams, pl.T)
}

// ParametersLiteral returns the ParametersLiteral of the target Parameters.
func (p Parameters) ParametersLiteral() ParametersLiteral {
	return ParametersLiteral{
		LogN:     p.LogN(),
		Q:        p.Q(),
		P:        p.P(),
		Pow2Base: p.Pow2Base(),
		Xe:       p.Xe(),
		Xs:       p.Xs(),
		T:        p.T(),
		RingType: p.RingType(),
	}
}

// PlaintextDimensions returns the maximum dimension of the matrix that can be SIMD packed in a single plaintext polynomial.
func (p Parameters) PlaintextDimensions() [2]int {
	switch p.RingType() {
	case ring.Standard:
		return [2]int{2, p.RingT().N() >> 1}
	case ring.ConjugateInvariant:
		return [2]int{1, p.RingT().N()}
	default:
		panic("cannot PlaintextDimensions: invalid ring type")
	}
}

// PlaintextLogDimensions returns the log2 of maximum dimension of the matrix that can be SIMD packed in a single plaintext polynomial.
func (p Parameters) PlaintextLogDimensions() [2]int {
	switch p.RingType() {
	case ring.Standard:
		return [2]int{1, p.RingT().LogN() - 1}
	case ring.ConjugateInvariant:
		return [2]int{0, p.RingT().LogN()}
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

// RingQMul returns a pointer to the ring of the extended basis for multiplication.
func (p Parameters) RingQMul() *ring.Ring {
	return p.ringQMul
}

// T returns the plaintext coefficient modulus t.
func (p Parameters) T() uint64 {
	return p.ringT.SubRings[0].Modulus
}

// LogT returns log2(plaintext coefficient modulus).
func (p Parameters) LogT() float64 {
	return math.Log2(float64(p.T()))
}

// RingT returns a pointer to the plaintext ring.
func (p Parameters) RingT() *ring.Ring {
	return p.ringT
}

// Equal compares two sets of parameters for equality.
func (p Parameters) Equal(other rlwe.ParametersInterface) bool {
	switch other := other.(type) {
	case Parameters:
		return p.Parameters.Equal(other.Parameters) && (p.T() == other.T())
	}

	panic(fmt.Errorf("cannot Equal: type do not match: %T != %T", p, other))
}

// MarshalBinary returns a []byte representation of the parameter set.
func (p Parameters) MarshalBinary() ([]byte, error) {
	return p.MarshalJSON()
}

// UnmarshalBinary decodes a []byte into a parameter set struct.
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
