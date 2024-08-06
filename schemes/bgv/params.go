package bgv

import (
	"encoding/json"
	"fmt"
	"math"
	"math/bits"
	"slices"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils"
)

const (
	NTTFlag = true
)

// ParametersLiteral is a literal representation of BGV parameters.  It has public
// fields and is used to express unchecked user-defined parameters literally into
// Go programs. The [NewParametersFromLiteral] function is used to generate the actual
// checked parameters from the literal representation.
//
// Users must set the polynomial degree (LogN) and the coefficient modulus, by either setting
// the Q and P fields to the desired moduli chain, or by setting the LogQ and LogP fields to
// the desired moduli sizes.
//
// Users must also specify the coefficient modulus in plaintext-space (T). This modulus must
// be an NTT-friendly prime in the plaintext space: it must be equal to 1 modulo 2n where
// n is the plaintext ring degree (i.e., the plaintext space has n slots).
//
// Optionally, users may specify the error variance (Sigma) and secrets' density (H). If left
// unset, standard default values for these field are substituted at parameter creation (see
// NewParametersFromLiteral).
type ParametersLiteral struct {
	LogN             int
	LogNthRoot       int
	Q                []uint64
	P                []uint64
	LogQ             []int `json:",omitempty"`
	LogP             []int `json:",omitempty"`
	Xe               ring.DistributionParameters
	Xs               ring.DistributionParameters
	PlaintextModulus uint64 // Plaintext modulus
}

// GetRLWEParametersLiteral returns the [rlwe.ParametersLiteral] from the target [bgv.ParametersLiteral].
// See the [ParametersLiteral] type for details on the BGV parameters.
func (p ParametersLiteral) GetRLWEParametersLiteral() rlwe.ParametersLiteral {
	return rlwe.ParametersLiteral{
		LogN:         p.LogN,
		LogNthRoot:   p.LogNthRoot,
		Q:            p.Q,
		P:            p.P,
		LogQ:         p.LogQ,
		LogP:         p.LogP,
		Xe:           p.Xe,
		Xs:           p.Xs,
		RingType:     ring.Standard,
		DefaultScale: rlwe.NewScaleModT(1, p.PlaintextModulus),
		NTTFlag:      NTTFlag,
	}
}

// Parameters represents a parameter set for the BGV cryptosystem. Its fields are private and
// immutable. See [ParametersLiteral] for user-specified parameters.
type Parameters struct {
	rlwe.Parameters
	ringQMul *ring.Ring
	ringT    *ring.Ring
}

// NewParameters instantiate a set of BGV parameters from the generic RLWE parameters and the BGV-specific ones.
// It returns the empty parameters [Parameters]{} and a non-nil error if the specified parameters are invalid.
// See the [ParametersLiteral] type for more details on the BGV parameters.
func NewParameters(rlweParams rlwe.Parameters, t uint64) (p Parameters, err error) {

	if !rlweParams.NTTFlag() {
		return Parameters{}, fmt.Errorf("provided RLWE parameters are invalid for BGV scheme (NTTFlag must be true)")
	}

	if t == 0 {
		return Parameters{}, fmt.Errorf("invalid parameters: t = 0")
	}

	if slices.Contains(rlweParams.Q(), t) {
		return Parameters{}, fmt.Errorf("insecure parameters: t|Q")
	}

	if rlweParams.Equal(&rlwe.Parameters{}) {
		return Parameters{}, fmt.Errorf("provided RLWE parameters are invalid")
	}

	if t > rlweParams.Q()[0] {
		return Parameters{}, fmt.Errorf("t=%d is larger than Q[0]=%d", t, rlweParams.Q()[0])
	}

	var ringQMul *ring.Ring
	nbQiMul := int(math.Ceil(float64(rlweParams.RingQ().ModulusAtLevel[rlweParams.MaxLevel()].BitLen()+rlweParams.LogN()) / 61.0))
	g := ring.NewNTTFriendlyPrimesGenerator(61, uint64(rlweParams.NthRoot()))
	primes, err := g.NextDownstreamPrimes(nbQiMul)
	if err != nil {
		return Parameters{}, err
	}
	if ringQMul, err = ring.NewRing(rlweParams.N(), primes); err != nil {
		return Parameters{}, err
	}

	// Find the largest cyclotomic order enabled by T
	var order uint64
	for order = uint64(1 << bits.Len64(t)); t&(order-1) != 1 && order != 0; order >>= 1 {
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

// NewParametersFromLiteral instantiate a set of BGV parameters from a [ParametersLiteral] specification.
// It returns the empty parameters [Parameters]{} and a non-nil error if the specified parameters are invalid.
//
// See [rlwe.NewParametersFromLiteral] for default values of the optional fields and other details on the BGV
// parameters.
func NewParametersFromLiteral(pl ParametersLiteral) (Parameters, error) {
	rlweParams, err := rlwe.NewParametersFromLiteral(pl.GetRLWEParametersLiteral())
	if err != nil {
		return Parameters{}, err
	}
	return NewParameters(rlweParams, pl.PlaintextModulus)
}

// ParametersLiteral returns the [ParametersLiteral] of the target Parameters.
func (p Parameters) ParametersLiteral() ParametersLiteral {
	return ParametersLiteral{
		LogN:             p.LogN(),
		LogNthRoot:       p.LogNthRoot(),
		Q:                p.Q(),
		P:                p.P(),
		Xe:               p.Xe(),
		Xs:               p.Xs(),
		PlaintextModulus: p.PlaintextModulus(),
	}
}

// GetRLWEParameters returns a pointer to the underlying RLWE parameters.
func (p Parameters) GetRLWEParameters() *rlwe.Parameters {
	return &p.Parameters
}

// MaxDimensions returns the maximum dimension of the matrix that can be SIMD packed in a single plaintext polynomial.
func (p Parameters) MaxDimensions() ring.Dimensions {
	switch p.RingType() {
	case ring.Standard:
		return ring.Dimensions{Rows: 2, Cols: p.RingT().N() >> 1}
	case ring.ConjugateInvariant:
		return ring.Dimensions{Rows: 1, Cols: p.RingT().N()}
	default:
		panic("cannot MaxDimensions: invalid ring type")
	}
}

// LogMaxDimensions returns the log2 of maximum dimension of the matrix that can be SIMD packed in a single plaintext polynomial.
func (p Parameters) LogMaxDimensions() ring.Dimensions {
	switch p.RingType() {
	case ring.Standard:
		return ring.Dimensions{Rows: 1, Cols: p.RingT().LogN() - 1}
	case ring.ConjugateInvariant:
		return ring.Dimensions{Rows: 0, Cols: p.RingT().LogN()}
	default:
		panic("cannot LogMaxDimensions: invalid ring type")
	}
}

// MaxSlots returns the total number of entries (slots) that a plaintext can store.
// This value is obtained by multiplying all dimensions from MaxDimensions.
func (p Parameters) MaxSlots() int {
	dims := p.MaxDimensions()
	return dims.Rows * dims.Cols
}

// LogMaxSlots returns the total number of entries (slots) that a plaintext can store.
// This value is obtained by summing all log dimensions from LogDimensions.
func (p Parameters) LogMaxSlots() int {
	dims := p.LogMaxDimensions()
	return dims.Rows + dims.Cols
}

// RingQMul returns a pointer to the ring of the extended basis for multiplication.
func (p Parameters) RingQMul() *ring.Ring {
	return p.ringQMul
}

// PlaintextModulus returns the plaintext coefficient modulus t.
func (p Parameters) PlaintextModulus() uint64 {
	return p.ringT.SubRings[0].Modulus
}

// LogT returns log2(plaintext coefficient modulus).
func (p Parameters) LogT() float64 {
	return math.Log2(float64(p.PlaintextModulus()))
}

// RingT returns a pointer to the plaintext ring.
func (p Parameters) RingT() *ring.Ring {
	return p.ringT
}

// GaloisElementForColRotation returns the Galois element for generating the
// automorphism phi(k): X -> X^{5^k mod 2N} mod (X^{N} + 1), which acts as a
// column-wise cyclic rotation by k position to the left on batched plaintexts.
//
// Example:
// Recall that batched plaintexts are 2xN/2 matrices, thus given the following
// plaintext matrix:
//
// [a, b, c, d][e, f, g, h]
//
// a rotation by k=3 will change the plaintext to:
//
// [d, a, b, d][h, e, f, g]
//
// Providing a negative k will change direction of the cyclic rotation do the right.
func (p Parameters) GaloisElementForColRotation(k int) uint64 {
	return p.Parameters.GaloisElement(k)
}

// GaloisElementForRowRotation returns the Galois element for generating the
// automorphism X -> X^{-1 mod NthRoot} mod (X^{N} + 1). This automorphism
// acts as a swapping the rows of the plaintext algebra when the plaintext
// is batched.
//
// Example:
// Recall that batched plaintexts are 2xN/2 matrices, thus given the following
// plaintext matrix:
//
// [a, b, c, d][e, f, g, h]
//
// a row rotation will change the plaintext to:
//
// [e, f, g, h][a, b, c, d]
func (p Parameters) GaloisElementForRowRotation() uint64 {
	return p.Parameters.GaloisElementOrderTwoOrthogonalSubgroup()
}

// GaloisElementsForInnerSum returns the list of Galois elements necessary to apply the method
// InnerSum operation with parameters batch and n.
func (p Parameters) GaloisElementsForInnerSum(batch, n int) (galEls []uint64) {
	galEls = rlwe.GaloisElementsForInnerSum(p, batch, n)
	if n > p.N()>>1 {
		galEls = append(galEls, p.GaloisElementForRowRotation())
	}
	return
}

// GaloisElementsForReplicate returns the list of Galois elements necessary to perform the
// Replicate operation with parameters batch and n.
func (p Parameters) GaloisElementsForReplicate(batch, n int) (galEls []uint64) {
	galEls = rlwe.GaloisElementsForReplicate(p, batch, n)
	if n > p.N()>>1 {
		galEls = append(galEls, p.GaloisElementForRowRotation())
	}
	return
}

// GaloisElementsForTrace returns the list of Galois elements required for the for the Trace operation.
// Trace maps X -> sum((-1)^i * X^{i*n+1}) for 2^{LogN} <= i < N.
func (p Parameters) GaloisElementsForTrace(logN int) []uint64 {
	return rlwe.GaloisElementsForTrace(p, logN)
}

// Equal compares two sets of parameters for equality.
func (p Parameters) Equal(other *Parameters) bool {
	return p.Parameters.Equal(&other.Parameters) && (p.PlaintextModulus() == other.PlaintextModulus())
}

// MarshalBinary returns a []byte representation of the parameter set.
// The representation corresponds to the JSON representation obtained
// from MarshalJSON.
func (p Parameters) MarshalBinary() ([]byte, error) {
	return p.MarshalJSON()
}

// UnmarshalBinary decodes a []byte into a parameter set struct.
func (p *Parameters) UnmarshalBinary(data []byte) (err error) {
	return p.UnmarshalJSON(data)
}

// MarshalJSON returns a JSON representation of this parameter set. See Marshal from the [encoding/json] package.
func (p Parameters) MarshalJSON() ([]byte, error) {
	return json.Marshal(p.ParametersLiteral())
}

// UnmarshalJSON reads a JSON representation of a parameter set into the receiver Parameter. See Unmarshal from the [encoding/json] package.
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
		LogN             int
		LogNthRoot       int
		Q                []uint64
		P                []uint64
		LogQ             []int
		LogP             []int
		Pow2Base         int
		Xe               map[string]interface{}
		Xs               map[string]interface{}
		RingType         ring.Type
		PlaintextModulus uint64
	}

	err = json.Unmarshal(b, &pl)
	if err != nil {
		return err
	}

	p.LogN = pl.LogN
	p.LogNthRoot = pl.LogNthRoot
	p.Q, p.P, p.LogQ, p.LogP = pl.Q, pl.P, pl.LogQ, pl.LogP
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
	p.PlaintextModulus = pl.PlaintextModulus
	return err
}
