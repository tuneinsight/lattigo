package ckks

import (
	"encoding/json"
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/ring/distribution"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

const (
	DefaultNTTFlag = true
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
	LogN     int
	Q        []uint64
	P        []uint64
	LogQ     []int `json:",omitempty"`
	LogP     []int `json:",omitempty"`
	Pow2Base int
	Xe       distribution.Distribution
	Xs       distribution.Distribution
	RingType ring.Type
	LogScale int
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
		DefaultNTTFlag: DefaultNTTFlag,
		DefaultScale:   rlwe.NewScale(math.Exp2(float64(p.LogScale))),
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

	if !rlweParams.DefaultNTTFlag() {
		return Parameters{}, fmt.Errorf("provided RLWE parameters are invalid for CKKS scheme (DefaultNTTFlag must be true)")
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
		LogN:     p.LogN(),
		Q:        p.Q(),
		P:        p.P(),
		Pow2Base: p.Pow2Base(),
		Xe:       p.Xe(),
		Xs:       p.Xs(),
		RingType: p.RingType(),
		LogScale: p.LogScale(),
	}
}

// DefaultPrecision returns the default precision in bits of the plaintext values which
// is max(53, log2(DefaultScale)).
func (p Parameters) DefaultPrecision() (prec uint) {
	if log2scale := math.Log2(p.DefaultScale().Float64()); log2scale <= 53 {
		prec = 53
	} else {
		prec = uint(log2scale)
	}

	return
}

// MaxDepth returns MaxLevel / DefaultScaleModuliRatio which is the maximum number of multiplicaitons
// followed by a rescaling that can be carried out with on a ciphertext with the DefaultScale.
func (p Parameters) MaxDepth() int {
	return p.MaxLevel() / p.DefaultScaleModuliRatio()
}

// DefaultScaleModuliRatio returns the default ratio between the scaling factor and moduli.
// This default ratio is computed as ceil(DefaultScalingFactor/2^{60}).
func (p Parameters) DefaultScaleModuliRatio() int {
	return int(math.Ceil(math.Log2(p.DefaultScale().Float64()) / 60.0))
}

// MaxLevel returns the maximum ciphertext level
func (p Parameters) MaxLevel() int {
	return p.QCount() - 1
}

// MaxSlots returns the theoretical maximum of plaintext slots allowed by the ring degree
func (p Parameters) MaxSlots() int {
	switch p.RingType() {
	case ring.Standard:
		return p.N() >> 1
	case ring.ConjugateInvariant:
		return p.N()
	default:
		panic("cannot MaxSlots: invalid ring type")
	}
}

// MaxLogSlots returns the log of the maximum number of slots enabled by the parameters
func (p Parameters) MaxLogSlots() int {
	switch p.RingType() {
	case ring.Standard:
		return p.LogN() - 1
	case ring.ConjugateInvariant:
		return p.LogN()
	default:
		panic("cannot MaxLogSlots: invalid ring type")
	}
}

// LogScale returns the log2 of the default scaling factor.
func (p Parameters) LogScale() int {
	return int(math.Round(math.Log2(p.DefaultScale().Float64())))
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

// GaloisElementsForLinearTransform generates the list of rotations needed for the evaluation of a linear transform
// with the provided list of non-zero diagonals, logSlots encoding and BSGSratio.
// If logBSGSRatio < 0, then provides the rotations needed for an evaluation without the BSGS approach.
func (p Parameters) GaloisElementsForLinearTransform(nonZeroDiags interface{}, logSlots, logBSGSRatio int) (galEls []uint64) {
	slots := 1 << logSlots
	if logBSGSRatio < 0 {
		_, _, rotN2 := rlwe.BSGSIndex(nonZeroDiags, slots, slots)

		galEls = make([]uint64, len(rotN2))

		for i := range rotN2 {
			galEls[i] = p.GaloisElementForColumnRotationBy(rotN2[i])
		}

		return
	}

	N1 := rlwe.FindBestBSGSRatio(nonZeroDiags, slots, logBSGSRatio)

	_, rotN1, rotN2 := rlwe.BSGSIndex(nonZeroDiags, slots, N1)

	rots := utils.GetDistincts(append(rotN1, rotN2...))

	galEls = make([]uint64, len(rots))
	for i, k := range rots {
		galEls[i] = p.GaloisElementForColumnRotationBy(k)
	}

	return
}

// Equal compares two sets of parameters for equality.
func (p Parameters) Equal(other Parameters) bool {
	return p.Parameters.Equal(other.Parameters)
}

// MarshalBinary returns a []byte representation of the parameter set.
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
