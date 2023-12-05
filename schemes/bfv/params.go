package bfv

import (
	"encoding/json"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/schemes/bgv"
)

// NewParameters instantiate a set of BFV parameters from the generic RLWE parameters and a plaintext modulus t.
// User must ensure that t = 1 mod 2n for 4 < n <= N, where N is the ring degree.
// It returns the empty parameters Parameters{} and a non-nil error if the specified parameters are invalid.
func NewParameters(rlweParams rlwe.Parameters, t uint64) (p Parameters, err error) {
	var pbgv bgv.Parameters
	pbgv, err = bgv.NewParameters(rlweParams, t)
	return Parameters{pbgv}, err
}

// NewParametersFromLiteral instantiate a set of BFV parameters from a ParametersLiteral specification.
// It returns the empty parameters Parameters{} and a non-nil error if the specified parameters are invalid.
//
// See `rlwe.NewParametersFromLiteral` for default values of the optional fields.
func NewParametersFromLiteral(pl ParametersLiteral) (p Parameters, err error) {
	var pbgv bgv.Parameters
	pbgv, err = bgv.NewParametersFromLiteral(bgv.ParametersLiteral(pl))
	return Parameters{pbgv}, err
}

// ParametersLiteral is a literal representation of BFV parameters.  It has public
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
type ParametersLiteral bgv.ParametersLiteral

// GetRLWEParametersLiteral returns the rlwe.ParametersLiteral from the target bfv.ParametersLiteral.
func (p ParametersLiteral) GetRLWEParametersLiteral() rlwe.ParametersLiteral {
	return bgv.ParametersLiteral(p).GetRLWEParametersLiteral()
}

// Parameters represents a parameter set for the BFV cryptosystem. Its fields are private and
// immutable. See ParametersLiteral for user-specified parameters.
type Parameters struct {
	bgv.Parameters
}

// Equal compares two sets of parameters for equality.
func (p Parameters) Equal(other *Parameters) bool {
	return p.Parameters.Equal(&other.Parameters)
}

// UnmarshalBinary decodes a []byte into a parameter set struct.
func (p *Parameters) UnmarshalBinary(data []byte) (err error) {
	return p.Parameters.UnmarshalJSON(data)
}

// UnmarshalJSON reads a JSON representation of a parameter set into the receiver Parameter. See `Unmarshal` from the `encoding/json` package.
func (p *Parameters) UnmarshalJSON(data []byte) (err error) {
	return p.Parameters.UnmarshalJSON(data)
}

func (p *ParametersLiteral) UnmarshalJSON(b []byte) (err error) {
	var pl struct {
		LogN             int
		LogNthRoot       int
		Q                []uint64
		P                []uint64
		LogQ             []int
		LogP             []int
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
