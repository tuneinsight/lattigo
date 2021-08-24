package drlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// CKGCRP is a type for the common reference polynomial of the collective key protocol.
type CKGCRP rlwe.PolyQP

// RKGCRP is a type for the common reference polynomial of the relinearization key protocol.
type RKGCRP []rlwe.PolyQP

// RTGCRP is a type for the common reference polynomial of the rotation key protocol.
type RTGCRP []rlwe.PolyQP

// RefreshCRP is a type for the common reference polynomial of the refresh protocol.
type RefreshCRP *ring.Poly

// UniformSampler is a struct storing samples for the drlwe package.
type UniformSampler struct {
	uniformSamplerQ *ring.UniformSampler
	uniformSamplerP *ring.UniformSampler
}

// NewUniformSampler creates a new UniformSampler.
func NewUniformSampler(key []byte, params rlwe.Parameters) (uniSampler UniformSampler, err error) {
	var prng utils.PRNG
	if prng, err = utils.NewKeyedPRNG(key); err != nil {
		return UniformSampler{}, nil
	}
	uniSampler.uniformSamplerQ = ring.NewUniformSampler(prng, params.RingQ())
	uniSampler.uniformSamplerP = ring.NewUniformSampler(prng, params.RingP())
	return
}

// Read samples new random polynomials on the input crp.
func (uniSampler *UniformSampler) Read(crp interface{}) {
	switch crp := crp.(type) {
	case *ring.Poly:
		uniSampler.uniformSamplerQ.Read(crp)
	case RefreshCRP:
		uniSampler.uniformSamplerQ.Read(crp)
	case rlwe.PolyQP:
		uniSampler.uniformSamplerQ.Read(crp.Q)
		uniSampler.uniformSamplerP.Read(crp.P)
	case CKGCRP:
		uniSampler.uniformSamplerQ.Read(crp.Q)
		uniSampler.uniformSamplerP.Read(crp.P)
	case []rlwe.PolyQP:
		for i := range crp {
			uniSampler.uniformSamplerQ.Read(crp[i].Q)
			uniSampler.uniformSamplerP.Read(crp[i].P)
		}
	case RKGCRP:
		for i := range crp {
			uniSampler.uniformSamplerQ.Read(crp[i].Q)
			uniSampler.uniformSamplerP.Read(crp[i].P)
		}
	case RTGCRP:
		for i := range crp {
			uniSampler.uniformSamplerQ.Read(crp[i].Q)
			uniSampler.uniformSamplerP.Read(crp[i].P)
		}
	default:
		panic("invalid crp struct, must be either *ring.Poly, [2]rlwe.PolyQP or [][2]rlwe.PolyQP")
	}
}
