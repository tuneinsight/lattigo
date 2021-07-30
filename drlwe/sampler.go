package drlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

type UniformSampler struct {
	beta            int
	uniformSamplerQ *ring.UniformSampler
	uniformSamplerP *ring.UniformSampler
}

func NewUniformSampler(key []byte, params rlwe.Parameters) (uniSampler UniformSampler, err error) {
	var prng utils.PRNG
	if prng, err = utils.NewKeyedPRNG(key); err != nil {
		return UniformSampler{}, nil
	}

	uniSampler.uniformSamplerQ = ring.NewUniformSampler(prng, params.RingQ())
	uniSampler.uniformSamplerP = ring.NewUniformSampler(prng, params.RingP())

	return
}

// ReadQP samples a pair of random polynomials on crp.
// The first polynomial will be in modulus Q and the second in modulus P.
func (uniSampler *UniformSampler) ReadQP(crp [2]*ring.Poly) {
	uniSampler.uniformSamplerQ.Read(crp[0])
	uniSampler.uniformSamplerP.Read(crp[1])
}

// ReadNew samples a new random polynomial in modulus Q.
func (uniSampler *UniformSampler) ReadQNew() *ring.Poly {
	return uniSampler.uniformSamplerQ.ReadNew()
}

// ReadLvlQNew samples a new random polynomial in modulus Q.
func (uniSampler *UniformSampler) ReadLvlQNew(level int) *ring.Poly {
	return uniSampler.uniformSamplerQ.ReadLvlNew(level)
}

// ReadQPNew samples a new pair of random polynomials.
// The first polynomial will be in modulus Q and the second in modulus P.
func (uniSampler *UniformSampler) ReadQPNew() [2]*ring.Poly {
	return [2]*ring.Poly{uniSampler.uniformSamplerQ.ReadNew(), uniSampler.uniformSamplerP.ReadNew()}
}

// ReadLvlQPNew samples a new pair of random polynomials at the specified levels.
// The first polynomial will be in modulus Q and the second in modulus P.
func (uniSampler *UniformSampler) ReadLvlQPNew(levelQ, levelP int) [2]*ring.Poly {
	return [2]*ring.Poly{uniSampler.uniformSamplerQ.ReadLvlNew(levelQ), uniSampler.uniformSamplerP.ReadLvlNew(levelP)}
}

// ReadQPVectorNew samples a new vector of beta random polynomials.
// The first polynomials will be in modulus Q and the second in modulus P.
func (uniSampler *UniformSampler) ReadQPVectorNew(beta int) [][2]*ring.Poly {
	crp := make([][2]*ring.Poly, beta)
	for i := range crp {
		crp[i][0] = uniSampler.uniformSamplerQ.ReadNew()
		crp[i][1] = uniSampler.uniformSamplerP.ReadNew()
	}
	return crp
}

// ReadQPVectorLvlNew samples a new vector of beta random polynomials at the specified levels.
// The first polynomials will be in modulus Q and the second in modulus P.
func (uniSampler *UniformSampler) ReadQPVectorLvlNew(beta int, levelQ, levelP int) [][2]*ring.Poly {
	crp := make([][2]*ring.Poly, beta)
	for i := range crp {
		crp[i][0] = uniSampler.uniformSamplerQ.ReadLvlNew(levelQ)
		crp[i][1] = uniSampler.uniformSamplerP.ReadLvlNew(levelP)
	}
	return crp
}

// ReadQPVector samples a new vector of beta random polynomials on crp.
// The first polynomials will be in modulus Q and the second in modulus P.
func (uniSampler *UniformSampler) ReadQPVector(crp [][2]*ring.Poly) {
	for i := range crp {
		uniSampler.uniformSamplerQ.Read(crp[i][0])
		uniSampler.uniformSamplerP.Read(crp[i][1])
	}
}
