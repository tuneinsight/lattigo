package drlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// CRS is an interface for Common Reference Strings.
// CRSs are PRNGs for which the read bits are the same for
// all parties.
type CRS interface {
	utils.PRNG
}

// CRP is an interface for Common References Polynomials.
// CRP are arrays of rlwe.PolyQP objects that can be retrieved
// using the Get method.
type CRP interface {
	Get(i int) rlwe.PolyQP
}

// crp is an implementation of the CRP interface that stores the CRPs
// in a slice pre-allocated slice.
type crp struct {
	crps []rlwe.PolyQP
}

// NewCRP creates a new CRP object with size common random polynomials
// at maximum level.
func NewCRP(params rlwe.Parameters, size int, crs CRS) CRP {
	return NewCRPAtLvl(params, size, params.QCount()-1, params.PCount()-1, crs)
}

// NewCRPAtLvl creates a new CRP object with size common random polynomials
// at the determined level.
func NewCRPAtLvl(params rlwe.Parameters, size, levelQ, levelP int, crs CRS) CRP {
	c := new(crp)
	c.crps = make([]rlwe.PolyQP, size)
	rQ := ring.NewUniformSampler(crs, params.RingQ())
	rP := ring.NewUniformSampler(crs, params.RingP())
	for i := range c.crps {
		c.crps[i] = params.RingQP().NewPolyLvl(levelQ, levelP)
		rQ.Read(c.crps[i].Q)
		rP.Read(c.crps[i].P)
	}
	return c
}

func (crp *crp) Get(i int) rlwe.PolyQP {
	if i < 0 || i >= len(crp.crps) {
		panic("invalid CRS index")
	}
	return crp.crps[i]
}
