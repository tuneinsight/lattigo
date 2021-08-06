package dckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/drlwe"
)

// Thresholdizer is the structure holding parameters for the Thresholdizing step
// of a Threshold MHE scheme.
type Thresholdizer struct {
	*drlwe.Thresholdizer
}

// Combiner is the structure holding parameters for the Combining step
// of a Threshold MHE scheme.
type Combiner struct {
	*drlwe.Combiner
}

// Combiner is the structure holding parameters for the Combining step
// of a Threshold MHE scheme, that caches the computed ring inverses.
type CombinerCache struct {
	*drlwe.CombinerCache
}

// NewThresholdizer creates a new Thresholdizer from parameters, that will be
// used to create threshold shares and threshold secret keys for every party.
func NewThresholdizer(params ckks.Parameters) *Thresholdizer {
	thresholdizer := new(Thresholdizer)
	thresholdizer.Thresholdizer = drlwe.NewThresholdizer(params.Parameters)
	return thresholdizer
}

// NewCombiner creates a new Combiner from parameters, that will be used to
// combine threshold secret and public keys into an additive share of the
// distributed secret key.
func NewCombiner(params ckks.Parameters, threshold int) *Combiner {
	combiner := new(Combiner)
	combiner.Combiner = drlwe.NewCombiner(params.Parameters, threshold)
	return combiner
}

// NewCombinerCache creates a new Combiner from an existing Combiner, that will
// be used to combine threshold secret and public keys into an additive share of
// the distributed secret key.
// When instantiated, a CombinerCache will precompute the inverses for the
// ThresholdPublicKeys it was given and store them in its cache.
// tpk is the party's public key.
// pks is a slice containing the public keys, of which we want the inverses of
// the difference with tpk precomputed.

func NewCombinerCache(comb *Combiner, tpk *drlwe.ThreshPublicKey, pks []*drlwe.ThreshPublicKey) *CombinerCache {
	combinercache := new(CombinerCache)
	combinercache.CombinerCache = drlwe.NewCombinerCache(comb.Combiner, tpk, pks)
	return combinercache
}
