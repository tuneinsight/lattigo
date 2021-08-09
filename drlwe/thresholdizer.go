package drlwe

import (
	"fmt"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// ShamirPublicKey is a type for Shamir public keys in a thresholdizer.
type ShamirPublicKey struct {
	*ring.Poly
}

// ShamirPolynomial is a type for a share-generating polynomial.
type ShamirPolynomial struct {
	coeffs []*ring.Poly
}

// ShamirSecretShare is a type for a share of a secret key in a thresholdizer.
type ShamirSecretShare struct {
	*ring.Poly
}

// ThresholdizerProtocol is an interface describing the local steps of a generic
// DRLWE thresholdizer.
type ThresholdizerProtocol interface {
	GenShamirPublicKey() *ShamirPublicKey
	GenShamirPolynomial(threshold int, sk *rlwe.SecretKey) (*ShamirPolynomial, error)

	AllocateShamirSecretShare() *ShamirSecretShare
	GenShamirSecretShare(pubPoint *ShamirPublicKey, secretPoly *ShamirPolynomial, shareOut *ShamirSecretShare)
	AggregateShares(share1, share2, outShare *ShamirSecretShare)
}

//--------------------------------THRESHOLDIZING--------------------------------

//Thresholdizer is the structure containing the parameters for a thresholdizer.
type Thresholdizer struct {
	ringQP    *ring.Ring
	samplerQP *ring.UniformSampler
}

// NewThresholdizer creates a new Thresholdizer instance from parameters.
func NewThresholdizer(params rlwe.Parameters) *Thresholdizer {

	thresholdizer := new(Thresholdizer)
	var err error
	thresholdizer.ringQP = params.RingQP()

	prng, err := utils.NewPRNG()
	if err != nil {
		panic("Error in Thresholdizer initalization : error in PRNG generation")
	}

	thresholdizer.samplerQP = ring.NewUniformSampler(prng, thresholdizer.ringQP)

	return thresholdizer
}

// GenShamirPublicKey generates a threshold public key from an id. Useful to avoid
// having to broadcast polynomials during the setup.
func (thresholdizer *Thresholdizer) GenShamirPublicKey() *ShamirPublicKey {
	idPRNG, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	gen := ring.NewUniformSampler(idPRNG, thresholdizer.ringQP)
	return &ShamirPublicKey{gen.ReadNew()}
}

// GenShamirPolynomial initiates a ShareGenPoly by sampling a random polynomial of
// degree threshold-1 with a constant term equal to the given secret key's value.
func (thresholdizer *Thresholdizer) GenShamirPolynomial(threshold int, sk *rlwe.SecretKey) (*ShamirPolynomial, error) {
	if threshold < 1 {
		return nil, fmt.Errorf("threshold should be >= 1")
	}
	gen := &ShamirPolynomial{coeffs: make([]*ring.Poly, int(threshold))}
	gen.coeffs[0] = sk.Value.CopyNew()
	for i := 1; i < threshold; i++ {
		gen.coeffs[i] = thresholdizer.samplerQP.ReadNew()
	}
	return gen, nil
}

// AllocateThresholdSecretShare allocates a Threshold secret share.
func (thresholdizer *Thresholdizer) AllocateThresholdSecretShare() *ShamirSecretShare {
	return &ShamirSecretShare{thresholdizer.ringQP.NewPoly()}
}

// GenShamirSecretShare generates a secret share for a given threshold public key.
// Stores the result in share_out. This result should be sent to the given
// threshold public key's owner.
func (thresholdizer *Thresholdizer) GenShamirSecretShare(pubPoint *ShamirPublicKey, gen *ShamirPolynomial, shareOut *ShamirSecretShare) {
	thresholdizer.ringQP.EvalPolMontgomeryNTT(gen.coeffs, pubPoint.Poly, shareOut.Poly)
}

// AggregateShares aggregates two secret shares(by adding them), and stores them
// in outShare.
func (thresholdizer *Thresholdizer) AggregateShares(share1, share2, outShare *ShamirSecretShare) {
	thresholdizer.ringQP.Add(share1.Poly, share2.Poly, outShare.Poly)
}

//-------------------------------COMBINING -------------------------------------

// Combiner is a structure that holds the parameters for the combining phase of
// a threshold secret sharing protocol.
type Combiner struct {
	ringQP    *ring.Ring
	threshold uint64
	tmp       *ring.Poly
}

//NewCombiner creates a new Combiner.
func NewCombiner(params rlwe.Parameters, threshold int) *Combiner {
	combiner := new(Combiner)
	combiner.ringQP = params.RingQP()
	combiner.tmp = combiner.ringQP.NewPoly()
	return combiner
}

// GenAdditiveShare generates an additive share of a cohort's secret key from a slice con-
// taining all active player's threshold public keys and a party's public and
// secret keys. Stores the result in out_key.
func (combiner *Combiner) GenAdditiveShare(activePoints []*ShamirPublicKey, tpk *ShamirPublicKey, tsks *ShamirSecretShare, tsk *rlwe.SecretKey) {

	if uint64(len(activePoints)) < combiner.threshold {
		panic("Not enough active players to combine threshold shares.")
	}

	tsk.Value.Copy(tsks.Poly)
	for _, activePoint := range activePoints {
		//Lagrange Interpolation with the public threshold key of other active players
		if !combiner.Equal(activePoint, tpk) {
			combiner.ringQP.MulCoeffsMontgomeryConstant(tsk.Value, activePoint.Poly, tsk.Value)
			combiner.ringQP.SubNoMod(activePoint.Poly, tpk.Poly, combiner.tmp)
			combiner.ringQP.InvMultPolyMontgomeryNTT(combiner.tmp, combiner.tmp)
			combiner.ringQP.MulCoeffsMontgomeryConstant(tsk.Value, combiner.tmp, tsk.Value)
		}
	}

	combiner.ringQP.Reduce(tsk.Value, tsk.Value)
}

//Equal compares two ThreshPublicKey for equality.
func (combiner *Combiner) Equal(tpk1, tpk2 *ShamirPublicKey) bool {
	return combiner.ringQP.Equal(tpk1.Poly, tpk2.Poly)
}

// CombinerCached is a structure that holds the parameters for the combining phase of
// a threshold secret sharing protocol, augmented with a stateful cache.
type CombinerCached struct {
	*Combiner
	inverses map[*ShamirPublicKey]*ring.Poly
}

// NewCombinerCache creates a new combiner with cache from parameters.
func NewCombinerCache(params rlwe.Parameters, threshold int) *CombinerCached {
	combinercache := new(CombinerCached)
	combinercache.Combiner = NewCombiner(params, threshold)
	combinercache.inverses = make(map[*ShamirPublicKey]*ring.Poly)
	return combinercache
}

// Precompute caches the inverses of the differences between tpk and each of
// pks (as needed for Lagrange interpolation)
func (combiner *CombinerCached) Precompute(tpk *ShamirPublicKey, pks []*ShamirPublicKey) {
	for _, key := range pks {
		if !combiner.Equal(tpk, key) {
			combiner.getLagrangeCoeff(tpk, key)
		}
	}
}

// getLagrangeCoeff computes the difference between the two given keys and stores
// its multiplicative inverse in pol_out, caching it as well.
func (combiner *CombinerCached) getLagrangeCoeff(thisKey *ShamirPublicKey, thatKey *ShamirPublicKey) (polOut *ring.Poly) {
	polOut, found := combiner.inverses[thatKey]
	if !found {
		//Inverse not in the cache, we have to compute it
		polOut = combiner.ringQP.NewPoly()
		combiner.ringQP.SubNoMod(thatKey.Poly, thisKey.Poly, polOut)
		if !combiner.ringQP.IsInvertible(polOut) {
			panic("keys yield a non-invertible difference")
		}
		combiner.ringQP.InvMultPolyMontgomeryNTT(polOut, polOut)
		combiner.ringQP.MulCoeffsMontgomeryConstant(polOut, thisKey.Poly, polOut)
		//Cache the result.
		combiner.inverses[thatKey] = polOut
	}
	return polOut
}

// GenFinalShare generates an additive share of a cohort's secret key from the values
// in the cache and a party's secret keys. Assumes the inverse corresponding
// to every active's party is in the cache. Stores the result in out_key.
func (combiner *CombinerCached) GenFinalShare(tsks *ShamirSecretShare, tsk *rlwe.SecretKey) {

	r := combiner.ringQP
	tsk.Value.Copy(tsks.Poly)
	for _, inv := range combiner.inverses {
		//Lagrange Interpolation with the threshold public key of other active players
		r.MulCoeffsMontgomeryConstant(tsk.Value, inv, tsk.Value)
	}

	r.Reduce(tsk.Value, tsk.Value)
}

// ClearCache replaces the cache of a combiner by an empty one.
func (combiner *CombinerCached) ClearCache() {
	combiner.inverses = make(map[*ShamirPublicKey]*ring.Poly)
}
