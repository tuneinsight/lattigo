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

//Equals compares two ThreshPublicKey for equality.
func (tpk *ShamirPublicKey) Equals(other *ShamirPublicKey) bool {
	return tpk.Poly == other.Poly || tpk.Poly.Equals(other.Poly)
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
// RLWE thresholdizer.
type ThresholdizerProtocol interface {
	GenShamirPublicKey() *ShamirPublicKey
	GenShamirPolynomial(threshold int, sk *rlwe.SecretKey) (*ShamirPolynomial, error)

	AllocateShamirSecretShare() *ShamirSecretShare
	GenShamirSecretShare(recipient *ShamirPublicKey, secretPoly *ShamirPolynomial, shareOut *ShamirSecretShare)
	AggregateShares(share1, share2, outShare *ShamirSecretShare)
}

// Combiner is an interface for the combining phase of a RLWE threshold secret sharing protocol.
type Combiner interface {
	GenAdditiveShare(activePk []*ShamirPublicKey, tpk *ShamirPublicKey, tsks *ShamirSecretShare, tsk *rlwe.SecretKey)
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
func (thresholdizer *Thresholdizer) GenShamirSecretShare(recipient *ShamirPublicKey, secretPoly *ShamirPolynomial, shareOut *ShamirSecretShare) {
	thresholdizer.ringQP.EvalPolMontgomeryNTT(secretPoly.coeffs, recipient.Poly, shareOut.Poly)
}

// AggregateShares aggregates two secret shares(by adding them), and stores them
// in outShare.
func (thresholdizer *Thresholdizer) AggregateShares(share1, share2, outShare *ShamirSecretShare) {
	thresholdizer.ringQP.Add(share1.Poly, share2.Poly, outShare.Poly)
}

// baseCombiner is a structure that holds the parameters for the combining phase of
// a threshold secret sharing protocol.
type baseCombiner struct {
	ringQP    *ring.Ring
	threshold int
	tmp       *ring.Poly
}

//NewCombiner creates a new Combiner.
func NewCombiner(params rlwe.Parameters, threshold int) Combiner {
	cmb := new(baseCombiner)
	cmb.ringQP = params.RingQP()
	cmb.threshold = threshold
	cmb.tmp = cmb.ringQP.NewPoly()
	return cmb
}

// GenAdditiveShare generates an additive share of a cohort's secret key from a slice con-
// taining all active player's threshold public keys and a party's public and
// secret keys. Stores the result in out_key.
func (cmb *baseCombiner) GenAdditiveShare(actives []*ShamirPublicKey, ownPublic *ShamirPublicKey, ownSecret *ShamirSecretShare, skOut *rlwe.SecretKey) {

	if len(actives) < cmb.threshold {
		panic("Not enough active players to combine threshold shares.")
	}

	skOut.Value.Copy(ownSecret.Poly)
	for _, active := range actives[:cmb.threshold] {
		//Lagrange Interpolation with the public threshold key of other active players
		if !active.Equals(ownPublic) {
			cmb.getLagrangeCoeff(ownPublic, active, cmb.tmp)
			cmb.ringQP.MulCoeffsMontgomeryConstant(skOut.Value, cmb.tmp, skOut.Value)
		}
	}

	cmb.ringQP.Reduce(skOut.Value, skOut.Value)
}

// getLagrangeCoeff computes the difference between the two given keys and stores
// its multiplicative inverse in pol_out, caching it as well.
func (cmb *baseCombiner) getLagrangeCoeff(thisKey *ShamirPublicKey, thatKey *ShamirPublicKey, polOut *ring.Poly) {
	//Inverse not in the cache, we have to compute it
	cmb.ringQP.SubNoMod(thatKey.Poly, thisKey.Poly, polOut)
	if !cmb.ringQP.IsInvertible(polOut) {
		panic("keys yield a non-invertible difference")
	}
	cmb.ringQP.InvMultPolyMontgomeryNTT(polOut, polOut)
	cmb.ringQP.MulCoeffsMontgomeryConstant(polOut, thatKey.Poly, polOut)
}

// CachedCombiner is a structure that holds the parameters for the combining phase of
// a threshold secret sharing protocol, augmented with a stateful cache.
type CachedCombiner struct {
	*baseCombiner
	lagrangeCoeffs map[*ShamirPublicKey]*ring.Poly
}

// NewCachedCombiner creates a new combiner with cache from parameters.
func NewCachedCombiner(params rlwe.Parameters, threshold int) *CachedCombiner {
	ccmb := new(CachedCombiner)
	ccmb.baseCombiner = NewCombiner(params, threshold).(*baseCombiner)
	ccmb.lagrangeCoeffs = make(map[*ShamirPublicKey]*ring.Poly)
	return ccmb
}

// GenAdditiveShare generates an additive share of a cohort's secret key from the values
// in the cache and a party's secret keys. Assumes the inverse corresponding
// to every active's party is in the cache. Stores the result in out_key.
func (cmb *CachedCombiner) GenAdditiveShare(actives []*ShamirPublicKey, ownPublic *ShamirPublicKey, ownSecret *ShamirSecretShare, skOut *rlwe.SecretKey) {

	r := cmb.ringQP
	skOut.Value.Copy(ownSecret.Poly)
	for _, active := range actives {
		//Lagrange Interpolation with the public threshold key of other active players
		if !active.Equals(ownPublic) {
			lagrangeCoeff := cmb.getLagrangeCoeff(ownPublic, active)
			cmb.ringQP.MulCoeffsMontgomeryConstant(skOut.Value, lagrangeCoeff, skOut.Value)
		}
	}

	r.Reduce(skOut.Value, skOut.Value)
}

// ClearCache replaces the cache of a combiner by an empty one.
func (cmb *CachedCombiner) ClearCache() {
	cmb.lagrangeCoeffs = make(map[*ShamirPublicKey]*ring.Poly)
}

// Precompute caches the inverses of the differences between tpk and each of
// pks (as needed for Lagrange interpolation)
func (cmb *CachedCombiner) Precompute(others []*ShamirPublicKey, own *ShamirPublicKey) {
	for _, key := range others {
		if !own.Equals(key) {
			_ = cmb.getLagrangeCoeff(own, key)
		}
	}
}

// getLagrangeCoeff computes the difference between the two given keys and stores
// its multiplicative inverse in pol_out, caching it as well.
func (cmb *CachedCombiner) getLagrangeCoeff(thisKey *ShamirPublicKey, thatKey *ShamirPublicKey) (polOut *ring.Poly) {
	_, found := cmb.lagrangeCoeffs[thatKey]
	if !found {
		//Inverse not in the cache, we have to compute it
		cmb.lagrangeCoeffs[thatKey] = cmb.ringQP.NewPoly()
		cmb.baseCombiner.getLagrangeCoeff(thisKey, thatKey, cmb.lagrangeCoeffs[thatKey])
	}
	return cmb.lagrangeCoeffs[thatKey]
}
