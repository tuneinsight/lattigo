package drlwe

import (
	"fmt"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// ShamirPublicKey is a type for Shamir public keys in a thresholdizer.
type ShamirPublicKey uint64

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
	GenShamirPolynomial(threshold int, sk *rlwe.SecretKey) (*ShamirPolynomial, error)

	AllocateShamirSecretShare() *ShamirSecretShare
	GenShamirSecretShare(recipient ShamirPublicKey, secretPoly *ShamirPolynomial, shareOut *ShamirSecretShare)
	AggregateShares(share1, share2, outShare *ShamirSecretShare)
}

// Combiner is an interface for the combining phase of a RLWE threshold secret sharing protocol.
type Combiner interface {
	GenAdditiveShare(activePk []ShamirPublicKey, tpk ShamirPublicKey, tsks *ShamirSecretShare, tsk *rlwe.SecretKey)
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

// GenShamirPolynomial initiates a ShareGenPoly by sampling a random polynomial of
// degree threshold-1 with a constant term equal to the given secret key's value.
func (thresholdizer *Thresholdizer) GenShamirPolynomial(threshold int, sk *rlwe.SecretKey) (*ShamirPolynomial, error) {
	if threshold < 1 {
		return nil, fmt.Errorf("threshold should be >= 1")
	}
	gen := &ShamirPolynomial{coeffs: make([]*ring.Poly, int(threshold))}
	gen.coeffs[0] = sk.Value // using the sk Poly directly since gen.coeffs is private and never modified internally.
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
func (thresholdizer *Thresholdizer) GenShamirSecretShare(recipient ShamirPublicKey, secretPoly *ShamirPolynomial, shareOut *ShamirSecretShare) {
	thresholdizer.ringQP.EvalPolMontgomeryScalarNTT(secretPoly.coeffs, uint64(recipient), shareOut.Poly)
}

// AggregateShares aggregates two secret shares(by adding them), and stores them
// in outShare.
func (thresholdizer *Thresholdizer) AggregateShares(share1, share2, outShare *ShamirSecretShare) {
	thresholdizer.ringQP.Add(share1.Poly, share2.Poly, outShare.Poly)
}

// baseCombiner is a structure that holds the parameters for the combining phase of
// a threshold secret sharing protocol.
type baseCombiner struct {
	ringQP     *ring.Ring
	threshold  int
	tmp1, tmp2 []uint64
}

//NewCombiner creates a new Combiner.
func NewCombiner(params rlwe.Parameters, threshold int) Combiner {
	cmb := new(baseCombiner)
	cmb.ringQP = params.RingQP()
	cmb.threshold = threshold
	cmb.tmp1, cmb.tmp2 = make([]uint64, params.QPCount()), make([]uint64, params.QPCount())
	return cmb
}

// GenAdditiveShare generates an additive share of a cohort's secret key from a slice con-
// taining all active player's threshold public keys and a party's public and
// secret keys. Stores the result in out_key.
func (cmb *baseCombiner) GenAdditiveShare(actives []ShamirPublicKey, ownPublic ShamirPublicKey, ownSecret *ShamirSecretShare, skOut *rlwe.SecretKey) {

	if len(actives) < cmb.threshold {
		panic("Not enough active players to combine threshold shares.")
	}

	prod := cmb.tmp2
	for i, qi := range cmb.ringQP.Modulus {
		prod[i] = ring.MForm(1, qi, cmb.ringQP.BredParams[i])
	}

	for _, active := range actives[:cmb.threshold] {
		//Lagrange Interpolation with the public threshold key of other active players
		if active != ownPublic {
			cmb.lagrangeCoeff(ownPublic, active, cmb.tmp1)
			for i, qi := range cmb.ringQP.Modulus {
				prod[i] = ring.MRedConstant(prod[i], cmb.tmp1[i], qi, cmb.ringQP.MredParams[i])
			}
		}
	}

	cmb.ringQP.MulScalarCRT(ownSecret.Poly, prod, skOut.Value)
}

// lagrangeCoeff computes the difference between the two given keys and stores
// its multiplicative inverse in pol_out, caching it as well.
func (cmb *baseCombiner) lagrangeCoeff(thisKey ShamirPublicKey, thatKey ShamirPublicKey, lagCoeff []uint64) {

	this := uint64(thisKey)
	that := uint64(thatKey)

	for i, qi := range cmb.ringQP.Modulus {
		if this > that {
			lagCoeff[i] = that + qi - this
		} else {
			lagCoeff[i] = that - this
		}
	}

	cmb.ringQP.InverseCRT(lagCoeff)

	for i, invi := range lagCoeff {
		lagCoeff[i] = ring.MRedConstant(invi, that, cmb.ringQP.Modulus[i], cmb.ringQP.MredParams[i])
	}
}

// CachedCombiner is a structure that holds the parameters for the combining phase of
// a threshold secret sharing protocol, augmented with a stateful cache.
type CachedCombiner struct {
	*baseCombiner
	lagrangeCoeffs map[ShamirPublicKey][]uint64
}

// NewCachedCombiner creates a new combiner with cache from parameters.
func NewCachedCombiner(params rlwe.Parameters, threshold int) *CachedCombiner {
	ccmb := new(CachedCombiner)
	ccmb.baseCombiner = NewCombiner(params, threshold).(*baseCombiner)
	ccmb.lagrangeCoeffs = make(map[ShamirPublicKey][]uint64)
	return ccmb
}

// GenAdditiveShare generates an additive share of a cohort's secret key from the values
// in the cache and a party's secret keys. Assumes the inverse corresponding
// to every active's party is in the cache. Stores the result in out_key.
func (cmb *CachedCombiner) GenAdditiveShare(actives []ShamirPublicKey, ownPublic ShamirPublicKey, ownSecret *ShamirSecretShare, skOut *rlwe.SecretKey) {

	r := cmb.ringQP
	skOut.Value.Copy(ownSecret.Poly)
	for _, active := range actives {
		//Lagrange Interpolation with the public threshold key of other active players
		if active != ownPublic {
			lagrangeCoeff := cmb.lagrangeCoeff(ownPublic, active)
			cmb.ringQP.MulScalarCRT(skOut.Value, lagrangeCoeff, skOut.Value)
		}
	}

	r.Reduce(skOut.Value, skOut.Value)
}

// ClearCache replaces the cache of a combiner by an empty one.
func (cmb *CachedCombiner) ClearCache() {
	cmb.lagrangeCoeffs = make(map[ShamirPublicKey][]uint64)
}

// Precompute caches the inverses of the differences between tpk and each of
// pks (as needed for Lagrange interpolation)
func (cmb *CachedCombiner) Precompute(others []ShamirPublicKey, own ShamirPublicKey) {
	for _, key := range others {
		if own != key {
			_ = cmb.lagrangeCoeff(own, key)
		}
	}
}

// lagrangeCoeff computes the difference between the two given keys and stores
// its multiplicative inverse in pol_out, caching it as well.
func (cmb *CachedCombiner) lagrangeCoeff(thisKey ShamirPublicKey, thatKey ShamirPublicKey) (lagCoeff []uint64) {
	_, found := cmb.lagrangeCoeffs[thatKey]
	if !found {
		//Inverse not in the cache, we have to compute it
		cmb.lagrangeCoeffs[thatKey] = make([]uint64, len(cmb.ringQP.Modulus))
		cmb.baseCombiner.lagrangeCoeff(thisKey, thatKey, cmb.lagrangeCoeffs[thatKey])
	}
	return cmb.lagrangeCoeffs[thatKey]
}
