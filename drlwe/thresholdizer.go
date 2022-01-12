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
	coeffs []rlwe.PolyQP
}

// ShamirSecretShare is a type for a share of a secret key in a thresholdizer.
type ShamirSecretShare struct {
	rlwe.PolyQP
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
	params   *rlwe.Parameters
	ringQP   *rlwe.RingQP
	samplerQ *ring.UniformSampler
	samplerP *ring.UniformSampler
}

// NewThresholdizer creates a new Thresholdizer instance from parameters.
func NewThresholdizer(params rlwe.Parameters) *Thresholdizer {

	thresholdizer := new(Thresholdizer)
	thresholdizer.params = &params
	thresholdizer.ringQP = params.RingQP()

	prng, err := utils.NewPRNG()
	if err != nil {
		panic("Error in Thresholdizer initalization : error in PRNG generation")
	}

	thresholdizer.samplerQ = ring.NewUniformSampler(prng, params.RingQ())
	thresholdizer.samplerP = ring.NewUniformSampler(prng, params.RingP())

	return thresholdizer
}

// GenShamirPolynomial initiates a ShareGenPoly by sampling a random polynomial of
// degree threshold-1 with a constant term equal to the given secret key's value.
func (thresholdizer *Thresholdizer) GenShamirPolynomial(threshold int, sk *rlwe.SecretKey) (*ShamirPolynomial, error) {
	if threshold < 1 {
		return nil, fmt.Errorf("threshold should be >= 1")
	}
	gen := &ShamirPolynomial{coeffs: make([]rlwe.PolyQP, int(threshold))}
	gen.coeffs[0] = sk.Value // using the sk Poly directly since gen.coeffs is private and never modified internally.
	for i := 1; i < threshold; i++ {
		gen.coeffs[i].Q = thresholdizer.samplerQ.ReadNew()
		gen.coeffs[i].P = thresholdizer.samplerP.ReadNew()
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
	thresholdizer.ringQP.EvalPolMontgomeryScalarNTT(secretPoly.coeffs, uint64(recipient), shareOut.PolyQP)
}

// AggregateShares aggregates two secret shares(by adding them), and stores them
// in outShare.
func (thresholdizer *Thresholdizer) AggregateShares(share1, share2, outShare *ShamirSecretShare) {
	lvlQ, lvlP := thresholdizer.params.QCount()-1, thresholdizer.params.PCount()-1
	thresholdizer.ringQP.AddLvl(lvlQ, lvlP, share1.PolyQP, share2.PolyQP, outShare.PolyQP)
}

// baseCombiner is a structure that holds the parameters for the combining phase of
// a threshold secret sharing protocol.
type baseCombiner struct {
	ringQP     *rlwe.RingQP
	threshold  int
	tmp1, tmp2 []uint64
	one        ring.Scalar
}

//NewCombiner creates a new Combiner.
func NewCombiner(params rlwe.Parameters, threshold int) Combiner {
	cmb := new(baseCombiner)
	cmb.ringQP = params.RingQP()
	cmb.threshold = threshold
	cmb.tmp1, cmb.tmp2 = make([]uint64, params.QPCount()), make([]uint64, params.QPCount())
	cmb.one = cmb.ringQP.NewScalarFromUInt64(1)

	qlen := len(cmb.ringQP.RingQ.Modulus)
	for i, qi := range cmb.ringQP.RingQ.Modulus {
		cmb.one[i] = ring.MForm(cmb.one[i], qi, cmb.ringQP.RingQ.BredParams[i])
	}
	for i, pi := range cmb.ringQP.RingP.Modulus {
		cmb.one[i+qlen] = ring.MForm(cmb.one[i+qlen], pi, cmb.ringQP.RingP.BredParams[i])
	}
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
	copy(prod, cmb.one)

	for _, active := range actives[:cmb.threshold] {
		//Lagrange Interpolation with the public threshold key of other active players
		if active != ownPublic {
			cmb.lagrangeCoeff(ownPublic, active, cmb.tmp1)
			cmb.ringQP.ScalarMulCRT(prod, cmb.tmp1, prod)
		}
	}

	cmb.ringQP.MulScalarCRT(ownSecret.PolyQP, prod, skOut.Value)
}

// lagrangeCoeff computes the difference between the two given keys and stores
// its multiplicative inverse in pol_out, caching it as well.
func (cmb *baseCombiner) lagrangeCoeff(thisKey ShamirPublicKey, thatKey ShamirPublicKey, lagCoeff []uint64) {

	this := cmb.ringQP.NewScalarFromUInt64(uint64(thisKey))
	that := cmb.ringQP.NewScalarFromUInt64(uint64(thatKey))

	cmb.ringQP.SubScalarCRT(that, this, lagCoeff)

	cmb.ringQP.InverseCRT(lagCoeff)

	cmb.ringQP.ScalarMulCRT(lagCoeff, that, lagCoeff)
}

// CachedCombiner is a structure that holds the parameters for the combining phase of
// a threshold secret sharing protocol, augmented with a stateful cache.
type CachedCombiner struct {
	*baseCombiner
	lagrangeCoeffs map[ShamirPublicKey]ring.Scalar
}

// NewCachedCombiner creates a new combiner with cache from parameters.
func NewCachedCombiner(params rlwe.Parameters, threshold int) *CachedCombiner {
	ccmb := new(CachedCombiner)
	ccmb.baseCombiner = NewCombiner(params, threshold).(*baseCombiner)
	ccmb.lagrangeCoeffs = make(map[ShamirPublicKey]ring.Scalar)
	return ccmb
}

// // GenAdditiveShare generates an additive share of a cohort's secret key from the values
// // in the cache and a party's secret keys. Assumes the inverse corresponding
// // to every active's party is in the cache. Stores the result in out_key.
// func (cmb *CachedCombiner) GenAdditiveShare(actives []ShamirPublicKey, ownPublic ShamirPublicKey, ownSecret *ShamirSecretShare, skOut *rlwe.SecretKey) {

// 	prod := cmb.tmp2
// 	for i, qi := range cmb.ringQP.Modulus {
// 		prod[i] = ring.MForm(1, qi, cmb.ringQP.BredParams[i])
// 	}

// 	for _, active := range actives[:cmb.threshold] {
// 		//Lagrange Interpolation with the public threshold key of other active players
// 		if active != ownPublic {
// 			lagCoeff := cmb.lagrangeCoeff(ownPublic, active)
// 			for i, qi := range cmb.ringQP.Modulus {
// 				prod[i] = ring.MRedConstant(prod[i], lagCoeff[i], qi, cmb.ringQP.MredParams[i])
// 			}
// 		}
// 	}

// 	cmb.ringQP.MulScalarCRT(ownSecret.Poly, prod, skOut.Value)
// }

// ClearCache replaces the cache of a combiner by an empty one.
func (cmb *CachedCombiner) ClearCache() {
	cmb.lagrangeCoeffs = make(map[ShamirPublicKey]ring.Scalar)
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
		cmb.lagrangeCoeffs[thatKey] = cmb.ringQP.NewScalar()
		cmb.baseCombiner.lagrangeCoeff(thisKey, thatKey, cmb.lagrangeCoeffs[thatKey])
	}
	return cmb.lagrangeCoeffs[thatKey]
}
