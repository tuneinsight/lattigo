package drlwe

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v3/utils"
)

// ShamirPublicKey is a type for Shamir public keys in a thresholdizer.
type ShamirPublicKey uint64

// ShamirPolynomial is a type for a share-generating polynomial.
type ShamirPolynomial struct {
	coeffs []ringqp.Poly
}

// ShamirSecretShare is a type for a share of a secret key in a thresholdizer.
type ShamirSecretShare struct {
	ringqp.Poly
}

//Thresholdizer is the structure containing the parameters for a thresholdizer.
type Thresholdizer struct {
	params   *rlwe.Parameters
	ringQP   *ringqp.Ring
	usampler ringqp.UniformSampler
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

	thresholdizer.usampler = ringqp.NewUniformSampler(prng, *params.RingQP())

	return thresholdizer
}

// GenShamirPolynomial initiates a ShareGenPoly by sampling a random polynomial of
// degree threshold-1 with a constant term equal to the given secret key's value.
func (thresholdizer *Thresholdizer) GenShamirPolynomial(threshold int, sk *rlwe.SecretKey) (*ShamirPolynomial, error) {
	if threshold < 1 {
		return nil, fmt.Errorf("threshold should be >= 1")
	}
	gen := &ShamirPolynomial{coeffs: make([]ringqp.Poly, int(threshold))}
	gen.coeffs[0] = sk.Value // using the sk Poly directly since gen.coeffs is private and never modified internally.
	for i := 1; i < threshold; i++ {
		gen.coeffs[i] = thresholdizer.ringQP.NewPoly()
		thresholdizer.usampler.Read(gen.coeffs[i])
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
	lvlQ, lvlP := thresholdizer.params.QCount()-1, thresholdizer.params.PCount()-1
	thresholdizer.ringQP.AddLvl(lvlQ, lvlP, share1.Poly, share2.Poly, outShare.Poly)
}

// Combiner is a structure that holds the parameters for the combining phase of
// a threshold secret sharing protocol.
type Combiner struct {
	ringQP         *ringqp.Ring
	threshold      int
	tmp1, tmp2     []uint64
	one            ring.RNSScalar
	lagrangeCoeffs map[ShamirPublicKey]ring.RNSScalar
}

//NewCombiner creates a new Combiner.
func NewCombiner(params rlwe.Parameters, own ShamirPublicKey, others []ShamirPublicKey, threshold int) *Combiner {
	cmb := new(Combiner)
	cmb.ringQP = params.RingQP()
	cmb.threshold = threshold
	cmb.tmp1, cmb.tmp2 = cmb.ringQP.NewScalar(), cmb.ringQP.NewScalar()
	cmb.one = cmb.ringQP.NewScalarFromUInt64(1)

	qlen := len(cmb.ringQP.RingQ.Modulus)
	for i, qi := range cmb.ringQP.RingQ.Modulus {
		cmb.one[i] = ring.MForm(cmb.one[i], qi, cmb.ringQP.RingQ.BredParams[i])
	}
	if cmb.ringQP.RingP != nil {
		for i, pi := range cmb.ringQP.RingP.Modulus {
			cmb.one[i+qlen] = ring.MForm(cmb.one[i+qlen], pi, cmb.ringQP.RingP.BredParams[i])
		}
	}

	// precomputes lagrange coefficient factors
	cmb.lagrangeCoeffs = make(map[ShamirPublicKey]ring.RNSScalar)
	for _, spk := range others {
		if spk != own {
			cmb.lagrangeCoeffs[spk] = cmb.ringQP.NewScalar()
			cmb.lagrangeCoeff(own, spk, cmb.lagrangeCoeffs[spk])
		}
	}

	return cmb
}

// GenAdditiveShare generates an additive share of a cohort's secret key from a slice con-
// taining all active player's threshold public keys and a party's public and
// secret keys. Stores the result in out_key.
func (cmb *Combiner) GenAdditiveShare(actives []ShamirPublicKey, ownPublic ShamirPublicKey, ownSecret *ShamirSecretShare, skOut *rlwe.SecretKey) {

	if len(actives) < cmb.threshold {
		panic("Not enough active players to combine threshold shares.")
	}

	prod := cmb.tmp2
	copy(prod, cmb.one)

	for _, active := range actives[:cmb.threshold] {
		//Lagrange Interpolation with the public threshold key of other active players
		if active != ownPublic {
			cmb.tmp1 = cmb.lagrangeCoeffs[active]
			cmb.ringQP.MulRNSScalar(prod, cmb.tmp1, prod)
		}
	}

	cmb.ringQP.MulScalarCRT(ownSecret.Poly, prod, skOut.Value)
}

// lagrangeCoeff computes the difference between the two given keys and stores
// its multiplicative inverse in pol_out, caching it as well.
func (cmb *Combiner) lagrangeCoeff(thisKey ShamirPublicKey, thatKey ShamirPublicKey, lagCoeff []uint64) {

	this := cmb.ringQP.NewScalarFromUInt64(uint64(thisKey))
	that := cmb.ringQP.NewScalarFromUInt64(uint64(thatKey))

	cmb.ringQP.SubRNSScalar(that, this, lagCoeff)

	cmb.ringQP.InverseCRT(lagCoeff)

	cmb.ringQP.MulRNSScalar(lagCoeff, that, lagCoeff)
}
