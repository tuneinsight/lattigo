package drlwe

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// Thresholdizer is a type for generating secret-shares of ringqp.Poly types such that
// the resulting sharing has a t-out-of-N-threshold access-structure. It implements the
// `Thresholdize` operation as presented in "An Efficient Threshold Access-Structure
// for RLWE-Based Multiparty Homomorphic Encryption" (2022) by Mouchet, C., Bertrand, E.,
// and Hubaux, J. P. (https://eprint.iacr.org/2022/780).
//
// See the `drlwe` package README.md.
type Thresholdizer struct {
	params   *rlwe.Parameters
	ringQP   *ringqp.Ring
	usampler ringqp.UniformSampler
}

// Combiner is a type for generating t-out-of-t additive shares from local t-out-of-N
// shares. It implements the `Combine` operation as presented in "An Efficient Threshold
// Access-Structure for RLWE-Based Multiparty Homomorphic Encryption" (2022) by Mouchet, C.,
// Bertrand, E., and Hubaux, J. P. (https://eprint.iacr.org/2022/780).
type Combiner struct {
	ringQP         *ringqp.Ring
	threshold      int
	tmp1, tmp2     []uint64
	one            ring.RNSScalar
	lagrangeCoeffs map[ShamirPublicPoint]ring.RNSScalar
}

// ShamirPublicPoint is a type for Shamir public point associated with a party identity within
// the t-out-of-N-threshold scheme.
//
// See Thresholdizer and Combiner types.
type ShamirPublicPoint uint64

// ShamirPolynomial represents a polynomial with ringqp.Poly coefficients. It is used by the
// Thresholdizer type to produce t-out-of-N-threshold shares of an ringqp.Poly.
//
// See Thresholdizer type.
type ShamirPolynomial struct {
	Coeffs []ringqp.Poly
}

// ShamirSecretShare represents a t-out-of-N-threshold secret-share.
//
// See Thresholdizer and Combiner types.
type ShamirSecretShare struct {
	ringqp.Poly
}

// NewThresholdizer creates a new Thresholdizer instance from parameters.
func NewThresholdizer(params rlwe.Parameters) *Thresholdizer {

	thr := new(Thresholdizer)
	thr.params = &params
	thr.ringQP = params.RingQP()

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(fmt.Errorf("could not initialize PRNG: %s", err))
	}

	thr.usampler = ringqp.NewUniformSampler(prng, *params.RingQP())

	return thr
}

// GenShamirPolynomial generates a new secret ShamirPolynomial to be used in the Thresholdizer.GenShamirSecretShare method.
// It does so by sampling a random polynomial of degree threshold - 1 and with its constant term equal to secret.
func (thr *Thresholdizer) GenShamirPolynomial(threshold int, secret *rlwe.SecretKey) (*ShamirPolynomial, error) {
	if threshold < 1 {
		return nil, fmt.Errorf("threshold should be >= 1")
	}
	gen := &ShamirPolynomial{Coeffs: make([]ringqp.Poly, int(threshold))}
	gen.Coeffs[0] = secret.Value.CopyNew()
	for i := 1; i < threshold; i++ {
		gen.Coeffs[i] = thr.ringQP.NewPoly()
		thr.usampler.Read(gen.Coeffs[i])
	}
	return gen, nil
}

// AllocateThresholdSecretShare allocates a ShamirSecretShare struct.
func (thr *Thresholdizer) AllocateThresholdSecretShare() *ShamirSecretShare {
	return &ShamirSecretShare{thr.ringQP.NewPoly()}
}

// GenShamirSecretShare generates a secret share for the given recipient, identified by its ShamirPublicPoint.
// The result is stored in ShareOut and should be sent to this party.
func (thr *Thresholdizer) GenShamirSecretShare(recipient ShamirPublicPoint, secretPoly *ShamirPolynomial, shareOut *ShamirSecretShare) {
	thr.ringQP.EvalPolyScalar(secretPoly.Coeffs, uint64(recipient), shareOut.Poly)
}

// AggregateShares aggregates two ShamirSecretShare and stores the result in outShare.
func (thr *Thresholdizer) AggregateShares(share1, share2, outShare *ShamirSecretShare) {
	if share1.LevelQ() != share2.LevelQ() || share1.LevelQ() != outShare.LevelQ() || share1.LevelP() != share2.LevelP() || share1.LevelP() != outShare.LevelP() {
		panic("shares level do not match")
	}
	thr.ringQP.AtLevel(share1.LevelQ(), share1.LevelP()).Add(share1.Poly, share2.Poly, outShare.Poly)
}

// NewCombiner creates a new Combiner struct from the parameters and the set of ShamirPublicPoints. Note that the other
// parameter may contain the instantiator's own ShamirPublicPoint.
func NewCombiner(params rlwe.Parameters, own ShamirPublicPoint, others []ShamirPublicPoint, threshold int) *Combiner {
	cmb := new(Combiner)
	cmb.ringQP = params.RingQP()
	cmb.threshold = threshold
	cmb.tmp1, cmb.tmp2 = cmb.ringQP.NewRNSScalar(), cmb.ringQP.NewRNSScalar()
	cmb.one = cmb.ringQP.NewRNSScalarFromUInt64(1)

	qlen := cmb.ringQP.RingQ.NbModuli()
	for i, table := range cmb.ringQP.RingQ.Tables {
		cmb.one[i] = ring.MForm(cmb.one[i], table.Modulus, table.BRedParams)
	}
	if cmb.ringQP.RingP != nil {
		for i, table := range cmb.ringQP.RingP.Tables {
			cmb.one[i+qlen] = ring.MForm(cmb.one[i+qlen], table.Modulus, table.BRedParams)
		}
	}

	// precomputes lagrange coefficient factors
	cmb.lagrangeCoeffs = make(map[ShamirPublicPoint]ring.RNSScalar)
	for _, spk := range others {
		if spk != own {
			cmb.lagrangeCoeffs[spk] = cmb.ringQP.NewRNSScalar()
			cmb.lagrangeCoeff(own, spk, cmb.lagrangeCoeffs[spk])
		}
	}

	return cmb
}

// GenAdditiveShare generates a t-out-of-t additive share of the secret from a local aggregated share ownSecret and the set of active identities, identified
// by their ShamirPublicPoint. It stores the resulting additive share in skOut.
func (cmb *Combiner) GenAdditiveShare(activesPoints []ShamirPublicPoint, ownPoint ShamirPublicPoint, ownShare *ShamirSecretShare, skOut *rlwe.SecretKey) {

	if len(activesPoints) < cmb.threshold {
		panic("cannot GenAdditiveShare: Not enough active players to combine threshold shares.")
	}

	prod := cmb.tmp2
	copy(prod, cmb.one)

	for _, active := range activesPoints[:cmb.threshold] {
		//Lagrange Interpolation with the public threshold key of other active players
		if active != ownPoint {
			cmb.tmp1 = cmb.lagrangeCoeffs[active]
			cmb.ringQP.MulRNSScalar(prod, cmb.tmp1, prod)
		}
	}

	cmb.ringQP.MulRNSScalarMontgomery(ownShare.Poly, prod, skOut.Value)
}

func (cmb *Combiner) lagrangeCoeff(thisKey ShamirPublicPoint, thatKey ShamirPublicPoint, lagCoeff []uint64) {

	this := cmb.ringQP.NewRNSScalarFromUInt64(uint64(thisKey))
	that := cmb.ringQP.NewRNSScalarFromUInt64(uint64(thatKey))

	cmb.ringQP.SubRNSScalar(that, this, lagCoeff)

	cmb.ringQP.Inverse(lagCoeff)

	cmb.ringQP.MulRNSScalar(lagCoeff, that, lagCoeff)
}

func (s *ShamirSecretShare) MarshalBinary() ([]byte, error) {
	return s.Poly.MarshalBinary()
}

func (s *ShamirSecretShare) UnmarshalBinary(b []byte) error {
	return s.Poly.UnmarshalBinary(b)
}
