// Package drlwe implements a distributed (or threshold) version of the CKKS scheme that enables secure multiparty computation solutions with secret-shared secret keys.
package drlwe

import (
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// CKGProtocol is the structure storing the parameters and and precomputations for the collective key generation protocol.
type CKGProtocol struct {
	params           rlwe.Parameters
	gaussianSamplerQ *ring.GaussianSampler
}

// ShallowCopy creates a shallow copy of CKGProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// CKGProtocol can be used concurrently.
func (ckg *CKGProtocol) ShallowCopy() *CKGProtocol {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	return &CKGProtocol{ckg.params, ring.NewGaussianSampler(prng, ckg.params.RingQ(), ckg.params.Sigma(), int(6*ckg.params.Sigma()))}
}

// CKGShare is a struct storing the CKG protocol's share.
type CKGShare struct {
	Value ringqp.Poly
}

// CKGCRP is a type for common reference polynomials in the CKG protocol.
type CKGCRP ringqp.Poly

// MarshalBinary encodes the target element on a slice of bytes.
func (share *CKGShare) MarshalBinary() (data []byte, err error) {
	data = make([]byte, share.Value.MarshalBinarySize64())
	if _, err = share.Value.Encode64(data); err != nil {
		return nil, err
	}
	return
}

// UnmarshalBinary decodes a slice of bytes on the target element.
func (share *CKGShare) UnmarshalBinary(data []byte) (err error) {
	_, err = share.Value.Decode64(data)
	return err
}

// NewCKGProtocol creates a new CKGProtocol instance
func NewCKGProtocol(params rlwe.Parameters) *CKGProtocol {
	ckg := new(CKGProtocol)
	ckg.params = params
	var err error
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ckg.gaussianSamplerQ = ring.NewGaussianSampler(prng, params.RingQ(), params.Sigma(), int(6*params.Sigma()))
	return ckg
}

// AllocateShare allocates the share of the CKG protocol.
func (ckg *CKGProtocol) AllocateShare() *CKGShare {
	return &CKGShare{ckg.params.RingQP().NewPoly()}
}

// SampleCRP samples a common random polynomial to be used in the CKG protocol from the provided
// common reference string.
func (ckg *CKGProtocol) SampleCRP(crs CRS) CKGCRP {
	crp := ckg.params.RingQP().NewPoly()
	ringqp.NewUniformSampler(crs, *ckg.params.RingQP()).Read(crp)
	return CKGCRP(crp)
}

// GenShare generates the party's public key share from its secret key as:
//
// crp*s_i + e_i
//
// for the receiver protocol. Has no effect is the share was already generated.
func (ckg *CKGProtocol) GenShare(sk *rlwe.SecretKey, crp CKGCRP, shareOut *CKGShare) {
	ringQP := ckg.params.RingQP()

	ckg.gaussianSamplerQ.Read(shareOut.Value.Q)

	if ringQP.RingP != nil {
		ringQP.ExtendBasisSmallNormAndCenter(shareOut.Value.Q, ckg.params.PCount()-1, nil, shareOut.Value.P)
	}

	levelQ, levelP := ckg.params.QCount()-1, ckg.params.PCount()-1
	ringQP.NTTLvl(levelQ, levelP, shareOut.Value, shareOut.Value)
	ringQP.MFormLvl(levelQ, levelP, shareOut.Value, shareOut.Value)

	ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, sk.Poly, ringqp.Poly(crp), shareOut.Value)
}

// AggregateShares aggregates a new share to the aggregate key
func (ckg *CKGProtocol) AggregateShares(share1, share2, shareOut *CKGShare) {
	ckg.params.RingQP().AddLvl(ckg.params.QCount()-1, ckg.params.PCount()-1, share1.Value, share2.Value, shareOut.Value)
}

// GenPublicKey return the current aggregation of the received shares as a bfv.PublicKey.
func (ckg *CKGProtocol) GenPublicKey(roundShare *CKGShare, crp CKGCRP, pubkey *rlwe.PublicKey) {
	pubkey.Value[0].Copy(roundShare.Value)
	pubkey.Value[1].Copy(ringqp.Poly(crp))
}
