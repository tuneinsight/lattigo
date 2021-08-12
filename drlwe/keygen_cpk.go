//Package drlwe implements a distributed (or threshold) version of the CKKS scheme that enables secure multiparty computation solutions with secret-shared secret keys.
package drlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// CollectivePublicKeyGenerator is an interface describing the local steps of a generic RLWE CKG protocol.
type CollectivePublicKeyGenerator interface {
	AllocateShares() *CKGShare
	GenShare(sk *rlwe.SecretKey, crs rlwe.PolyQP, shareOut *CKGShare)
	AggregateShares(share1, share2, shareOut *CKGShare)
	GenPublicKey(aggregatedShare *CKGShare, crs rlwe.PolyQP, pubkey *rlwe.PublicKey)
}

// CKGProtocol is the structure storing the parameters and and precomputations for the collective key generation protocol.
type CKGProtocol struct {
	params           rlwe.Parameters
	gaussianSamplerQ *ring.GaussianSampler
}

// CKGShare is a struct storing the CKG protocol's share.
type CKGShare struct {
	Value rlwe.PolyQP
}

// MarshalBinary encodes the target element on a slice of bytes.
func (share *CKGShare) MarshalBinary() (data []byte, err error) {
	data = make([]byte, share.Value.GetDataLen(true))
	if _, err = share.Value.WriteTo(data); err != nil {
		return nil, err
	}
	return
}

// UnmarshalBinary decodes a slice of bytes on the target element.
func (share *CKGShare) UnmarshalBinary(data []byte) (err error) {
	_, err = share.Value.DecodePolyNew(data)
	return err
}

// NewCKGProtocol creates a new CKGProtocol instance
func NewCKGProtocol(params rlwe.Parameters) *CKGProtocol { // TODO drlwe.Params
	ckg := new(CKGProtocol)
	ckg.params = params
	var err error
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ckg.gaussianSamplerQ = ring.NewGaussianSampler(prng, ckg.params.RingQ(), params.Sigma(), int(6*params.Sigma()))
	return ckg
}

// AllocateShares allocates the share of the CKG protocol.
func (ckg *CKGProtocol) AllocateShares() *CKGShare {
	return &CKGShare{rlwe.PolyQP{ckg.params.RingQ().NewPoly(), ckg.params.RingP().NewPoly()}}
}

// GenShare generates the party's public key share from its secret key as:
//
// crs*s_i + e_i
//
// for the receiver protocol. Has no effect is the share was already generated.
func (ckg *CKGProtocol) GenShare(sk *rlwe.SecretKey, crs rlwe.PolyQP, shareOut *CKGShare) {
	ringQ := ckg.params.RingQ()
	ringP := ckg.params.RingP()
	ckg.gaussianSamplerQ.Read(shareOut.Value[0])
	extendBasisSmallNormAndCenter(ringQ, ringP, shareOut.Value[0], shareOut.Value[1])
	ringQ.NTT(shareOut.Value[0], shareOut.Value[0])
	ringP.NTT(shareOut.Value[1], shareOut.Value[1])
	ringQ.MulCoeffsMontgomeryAndSub(sk.Value[0], crs[0], shareOut.Value[0])
	ringP.MulCoeffsMontgomeryAndSub(sk.Value[1], crs[1], shareOut.Value[1])
}

// AggregateShares aggregates a new share to the aggregate key
func (ckg *CKGProtocol) AggregateShares(share1, share2, shareOut *CKGShare) {
	ckg.params.RingQ().Add(share1.Value[0], share2.Value[0], shareOut.Value[0])
	ckg.params.RingP().Add(share1.Value[1], share2.Value[1], shareOut.Value[1])
}

// GenPublicKey return the current aggregation of the received shares as a bfv.PublicKey.
func (ckg *CKGProtocol) GenPublicKey(roundShare *CKGShare, crs rlwe.PolyQP, pubkey *rlwe.PublicKey) {
	pubkey.Value[0].Copy(roundShare.Value)
	pubkey.Value[1].Copy(crs)
}
