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
	GenShare(sk *rlwe.SecretKey, crs [2]*ring.Poly, shareOut *CKGShare)
	AggregateShares(share1, share2, shareOut *CKGShare)
	GenPublicKey(aggregatedShare *CKGShare, crs [2]*ring.Poly, pubkey *rlwe.PublicKey)
}

// CKGProtocol is the structure storing the parameters and and precomputations for the collective key generation protocol.
type CKGProtocol struct {
	params rlwe.Parameters

	ringQ            *ring.Ring
	ringP            *ring.Ring
	gaussianSamplerQ *ring.GaussianSampler
}

// CKGShare is a struct storing the CKG protocol's share.
type CKGShare struct {
	Value [2]*ring.Poly
}

// MarshalBinary encodes the target element on a slice of bytes.
func (share *CKGShare) MarshalBinary() (data []byte, err error) {
	data = make([]byte, share.Value[0].GetDataLen(true)+share.Value[1].GetDataLen(true))
	var inc int
	ptr := 0
	if inc, err = share.Value[0].WriteTo(data[ptr:]); err != nil {
		return []byte{}, err
	}
	ptr += inc

	if _, err = share.Value[1].WriteTo(data[ptr:]); err != nil {
		return []byte{}, err
	}
	return
}

// UnmarshalBinary decodes a slice of bytes on the target element.
func (share *CKGShare) UnmarshalBinary(data []byte) (err error) {

	var pt, inc int
	share.Value[0] = new(ring.Poly)
	if inc, err = share.Value[0].DecodePolyNew(data[pt:]); err != nil {
		return err
	}

	pt += inc

	share.Value[1] = new(ring.Poly)
	if _, err = share.Value[1].DecodePolyNew(data[pt:]); err != nil {
		return err
	}
	return err
}

// NewCKGProtocol creates a new CKGProtocol instance
func NewCKGProtocol(params rlwe.Parameters) *CKGProtocol { // TODO drlwe.Params

	ckg := new(CKGProtocol)
	ckg.params = params
	ckg.ringQ = params.RingQ()
	ckg.ringP = params.RingP()

	var err error
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ckg.gaussianSamplerQ = ring.NewGaussianSampler(prng, ckg.ringQ, params.Sigma(), int(6*params.Sigma()))
	return ckg
}

// AllocateShares allocates the share of the CKG protocol.
func (ckg *CKGProtocol) AllocateShares() *CKGShare {
	return &CKGShare{[2]*ring.Poly{ckg.ringQ.NewPoly(), ckg.ringP.NewPoly()}}
}

// GenShare generates the party's public key share from its secret key as:
//
// crs*s_i + e_i
//
// for the receiver protocol. Has no effect is the share was already generated.
func (ckg *CKGProtocol) GenShare(sk *rlwe.SecretKey, crs [2]*ring.Poly, shareOut *CKGShare) {
	ringQ := ckg.ringQ
	ringP := ckg.ringP
	ckg.gaussianSamplerQ.Read(shareOut.Value[0])
	extendBasisSmallNormAndCenter(ringQ, ringP, shareOut.Value[0], shareOut.Value[1])
	ringQ.NTT(shareOut.Value[0], shareOut.Value[0])
	ringP.NTT(shareOut.Value[1], shareOut.Value[1])
	ringQ.MulCoeffsMontgomeryAndSub(sk.Value[0], crs[0], shareOut.Value[0])
	ringP.MulCoeffsMontgomeryAndSub(sk.Value[1], crs[1], shareOut.Value[1])
}

// AggregateShares aggregates a new share to the aggregate key
func (ckg *CKGProtocol) AggregateShares(share1, share2, shareOut *CKGShare) {
	ckg.ringQ.Add(share1.Value[0], share2.Value[0], shareOut.Value[0])
	ckg.ringP.Add(share1.Value[1], share2.Value[1], shareOut.Value[1])
}

// GenPublicKey return the current aggregation of the received shares as a bfv.PublicKey.
func (ckg *CKGProtocol) GenPublicKey(roundShare *CKGShare, crs [2]*ring.Poly, pubkey *rlwe.PublicKey) {
	pubkey.Value[0][0].Copy(roundShare.Value[0])
	pubkey.Value[0][1].Copy(roundShare.Value[1])
	pubkey.Value[1][0].Copy(crs[0])
	pubkey.Value[1][1].Copy(crs[1])
}
