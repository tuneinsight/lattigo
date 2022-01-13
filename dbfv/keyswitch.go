package dbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/drlwe"
)

// CKSProtocol is a structure storing the parameters for the collective key-switching protocol.
type CKSProtocol struct {
	drlwe.CKSProtocol
	maxLevel int
}

// NewCKSProtocol creates a new CKSProtocol that will be used to operate a collective key-switching on a ciphertext encrypted under a collective public-key, whose
// secret-shares are distributed among j parties, re-encrypting the ciphertext under another public-key, whose secret-shares are also known to the
// parties.
func NewCKSProtocol(params bfv.Parameters, sigmaSmudging float64) *CKSProtocol {
	return &CKSProtocol{*drlwe.NewCKSProtocol(params.Parameters, sigmaSmudging), params.MaxLevel()}
}

// AllocateShareBFV allocates the shares of one party in the CKS protocol for BFV.
func (cks *CKSProtocol) AllocateShareBFV() *drlwe.CKSShare {
	return cks.AllocateShare(cks.maxLevel)
}

// ShallowCopy creates a shallow copy of CKSProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// CKSProtocol can be used concurrently.
func (cks *CKSProtocol) ShallowCopy() *CKSProtocol {
	return &CKSProtocol{*cks.CKSProtocol.ShallowCopy(), cks.maxLevel}
}

// PCKSProtocol is the structure storing the parameters for the collective public key-switching.
type PCKSProtocol struct {
	drlwe.PCKSProtocol
	maxLevel int
}

// NewPCKSProtocol creates a new PCKSProtocol object and will be used to re-encrypt a ciphertext ctx encrypted under a secret-shared key among j parties under a new
// collective public-key.
func NewPCKSProtocol(params bfv.Parameters, sigmaSmudging float64) *PCKSProtocol {
	return &PCKSProtocol{*drlwe.NewPCKSProtocol(params.Parameters, sigmaSmudging), params.MaxLevel()}
}

// AllocateShareBFV allocates the shares of one party in the PCKS protocol for BFV.
func (pcks *PCKSProtocol) AllocateShareBFV() *drlwe.PCKSShare {
	return pcks.AllocateShare(pcks.maxLevel)
}

// ShallowCopy creates a shallow copy of PCKSProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// PCKSProtocol can be used concurrently.
func (pcks *PCKSProtocol) ShallowCopy() *PCKSProtocol {
	return &PCKSProtocol{*pcks.PCKSProtocol.ShallowCopy(), pcks.maxLevel}
}
