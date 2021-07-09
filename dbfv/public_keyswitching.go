package dbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/drlwe"
)

// PCKSProtocol is the structure storing the parameters for the collective public key-switching.
type PCKSProtocol struct {
	drlwe.PCKSProtocol

	params bfv.Parameters
}

// NewPCKSProtocol creates a new PCKSProtocol object and will be used to re-encrypt a ciphertext ctx encrypted under a secret-shared key among j parties under a new
// collective public-key.
func NewPCKSProtocol(params bfv.Parameters, sigmaSmudging float64) *PCKSProtocol {
	return &PCKSProtocol{*drlwe.NewPCKSProtocol(params.Parameters, sigmaSmudging), params}
}

// AllocateShareBFV allocates the shares of one party in the PCKS protocol for BFV.
func (pcks *PCKSProtocol) AllocateShareBFV() *drlwe.PCKSShare {
	return pcks.AllocateShare(pcks.params.MaxLevel())
}
