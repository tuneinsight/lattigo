package dckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/drlwe"
)

// CKSProtocol is a structure storing the parameters for the collective key-switching protocol.
type CKSProtocol struct {
	drlwe.CKSProtocol
}

// NewCKSProtocol creates a new CKSProtocol that will be used to operate a collective key-switching on a ciphertext encrypted under a collective public-key, whose
// secret-shares are distributed among j parties, re-encrypting the ciphertext under another public-key, whose secret-shares are also known to the
// parties.
func NewCKSProtocol(params ckks.Parameters, sigmaSmudging float64) (cks *CKSProtocol) {
	return &CKSProtocol{*drlwe.NewCKSProtocol(params.Parameters, sigmaSmudging)}
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (cks *CKSProtocol) KeySwitch(ctIn *ckks.Ciphertext, combined *drlwe.CKSShare, ctOut *ckks.Ciphertext) {
	cks.CKSProtocol.KeySwitch(ctIn.Ciphertext, combined, ctOut.Ciphertext)
	ctOut.Scale = ctIn.Scale
}

// ShallowCopy creates a shallow copy of CKSProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// CKSProtocol can be used concurrently.
func (cks *CKSProtocol) ShallowCopy() *CKSProtocol {
	return &CKSProtocol{*cks.CKSProtocol.ShallowCopy()}
}

// PCKSProtocol is the structure storing the parameters for the collective public key-switching.
type PCKSProtocol struct {
	drlwe.PCKSProtocol
}

// NewPCKSProtocol creates a new PCKSProtocol object and will be used to re-encrypt a ciphertext ctx encrypted under a secret-shared key mong j parties under a new
// collective public-key.
func NewPCKSProtocol(params ckks.Parameters, sigmaSmudging float64) *PCKSProtocol {
	return &PCKSProtocol{*drlwe.NewPCKSProtocol(params.Parameters, sigmaSmudging)}
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut.
func (pcks *PCKSProtocol) KeySwitch(ctIn *ckks.Ciphertext, combined *drlwe.PCKSShare, ctOut *ckks.Ciphertext) {
	pcks.PCKSProtocol.KeySwitch(ctIn.Ciphertext, combined, ctOut.Ciphertext)
	ctOut.Scale = ctIn.Scale
}

// ShallowCopy creates a shallow copy of PCKSProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// PCKSProtocol can be used concurrently.
func (pcks *PCKSProtocol) ShallowCopy() *PCKSProtocol {
	return &PCKSProtocol{*pcks.PCKSProtocol.ShallowCopy()}
}
