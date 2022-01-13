//Package dbfv implements a distributed (or threshold) version of the BFV scheme that enables secure multiparty computation solutions with secret-shared secret keys.
package dbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/drlwe"
)

// CKGProtocol is the structure storing the parameters and state for a party in the collective key generation protocol.
type CKGProtocol struct {
	drlwe.CKGProtocol
}

// NewCKGProtocol creates a new CKGProtocol instance
func NewCKGProtocol(params bfv.Parameters) *CKGProtocol {
	return &CKGProtocol{*drlwe.NewCKGProtocol(params.Parameters)}
}

// ShallowCopy creates a shallow copy of CKGProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// CKGProtocol can be used concurrently.
func (ckg *CKGProtocol) ShallowCopy() *CKGProtocol {
	return &CKGProtocol{*ckg.CKGProtocol.ShallowCopy()}
}

// RKGProtocol is the structure storing the parameters and state for a party in the collective relinearization key
// generation protocol.
type RKGProtocol struct {
	drlwe.RKGProtocol
}

// NewRKGProtocol creates a new RKGProtocol object that will be used to generate a collective evaluation-key
// among j parties in the given context with the given bit-decomposition.
func NewRKGProtocol(params bfv.Parameters) *RKGProtocol {
	return &RKGProtocol{*drlwe.NewRKGProtocol(params.Parameters, 0.5)}
}

// ShallowCopy creates a shallow copy of RKGProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// RKGProtocol can be used concurrently.
func (rkg *RKGProtocol) ShallowCopy() *RKGProtocol {
	return &RKGProtocol{*rkg.RKGProtocol.ShallowCopy()}
}

// RTGProtocol is the structure storing the parameters for the collective rotation-keys generation.
type RTGProtocol struct {
	drlwe.RTGProtocol
}

// NewRotKGProtocol creates a new rotkg object and will be used to generate collective rotation-keys from a shared secret-key among j parties.
func NewRotKGProtocol(params bfv.Parameters) (rtg *RTGProtocol) {
	return &RTGProtocol{*drlwe.NewRTGProtocol(params.Parameters)}
}

// ShallowCopy creates a shallow copy of RTGProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// RTGProtocol can be used concurrently.
func (rtg *RTGProtocol) ShallowCopy() *RTGProtocol {
	return &RTGProtocol{*rtg.RTGProtocol.ShallowCopy()}
}
