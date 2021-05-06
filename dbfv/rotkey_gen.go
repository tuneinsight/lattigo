package dbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/drlwe"
)

// RTGProtocol is the structure storing the parameters for the collective rotation-keys generation.
// TODO: extract galois parameters type and remove the rotation type from the interface
type RTGProtocol struct {
	drlwe.RTGProtocol
}

// NewRotKGProtocol creates a new rotkg object and will be used to generate collective rotation-keys from a shared secret-key among j parties.
func NewRotKGProtocol(params bfv.Parameters) (rtg *RTGProtocol) {
	return &RTGProtocol{*drlwe.NewRTGProtocol(params.Parameters)}
}
