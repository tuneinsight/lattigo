package dbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/drlwe"
)

// RKGProtocol is the structure storing the parameters and state for a party in the collective relinearization key
// generation protocol.
type RKGProtocol struct {
	drlwe.RKGProtocol
}

// NewRKGProtocol creates a new RKGProtocol object that will be used to generate a collective evaluation-key
// among j parties in the given context with the given bit-decomposition.
func NewRKGProtocol(params bfv.Parameters) *RKGProtocol {
	return &RKGProtocol{*drlwe.NewRKGProtocol(params, 0.5)}
}
