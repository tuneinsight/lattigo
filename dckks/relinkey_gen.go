package dckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/drlwe"
)

// RKGProtocol is the structure storing the parameters and state for a party in the collective relinearization key
// generation protocol.
type RKGProtocol struct {
	drlwe.RKGProtocol
}

// NewRKGProtocol creates a new RKGProtocol object that will be used to generate a collective evaluation-key
// among j parties in the given context with the given bit-decomposition.
func NewRKGProtocol(params *ckks.Parameters) *RKGProtocol {
	return &RKGProtocol{*drlwe.NewRKGProtocol(params.N(), params.Qi(), params.Pi(), 0.5, params.Sigma())}
}

// GenCKKSRelinearizationKey finalizes the protocol and returns the common EvaluationKey.
func (ekg *RKGProtocol) GenCKKSRelinearizationKey(round1 *drlwe.RKGShare, round2 *drlwe.RKGShare, evalKeyOut *ckks.RelinearizationKey) {
	ekg.GenRelinearizationKey(round1, round2, &evalKeyOut.RelinearizationKey)
}
