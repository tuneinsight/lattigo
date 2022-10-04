package dbfv

import (
	"github.com/tuneinsight/lattigo/v4/bfv"
	"github.com/tuneinsight/lattigo/v4/drlwe"
)

// NewCKGProtocol creates a new drlwe.CKGProtocol instance from the BFV parameters
func NewCKGProtocol(params bfv.Parameters) *drlwe.CKGProtocol {
	return drlwe.NewCKGProtocol(params.Parameters)
}

// NewRKGProtocol creates a new drlwe.RKGProtocol instance from the BFV parameters
func NewRKGProtocol(params bfv.Parameters) *drlwe.RKGProtocol {
	return drlwe.NewRKGProtocol(params.Parameters)
}

// NewRTGProtocol creates a new drlwe.RTGProtocol instance from the BFV parameters
func NewRTGProtocol(params bfv.Parameters) *drlwe.RTGProtocol {
	return drlwe.NewRTGProtocol(params.Parameters)
}
