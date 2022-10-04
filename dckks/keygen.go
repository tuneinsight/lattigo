package dckks

import (
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/drlwe"
)

// NewCKGProtocol creates a new drlwe.CKGProtocol instance from the CKKS parameters
func NewCKGProtocol(params ckks.Parameters) *drlwe.CKGProtocol {
	return drlwe.NewCKGProtocol(params.Parameters)
}

// NewRKGProtocol creates a new drlwe.RKGProtocol instance from the CKKS parameters
func NewRKGProtocol(params ckks.Parameters) *drlwe.RKGProtocol {
	return drlwe.NewRKGProtocol(params.Parameters)
}

// NewRTGProtocol creates a new drlwe.RTGProtocol instance from the CKKS parameters
func NewRTGProtocol(params ckks.Parameters) *drlwe.RTGProtocol {
	return drlwe.NewRTGProtocol(params.Parameters)
}
