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
	ckg := new(CKGProtocol)
	ckg.CKGProtocol = *drlwe.NewCKGProtocol(params)
	return ckg
}
