package dbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// CKSProtocol is a structure storing the parameters for the collective key-switching protocol.
type CKSProtocol struct {
	drlwe.CKSProtocol
}

// NewCKSProtocol creates a new CKSProtocol that will be used to operate a collective key-switching on a ciphertext encrypted under a collective public-key, whose
// secret-shares are distributed among j parties, re-encrypting the ciphertext under another public-key, whose secret-shares are also known to the
// parties.
func NewCKSProtocol(params bfv.Parameters, sigmaSmudging float64) *CKSProtocol {
	return &CKSProtocol{*drlwe.NewCKSProtocol(params.Parameters, sigmaSmudging)}
}

// GenShare is the first and unique round of the CKSProtocol protocol. Each party holding a ciphertext ctx encrypted under a collective publick-key musth
// compute the following :
//
// [(skInput_i - skOutput_i) * ctx[0] + e_i]
//
// Each party then broadcast the result of this computation to the other j-1 parties.
func (cks *CKSProtocol) GenShare(skInput, skOutput *rlwe.SecretKey, ct *bfv.Ciphertext, shareOut drlwe.CKSShare) {
	cks.CKSProtocol.GenShare(skInput, skOutput, ct, shareOut)
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (cks *CKSProtocol) KeySwitch(combined drlwe.CKSShare, ct *bfv.Ciphertext, ctOut *bfv.Ciphertext) {
	cks.CKSProtocol.KeySwitch(combined, ct, ctOut)
}
