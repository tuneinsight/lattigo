package rgsw

import (
	"github.com/tuneinsight/lattigo/v3/rlwe/gadget"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
)

// Plaintext stores an RGSW plaintext value.
type Plaintext gadget.Plaintext

// NewPlaintext creates a new RGSW plaintext from value, which can be either uint64, int64 or *ring.Poly.
// Plaintext is returned in the NTT and Mongtomery domain.
func NewPlaintext(value interface{}, levelQ, levelP, logBase2, decompBIT int, ringQP ringqp.Ring) (pt *Plaintext) {
	return &Plaintext{Value: gadget.NewPlaintext(value, levelQ, levelP, logBase2, decompBIT, ringQP).Value}
}
