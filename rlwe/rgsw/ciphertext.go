// Package rgsw implements R-LWE based RGSW ciphertexts. An RGSW ciphertext is a tuple of gadget ciphertexts (see package gadget), where
// the first element is a gadget ciphertext encrypting the message and the second element a gadget ciphertext encryping the message times
// the secret.
package rgsw

import (
	"github.com/tuneinsight/lattigo/v3/rlwe/gadget"
	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
)

// Ciphertext is a generic type for RGSW ciphertext.
type Ciphertext struct {
	Value [2]gadget.GadgetCiphertext
}

// LevelQ returns the level of the modulus Q of the target.
func (ct *Ciphertext) LevelQ() int {
	return ct.Value[0].LevelQ()
}

// LevelP returns the level of the modulus P of the target.
func (ct *Ciphertext) LevelP() int {
	return ct.Value[0].LevelP()
}

// NewCiphertext allocates a new RGSW ciphertext in the NTT domain.
func NewCiphertext(levelQ, levelP, decompRNS, decompBit int, ringQP ringqp.Ring) (ct *Ciphertext) {
	return &Ciphertext{
		Value: [2]gadget.GadgetCiphertext{
			*gadget.NewGadgetCiphertext(levelQ, levelP, decompRNS, decompBit, ringQP),
			*gadget.NewGadgetCiphertext(levelQ, levelP, decompRNS, decompBit, ringQP),
		},
	}
}
