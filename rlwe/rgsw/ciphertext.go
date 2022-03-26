package rgsw

import (
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/rlwe/gadget"
)

// Ciphertext is a generic type for RGSW ciphertext.
type Ciphertext struct {
	Value [2]gadget.Ciphertext
}

// LevelQ returns the level of the modulus Q of the target.
func (ct *Ciphertext) LevelQ() int {
	return ct.Value[0].LevelQ()
}

// LevelP returns the level of the modulus P of the target.
func (ct *Ciphertext) LevelP() int {
	return ct.Value[0].LevelP()
}

// NewCiphertextNTT allocates a new RGSW ciphertext in the NTT domain.
func NewCiphertextNTT(params rlwe.Parameters, levelQ int) (ct *Ciphertext) {
	levelP := params.PCount() - 1
	return &Ciphertext{
		Value: [2]gadget.Ciphertext{
			*gadget.NewCiphertextNTT(levelQ, levelP, params.DecompRNS(levelQ, levelP), params.DecompBIT(levelQ, levelP), *params.RingQP()),
			*gadget.NewCiphertextNTT(levelQ, levelP, params.DecompRNS(levelQ, levelP), params.DecompBIT(levelQ, levelP), *params.RingQP()),
		},
	}
}
