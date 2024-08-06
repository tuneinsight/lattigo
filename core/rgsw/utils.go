package rgsw

import (
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

// NoiseRGSWCiphertext returns the log2 of the standard deviation of the noise of each component of the RGSW ciphertext.
// pt must be in the NTT and Montgomery domain
func NoiseRGSWCiphertext(ct *Ciphertext, pt ring.Poly, sk *rlwe.SecretKey, params rlwe.Parameters) (float64, float64) {
	ptsk := *pt.CopyNew()
	params.RingQ().AtLevel(ct.LevelQ()).MulCoeffsMontgomery(ptsk, sk.Value.Q, ptsk)
	return rlwe.NoiseGadgetCiphertext(&ct.Value[0], pt, sk, params), rlwe.NoiseGadgetCiphertext(&ct.Value[1], ptsk, sk, params)
}
