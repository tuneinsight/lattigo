package lwe

import (
	"github.com/tuneinsight/lattigo/v3/rlwe"
)

// LWEToRLWE transforms a set of LWE samples into their respective RLWE ciphertext such that decrypt(RLWE)[0] = decrypt(LWE)
func (h *Handler) LWEToRLWE(ctLWE map[int]*Ciphertext) (ctRLWE map[int]*rlwe.Ciphertext) {

	var level int
	for i := range ctLWE {
		level = ctLWE[i].Level()
		break
	}

	ringQ := h.paramsLUT.RingQ()
	acc := ringQ.NewPolyLvl(level)
	ctRLWE = make(map[int]*rlwe.Ciphertext)
	for i := range ctLWE {

		// Alocates ciphertext
		ctRLWE[i] = rlwe.NewCiphertextNTT(h.paramsLUT, 1, level)

		for u := 0; u < level+1; u++ {

			ctRLWE[i].Value[0].Coeffs[u][0] = ctLWE[i].Value[u][0]

			// Copy coefficients multiplied by X^{N-1} in reverse order:
			// a_{0} -a_{N-1} -a2_{N-2} ... -a_{1}
			tmp0, tmp1 := acc.Coeffs[u], ctLWE[i].Value[u][1:]
			tmp0[0] = tmp1[0]
			for k := 1; k < ringQ.N; k++ {
				tmp0[k] = ringQ.Modulus[u] - tmp1[ringQ.N-k]
			}

			copy(ctRLWE[i].Value[1].Coeffs[u], acc.Coeffs[u])
		}

		// Switches to NTT domain
		ringQ.NTTLvl(level, ctRLWE[i].Value[0], ctRLWE[i].Value[0])
		ringQ.NTTLvl(level, ctRLWE[i].Value[1], ctRLWE[i].Value[1])
	}

	return
}
