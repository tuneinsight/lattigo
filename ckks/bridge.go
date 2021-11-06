package ckks

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// UnfoldConjugateInvariantNTTLvl maps the compressed representation of Z_Q[X+X^-1]/(X^2N + 1) to full representation in Z_Q[X]/(X^2N+1)
func UnfoldConjugateInvariantNTTLvl(level int, p1, p2 *ring.Poly) {

	if 2*len(p1.Coeffs[0]) != len(p2.Coeffs[0]) {
		panic("Ring degree of p2 must be twice the ring degree of p1")
	}

	N := len(p1.Coeffs[0])

	for i := 0; i < level+1; i++ {
		tmp2, tmp1 := p2.Coeffs[i], p1.Coeffs[i]
		copy(tmp2, tmp1)
		for idx, jdx := N-1, N; jdx < 2*N; idx, jdx = idx-1, jdx+1 {
			tmp2[jdx] = tmp1[idx]
		}
	}

	return
}

// FoldConjugateInvariantNTTLvl folds [X] to [X+X^-1] in compressed form.
// Requires ringQ in ConjugateInvariantRing and p1 in (X^2N+1) and p2 in (X^N+1)
func FoldConjugateInvariantNTTLvl(level int, p1 *ring.Poly, ringQ *ring.Ring, permuteNTTIndexInv []uint64, p2 *ring.Poly) {

	if len(p1.Coeffs[0]) != 2*len(p2.Coeffs[0]) {
		panic("Ring degree of p2 must be 2N and ring degree of p1 must be N")
	}

	N := ringQ.N
	ringQ.N = len(p2.Coeffs[0])
	ringQ.PermuteNTTWithIndexLvl(level, p1, permuteNTTIndexInv, p2)
	ringQ.Add(p2, p1, p2)
	ringQ.N = N
}

func isConjugateInvariant(r *ring.Ring) bool {
	switch r.NumberTheoreticTransformer.(type) {
	case ring.NumberTheoreticTransformerConjugateInvariant:
		return true
	default:
		return false
	}
}

// SwkComplexToReal is a SwitchingKey to switch from CKKS to RCKKS.
type SwkComplexToReal struct {
	rlwe.SwitchingKey
}

// SwkRealToComplex is a Switchingkey to switch from RCKKS to CKKS
type SwkRealToComplex struct {
	rlwe.SwitchingKey
}

// GenSwitchingKeysForRingSwap generates the necessary switching keys to switch from CKKS to RCKKS and vice-versa.
func GenSwitchingKeysForRingSwap(params Parameters, skCKKS, skRCKKS *rlwe.SecretKey) (SwkComplexToReal, SwkRealToComplex) {

	if params.RingType() != ring.Standard {
		panic("ringType must be ring.Standard")
	}

	kgen := NewKeyGenerator(params)

	skRCKKSMappedToCKKS := NewSecretKey(params)
	UnfoldConjugateInvariantNTTLvl(skRCKKS.Value.Q.Level(), skRCKKS.Value.Q, skRCKKSMappedToCKKS.Value.Q)
	UnfoldConjugateInvariantNTTLvl(skRCKKS.Value.P.Level(), skRCKKS.Value.P, skRCKKSMappedToCKKS.Value.P)

	swkRtoC := kgen.GenSwitchingKey(skRCKKSMappedToCKKS, skCKKS)
	swkCtoR := kgen.GenSwitchingKey(skCKKS, skRCKKSMappedToCKKS)

	return SwkComplexToReal{*swkCtoR}, SwkRealToComplex{*swkRtoC}

}

// Bridge switches a CKKS to an RCKKS ciphertext or vice-versa.
// The bridge from CKKS to RCKKS multiplies the ciphertext scale by 2.
func (eval *evaluator) Bridge(ctIn *Ciphertext, swk interface{}, ctOut *Ciphertext) {

	ks := eval.GetKeySwitcher()

	ringQ := eval.params.RingQ()

	if isConjugateInvariant(ringQ) {
		panic("Evaluator must be instantiated with Default Ring.")
	}

	level := utils.MinInt(ctIn.Level(), ctOut.Level())

	switch swk := swk.(type) {
	case SwkComplexToReal: // Convert CKKS to RCKKS

		if len(ctIn.Value[0].Coeffs[0]) != 2*len(ctOut.Value[0].Coeffs[0]) {
			panic("ctIn ring degree must be twice ctOut ring degree")
		}

		ks.SwitchKeysInPlace(level, ctIn.Value[1], &swk.SwitchingKey, ks.Pool[1].Q, ks.Pool[2].Q)
		ringQ.Add(ks.Pool[1].Q, ctIn.Value[0], ks.Pool[1].Q)

		galEl := uint64(2*ringQ.N - 1)
		index, ok := eval.permuteNTTIndex[galEl]
		if !ok {
			if eval.permuteNTTIndex == nil {
				eval.permuteNTTIndex = make(map[uint64][]uint64)
			}
			index = ringQ.PermuteNTTIndex(galEl)
			eval.permuteNTTIndex[galEl] = index
		}
		FoldConjugateInvariantNTTLvl(level, ks.Pool[1].Q, ringQ, index, ctOut.Value[0])
		FoldConjugateInvariantNTTLvl(level, ks.Pool[2].Q, ringQ, index, ctOut.Value[1])
		ctOut.Scale = 2 * ctIn.Scale
	case SwkRealToComplex: // Convert RCKKS to CKKS

		if 2*len(ctIn.Value[0].Coeffs[0]) != len(ctOut.Value[0].Coeffs[0]) {
			panic("ctOut ring degree must be twice ctIn ring degree")
		}

		UnfoldConjugateInvariantNTTLvl(level, ctIn.Value[0], ctOut.Value[0])
		UnfoldConjugateInvariantNTTLvl(level, ctIn.Value[1], ctOut.Value[1])

		// Switches the RCKKS key [X+X^-1] to a CKKS key [X]
		ks.SwitchKeysInPlace(level, ctOut.Value[1], &swk.SwitchingKey, ks.Pool[1].Q, ks.Pool[2].Q)
		ringQ.Add(ctOut.Value[0], ks.Pool[1].Q, ctOut.Value[0])
		ring.CopyValues(ks.Pool[2].Q, ctOut.Value[1])

	default:
		panic("invalud swk, must be either SwkComplexToReal or SwkRealToComplex")
	}
}
