package ckks

import (
<<<<<<< HEAD
=======
	"fmt"

>>>>>>> dev_rckks
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

<<<<<<< HEAD
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

	if params.RingType() != rlwe.RingStandard {
		panic("ringType must be rlwe.RingStandard")
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
=======
// DomainSwitcher is a type for switching between the standard CKKS domain (which encrypts vectors of complex numbers)
// and the conjugate invariant variant of CKKS (which encrypts vectors of real numbers).
type DomainSwitcher struct {
	rlwe.KeySwitcher

	stdRingQ, conjugateRingQ *ring.Ring

	*SwkComplexToReal
	*SwkRealToComplex

	permuteNTTIndex []uint64
}

// NewDomainSwitcher instantiate a new DomainSwitcher type. It may be instantiated from parameters from either RingType.
// The method returns an error if the parameters cannot support the switching (e.g., the NTT transforms are undefined for
// either of the two ring type).
func NewDomainSwitcher(params Parameters, comlexToRealSwk *SwkComplexToReal, RealToComplexSwk *SwkRealToComplex) (DomainSwitcher, error) {

	s := DomainSwitcher{
		SwkComplexToReal: comlexToRealSwk,
		SwkRealToComplex: RealToComplexSwk,
	}
	var err error
	if s.stdRingQ, err = params.RingQ().StandardRing(); err != nil {
		return DomainSwitcher{}, fmt.Errorf("cannot switch between schemes because the standard NTT is undefined for params: %f", err)
	}
	if s.conjugateRingQ, err = params.RingQ().ConjugateInvariantRing(); err != nil {
		return DomainSwitcher{}, fmt.Errorf("cannot switch between schemes because the standard NTT is undefined for params: %f", err)
	}

	stdParams, err := params.StandardParameters()
	s.KeySwitcher = *rlwe.NewKeySwitcher(stdParams.Parameters)

	s.permuteNTTIndex = s.stdRingQ.PermuteNTTIndex((uint64(s.stdRingQ.N) << 1) - 1)
	return s, nil
}

// ComplexToReal switches the provided ciphertext `ctIn` from the standard domain to the conjugate
// invariant domain and write the result into `ctOut`.
// Given ctInCKKS = enc(real(m) + imag(m)) in Z[X](X^N + 1), returns ctOutCI = enc(real(m))
// in Z[X+X^-1]/(X^N + 1) in compressed form (N/2 coefficients).
// The scale of the output ciphertext is 2 times the scale of the input one.
// Requires the ring degree of ctOut to be half the ring degree of ctIn.
// The security is changed from Z[X]/(X^N+1) to Z[X]/(X^N/2+1).
// The method panics if the DomainSwitcher was not initialized with a SwkComplexToReal key.
func (switcher *DomainSwitcher) ComplexToReal(ctIn, ctOut *Ciphertext) {

	level := utils.MinInt(ctIn.Level(), ctOut.Level())

	if len(ctIn.Value[0].Coeffs[0]) != 2*len(ctOut.Value[0].Coeffs[0]) {
		panic("ctIn ring degree must be twice ctOut ring degree")
	}

	if switcher.SwkComplexToReal == nil {
		panic("no SwkComplexToReal provided to this DomainSwitcher")
	}

	switcher.SwitchKeysInPlace(level, ctIn.Value[1], &switcher.SwkComplexToReal.SwitchingKey, switcher.Pool[1].Q, switcher.Pool[2].Q)
	switcher.stdRingQ.Add(switcher.Pool[1].Q, ctIn.Value[0], switcher.Pool[1].Q)

	switcher.conjugateRingQ.FoldStandardToConjugateInvariant(level, switcher.Pool[1].Q, switcher.permuteNTTIndex, ctOut.Value[0])
	switcher.conjugateRingQ.FoldStandardToConjugateInvariant(level, switcher.Pool[2].Q, switcher.permuteNTTIndex, ctOut.Value[1])
	ctOut.Scale = 2 * ctIn.Scale
}

// RealToComplex switches the provided ciphertext `ctIn` from the conjugate invariant domain to the .
// standard domain and write the result into `ctOut`.
// Given ctInCI = enc(real(m)) in Z[X+X^-1]/(X^2N+1) in compressed form (N coefficients), returns
// ctOutCKKS = enc(real(m) + imag(0)) in Z[X]/(X^2N+1).
// Requires the ring degree of ctOut to be twice the ring degree of ctIn.
// The security is changed from Z[X]/(X^N+1) to Z[X]/(X^2N+1).
// The method panics if the DomainSwitcher was not initialized with a SwkRealToComplex key.
func (switcher *DomainSwitcher) RealToComplex(ctIn, ctOut *Ciphertext) {

	level := utils.MinInt(ctIn.Level(), ctOut.Level())

	if 2*len(ctIn.Value[0].Coeffs[0]) != len(ctOut.Value[0].Coeffs[0]) {
		panic("ctOut ring degree must be twice ctIn ring degree")
	}

	if switcher.SwkRealToComplex == nil {
		panic("no SwkRealToComplex provided to this DomainSwitcher")
	}

	switcher.stdRingQ.UnfoldConjugateInvariantToStandard(level, ctIn.Value[0], ctOut.Value[0])
	switcher.stdRingQ.UnfoldConjugateInvariantToStandard(level, ctIn.Value[1], ctOut.Value[1])

	// Switches the RCKswitcher key [X+X^-1] to a CKswitcher key [X]
	switcher.SwitchKeysInPlace(level, ctOut.Value[1], &switcher.SwkRealToComplex.SwitchingKey, switcher.Pool[1].Q, switcher.Pool[2].Q)
	switcher.stdRingQ.Add(ctOut.Value[0], switcher.Pool[1].Q, ctOut.Value[0])
	ring.CopyValues(switcher.Pool[2].Q, ctOut.Value[1])
>>>>>>> dev_rckks
}
