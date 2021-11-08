package ckks

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// SchemeSwitcher is a type for switching between the CKKS and RCKKS schemes.
type SchemeSwitcher struct {
	rlwe.KeySwitcher

	stdRingQ, conjugateRingQ *ring.Ring

	*rlwe.SwkComplexToReal
	*rlwe.SwkRealToComplex

	permuteNTTIndex []uint64
}

// NewSchemeSwitcher instantiate a new SchemeSwitcher type. It may be instantiated from parameters from either RingType.
// The method returns an error if the parameters cannot support the switching (e.g., the NTT transforms are undefined for
// either of the two ring type).
func NewSchemeSwitcher(params Parameters, comlexToRealSwk *rlwe.SwkComplexToReal, RealToComplexSwk *rlwe.SwkRealToComplex) (SchemeSwitcher, error) {

	s := SchemeSwitcher{
		KeySwitcher:      *rlwe.NewKeySwitcher(params.Parameters),
		SwkComplexToReal: comlexToRealSwk,
		SwkRealToComplex: RealToComplexSwk,
	}
	var err error
	if s.stdRingQ, err = params.RingQ().StandardRing(); err != nil {
		return SchemeSwitcher{}, err
	}
	if s.conjugateRingQ, err = params.RingQ().ConjugateInvariantRing(); err != nil {
		return SchemeSwitcher{}, err // TODO: better error handline in case no ConjugateInvariant ring exist
	}
	s.permuteNTTIndex = s.stdRingQ.PermuteNTTIndex((uint64(s.stdRingQ.N) << 1) - 1)
	return s, nil
}

// ComplexToReal switches the provided ciphertext `ctIn` from the  CKKS scheme to the RCKKS scheme
// and write the result into `ctOut`.
// The method panics if the SchemeSwitcher was not initialized with a SwkComplexToReal key.
func (switcher *SchemeSwitcher) ComplexToReal(ctIn, ctOut *Ciphertext) {

	level := utils.MinInt(ctIn.Level(), ctOut.Level())

	if len(ctIn.Value[0].Coeffs[0]) != 2*len(ctOut.Value[0].Coeffs[0]) {
		panic("ctIn ring degree must be twice ctOut ring degree")
	}

	if switcher.SwkComplexToReal == nil {
		panic("no SwkComplexToReal provided to this SchemeSwitcher")
	}

	switcher.SwitchKeysInPlace(level, ctIn.Value[1], &switcher.SwkComplexToReal.SwitchingKey, switcher.Pool[1].Q, switcher.Pool[2].Q)
	switcher.stdRingQ.Add(switcher.Pool[1].Q, ctIn.Value[0], switcher.Pool[1].Q)

	switcher.conjugateRingQ.FoldStandardToConjugateInvariant(level, switcher.Pool[1].Q, switcher.permuteNTTIndex, ctOut.Value[0])
	switcher.conjugateRingQ.FoldStandardToConjugateInvariant(level, switcher.Pool[2].Q, switcher.permuteNTTIndex, ctOut.Value[1])
	ctOut.Scale = 2 * ctIn.Scale
}

// RealToComplex switches the provided ciphertext `ctIn` from the  RCKKS scheme to the CKKS scheme.
// and write the result into `ctOut`.
// The method panics if the SchemeSwitcher was not initialized with a SwkRealToComplex key.
func (switcher *SchemeSwitcher) RealToComplex(ctIn, ctOut *Ciphertext) {

	level := utils.MinInt(ctIn.Level(), ctOut.Level())

	if 2*len(ctIn.Value[0].Coeffs[0]) != len(ctOut.Value[0].Coeffs[0]) {
		panic("ctOut ring degree must be twice ctIn ring degree")
	}

	if switcher.SwkRealToComplex == nil {
		panic("no SwkRealToComplex provided to this SchemeSwitcher")
	}

	switcher.stdRingQ.UnfoldConjugateInvariantToStandard(level, ctIn.Value[0], ctOut.Value[0])
	switcher.stdRingQ.UnfoldConjugateInvariantToStandard(level, ctIn.Value[1], ctOut.Value[1])

	// Switches the RCKswitcher key [X+X^-1] to a CKswitcher key [X]
	switcher.SwitchKeysInPlace(level, ctOut.Value[1], &switcher.SwkRealToComplex.SwitchingKey, switcher.Pool[1].Q, switcher.Pool[2].Q)
	switcher.stdRingQ.Add(ctOut.Value[0], switcher.Pool[1].Q, ctOut.Value[0])
	ring.CopyValues(switcher.Pool[2].Q, ctOut.Value[1])
}
