package ckks

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// DomainSwitcher is a type for switching between the standard CKKS domain (which encrypts vectors of complex numbers)
// and the conjugate invariant variant of CKKS (which encrypts vectors of real numbers).
type DomainSwitcher struct {
	stdRingQ, conjugateRingQ *ring.Ring

	stdToci, ciToStd  *rlwe.EvaluationKey
	automorphismIndex []uint64
}

// NewDomainSwitcher instantiate a new DomainSwitcher type. It may be instantiated from parameters from either RingType.
// The method returns an error if the parameters cannot support the switching (e.g., the NTTs are undefined for
// either of the two ring types).
// The comlexToRealEvk and comlexToRealEvk EvaluationKeys can be generated using the rlwe.KeyGenerator.GenEvaluationKeysForRingSwap(*).
func NewDomainSwitcher(params Parameters, comlexToRealEvk, realToComplexEvk *rlwe.EvaluationKey) (DomainSwitcher, error) {

	s := DomainSwitcher{
		stdToci: comlexToRealEvk,
		ciToStd: realToComplexEvk,
	}
	var err error
	if s.stdRingQ, err = params.RingQ().StandardRing(); err != nil {
		return DomainSwitcher{}, fmt.Errorf("cannot switch between schemes because the standard NTT is undefined for params: %s", err)
	}
	if s.conjugateRingQ, err = params.RingQ().ConjugateInvariantRing(); err != nil {
		return DomainSwitcher{}, fmt.Errorf("cannot switch between schemes because the standard NTT is undefined for params: %s", err)
	}

	s.automorphismIndex = ring.AutomorphismNTTIndex(s.stdRingQ.N(), s.stdRingQ.NthRoot(), s.stdRingQ.NthRoot()-1)

	return s, nil
}

// ComplexToReal switches the provided ciphertext `ctIn` from the standard domain to the conjugate
// invariant domain and writes the result into `ctOut`.
// Given ctInCKKS = enc(real(m) + imag(m)) in Z[X](X^N + 1), returns ctOutCI = enc(real(m))
// in Z[X+X^-1]/(X^N + 1) in compressed form (N/2 coefficients).
// The scale of the output ciphertext is twice the scale of the input one.
// Requires the ring degree of ctOut to be half the ring degree of ctIn.
// The security is changed from Z[X]/(X^N+1) to Z[X]/(X^N/2+1).
// The method panics if the DomainSwitcher was not initialized with a the appropriate EvaluationKeys.
func (switcher *DomainSwitcher) ComplexToReal(eval Evaluator, ctIn, ctOut *rlwe.Ciphertext) {

	evalRLWE := eval.GetRLWEEvaluator()

	if evalRLWE.Parameters().RingType() != ring.Standard {
		panic("cannot ComplexToReal: provided evaluator is not instantiated with RingType ring.Standard")
	}

	level := utils.MinInt(ctIn.Level(), ctOut.Level())

	if len(ctIn.Value[0].Coeffs[0]) != 2*len(ctOut.Value[0].Coeffs[0]) {
		panic("cannot ComplexToReal: ctIn ring degree must be twice ctOut ring degree")
	}

	ctOut.Resize(1, level)

	if switcher.stdToci == nil {
		panic("cannot ComplexToReal: no realToComplexEvk provided to this DomainSwitcher")
	}

	ctTmp := &rlwe.Ciphertext{Value: []*ring.Poly{evalRLWE.BuffQP[1].Q, evalRLWE.BuffQP[2].Q}}
	ctTmp.MetaData = ctIn.MetaData

	evalRLWE.GadgetProduct(level, ctIn.Value[1], switcher.stdToci.GadgetCiphertext, ctTmp)
	switcher.stdRingQ.AtLevel(level).Add(evalRLWE.BuffQP[1].Q, ctIn.Value[0], evalRLWE.BuffQP[1].Q)

	switcher.conjugateRingQ.AtLevel(level).FoldStandardToConjugateInvariant(evalRLWE.BuffQP[1].Q, switcher.automorphismIndex, ctOut.Value[0])
	switcher.conjugateRingQ.AtLevel(level).FoldStandardToConjugateInvariant(evalRLWE.BuffQP[2].Q, switcher.automorphismIndex, ctOut.Value[1])
	ctOut.MetaData = ctIn.MetaData
	ctOut.Scale = ctIn.Scale.Mul(rlwe.NewScale(2))
}

// RealToComplex switches the provided ciphertext `ctIn` from the conjugate invariant domain to the
// standard domain and writes the result into `ctOut`.
// Given ctInCI = enc(real(m)) in Z[X+X^-1]/(X^2N+1) in compressed form (N coefficients), returns
// ctOutCKKS = enc(real(m) + imag(0)) in Z[X]/(X^2N+1).
// Requires the ring degree of ctOut to be twice the ring degree of ctIn.
// The security is changed from Z[X]/(X^N+1) to Z[X]/(X^2N+1).
// The method panics if the DomainSwitcher was not initialized with a the appropriate EvaluationKeys.
func (switcher *DomainSwitcher) RealToComplex(eval Evaluator, ctIn, ctOut *rlwe.Ciphertext) {

	evalRLWE := eval.GetRLWEEvaluator()

	if evalRLWE.Parameters().RingType() != ring.Standard {
		panic("cannot RealToComplex: provided evaluator is not instantiated with RingType ring.Standard")
	}

	level := utils.MinInt(ctIn.Level(), ctOut.Level())

	if 2*len(ctIn.Value[0].Coeffs[0]) != len(ctOut.Value[0].Coeffs[0]) {
		panic("cannot RealToComplex: ctOut ring degree must be twice ctIn ring degree")
	}

	ctOut.Resize(1, level)

	if switcher.ciToStd == nil {
		panic("cannot RealToComplex: no realToComplexEvk provided to this DomainSwitcher")
	}

	switcher.stdRingQ.AtLevel(level).UnfoldConjugateInvariantToStandard(ctIn.Value[0], ctOut.Value[0])
	switcher.stdRingQ.AtLevel(level).UnfoldConjugateInvariantToStandard(ctIn.Value[1], ctOut.Value[1])

	ctTmp := &rlwe.Ciphertext{Value: []*ring.Poly{evalRLWE.BuffQP[1].Q, evalRLWE.BuffQP[2].Q}}
	ctTmp.MetaData = ctIn.MetaData

	// Switches the RCKswitcher key [X+X^-1] to a CKswitcher key [X]
	evalRLWE.GadgetProduct(level, ctOut.Value[1], switcher.ciToStd.GadgetCiphertext, ctTmp)
	switcher.stdRingQ.AtLevel(level).Add(ctOut.Value[0], evalRLWE.BuffQP[1].Q, ctOut.Value[0])
	ring.CopyLvl(level, evalRLWE.BuffQP[2].Q, ctOut.Value[1])
	ctOut.MetaData = ctIn.MetaData
}
