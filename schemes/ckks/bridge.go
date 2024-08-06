package ckks

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils"
)

// DomainSwitcher is a type for switching between the standard CKKS domain (which encrypts vectors of complex numbers)
// and the conjugate invariant variant of CKKS (which encrypts vectors of real numbers).
type DomainSwitcher struct {
	stdRingQ, conjugateRingQ *ring.Ring

	stdToci, ciToStd  *rlwe.EvaluationKey
	automorphismIndex []uint64
}

// NewDomainSwitcher instantiate a new [DomainSwitcher] type. It may be instantiated from parameters from either RingType.
// The method returns an error if the parameters cannot support the switching (e.g., the NTTs are undefined for
// either of the two ring types).
// The comlexToRealEvk and comlexToRealEvk EvaluationKeys can be generated using [rlwe.KeyGenerator.GenEvaluationKeysForRingSwap].
func NewDomainSwitcher(params Parameters, comlexToRealEvk, realToComplexEvk *rlwe.EvaluationKey) (DomainSwitcher, error) {

	s := DomainSwitcher{
		stdToci: comlexToRealEvk,
		ciToStd: realToComplexEvk,
	}
	var err error
	if s.stdRingQ, err = params.RingQ().StandardRing(); err != nil {
		return DomainSwitcher{}, fmt.Errorf("cannot NewDomainSwitcher because the standard NTT is undefined for params: %s", err)
	}
	if s.conjugateRingQ, err = params.RingQ().ConjugateInvariantRing(); err != nil {
		return DomainSwitcher{}, fmt.Errorf("cannot NewDomainSwitcher because the standard NTT is undefined for params: %s", err)
	}

	// Sanity check, this error should not happen unless the
	// algorithm has been modified to provide invalid inputs.
	if s.automorphismIndex, err = ring.AutomorphismNTTIndex(s.stdRingQ.N(), s.stdRingQ.NthRoot(), s.stdRingQ.NthRoot()-1); err != nil {
		panic(err)
	}

	return s, nil
}

// ComplexToReal switches the provided ciphertext ctIn from the standard domain to the conjugate
// invariant domain and writes the result into opOut.
// Given ctInCKKS = enc(real(m) + imag(m)) in Z[X](X^N + 1), returns opOutCI = enc(real(m))
// in Z[X+X^-1]/(X^N + 1) in compressed form (N/2 coefficients).
// The scale of the output ciphertext is twice the scale of the input one.
// Requires the ring degree of opOut to be half the ring degree of ctIn.
// The security is changed from Z[X]/(X^N+1) to Z[X]/(X^N/2+1).
// The method will return an error if the DomainSwitcher was not initialized with a the appropriate EvaluationKeys.
func (switcher DomainSwitcher) ComplexToReal(eval *Evaluator, ctIn, opOut *rlwe.Ciphertext) (err error) {

	evalRLWE := eval.Evaluator

	if evalRLWE.GetRLWEParameters().RingType() != ring.Standard {
		return fmt.Errorf("cannot ComplexToReal: provided evaluator is not instantiated with RingType ring.Standard")
	}

	level := utils.Min(ctIn.Level(), opOut.Level())

	if ctIn.Value[0].N() != 2*opOut.Value[0].N() {
		return fmt.Errorf("cannot ComplexToReal: ctIn ring degree must be twice opOut ring degree")
	}

	opOut.Resize(1, level)

	if switcher.stdToci == nil {
		return fmt.Errorf("cannot ComplexToReal: no realToComplexEvk provided to this DomainSwitcher")
	}

	ctTmp := &rlwe.Ciphertext{}
	ctTmp.Value = []ring.Poly{evalRLWE.BuffQP[1].Q, evalRLWE.BuffQP[2].Q}
	ctTmp.MetaData = ctIn.MetaData

	evalRLWE.GadgetProduct(level, ctIn.Value[1], &switcher.stdToci.GadgetCiphertext, ctTmp)
	switcher.stdRingQ.AtLevel(level).Add(evalRLWE.BuffQP[1].Q, ctIn.Value[0], evalRLWE.BuffQP[1].Q)

	switcher.conjugateRingQ.AtLevel(level).FoldStandardToConjugateInvariant(evalRLWE.BuffQP[1].Q, switcher.automorphismIndex, opOut.Value[0])
	switcher.conjugateRingQ.AtLevel(level).FoldStandardToConjugateInvariant(evalRLWE.BuffQP[2].Q, switcher.automorphismIndex, opOut.Value[1])
	*opOut.MetaData = *ctIn.MetaData
	opOut.Scale = ctIn.Scale.Mul(rlwe.NewScale(2))
	return
}

// RealToComplex switches the provided ciphertext ctIn from the conjugate invariant domain to the
// standard domain and writes the result into opOut.
// Given ctInCI = enc(real(m)) in Z[X+X^-1]/(X^2N+1) in compressed form (N coefficients), returns
// opOutCKKS = enc(real(m) + imag(0)) in Z[X]/(X^2N+1).
// Requires the ring degree of opOut to be twice the ring degree of ctIn.
// The security is changed from Z[X]/(X^N+1) to Z[X]/(X^2N+1).
// The method will return an error if the [DomainSwitcher] was not initialized with a the appropriate EvaluationKeys.
func (switcher DomainSwitcher) RealToComplex(eval *Evaluator, ctIn, opOut *rlwe.Ciphertext) (err error) {

	evalRLWE := eval.Evaluator

	if evalRLWE.GetRLWEParameters().RingType() != ring.Standard {
		return fmt.Errorf("cannot RealToComplex: provided evaluator is not instantiated with RingType ring.Standard")
	}

	level := utils.Min(ctIn.Level(), opOut.Level())

	if 2*ctIn.Value[0].N() != opOut.Value[0].N() {
		return fmt.Errorf("cannot RealToComplex: opOut ring degree must be twice ctIn ring degree")
	}

	opOut.Resize(1, level)

	if switcher.ciToStd == nil {
		return fmt.Errorf("cannot RealToComplex: no realToComplexEvk provided to this DomainSwitcher")
	}

	switcher.stdRingQ.AtLevel(level).UnfoldConjugateInvariantToStandard(ctIn.Value[0], opOut.Value[0])
	switcher.stdRingQ.AtLevel(level).UnfoldConjugateInvariantToStandard(ctIn.Value[1], opOut.Value[1])

	ctTmp := &rlwe.Ciphertext{}
	ctTmp.Value = []ring.Poly{evalRLWE.BuffQP[1].Q, evalRLWE.BuffQP[2].Q}
	ctTmp.MetaData = ctIn.MetaData

	// Switches the RCKswitcher key [X+X^-1] to a CKswitcher key [X]
	evalRLWE.GadgetProduct(level, opOut.Value[1], &switcher.ciToStd.GadgetCiphertext, ctTmp)
	switcher.stdRingQ.AtLevel(level).Add(opOut.Value[0], evalRLWE.BuffQP[1].Q, opOut.Value[0])
	opOut.Value[1].CopyLvl(level, evalRLWE.BuffQP[2].Q)
	*opOut.MetaData = *ctIn.MetaData
	return
}
