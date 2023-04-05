// Package dckks implements a distributed (or threshold) version of the CKKS scheme that enables secure multiparty computation solutions with secret-shared secret keys.
package dckks

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// E2SProtocol is the structure storing the parameters and temporary buffers
// required by the encryption-to-shares protocol.
type E2SProtocol struct {
	*drlwe.CKSProtocol

	params     ckks.Parameters
	zero       *rlwe.SecretKey
	maskBigint []*big.Int
	buff       *ring.Poly
}

// ShallowCopy creates a shallow copy of E2SProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// E2SProtocol can be used concurrently.
func (e2s *E2SProtocol) ShallowCopy() *E2SProtocol {

	maskBigint := make([]*big.Int, len(e2s.maskBigint))
	for i := range maskBigint {
		maskBigint[i] = new(big.Int)
	}

	return &E2SProtocol{
		CKSProtocol: e2s.CKSProtocol.ShallowCopy(),
		params:      e2s.params,
		zero:        e2s.zero,
		maskBigint:  maskBigint,
		buff:        e2s.params.RingQ().NewPoly(),
	}
}

// NewE2SProtocol creates a new E2SProtocol struct from the passed CKKS parameters.
func NewE2SProtocol(params ckks.Parameters, sigmaSmudging float64) *E2SProtocol {
	e2s := new(E2SProtocol)
	e2s.CKSProtocol = drlwe.NewCKSProtocol(params.Parameters, sigmaSmudging)
	e2s.params = params
	e2s.zero = rlwe.NewSecretKey(params.Parameters)
	e2s.maskBigint = make([]*big.Int, params.N())
	for i := range e2s.maskBigint {
		e2s.maskBigint[i] = new(big.Int)
	}
	e2s.buff = e2s.params.RingQ().NewPoly()
	return e2s
}

// AllocateShare allocates a share of the E2S protocol
func (e2s *E2SProtocol) AllocateShare(level int) (share *drlwe.CKSShare) {
	return e2s.CKSProtocol.AllocateShare(level)
}

// GenShare generates a party's share in the encryption-to-shares protocol. This share consist in the additive secret-share of the party
// which is written in secretShareOut and in the public masked-decryption share written in publicShareOut.
// This protocol requires additional inputs which are :
// logBound : the bit length of the masks
// logSlots : the bit length of the number of slots
// ct1      : the degree 1 element the ciphertext to share, i.e. ct1 = ckk.Ciphertext.Value[1].
// The method "GetMinimumLevelForBootstrapping" should be used to get the minimum level at which E2S can be called while still ensure 128-bits of security, as well as the
// value for logBound.
func (e2s *E2SProtocol) GenShare(sk *rlwe.SecretKey, logBound uint, logSlots int, ct *rlwe.Ciphertext, secretShareOut *drlwe.AdditiveShareBigint, publicShareOut *drlwe.CKSShare) {

	ct1 := ct.Value[1]

	levelQ := utils.Min(ct1.Level(), publicShareOut.Value.Level())

	ringQ := e2s.params.RingQ().AtLevel(levelQ)

	// Get the upperbound on the norm
	// Ensures that bound >= 2^{128+logbound}
	bound := ring.NewUint(1)
	bound.Lsh(bound, uint(logBound))

	boundMax := new(big.Int).Set(ringQ.ModulusAtLevel[levelQ])

	var sign int

	sign = bound.Cmp(boundMax)

	if sign == 1 || bound.Cmp(boundMax) == 1 {
		panic("cannot GenShare: ciphertext level is not large enough for refresh correctness")
	}

	boundHalf := new(big.Int).Rsh(bound, 1)

	dslots := 1 << logSlots
	if ringQ.Type() == ring.Standard {
		dslots *= 2
	}

	// Generate the mask in Z[Y] for Y = X^{N/(2*slots)}
	for i := 0; i < dslots; i++ {
		e2s.maskBigint[i] = ring.RandInt(bound)
		sign = e2s.maskBigint[i].Cmp(boundHalf)
		if sign == 1 || sign == 0 {
			e2s.maskBigint[i].Sub(e2s.maskBigint[i], bound)
		}

		secretShareOut.Value[i].Set(e2s.maskBigint[i])
	}

	// Encrypt the mask
	// Generates an encryption of zero and subtracts the mask
	e2s.CKSProtocol.GenShare(sk, e2s.zero, ct, publicShareOut)

	ringQ.SetCoefficientsBigint(secretShareOut.Value[:dslots], e2s.buff)

	// Maps Y^{N/n} -> X^{N} in Montgomery and NTT
	ckks.NttSparseAndMontgomery(ringQ, logSlots, false, e2s.buff)

	// Subtracts the mask to the encryption of zero
	ringQ.Sub(publicShareOut.Value, e2s.buff, publicShareOut.Value)
}

// GetShare is the final step of the encryption-to-share protocol. It performs the masked decryption of the target ciphertext followed by a
// the removal of the caller's secretShare as generated in the GenShare method.
// If the caller is not secret-key-share holder (i.e., didn't generate a decryption share), `secretShare` can be set to nil.
// Therefore, in order to obtain an additive sharing of the message, only one party should call this method, and the other parties should use
// the secretShareOut output of the GenShare method.
func (e2s *E2SProtocol) GetShare(secretShare *drlwe.AdditiveShareBigint, aggregatePublicShare *drlwe.CKSShare, logSlots int, ct *rlwe.Ciphertext, secretShareOut *drlwe.AdditiveShareBigint) {

	levelQ := utils.Min(ct.Level(), aggregatePublicShare.Value.Level())

	ringQ := e2s.params.RingQ().AtLevel(levelQ)

	// Adds the decryption share on the ciphertext and stores the result in a buff
	ringQ.Add(aggregatePublicShare.Value, ct.Value[0], e2s.buff)

	// Switches the LSSS RNS NTT ciphertext outside of the NTT domain
	ringQ.INTT(e2s.buff, e2s.buff)

	dslots := 1 << logSlots
	if ringQ.Type() == ring.Standard {
		dslots *= 2
	}

	gap := ringQ.N() / dslots

	// Switches the LSSS RNS ciphertext outside of the RNS domain
	ringQ.PolyToBigintCentered(e2s.buff, gap, e2s.maskBigint)

	// Subtracts the last mask
	if secretShare != nil {
		a := secretShareOut.Value
		b := e2s.maskBigint
		c := secretShare.Value
		for i := range secretShareOut.Value[:dslots] {
			a[i].Add(c[i], b[i])
		}
	} else {
		a := secretShareOut.Value
		b := e2s.maskBigint
		for i := range secretShareOut.Value[:dslots] {
			a[i].Set(b[i])
		}
	}
}

// S2EProtocol is the structure storing the parameters and temporary buffers
// required by the shares-to-encryption protocol.
type S2EProtocol struct {
	*drlwe.CKSProtocol
	params   ckks.Parameters
	tmp      *ring.Poly
	ssBigint []*big.Int
	zero     *rlwe.SecretKey
}

// ShallowCopy creates a shallow copy of S2EProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// S2EProtocol can be used concurrently.
func (s2e *S2EProtocol) ShallowCopy() *S2EProtocol {
	return &S2EProtocol{
		CKSProtocol: s2e.CKSProtocol.ShallowCopy(),
		params:      s2e.params,
		tmp:         s2e.params.RingQ().NewPoly(),
		ssBigint:    make([]*big.Int, s2e.params.N()),
		zero:        s2e.zero,
	}
}

// NewS2EProtocol creates a new S2EProtocol struct from the passed CKKS parameters.
func NewS2EProtocol(params ckks.Parameters, sigmaSmudging float64) *S2EProtocol {
	s2e := new(S2EProtocol)
	s2e.CKSProtocol = drlwe.NewCKSProtocol(params.Parameters, sigmaSmudging)
	s2e.params = params
	s2e.tmp = s2e.params.RingQ().NewPoly()
	s2e.ssBigint = make([]*big.Int, s2e.params.N())
	s2e.zero = rlwe.NewSecretKey(params.Parameters)
	return s2e
}

// AllocateShare allocates a share of the S2E protocol
func (s2e S2EProtocol) AllocateShare(level int) (share *drlwe.CKSShare) {
	return s2e.CKSProtocol.AllocateShare(level)
}

// GenShare generates a party's in the shares-to-encryption protocol given the party's secret-key share `sk`, a common
// polynomial sampled from the CRS `crs` and the party's secret share of the message.
func (s2e *S2EProtocol) GenShare(sk *rlwe.SecretKey, crs drlwe.CKSCRP, logSlots int, secretShare *drlwe.AdditiveShareBigint, c0ShareOut *drlwe.CKSShare) {

	if crs.Value.Level() != c0ShareOut.Value.Level() {
		panic("cannot GenShare: crs and c0ShareOut level must be equal")
	}

	ringQ := s2e.params.RingQ().AtLevel(crs.Value.Level())

	// Generates an encryption share
	ct := &rlwe.Ciphertext{}
	ct.Value = []*ring.Poly{nil, &crs.Value}
	ct.MetaData.IsNTT = true
	s2e.CKSProtocol.GenShare(s2e.zero, sk, ct, c0ShareOut)

	dslots := 1 << logSlots
	if ringQ.Type() == ring.Standard {
		dslots *= 2
	}

	ringQ.SetCoefficientsBigint(secretShare.Value[:dslots], s2e.tmp)

	// Maps Y^{N/n} -> X^{N} in Montgomery and NTT
	ckks.NttSparseAndMontgomery(ringQ, logSlots, false, s2e.tmp)

	ringQ.Add(c0ShareOut.Value, s2e.tmp, c0ShareOut.Value)
}

// GetEncryption computes the final encryption of the secret-shared message when provided with the aggregation `c0Agg` of the parties'
// share in the protocol and with the common, CRS-sampled polynomial `crs`.
func (s2e *S2EProtocol) GetEncryption(c0Agg *drlwe.CKSShare, crs drlwe.CKSCRP, ctOut *rlwe.Ciphertext) {

	if ctOut.Degree() != 1 {
		panic("cannot GetEncryption: ctOut must have degree 1.")
	}

	if c0Agg.Value.Level() != crs.Value.Level() {
		panic("cannot GetEncryption: c0Agg level must be equal to crs level")
	}

	if ctOut.Level() != crs.Value.Level() {
		panic("cannot GetEncryption: ctOut level must be equal to crs level")
	}

	ctOut.Value[0].Copy(c0Agg.Value)
	ctOut.Value[1].Copy(&crs.Value)
}
