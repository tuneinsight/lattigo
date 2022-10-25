package dbgv

import (
	"github.com/tuneinsight/lattigo/v4/bgv"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// E2SProtocol is the structure storing the parameters and temporary buffers
// required by the encryption-to-shares protocol.
type E2SProtocol struct {
	CKSProtocol
	params bgv.Parameters

	maskSampler *ring.UniformSampler
	encoder     bgv.Encoder

	zero              *rlwe.SecretKey
	tmpPlaintextRingT *ring.Poly
	tmpPlaintextRingQ *ring.Poly
}

// ShallowCopy creates a shallow copy of E2SProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// E2SProtocol can be used concurrently.
func (e2s *E2SProtocol) ShallowCopy() *E2SProtocol {

	params := e2s.params

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	return &E2SProtocol{
		CKSProtocol:       *e2s.CKSProtocol.ShallowCopy(),
		params:            e2s.params,
		maskSampler:       ring.NewUniformSampler(prng, params.RingT()),
		encoder:           e2s.encoder.ShallowCopy(),
		zero:              e2s.zero,
		tmpPlaintextRingT: params.RingT().NewPoly(),
		tmpPlaintextRingQ: params.RingQ().NewPoly(),
	}
}

// NewE2SProtocol creates a new E2SProtocol struct from the passed bgv parameters.
func NewE2SProtocol(params bgv.Parameters, sigmaSmudging float64) *E2SProtocol {
	e2s := new(E2SProtocol)
	e2s.CKSProtocol = *NewCKSProtocol(params, sigmaSmudging)
	e2s.params = params
	e2s.encoder = bgv.NewEncoder(params)
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	e2s.maskSampler = ring.NewUniformSampler(prng, params.RingT())
	e2s.zero = rlwe.NewSecretKey(params.Parameters)
	e2s.tmpPlaintextRingQ = params.RingQ().NewPoly()
	e2s.tmpPlaintextRingT = params.RingT().NewPoly()
	return e2s
}

// AllocateShare allocates a share of the E2S protocol
func (e2s *E2SProtocol) AllocateShare(level int) (share *drlwe.CKSShare) {
	share = e2s.CKSProtocol.AllocateShare(level)
	share.Value.IsNTT = true
	return
}

// GenShare generates a party's share in the encryption-to-shares protocol. This share consist in the additive secret-share of the party
// which is written in secretShareOut and in the public masked-decryption share written in publicShareOut.
// ct1 is degree 1 element of a bgv.Ciphertext, i.e. bgv.Ciphertext.Value[1].
func (e2s *E2SProtocol) GenShare(sk *rlwe.SecretKey, ct1 *ring.Poly, secretShareOut *rlwe.AdditiveShare, publicShareOut *drlwe.CKSShare) {
	level := utils.MinInt(ct1.Level(), publicShareOut.Value.Level())
	e2s.CKSProtocol.GenShare(sk, e2s.zero, ct1, publicShareOut)
	e2s.maskSampler.Read(&secretShareOut.Value)
	e2s.encoder.RingT2Q(level, &secretShareOut.Value, e2s.tmpPlaintextRingQ)
	e2s.params.RingQ().NTTLvl(level, e2s.tmpPlaintextRingQ, e2s.tmpPlaintextRingQ)
	e2s.params.RingQ().SubLvl(level, publicShareOut.Value, e2s.tmpPlaintextRingQ, publicShareOut.Value)
}

// GetShare is the final step of the encryption-to-share protocol. It performs the masked decryption of the target ciphertext followed by a
// the removal of the caller's secretShare as generated in the GenShare method.
// If the caller is not secret-key-share holder (i.e., didn't generate a decryption share), `secretShare` can be set to nil.
// Therefore, in order to obtain an additive sharing of the message, only one party should call this method, and the other parties should use
// the secretShareOut output of the GenShare method.
func (e2s *E2SProtocol) GetShare(secretShare *rlwe.AdditiveShare, aggregatePublicShare *drlwe.CKSShare, ct *rlwe.Ciphertext, secretShareOut *rlwe.AdditiveShare) {
	level := utils.MinInt(ct.Level(), aggregatePublicShare.Value.Level())
	e2s.params.RingQ().AddLvl(level, aggregatePublicShare.Value, ct.Value[0], e2s.tmpPlaintextRingQ)
	e2s.params.RingQ().InvNTTLvl(level, e2s.tmpPlaintextRingQ, e2s.tmpPlaintextRingQ)
	e2s.encoder.RingQ2T(level, e2s.tmpPlaintextRingQ, e2s.tmpPlaintextRingT)
	if secretShare != nil {
		e2s.params.RingT().Add(&secretShare.Value, e2s.tmpPlaintextRingT, &secretShareOut.Value)
	} else {
		secretShareOut.Value.Copy(e2s.tmpPlaintextRingT)
	}
}

// S2EProtocol is the structure storing the parameters and temporary buffers
// required by the shares-to-encryption protocol.
type S2EProtocol struct {
	CKSProtocol
	params bgv.Parameters

	encoder bgv.Encoder

	zero              *rlwe.SecretKey
	tmpPlaintextRingQ *ring.Poly
}

// NewS2EProtocol creates a new S2EProtocol struct from the passed bgv parameters.
func NewS2EProtocol(params bgv.Parameters, sigmaSmudging float64) *S2EProtocol {
	s2e := new(S2EProtocol)
	s2e.CKSProtocol = *NewCKSProtocol(params, sigmaSmudging)
	s2e.params = params
	s2e.encoder = bgv.NewEncoder(params)
	s2e.zero = rlwe.NewSecretKey(params.Parameters)
	s2e.tmpPlaintextRingQ = params.RingQ().NewPoly()
	return s2e
}

// AllocateShare allocates a share of the S2E protocol
func (s2e S2EProtocol) AllocateShare(level int) (share *drlwe.CKSShare) {
	share = s2e.CKSProtocol.AllocateShare(level)
	share.Value.IsNTT = true
	return
}

// ShallowCopy creates a shallow copy of S2EProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// S2EProtocol can be used concurrently.
func (s2e *S2EProtocol) ShallowCopy() *S2EProtocol {
	params := s2e.params
	return &S2EProtocol{
		CKSProtocol:       *s2e.CKSProtocol.ShallowCopy(),
		encoder:           s2e.encoder.ShallowCopy(),
		params:            params,
		zero:              s2e.zero,
		tmpPlaintextRingQ: params.RingQ().NewPoly(),
	}
}

// GenShare generates a party's in the shares-to-encryption protocol given the party's secret-key share `sk`, a common
// polynomial sampled from the CRS `crp` and the party's secret share of the message.
func (s2e *S2EProtocol) GenShare(sk *rlwe.SecretKey, crp drlwe.CKSCRP, secretShare *rlwe.AdditiveShare, c0ShareOut *drlwe.CKSShare) {

	c1 := ring.Poly(crp)

	if c1.Level() != c0ShareOut.Value.Level() {
		panic("cannot GenShare: c1 and c0ShareOut level must be equal")
	}

	s2e.CKSProtocol.GenShare(s2e.zero, sk, &c1, c0ShareOut)
	s2e.encoder.RingT2Q(c1.Level(), &secretShare.Value, s2e.tmpPlaintextRingQ)
	s2e.params.RingQ().NTTLvl(c1.Level(), s2e.tmpPlaintextRingQ, s2e.tmpPlaintextRingQ)
	s2e.params.RingQ().AddLvl(c1.Level(), c0ShareOut.Value, s2e.tmpPlaintextRingQ, c0ShareOut.Value)
}

// GetEncryption computes the final encryption of the secret-shared message when provided with the aggregation `c0Agg` of the parties'
// shares in the protocol and with the common, CRS-sampled polynomial `crp`.
func (s2e *S2EProtocol) GetEncryption(c0Agg *drlwe.CKSShare, crp drlwe.CKSCRP, ctOut *rlwe.Ciphertext) {
	if ctOut.Degree() != 1 {
		panic("cannot GetEncryption: ctOut must have degree 1.")
	}
	ctOut.Value[0].Copy(c0Agg.Value)
	ctOut.Value[1].Copy((*ring.Poly)(&crp))
}
