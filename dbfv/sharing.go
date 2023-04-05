package dbfv

import (
	"github.com/tuneinsight/lattigo/v4/bfv"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

// E2SProtocol is the structure storing the parameters and temporary buffers
// required by the encryption-to-shares protocol.
type E2SProtocol struct {
	*drlwe.CKSProtocol
	params bfv.Parameters

	maskSampler *ring.UniformSampler
	encoder     bfv.Encoder

	zero              *rlwe.SecretKey
	tmpPlaintextRingT *bfv.PlaintextRingT
	tmpPlaintext      *rlwe.Plaintext
}

// ShallowCopy creates a shallow copy of E2SProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// E2SProtocol can be used concurrently.
func (e2s *E2SProtocol) ShallowCopy() *E2SProtocol {

	params := e2s.params

	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}

	return &E2SProtocol{
		CKSProtocol:       e2s.CKSProtocol.ShallowCopy(),
		params:            e2s.params,
		maskSampler:       ring.NewUniformSampler(prng, params.RingT()),
		encoder:           e2s.encoder.ShallowCopy(),
		zero:              e2s.zero,
		tmpPlaintextRingT: bfv.NewPlaintextRingT(params),
		tmpPlaintext:      bfv.NewPlaintext(params, params.MaxLevel()),
	}
}

// NewE2SProtocol creates a new E2SProtocol struct from the passed BFV parameters.
func NewE2SProtocol(params bfv.Parameters, sigmaSmudging float64) *E2SProtocol {
	e2s := new(E2SProtocol)
	e2s.CKSProtocol = drlwe.NewCKSProtocol(params.Parameters, sigmaSmudging)
	e2s.params = params
	e2s.encoder = bfv.NewEncoder(params)
	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}
	e2s.maskSampler = ring.NewUniformSampler(prng, params.RingT())
	e2s.zero = rlwe.NewSecretKey(params.Parameters)
	e2s.tmpPlaintext = bfv.NewPlaintext(params, params.MaxLevel())
	e2s.tmpPlaintextRingT = bfv.NewPlaintextRingT(params)
	return e2s
}

// GenShare generates a party's share in the encryption-to-shares protocol. This share consist in the additive secret-share of the party
// which is written in secretShareOut and in the public masked-decryption share written in publicShareOut.
// ct1 is degree 1 element of a bfv.Ciphertext, i.e. bfv.Ciphertext.Value[1].
func (e2s *E2SProtocol) GenShare(sk *rlwe.SecretKey, ct *rlwe.Ciphertext, secretShareOut *drlwe.AdditiveShare, publicShareOut *drlwe.CKSShare) {
	e2s.CKSProtocol.GenShare(sk, e2s.zero, ct, publicShareOut)
	e2s.maskSampler.Read(&secretShareOut.Value)
	e2s.encoder.ScaleUp(&bfv.PlaintextRingT{Plaintext: &rlwe.Plaintext{Value: &secretShareOut.Value}}, e2s.tmpPlaintext)
	e2s.params.RingQ().Sub(publicShareOut.Value, e2s.tmpPlaintext.Value, publicShareOut.Value)
}

// GetShare is the final step of the encryption-to-share protocol. It performs the masked decryption of the target ciphertext followed by a
// the removal of the caller's secretShare as generated in the GenShare method.
// If the caller is not secret-key-share holder (i.e., didn't generate a decryption share), `secretShare` can be set to nil.
// Therefore, in order to obtain an additive sharing of the message, only one party should call this method, and the other parties should use
// the secretShareOut output of the GenShare method.
func (e2s *E2SProtocol) GetShare(secretShare *drlwe.AdditiveShare, aggregatePublicShare *drlwe.CKSShare, ct *rlwe.Ciphertext, secretShareOut *drlwe.AdditiveShare) {
	e2s.params.RingQ().Add(aggregatePublicShare.Value, ct.Value[0], e2s.tmpPlaintext.Value)
	e2s.encoder.ScaleDown(e2s.tmpPlaintext, e2s.tmpPlaintextRingT)
	if secretShare != nil {
		e2s.params.RingT().Add(&secretShare.Value, e2s.tmpPlaintextRingT.Value, &secretShareOut.Value)
	} else {
		secretShareOut.Value.Copy(e2s.tmpPlaintextRingT.Value)
	}
}

// S2EProtocol is the structure storing the parameters and temporary buffers
// required by the shares-to-encryption protocol.
type S2EProtocol struct {
	*drlwe.CKSProtocol
	params bfv.Parameters

	encoder bfv.Encoder

	zero         *rlwe.SecretKey
	tmpPlaintext *rlwe.Plaintext
}

// NewS2EProtocol creates a new S2EProtocol struct from the passed BFV parameters.
func NewS2EProtocol(params bfv.Parameters, sigmaSmudging float64) *S2EProtocol {
	s2e := new(S2EProtocol)
	s2e.CKSProtocol = drlwe.NewCKSProtocol(params.Parameters, sigmaSmudging)
	s2e.params = params
	s2e.encoder = bfv.NewEncoder(params)
	s2e.zero = rlwe.NewSecretKey(params.Parameters)
	s2e.tmpPlaintext = bfv.NewPlaintext(params, params.MaxLevel())
	return s2e
}

// ShallowCopy creates a shallow copy of S2EProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// S2EProtocol can be used concurrently.
func (s2e *S2EProtocol) ShallowCopy() *S2EProtocol {
	params := s2e.params
	return &S2EProtocol{
		CKSProtocol:  s2e.CKSProtocol.ShallowCopy(),
		encoder:      s2e.encoder.ShallowCopy(),
		params:       params,
		zero:         s2e.zero,
		tmpPlaintext: bfv.NewPlaintext(params, params.MaxLevel()),
	}
}

// GenShare generates a party's in the shares-to-encryption protocol given the party's secret-key share `sk`, a common
// polynomial sampled from the CRS `crp` and the party's secret share of the message.
func (s2e *S2EProtocol) GenShare(sk *rlwe.SecretKey, crp drlwe.CKSCRP, secretShare *drlwe.AdditiveShare, c0ShareOut *drlwe.CKSShare) {
	s2e.encoder.ScaleUp(&bfv.PlaintextRingT{Plaintext: &rlwe.Plaintext{Value: &secretShare.Value}}, s2e.tmpPlaintext)
	ct := &rlwe.Ciphertext{}
	ct.Value = []*ring.Poly{nil, &crp.Value}
	ct.MetaData.IsNTT = false
	s2e.CKSProtocol.GenShare(s2e.zero, sk, ct, c0ShareOut)
	s2e.params.RingQ().Add(c0ShareOut.Value, s2e.tmpPlaintext.Value, c0ShareOut.Value)
}

// GetEncryption computes the final encryption of the secret-shared message when provided with the aggregation `c0Agg` of the parties'
// shares in the protocol and with the common, CRS-sampled polynomial `crp`.
func (s2e *S2EProtocol) GetEncryption(c0Agg *drlwe.CKSShare, crp drlwe.CKSCRP, ctOut *rlwe.Ciphertext) {
	if ctOut.Degree() != 1 {
		panic("cannot GetEncryption: ctOut must have degree 1.")
	}
	ctOut.Value[0].Copy(c0Agg.Value)
	ctOut.Value[1].Copy(&crp.Value)
}
