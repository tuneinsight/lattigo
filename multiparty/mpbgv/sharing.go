package mpbgv

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/multiparty"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/bgv"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

// EncToShareProtocol is the structure storing the parameters and temporary buffers
// required by the encryption-to-shares protocol.
type EncToShareProtocol struct {
	multiparty.KeySwitchProtocol
	params bgv.Parameters

	maskSampler *ring.UniformSampler
	encoder     *bgv.Encoder

	zero *rlwe.SecretKey
}

func NewAdditiveShare(params bgv.Parameters) multiparty.AdditiveShare {
	return multiparty.NewAdditiveShare(params.RingT())
}

// ShallowCopy creates a shallow copy of [EncToShareProtocol] in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// [EncToShareProtocol] can be used concurrently.
func (e2s EncToShareProtocol) ShallowCopy() EncToShareProtocol {

	params := e2s.params

	prng, err := sampling.NewPRNG()

	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}

	return EncToShareProtocol{
		KeySwitchProtocol: e2s.KeySwitchProtocol.ShallowCopy(),
		params:            e2s.params,
		maskSampler:       ring.NewUniformSampler(prng, params.RingT()),
		encoder:           e2s.encoder.ShallowCopy(),
		zero:              e2s.zero,
	}
}

// NewEncToShareProtocol creates a new EncToShareProtocol struct from the passed bgv.Parameters.
func NewEncToShareProtocol(params bgv.Parameters, noiseFlooding ring.DistributionParameters) (EncToShareProtocol, error) {
	e2s := EncToShareProtocol{}

	var err error
	if e2s.KeySwitchProtocol, err = multiparty.NewKeySwitchProtocol(params.Parameters, noiseFlooding); err != nil {
		return EncToShareProtocol{}, err
	}

	e2s.params = params
	e2s.encoder = bgv.NewEncoder(params)
	prng, err := sampling.NewPRNG()

	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}

	e2s.maskSampler = ring.NewUniformSampler(prng, params.RingT())

	e2s.zero = rlwe.NewSecretKey(params.Parameters)
	return e2s, nil
}

// AllocateShare allocates a share of the EncToShare protocol
func (e2s EncToShareProtocol) AllocateShare(level int) (share multiparty.KeySwitchShare) {
	return e2s.KeySwitchProtocol.AllocateShare(level)
}

// GenShare generates a party's share in the encryption-to-shares protocol. This share consist in the additive secret-share of the party
// which is written in secretShareOut and in the public masked-decryption share written in publicShareOut.
// ct1 is degree 1 element of a rlwe.Ciphertext, i.e. rlwe.Ciphertext.Value[1].
func (e2s EncToShareProtocol) GenShare(sk *rlwe.SecretKey, ct *rlwe.Ciphertext, secretShareOut *multiparty.AdditiveShare, publicShareOut *multiparty.KeySwitchShare) {
	level := utils.Min(ct.Level(), publicShareOut.Value.Level())
	e2s.KeySwitchProtocol.GenShare(sk, e2s.zero, ct, publicShareOut)
	e2s.maskSampler.Read(secretShareOut.Value)
	ringQ := e2s.params.RingQ().AtLevel(level)
	buffQ := ringQ.NewPoly()
	e2s.encoder.RingT2Q(level, true, secretShareOut.Value, buffQ)
	ringQ.NTT(buffQ, buffQ)
	ringQ.Sub(publicShareOut.Value, buffQ, publicShareOut.Value)
}

// GetShare is the final step of the encryption-to-share protocol. It performs the masked decryption of the target ciphertext followed by a
// the removal of the caller's secretShare as generated in the [EncToShareProtocol.GenShare] method.
// If the caller is not secret-key-share holder (i.e., didn't generate a decryption share), secretShare can be set to nil.
// Therefore, in order to obtain an additive sharing of the message, only one party should call this method, and the other parties should use
// the secretShareOut output of the GenShare method.
func (e2s EncToShareProtocol) GetShare(secretShare *multiparty.AdditiveShare, aggregatePublicShare multiparty.KeySwitchShare, ct *rlwe.Ciphertext, secretShareOut *multiparty.AdditiveShare) {
	level := utils.Min(ct.Level(), aggregatePublicShare.Value.Level())
	ringQ := e2s.params.RingQ().AtLevel(level)
	buffQ := ringQ.NewPoly()
	buffT := e2s.params.RingT().NewPoly()
	ringQ.Add(aggregatePublicShare.Value, ct.Value[0], buffQ)
	ringQ.INTT(buffQ, buffQ)
	e2s.encoder.RingQ2T(level, true, buffQ, buffT)
	if secretShare != nil {
		e2s.params.RingT().Add(secretShare.Value, buffT, secretShareOut.Value)
	} else {
		secretShareOut.Value.Copy(buffT)
	}
}

// ShareToEncProtocol is the structure storing the parameters and temporary buffers
// required by the shares-to-encryption protocol.
type ShareToEncProtocol struct {
	multiparty.KeySwitchProtocol
	params bgv.Parameters

	encoder *bgv.Encoder

	zero *rlwe.SecretKey
}

// NewShareToEncProtocol creates a new ShareToEncProtocol struct from the passed integer parameters.
func NewShareToEncProtocol(params bgv.Parameters, noiseFlooding ring.DistributionParameters) (ShareToEncProtocol, error) {
	s2e := ShareToEncProtocol{}

	var err error
	if s2e.KeySwitchProtocol, err = multiparty.NewKeySwitchProtocol(params.Parameters, noiseFlooding); err != nil {
		return ShareToEncProtocol{}, err
	}

	s2e.params = params
	s2e.encoder = bgv.NewEncoder(params)
	s2e.zero = rlwe.NewSecretKey(params.Parameters)
	return s2e, nil
}

// AllocateShare allocates a share of the ShareToEnc protocol
func (s2e ShareToEncProtocol) AllocateShare(level int) (share multiparty.KeySwitchShare) {
	return s2e.KeySwitchProtocol.AllocateShare(level)
}

// ShallowCopy creates a shallow copy of [ShareToEncProtocol] in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// ShareToEncProtocol can be used concurrently.
func (s2e ShareToEncProtocol) ShallowCopy() ShareToEncProtocol {
	params := s2e.params
	return ShareToEncProtocol{
		KeySwitchProtocol: s2e.KeySwitchProtocol.ShallowCopy(),
		encoder:           s2e.encoder.ShallowCopy(),
		params:            params,
		zero:              s2e.zero,
	}
}

// GenShare generates a party's in the shares-to-encryption protocol given the party's secret-key share `sk`, a common
// polynomial sampled from the CRS `crp` and the party's secret share of the message.
func (s2e ShareToEncProtocol) GenShare(sk *rlwe.SecretKey, crp multiparty.KeySwitchCRP, secretShare multiparty.AdditiveShare, c0ShareOut *multiparty.KeySwitchShare) (err error) {

	if crp.Value.Level() != c0ShareOut.Value.Level() {
		return fmt.Errorf("cannot GenShare: crp and c0ShareOut level must be equal")
	}

	ct := &rlwe.Ciphertext{}
	ct.Value = []ring.Poly{{}, crp.Value}
	ct.MetaData = &rlwe.MetaData{}
	ct.MetaData.IsNTT = true
	s2e.KeySwitchProtocol.GenShare(s2e.zero, sk, ct, c0ShareOut)
	ringQ := s2e.params.RingQ().AtLevel(crp.Value.Level())
	buffQ := ringQ.NewPoly()
	s2e.encoder.RingT2Q(crp.Value.Level(), true, secretShare.Value, buffQ)
	ringQ.NTT(buffQ, buffQ)
	ringQ.Add(c0ShareOut.Value, buffQ, c0ShareOut.Value)
	return
}

// GetEncryption computes the final encryption of the secret-shared message when provided with the aggregation `c0Agg` of the parties'
// shares in the protocol and with the common, CRS-sampled polynomial `crp`.
func (s2e ShareToEncProtocol) GetEncryption(c0Agg multiparty.KeySwitchShare, crp multiparty.KeySwitchCRP, opOut *rlwe.Ciphertext) (err error) {
	if opOut.Degree() != 1 {
		return fmt.Errorf("cannot GetEncryption: opOut must have degree 1")
	}
	opOut.Value[0].Copy(c0Agg.Value)
	opOut.Value[1].Copy(crp.Value)
	return
}
