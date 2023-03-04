package drlwe

import (
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// PCKSShare represents a party's share in the PCKS protocol.
type PCKSShare struct {
	Value [2]*ring.Poly
}

// PCKSProtocol is the structure storing the parameters for the collective public key-switching.
type PCKSProtocol struct {
	params        rlwe.Parameters
	sigmaSmudging float64

	buff *ring.Poly

	rlwe.Encryptor
	gaussianSampler *ring.GaussianSampler
}

// ShallowCopy creates a shallow copy of PCKSProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// PCKSProtocol can be used concurrently.
func (pcks *PCKSProtocol) ShallowCopy() *PCKSProtocol {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	params := pcks.params

	return &PCKSProtocol{
		params:          params,
		Encryptor:       rlwe.NewEncryptor(params, nil),
		sigmaSmudging:   pcks.sigmaSmudging,
		buff:            params.RingQ().NewPoly(),
		gaussianSampler: ring.NewGaussianSampler(prng, params.RingQ(), pcks.sigmaSmudging, int(6*pcks.sigmaSmudging)),
	}
}

// NewPCKSProtocol creates a new PCKSProtocol object and will be used to re-encrypt a ciphertext ctx encrypted under a secret-shared key among j parties under a new
// collective public-key.
func NewPCKSProtocol(params rlwe.Parameters, sigmaSmudging float64) (pcks *PCKSProtocol) {
	pcks = new(PCKSProtocol)
	pcks.params = params
	pcks.sigmaSmudging = sigmaSmudging

	pcks.buff = params.RingQ().NewPoly()

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	pcks.Encryptor = rlwe.NewEncryptor(params, nil)

	pcks.gaussianSampler = ring.NewGaussianSampler(prng, params.RingQ(), sigmaSmudging, int(6*sigmaSmudging))

	return pcks
}

// AllocateShare allocates the shares of the PCKS protocol.
func (pcks *PCKSProtocol) AllocateShare(levelQ int) (s *PCKSShare) {
	return &PCKSShare{[2]*ring.Poly{pcks.params.RingQ().AtLevel(levelQ).NewPoly(), pcks.params.RingQ().AtLevel(levelQ).NewPoly()}}
}

// GenShare computes a party's share in the PCKS protocol from secret-key sk to public-key pk.
// ct is the rlwe.Ciphertext to keyswitch. Note that ct.Value[0] is not used by the function and can be nil/zero.
//
// Expected noise: ctNoise + encFreshPk + smudging
func (pcks *PCKSProtocol) GenShare(sk *rlwe.SecretKey, pk *rlwe.PublicKey, ct *rlwe.Ciphertext, shareOut *PCKSShare) {

	levelQ := utils.MinInt(shareOut.Value[0].Level(), ct.Value[1].Level())

	ringQ := pcks.params.RingQ().AtLevel(levelQ)

	// Encrypt zero
	pcks.Encryptor.WithKey(pk).EncryptZero(&rlwe.Ciphertext{
		Value: []*ring.Poly{
			shareOut.Value[0],
			shareOut.Value[1],
		},
		MetaData: ct.MetaData,
	})

	// Add ct[1] * s and noise
	if ct.IsNTT {
		ringQ.MulCoeffsMontgomeryThenAdd(ct.Value[1], sk.Value.Q, shareOut.Value[0])
		pcks.gaussianSampler.Read(pcks.buff)
		ringQ.NTT(pcks.buff, pcks.buff)
		ringQ.Add(shareOut.Value[0], pcks.buff, shareOut.Value[0])
	} else {
		ringQ.NTTLazy(ct.Value[1], pcks.buff)
		ringQ.MulCoeffsMontgomeryLazy(pcks.buff, sk.Value.Q, pcks.buff)
		ringQ.INTT(pcks.buff, pcks.buff)
		pcks.gaussianSampler.ReadAndAdd(pcks.buff)
		ringQ.Add(shareOut.Value[0], pcks.buff, shareOut.Value[0])
	}
}

// AggregateShares is the second part of the first and unique round of the PCKSProtocol protocol. Each party upon receiving the j-1 elements from the
// other parties computes :
//
// [ctx[0] + sum(s_i * ctx[0] + u_i * pk[0] + e_0i), sum(u_i * pk[1] + e_1i)]
func (pcks *PCKSProtocol) AggregateShares(share1, share2, shareOut *PCKSShare) {
	levelQ1, levelQ2 := share1.Value[0].Level(), share1.Value[1].Level()
	if levelQ1 != levelQ2 {
		panic("cannot AggregateShares: the two shares are at different levelQ.")
	}
	pcks.params.RingQ().AtLevel(levelQ1).Add(share1.Value[0], share2.Value[0], shareOut.Value[0])
	pcks.params.RingQ().AtLevel(levelQ1).Add(share1.Value[1], share2.Value[1], shareOut.Value[1])

}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (pcks *PCKSProtocol) KeySwitch(ctIn *rlwe.Ciphertext, combined *PCKSShare, ctOut *rlwe.Ciphertext) {

	level := ctIn.Level()

	if ctIn != ctOut {
		ctOut.Resize(ctIn.Degree(), level)
		ctOut.MetaData = ctIn.MetaData
	}

	pcks.params.RingQ().AtLevel(level).Add(ctIn.Value[0], combined.Value[0], ctOut.Value[0])

	ring.CopyLvl(level, combined.Value[1], ctOut.Value[1])
}

// MarshalBinary encodes a PCKS share on a slice of bytes.
func (share *PCKSShare) MarshalBinary() (data []byte, err error) {
	data = make([]byte, share.Value[0].MarshalBinarySize64()+share.Value[1].MarshalBinarySize64())
	var inc, pt int
	if inc, err = share.Value[0].Encode64(data[pt:]); err != nil {
		return nil, err
	}
	pt += inc

	if _, err = share.Value[1].Encode64(data[pt:]); err != nil {
		return nil, err
	}
	return
}

// UnmarshalBinary decodes marshaled PCKS share on the target PCKS share.
func (share *PCKSShare) UnmarshalBinary(data []byte) (err error) {
	var pt, inc int
	share.Value[0] = new(ring.Poly)
	if inc, err = share.Value[0].Decode64(data[pt:]); err != nil {
		return
	}
	pt += inc

	share.Value[1] = new(ring.Poly)
	if _, err = share.Value[1].Decode64(data[pt:]); err != nil {
		return
	}
	return
}
