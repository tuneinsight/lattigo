package multiparty

import (
	"fmt"
	"io"

	"github.com/tuneinsight/lattigo/v6/ring"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

// PublicKeySwitchProtocol is the structure storing the parameters for the collective public key-switching.
type PublicKeySwitchProtocol struct {
	params rlwe.Parameters
	noise  ring.DistributionParameters

	buf ring.Poly

	*rlwe.Encryptor
	noiseSampler ring.Sampler
}

// PublicKeySwitchShare represents a party's share in the PublicKeySwitch protocol.
type PublicKeySwitchShare struct {
	rlwe.Element[ring.Poly]
}

// NewPublicKeySwitchProtocol creates a new PublicKeySwitchProtocol object and will be used to re-encrypt a ciphertext ctx encrypted under a secret-shared key among j parties under a new
// collective public-key.
func NewPublicKeySwitchProtocol(params rlwe.ParameterProvider, noiseFlooding ring.DistributionParameters) (pcks PublicKeySwitchProtocol, err error) {
	pcks = PublicKeySwitchProtocol{}
	pcks.params = *params.GetRLWEParameters()
	pcks.noise = noiseFlooding

	pcks.buf = pcks.params.RingQ().NewPoly()

	prng, err := sampling.NewPRNG()

	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}

	pcks.Encryptor = rlwe.NewEncryptor(pcks.params, nil)

	switch noiseFlooding.(type) {
	case ring.DiscreteGaussian:
	default:
		return pcks, fmt.Errorf("invalid distribution type, expected %T but got %T", ring.DiscreteGaussian{}, noiseFlooding)
	}

	pcks.noiseSampler, err = ring.NewSampler(prng, pcks.params.RingQ(), noiseFlooding, false)

	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}

	return pcks, nil
}

// AllocateShare allocates the shares of the PublicKeySwitch protocol.
func (pcks PublicKeySwitchProtocol) AllocateShare(levelQ int) (s PublicKeySwitchShare) {
	return PublicKeySwitchShare{*rlwe.NewElement(pcks.params, 1, levelQ)}
}

// GenShare computes a party's share in the PublicKeySwitch protocol from secret-key sk to public-key pk.
// ct is the rlwe.Ciphertext to keyswitch. Note that ct.Value[0] is not used by the function and can be nil/zero.
//
// Expected noise: ctNoise + encFreshPk + smudging
func (pcks PublicKeySwitchProtocol) GenShare(sk *rlwe.SecretKey, pk *rlwe.PublicKey, ct *rlwe.Ciphertext, shareOut *PublicKeySwitchShare) {

	levelQ := utils.Min(shareOut.Level(), ct.Value[1].Level())

	ringQ := pcks.params.RingQ().AtLevel(levelQ)

	// Encrypt zero
	enc := pcks.Encryptor.WithKey(pk)

	if err := enc.EncryptZero(&rlwe.Ciphertext{
		Element: rlwe.Element[ring.Poly]{
			Value: []ring.Poly{
				shareOut.Value[0],
				shareOut.Value[1],
			},
			MetaData: ct.MetaData,
		},
	}); err != nil {
		// Sanity check, this error should not happen.
		panic(err)
	}

	// Add ct[1] * s and noise
	if ct.IsNTT {
		ringQ.MulCoeffsMontgomeryThenAdd(ct.Value[1], sk.Value.Q, shareOut.Value[0])
		pcks.noiseSampler.Read(pcks.buf)
		ringQ.NTT(pcks.buf, pcks.buf)
		ringQ.Add(shareOut.Value[0], pcks.buf, shareOut.Value[0])
	} else {
		ringQ.NTTLazy(ct.Value[1], pcks.buf)
		ringQ.MulCoeffsMontgomeryLazy(pcks.buf, sk.Value.Q, pcks.buf)
		ringQ.INTT(pcks.buf, pcks.buf)
		pcks.noiseSampler.ReadAndAdd(pcks.buf)
		ringQ.Add(shareOut.Value[0], pcks.buf, shareOut.Value[0])
	}
}

// AggregateShares is the second part of the first and unique round of the PublicKeySwitchProtocol protocol. Each party upon receiving the j-1 elements from the
// other parties computes :
//
// [ctx[0] + sum(s_i * ctx[0] + u_i * pk[0] + e_0i), sum(u_i * pk[1] + e_1i)]
func (pcks PublicKeySwitchProtocol) AggregateShares(share1, share2 PublicKeySwitchShare, shareOut *PublicKeySwitchShare) (err error) {
	levelQ1, levelQ2 := share1.Value[0].Level(), share1.Value[1].Level()
	if levelQ1 != levelQ2 {
		return fmt.Errorf("cannot AggregateShares: the two shares are at different levelQ")
	}
	pcks.params.RingQ().AtLevel(levelQ1).Add(share1.Value[0], share2.Value[0], shareOut.Value[0])
	pcks.params.RingQ().AtLevel(levelQ1).Add(share1.Value[1], share2.Value[1], shareOut.Value[1])
	return
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in opOut
func (pcks PublicKeySwitchProtocol) KeySwitch(ctIn *rlwe.Ciphertext, combined PublicKeySwitchShare, opOut *rlwe.Ciphertext) {

	level := ctIn.Level()

	if ctIn != opOut {
		opOut.Resize(ctIn.Degree(), level)
		*opOut.MetaData = *ctIn.MetaData
	}

	pcks.params.RingQ().AtLevel(level).Add(ctIn.Value[0], combined.Value[0], opOut.Value[0])

	opOut.Value[1].CopyLvl(level, combined.Value[1])
}

// ShallowCopy creates a shallow copy of [PublicKeySwitchProtocol] in which all the read-only data-structures are
// shared with the receiver and the temporary bufers are reallocated. The receiver and the returned
// [PublicKeySwitchProtocol] can be used concurrently.
func (pcks PublicKeySwitchProtocol) ShallowCopy() PublicKeySwitchProtocol {
	prng, err := sampling.NewPRNG()

	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}

	params := pcks.params

	Xe, err := ring.NewSampler(prng, params.RingQ(), pcks.noise, false)

	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}

	return PublicKeySwitchProtocol{
		noiseSampler: Xe,
		noise:        pcks.noise,
		Encryptor:    pcks.Encryptor.ShallowCopy(),
		params:       params,
		buf:          params.RingQ().NewPoly(),
	}
}

// BinarySize returns the serialized size of the object in bytes.
func (share PublicKeySwitchShare) BinarySize() int {
	return share.Element.BinarySize()
}

// WriteTo writes the object on an [io.Writer]. It implements the [io.WriterTo]
// interface, and will write exactly object.BinarySize() bytes on w.
//
// Unless w implements the [buffer.Writer] interface (see lattigo/utils/buffer/writer.go),
// it will be wrapped into a [bufio.Writer]. Since this requires allocations, it
// is preferable to pass a [buffer.Writer] directly:
//
//   - When writing multiple times to a [io.Writer], it is preferable to first wrap the
//     [io.Writer] in a pre-allocated [bufio.Writer].
//   - When writing to a pre-allocated var b []byte, it is preferable to pass
//     buffer.NewBuffer(b) as w (see lattigo/utils/buffer/buffer.go).
func (share PublicKeySwitchShare) WriteTo(w io.Writer) (n int64, err error) {
	return share.Element.WriteTo(w)
}

// ReadFrom reads on the object from an [io.Writer]. It implements the
// [io.ReaderFrom] interface.
//
// Unless r implements the [buffer.Reader] interface (see see lattigo/utils/buffer/reader.go),
// it will be wrapped into a [bufio.Reader]. Since this requires allocation, it
// is preferable to pass a [buffer.Reader] directly:
//
//   - When reading multiple values from a [io.Reader], it is preferable to first
//     first wrap [io.Reader] in a pre-allocated [bufio.Reader].
//   - When reading from a var b []byte, it is preferable to pass a buffer.NewBuffer(b)
//     as w (see lattigo/utils/buffer/buffer.go).
func (share *PublicKeySwitchShare) ReadFrom(r io.Reader) (n int64, err error) {
	return share.Element.ReadFrom(r)
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (share PublicKeySwitchShare) MarshalBinary() (p []byte, err error) {
	return share.Element.MarshalBinary()
}

// UnmarshalBinary decodes a slice of bytes generated by
// [PublicKeySwitchShare.MarshalBinary] or [PublicKeySwitchShare.WriteTo] on the object.
func (share *PublicKeySwitchShare) UnmarshalBinary(p []byte) (err error) {
	return share.Element.UnmarshalBinary(p)
}
