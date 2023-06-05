package drlwe

import (
	"fmt"
	"io"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/ring/distribution"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

// PublicKeySwitchProtocol is the structure storing the parameters for the collective public key-switching.
type PublicKeySwitchProtocol struct {
	params rlwe.Parameters
	noise  distribution.Distribution

	buf *ring.Poly

	rlwe.EncryptorInterface
	noiseSampler ring.Sampler
}

// PublicKeySwitchShare represents a party's share in the PublicKeySwitch protocol.
type PublicKeySwitchShare struct {
	rlwe.OperandQ
}

// ShallowCopy creates a shallow copy of PublicKeySwitchProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary bufers are reallocated. The receiver and the returned
// PublicKeySwitchProtocol can be used concurrently.
func (pcks *PublicKeySwitchProtocol) ShallowCopy() *PublicKeySwitchProtocol {
	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}

	params := pcks.params

		noiseSampler:       ring.NewSampler(prng, params.RingQ(), pcks.noise, false),
		noise:              pcks.noise,
		EncryptorInterface: rlwe.NewEncryptor(params, nil),
		params:             params,
		buf:                params.RingQ().NewPoly(),
	return &PublicKeySwitchProtocol{
	}
}

// NewPublicKeySwitchProtocol creates a new PublicKeySwitchProtocol object and will be used to re-encrypt a ciphertext ctx encrypted under a secret-shared key among j parties under a new
// collective public-key.
func NewPublicKeySwitchProtocol(params rlwe.Parameters, noise distribution.Distribution) (pcks *PublicKeySwitchProtocol) {
	pcks = new(PublicKeySwitchProtocol)
	pcks.params = params
	pcks.noise = noise.CopyNew()

	pcks.buf = params.RingQ().NewPoly()

	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}

	pcks.EncryptorInterface = rlwe.NewEncryptor(params, nil)

	switch noise.(type) {
	case *distribution.DiscreteGaussian:
	default:
		panic(fmt.Sprintf("invalid distribution type, expected %T but got %T", &distribution.DiscreteGaussian{}, noise))
	}

	pcks.noiseSampler = ring.NewSampler(prng, params.RingQ(), noise, false)

	return pcks
}

// AllocateShare allocates the shares of the PublicKeySwitch protocol.
func (pcks *PublicKeySwitchProtocol) AllocateShare(levelQ int) (s *PublicKeySwitchShare) {
	return &PublicKeySwitchShare{*rlwe.NewOperandQ(pcks.params, 1, levelQ)}
}

// GenShare computes a party's share in the PublicKeySwitch protocol from secret-key sk to public-key pk.
// ct is the rlwe.Ciphertext to keyswitch. Note that ct.Value[0] is not used by the function and can be nil/zero.
//
// Expected noise: ctNoise + encFreshPk + smudging
func (pcks *PublicKeySwitchProtocol) GenShare(sk *rlwe.SecretKey, pk *rlwe.PublicKey, ct *rlwe.Ciphertext, shareOut *PublicKeySwitchShare) {

	levelQ := utils.Min(shareOut.Level(), ct.Value[1].Level())

	ringQ := pcks.params.RingQ().AtLevel(levelQ)

	// Encrypt zero
	pcks.EncryptorInterface.WithKey(pk).EncryptZero(&rlwe.Ciphertext{
		OperandQ: rlwe.OperandQ{
			Value: []*ring.Poly{
				shareOut.Value[0],
				shareOut.Value[1],
			},
			MetaData: ct.MetaData,
		},
	})

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
func (pcks *PublicKeySwitchProtocol) AggregateShares(share1, share2, shareOut *PublicKeySwitchShare) {
	levelQ1, levelQ2 := share1.Value[0].Level(), share1.Value[1].Level()
	if levelQ1 != levelQ2 {
		panic("cannot AggregateShares: the two shares are at different levelQ.")
	}
	pcks.params.RingQ().AtLevel(levelQ1).Add(share1.Value[0], share2.Value[0], shareOut.Value[0])
	pcks.params.RingQ().AtLevel(levelQ1).Add(share1.Value[1], share2.Value[1], shareOut.Value[1])

}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (pcks *PublicKeySwitchProtocol) KeySwitch(ctIn *rlwe.Ciphertext, combined *PublicKeySwitchShare, ctOut *rlwe.Ciphertext) {

	level := ctIn.Level()

	if ctIn != ctOut {
		ctOut.Resize(ctIn.Degree(), level)
		ctOut.MetaData = ctIn.MetaData
	}

	pcks.params.RingQ().AtLevel(level).Add(ctIn.Value[0], combined.Value[0], ctOut.Value[0])

	ring.CopyLvl(level, combined.Value[1], ctOut.Value[1])
}

// BinarySize returns the size in bytes of the object
// when encoded using Encode.
func (share *PublicKeySwitchShare) BinarySize() int {
	return share.OperandQ.BinarySize()
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (share *PublicKeySwitchShare) MarshalBinary() (p []byte, err error) {
	return share.OperandQ.MarshalBinary()
}

// Encode encodes the object into a binary form on a preallocated slice of bytes
// and returns the number of bytes written.
func (share *PublicKeySwitchShare) Encode(p []byte) (n int, err error) {
	return share.OperandQ.Encode(p)
}

// WriteTo writes the object on an io.Writer.
// To ensure optimal efficiency and minimal allocations, the user is encouraged
// to provide a struct implementing the interface bufer.Writer, which defines
// a subset of the method of the bufio.Writer.
// If w is not compliant to the bufer.Writer interface, it will be wrapped in
// a new bufio.Writer.
// For additional information, see lattigo/utils/bufer/writer.go.
func (share *PublicKeySwitchShare) WriteTo(w io.Writer) (n int64, err error) {
	return share.OperandQ.WriteTo(w)
}

// UnmarshalBinary decodes a slice of bytes generated by
// MarshalBinary or WriteTo on the object.
func (share *PublicKeySwitchShare) UnmarshalBinary(p []byte) (err error) {
	return share.OperandQ.UnmarshalBinary(p)
}

// Decode decodes a slice of bytes generated by Encode
// on the object and returns the number of bytes read.
func (share *PublicKeySwitchShare) Decode(p []byte) (n int, err error) {
	return share.OperandQ.Decode(p)
}

// ReadFrom reads on the object from an io.Writer.
// To ensure optimal efficiency and minimal allocations, the user is encouraged
// to provide a struct implementing the interface bufer.Reader, which defines
// a subset of the method of the bufio.Reader.
// If r is not compliant to the bufer.Reader interface, it will be wrapped in
// a new bufio.Reader.
// For additional information, see lattigo/utils/bufer/reader.go.
func (share *PublicKeySwitchShare) ReadFrom(r io.Reader) (n int64, err error) {
	return share.OperandQ.ReadFrom(r)
}
