package multiparty

import (
	"io"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

// PublicKeyGenProtocol is the structure storing the parameters and and precomputations for
// the collective encryption key generation protocol.
type PublicKeyGenProtocol struct {
	params           rlwe.Parameters
	gaussianSamplerQ ring.Sampler
}

// PublicKeyGenShare is a struct storing the PublicKeyGen protocol's share.
type PublicKeyGenShare struct {
	Value ringqp.Poly
}

// PublicKeyGenCRP is a type for common reference polynomials in the PublicKeyGen protocol.
type PublicKeyGenCRP struct {
	Value ringqp.Poly
}

// NewPublicKeyGenProtocol creates a new [PublicKeyGenProtocol] instance
func NewPublicKeyGenProtocol(params rlwe.ParameterProvider) PublicKeyGenProtocol {
	ckg := PublicKeyGenProtocol{}
	ckg.params = *params.GetRLWEParameters()

	var err error
	prng, err := sampling.NewPRNG()

	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}

	ckg.gaussianSamplerQ, err = ring.NewSampler(prng, ckg.params.RingQ(), ckg.params.Xe(), false)

	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}

	return ckg
}

// AllocateShare allocates the share of the PublicKeyGen protocol.
func (ckg PublicKeyGenProtocol) AllocateShare() PublicKeyGenShare {
	return PublicKeyGenShare{ckg.params.RingQP().NewPoly()}
}

// SampleCRP samples a common random polynomial to be used in the PublicKeyGen protocol from the provided
// common reference string.
func (ckg PublicKeyGenProtocol) SampleCRP(crs CRS) PublicKeyGenCRP {
	crp := ckg.params.RingQP().NewPoly()
	ringqp.NewUniformSampler(crs, *ckg.params.RingQP()).Read(crp)
	return PublicKeyGenCRP{crp}
}

// GenShare generates the party's public key share from its secret key as:
//
// crp*s_i + e_i
//
// for the receiver protocol. Has no effect is the share was already generated.
func (ckg PublicKeyGenProtocol) GenShare(sk *rlwe.SecretKey, crp PublicKeyGenCRP, shareOut *PublicKeyGenShare) {
	ringQP := ckg.params.RingQP()

	ckg.gaussianSamplerQ.Read(shareOut.Value.Q)

	if ringQP.RingP != nil {
		ringQP.ExtendBasisSmallNormAndCenter(shareOut.Value.Q, ckg.params.MaxLevelP(), shareOut.Value.Q, shareOut.Value.P)
	}

	ringQP.NTT(shareOut.Value, shareOut.Value)
	ringQP.MForm(shareOut.Value, shareOut.Value)

	ringQP.MulCoeffsMontgomeryThenSub(sk.Value, crp.Value, shareOut.Value)
}

// AggregateShares aggregates a new share to the aggregate key
func (ckg PublicKeyGenProtocol) AggregateShares(share1, share2 PublicKeyGenShare, shareOut *PublicKeyGenShare) {
	ckg.params.RingQP().Add(share1.Value, share2.Value, shareOut.Value)
}

// GenPublicKey return the current aggregation of the received shares as a [rlwe.PublicKey].
func (ckg PublicKeyGenProtocol) GenPublicKey(roundShare PublicKeyGenShare, crp PublicKeyGenCRP, pubkey *rlwe.PublicKey) {
	pubkey.Value[0].Copy(roundShare.Value)
	pubkey.Value[1].Copy(crp.Value)
}

// ShallowCopy creates a shallow copy of [PublicKeyGenProtocol] in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// [PublicKeyGenProtocol] can be used concurrently.
func (ckg PublicKeyGenProtocol) ShallowCopy() PublicKeyGenProtocol {
	prng, err := sampling.NewPRNG()

	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}

	sampler, err := ring.NewSampler(prng, ckg.params.RingQ(), ckg.params.Xe(), false)

	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}

	return PublicKeyGenProtocol{ckg.params, sampler}
}

// BinarySize returns the serialized size of the object in bytes.
func (share PublicKeyGenShare) BinarySize() int {
	return share.Value.BinarySize()
}

// WriteTo writes the object on an [io.Writer]. It implements the [io.WriterTo]
// interface, and will write exactly object.BinarySize() bytes on w.
//
// Unless w implements the [buffer.Writer] interface (see lattigo/utils/buffer/writer.go),
// it will be wrapped into a [bufio.Writer]. Since this requires allocations, it
// is preferable to pass a [buffer.Writer] directly:
//
//   - When writing multiple times to a io.Writer, it is preferable to first wrap the
//     io.Writer in a pre-allocated [bufio.Writer].
//   - When writing to a pre-allocated var b []byte, it is preferable to pass
//     buffer.NewBuffer(b) as w (see lattigo/utils/buffer/buffer.go).
func (share PublicKeyGenShare) WriteTo(w io.Writer) (n int64, err error) {
	return share.Value.WriteTo(w)
}

// ReadFrom reads on the object from an [io.Writer]. It implements the
// [io.ReaderFrom] interface.
//
// Unless r implements the [buffer.Reader] interface (see see lattigo/utils/buffer/reader.go),
// it will be wrapped into a [bufio.Reader]. Since this requires allocation, it
// is preferable to pass a [buffer.Reader] directly:
//
//   - When reading multiple values from a io.Reader, it is preferable to first
//     first wrap [io.Reader] in a pre-allocated bufio.Reader.
//   - When reading from a var b []byte, it is preferable to pass a buffer.NewBuffer(b)
//     as w (see lattigo/utils/buffer/buffer.go).
func (share *PublicKeyGenShare) ReadFrom(r io.Reader) (n int64, err error) {
	return share.Value.ReadFrom(r)
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (share PublicKeyGenShare) MarshalBinary() (p []byte, err error) {
	return share.Value.MarshalBinary()
}

// UnmarshalBinary decodes a slice of bytes generated by
// [PublicKeyGenShare.MarshalBinary] or [PublicKeyGenShare.WriteTo] on the object.
func (share *PublicKeyGenShare) UnmarshalBinary(p []byte) (err error) {
	return share.Value.UnmarshalBinary(p)
}
