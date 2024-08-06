package multiparty

import (
	"fmt"
	"io"
	"math"

	"github.com/tuneinsight/lattigo/v6/ring"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

// KeySwitchProtocol is the structure storing the parameters and and precomputations for the collective key-switching protocol.
type KeySwitchProtocol struct {
	params       rlwe.Parameters
	noise        ring.DistributionParameters
	noiseSampler ring.Sampler
	buf          ring.Poly
	bufDelta     ring.Poly
}

// KeySwitchShare is a type for the KeySwitch protocol shares.
type KeySwitchShare struct {
	Value ring.Poly
}

// ShallowCopy creates a shallow copy of [KeySwitchProtocol] in which all the read-only data-structures are
// shared with the receiver and the temporary bufers are reallocated. The receiver and the returned
// [KeySwitchProtocol] can be used concurrently.
func (cks KeySwitchProtocol) ShallowCopy() KeySwitchProtocol {
	prng, err := sampling.NewPRNG()

	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}

	params := cks.params

	Xe, err := ring.NewSampler(prng, cks.params.RingQ(), cks.noise, false)

	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}

	return KeySwitchProtocol{
		params:       params,
		noiseSampler: Xe,
		buf:          params.RingQ().NewPoly(),
		bufDelta:     params.RingQ().NewPoly(),
		noise:        cks.noise,
	}
}

// KeySwitchCRP is a type for common reference polynomials in the KeySwitch protocol.
type KeySwitchCRP struct {
	Value ring.Poly
}

// NewKeySwitchProtocol creates a new [KeySwitchProtocol] that will be used to perform a collective key-switching on a ciphertext encrypted under a collective public-key, whose
// secret-shares are distributed among j parties, re-encrypting the ciphertext under another public-key, whose secret-shares are also known to the
// parties.
func NewKeySwitchProtocol(params rlwe.ParameterProvider, noiseFlooding ring.DistributionParameters) (KeySwitchProtocol, error) {
	cks := KeySwitchProtocol{}
	cks.params = *params.GetRLWEParameters()
	prng, err := sampling.NewPRNG()

	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}

	// EncFreshSK + sigmaSmudging

	switch noise := noiseFlooding.(type) {
	case ring.DiscreteGaussian:
		eFresh := cks.params.NoiseFreshSK()
		eNoise := noise.Sigma
		eSigma := math.Sqrt(eFresh*eFresh + eNoise*eNoise)
		cks.noise = ring.DiscreteGaussian{Sigma: eSigma, Bound: 6 * eSigma}
	default:
		return cks, fmt.Errorf("invalid distribution type, expected %T but got %T", ring.DiscreteGaussian{}, noise)
	}

	cks.noiseSampler, err = ring.NewSampler(prng, cks.params.RingQ(), cks.noise, false)

	// Sanity check, this error should not happen.
	if err != nil {
		panic(err)
	}

	cks.buf = cks.params.RingQ().NewPoly()
	cks.bufDelta = cks.params.RingQ().NewPoly()
	return cks, nil
}

// AllocateShare allocates the shares of the KeySwitchProtocol
func (cks KeySwitchProtocol) AllocateShare(level int) KeySwitchShare {
	return KeySwitchShare{cks.params.RingQ().AtLevel(level).NewPoly()}
}

// SampleCRP samples a common random polynomial to be used in the KeySwitch protocol from the provided
// common reference string.
func (cks KeySwitchProtocol) SampleCRP(level int, crs CRS) KeySwitchCRP {
	ringQ := cks.params.RingQ().AtLevel(level)
	crp := ringQ.NewPoly()
	ring.NewUniformSampler(crs, ringQ).Read(crp)
	return KeySwitchCRP{Value: crp}
}

// GenShare computes a party's share in the KeySwitchcol from secret-key skInput to secret-key skOutput.
// ct is the [rlwe.Ciphertext] to keyswitch. Note that ct.Value[0] is not used by the function and can be nil/zero.
//
// Expected noise: ctNoise + encFreshSk + smudging
func (cks KeySwitchProtocol) GenShare(skInput, skOutput *rlwe.SecretKey, ct *rlwe.Ciphertext, shareOut *KeySwitchShare) {

	levelQ := utils.Min(shareOut.Value.Level(), ct.Value[1].Level())

	shareOut.Value.Resize(levelQ)

	ringQ := cks.params.RingQ().AtLevel(levelQ)

	ringQ.Sub(skInput.Value.Q, skOutput.Value.Q, cks.bufDelta)

	var c1NTT ring.Poly
	if !ct.IsNTT {
		ringQ.NTTLazy(ct.Value[1], cks.buf)
		c1NTT = cks.buf
	} else {
		c1NTT = ct.Value[1]
	}

	// c1NTT * (skIn - skOut)
	ringQ.MulCoeffsMontgomeryLazy(c1NTT, cks.bufDelta, shareOut.Value)

	if !ct.IsNTT {
		// InvNTT(c1NTT * (skIn - skOut)) + e
		ringQ.INTTLazy(shareOut.Value, shareOut.Value)
		cks.noiseSampler.AtLevel(levelQ).ReadAndAdd(shareOut.Value)
	} else {
		// c1NTT * (skIn - skOut) + e
		cks.noiseSampler.AtLevel(levelQ).Read(cks.buf)
		ringQ.NTT(cks.buf, cks.buf)
		ringQ.Add(shareOut.Value, cks.buf, shareOut.Value)
	}
}

// AggregateShares is the second part of the unique round of the [KeySwitchProtocol] protocol. Upon receiving the j-1 elements each party computes:
//
// [ctx[0] + sum((skInput_i - skOutput_i) * ctx[0] + e_i), ctx[1]]
func (cks KeySwitchProtocol) AggregateShares(share1, share2 KeySwitchShare, shareOut *KeySwitchShare) (err error) {
	if share1.Level() != share2.Level() || share1.Level() != shareOut.Level() {
		return fmt.Errorf("cannot AggregateShares: shares levels do not match")
	}

	cks.params.RingQ().AtLevel(share1.Level()).Add(share1.Value, share2.Value, shareOut.Value)
	return
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in opOut
func (cks KeySwitchProtocol) KeySwitch(ctIn *rlwe.Ciphertext, combined KeySwitchShare, opOut *rlwe.Ciphertext) {

	level := ctIn.Level()

	if ctIn != opOut {

		opOut.Resize(ctIn.Degree(), level)

		opOut.Value[1].CopyLvl(level, ctIn.Value[1])

		*opOut.MetaData = *ctIn.MetaData
	}

	cks.params.RingQ().AtLevel(level).Add(ctIn.Value[0], combined.Value, opOut.Value[0])
}

// Level returns the level of the target share.
func (ckss KeySwitchShare) Level() int {
	return ckss.Value.Level()
}

// BinarySize returns the serialized size of the object in bytes.
func (ckss KeySwitchShare) BinarySize() int {
	return ckss.Value.BinarySize()
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
func (ckss KeySwitchShare) WriteTo(w io.Writer) (n int64, err error) {
	return ckss.Value.WriteTo(w)
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
func (ckss *KeySwitchShare) ReadFrom(r io.Reader) (n int64, err error) {
	return ckss.Value.ReadFrom(r)
}

// MarshalBinary encodes a KeySwitch share on a slice of bytes.
func (ckss KeySwitchShare) MarshalBinary() (p []byte, err error) {
	return ckss.Value.MarshalBinary()
}

// UnmarshalBinary decodes a slice of bytes generated by
// [KeySwitchShare.MarshalBinary] or[KeySwitchShare.WriteTo] on the object.
func (ckss *KeySwitchShare) UnmarshalBinary(p []byte) (err error) {
	return ckss.Value.UnmarshalBinary(p)
}
