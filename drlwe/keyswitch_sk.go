package drlwe

import (
	"bytes"
	"fmt"
	"io"
	"math"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/ring/distribution"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

// KeySwitchProtocol is the structure storing the parameters and and precomputations for the collective key-switching protocol.
type KeySwitchProtocol struct {
	params       rlwe.Parameters
	noise        distribution.Distribution
	noiseSampler ring.Sampler
	buf          *ring.Poly
	bufDelta     *ring.Poly
}

// KeySwitchShare is a type for the KeySwitch protocol shares.
type KeySwitchShare struct {
	Value *ring.Poly
}

// ShallowCopy creates a shallow copy of KeySwitchProtocol in which all the read-only data-structures are
// shared with the receiver and the temporary bufers are reallocated. The receiver and the returned
// KeySwitchProtocol can be used concurrently.
func (cks *KeySwitchProtocol) ShallowCopy() *KeySwitchProtocol {
	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}

	params := cks.params

	return &KeySwitchProtocol{
		params:       params,
		noiseSampler: ring.NewSampler(prng, cks.params.RingQ(), cks.noise, false),
		buf:          params.RingQ().NewPoly(),
		bufDelta:     params.RingQ().NewPoly(),
	}
}

// KeySwitchCRP is a type for common reference polynomials in the KeySwitch protocol.
type KeySwitchCRP struct {
	Value ring.Poly
}

// NewKeySwitchProtocol creates a new KeySwitchProtocol that will be used to perform a collective key-switching on a ciphertext encrypted under a collective public-key, whose
// secret-shares are distributed among j parties, re-encrypting the ciphertext under another public-key, whose secret-shares are also known to the
// parties.
func NewKeySwitchProtocol(params rlwe.Parameters, noise distribution.Distribution) *KeySwitchProtocol {
	cks := new(KeySwitchProtocol)
	cks.params = params
	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}

	// EncFreshSK + sigmaSmudging

	switch noise.(type) {
	case *distribution.DiscreteGaussian:
		eFresh := params.NoiseFreshSK()
		eNoise := noise.StandardDeviation(0, 0)
		eSigma := math.Sqrt(eFresh*eFresh + eNoise*eNoise)
		cks.noise = &distribution.DiscreteGaussian{Sigma: eSigma, Bound: 6 * eSigma}
	default:
		panic(fmt.Sprintf("invalid distribution type, expected %T but got %T", &distribution.DiscreteGaussian{}, noise))
	}

	cks.noiseSampler = ring.NewSampler(prng, params.RingQ(), cks.noise, false)

	cks.buf = params.RingQ().NewPoly()
	cks.bufDelta = params.RingQ().NewPoly()
	return cks
}

// AllocateShare allocates the shares of the KeySwitchProtocol
func (cks *KeySwitchProtocol) AllocateShare(level int) *KeySwitchShare {
	return &KeySwitchShare{cks.params.RingQ().AtLevel(level).NewPoly()}
}

// SampleCRP samples a common random polynomial to be used in the KeySwitch protocol from the provided
// common reference string.
func (cks *KeySwitchProtocol) SampleCRP(level int, crs CRS) KeySwitchCRP {
	ringQ := cks.params.RingQ().AtLevel(level)
	crp := ringQ.NewPoly()
	ring.NewUniformSampler(crs, ringQ).Read(crp)
	return KeySwitchCRP{Value: *crp}
}

// GenShare computes a party's share in the KeySwitchcol from secret-key skInput to secret-key skOutput.
// ct is the rlwe.Ciphertext to keyswitch. Note that ct.Value[0] is not used by the function and can be nil/zero.
//
// Expected noise: ctNoise + encFreshSk + smudging
func (cks *KeySwitchProtocol) GenShare(skInput, skOutput *rlwe.SecretKey, ct *rlwe.Ciphertext, shareOut *KeySwitchShare) {

	levelQ := utils.Min(shareOut.Value.Level(), ct.Value[1].Level())

	shareOut.Value.Resize(levelQ)

	ringQ := cks.params.RingQ().AtLevel(levelQ)

	ringQ.Sub(skInput.Value.Q, skOutput.Value.Q, cks.bufDelta)

	var c1NTT *ring.Poly
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

// AggregateShares is the second part of the unique round of the KeySwitchProtocol protocol. Upon receiving the j-1 elements each party computes :
//
// [ctx[0] + sum((skInput_i - skOutput_i) * ctx[0] + e_i), ctx[1]]
func (cks *KeySwitchProtocol) AggregateShares(share1, share2, shareOut *KeySwitchShare) {
	if share1.Level() != share2.Level() || share1.Level() != shareOut.Level() {
		panic("shares levels do not match")
	}

	cks.params.RingQ().AtLevel(share1.Level()).Add(share1.Value, share2.Value, shareOut.Value)
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (cks *KeySwitchProtocol) KeySwitch(ctIn *rlwe.Ciphertext, combined *KeySwitchShare, ctOut *rlwe.Ciphertext) {

	level := ctIn.Level()

	if ctIn != ctOut {

		ctOut.Resize(ctIn.Degree(), level)

		ring.CopyLvl(level, ctIn.Value[1], ctOut.Value[1])

		ctOut.MetaData = ctIn.MetaData
	}

	cks.params.RingQ().AtLevel(level).Add(ctIn.Value[0], combined.Value, ctOut.Value[0])
}

// Level returns the level of the target share.
func (ckss *KeySwitchShare) Level() int {
	return ckss.Value.Level()
}

// BinarySize returns the size in bytes of the object
// when encoded using Encode.
func (ckss *KeySwitchShare) BinarySize() int {
	return ckss.Value.BinarySize()
}

// MarshalBinary encodes a KeySwitch share on a slice of bytes.
func (ckss *KeySwitchShare) MarshalBinary() (p []byte, err error) {
	buf := bytes.NewBuffer([]byte{})
	_, err = ckss.WriteTo(buf)
	return buf.Bytes(), err
}

// Encode encodes the object into a binary form on a preallocated slice of bytes
// and returns the number of bytes written.
func (ckss *KeySwitchShare) Encode(p []byte) (ptr int, err error) {
	return ckss.Value.Encode(p)
}

// WriteTo writes the object on an io.Writer.
// To ensure optimal efficiency and minimal allocations, the user is encouraged
// to provide a struct implementing the interface bufer.Writer, which defines
// a subset of the method of the bufio.Writer.
// If w is not compliant to the bufer.Writer interface, it will be wrapped in
// a new bufio.Writer.
// For additional information, see lattigo/utils/bufer/writer.go.
func (ckss *KeySwitchShare) WriteTo(w io.Writer) (n int64, err error) {
	return ckss.Value.WriteTo(w)
}

// UnmarshalBinary decodes a slice of bytes generated by
// MarshalBinary or WriteTo on the object.
func (ckss *KeySwitchShare) UnmarshalBinary(p []byte) (err error) {
	_, err = ckss.ReadFrom(bytes.NewBuffer(p))
	return
}

// Decode decodes a slice of bytes generated by Encode
// on the object and returns the number of bytes read.
func (ckss *KeySwitchShare) Decode(p []byte) (ptr int, err error) {
	if ckss.Value == nil {
		ckss.Value = new(ring.Poly)
	}

	return ckss.Value.Decode(p)
}

// ReadFrom reads on the object from an io.Writer.
// To ensure optimal efficiency and minimal allocations, the user is encouraged
// to provide a struct implementing the interface bufer.Reader, which defines
// a subset of the method of the bufio.Reader.
// If r is not compliant to the bufer.Reader interface, it will be wrapped in
// a new bufio.Reader.
// For additional information, see lattigo/utils/bufer/reader.go.
func (ckss *KeySwitchShare) ReadFrom(r io.Reader) (n int64, err error) {
	if ckss.Value == nil {
		ckss.Value = new(ring.Poly)
	}

	return ckss.Value.ReadFrom(r)
}
